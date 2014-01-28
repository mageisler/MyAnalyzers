import FWCore.ParameterSet.Config as cms
import sys, os


#_____________________
#
# READ INPUT PARAMETERS
#_____________________


# maximum number of events
MaxNumberOfEvents = cms.untracked.int32(-1) # reduce for testing

#if len(sys.argv) > 2: 
    #MaxNumberOfEvents = cms.untracked.int32(int(sys.argv[2]))
	

Outfile = "metvalidator.root"

process = cms.Process("METVALIDATOR")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:/user/geisler/QCD_Pt-0to5_TuneZ2star_Summer12_PU_S10_START53_AODSIM.root"),
)
		
process.maxEvents = cms.untracked.PSet(
    input = MaxNumberOfEvents
)	

print " Outfile set to " + Outfile

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START53_V18PR::All'
		
### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(Outfile),
    closeFileFast = cms.untracked.bool(True),
)
	
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
   	src = cms.InputTag('offlinePrimaryVertices'),
	cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
	filter = cms.bool(True),
)

### PFMET-Type0-specific includes
from JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi import *
	
process.myfirstPrimaryVertexQuality = selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0.clone(
          vertices = cms.InputTag('selectedPrimaryVertexQuality'),
)
		
### IVF-specific includes
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")

###
### AssociationMap configuration
###
		
### PFCandidate AssociationMap-specific includes
from CommonTools.RecoUtils.pfcand_assomap_cfi import *
		
process.PFCand2VertexDefault = PFCandAssoMap.clone(
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		AssociationType = cms.InputTag('PFCandsToVertex'),	    
		MaxNumberOfAssociations = cms.int32(2),
		IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'), 
)
		
process.PFCand2VertexJetMet = PFCandAssoMapJetMet.clone(
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		AssociationType = cms.InputTag('PFCandsToVertex'),
		MaxNumberOfAssociations = cms.int32(1),
		IVFVertexCollection = cms.InputTag(''), 
)

process.CorrPFMETDefaultQ1 = pfMETcorrType0.clone(
	srcPFCandidateToVertexAssociations = cms.InputTag('PFCand2VertexDefault'),	
	srcHardScatterVertex = cms.InputTag('myfirstPrimaryVertexQuality'),
	quality = cms.int32(1),
)

process.CorrPFMETDefaultQ2 = pfMETcorrType0.clone(
	srcPFCandidateToVertexAssociations = cms.InputTag('PFCand2VertexDefault'),	
	srcHardScatterVertex = cms.InputTag('myfirstPrimaryVertexQuality'),
	quality = cms.int32(2),
)

process.CorrPFMETDefaultQ3 = pfMETcorrType0.clone(
	srcPFCandidateToVertexAssociations = cms.InputTag('PFCand2VertexDefault'),	
	srcHardScatterVertex = cms.InputTag('myfirstPrimaryVertexQuality'),
	quality = cms.int32(3),
)
		
process.CorrPFMETJetMet = pfMETcorrType0.clone(
	srcPFCandidateToVertexAssociations = cms.InputTag('PFCand2VertexJetMet'),	
	srcHardScatterVertex = cms.InputTag('myfirstPrimaryVertexQuality'),
	quality = cms.int32(1),
)
		
### validation
				
process.metvalidator = cms.EDAnalyzer('METValidator',
	recoMETLabels = cms.VInputTag(cms.InputTag("CorrPFMETDefaultQ1"), cms.InputTag("CorrPFMETDefaultQ2"), cms.InputTag("CorrPFMETDefaultQ3"), cms.InputTag("CorrPFMETJetMet")),
    PULabel = cms.InputTag("addPileupInfo"),
)

# Output definition

process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands =  cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('mv_output.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    )
)

### sequences and paths

process.startingstuff = cms.Sequence(
	process.selectedPrimaryVertexQuality *
	process.myfirstPrimaryVertexQuality *
	process.inclusiveVertexing
)

process.assomaps = cms.Sequence(
    	process.PFCand2VertexDefault *
    	process.PFCand2VertexJetMet 
)

process.corrMETs = cms.Sequence(
    	process.CorrPFMETDefaultQ1 *
    	process.CorrPFMETDefaultQ2 *
    	process.CorrPFMETDefaultQ3 *
    	process.CorrPFMETJetMet 
)

process.validation = cms.Sequence(
      	process.metvalidator 
)

process.p = cms.Path(
	process.startingstuff *
	process.assomaps *
	process.corrMETs *
	process.validation
)

process.pout = cms.EndPath(
	#process.output
)
