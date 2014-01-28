import FWCore.ParameterSet.Config as cms
import sys, os
	
datasample = ['QCD','DYToMuMu','GluGluToHToGG','RelValTTbar','RelValZmumuJets']
	
spectra = ['15to30','30to50','50to80','80to120','120to170','170to300','300to470','470to600','600to800','800to1000']
	
	
OutfileAnlzr = "FilterFinder"	

process = cms.Process("FILTERFINDER")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-RECO.root'),
    secondaryFileNames =  cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-DIGI-RAW-HLTDEBUG.root'),
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

for i in range(len(datasample)):
	if datasample[i] in process.source.fileNames[0]:
		OutfileAnlzr+= "_" + str(datasample[i])	
		
for i in range(len(spectra)):
	if spectra[i] in process.source.fileNames[0]:
		OutfileAnlzr+= "_Pt-" + str(spectra[i])

Outfile= OutfileAnlzr+".root"	

print " Outfile set to " + Outfile

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START53_V7E::All'

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(Outfile),
    closeFileFast = cms.untracked.bool(True),
)
		
### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

### validation-specific includes
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')	
	
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
   	src = cms.InputTag('offlinePrimaryVertices'),
	cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
	filter = cms.bool(True),
)

process.demo = cms.EDAnalyzer('FilterFinder',
	TrackingParticles = cms.InputTag("mergedtruth","MergedTrackTruth"),
	TrackCollection = cms.InputTag("generalTracks"),
	VertexCollection = cms.InputTag("selectedPrimaryVertexQuality"),
	BeamSpot = cms.InputTag('offlineBeamSpot'),	
	ConversionsCollection = cms.InputTag('allConversions'),
	KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
)

process.p = cms.Path(
      process.selectedPrimaryVertexQuality *
      process.demo
)