import FWCore.ParameterSet.Config as cms
import sys, os

process = cms.Process("FINALASSOCTEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-RECO.root'),
    secondaryFileNames =  cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-DIGI-RAW-HLTDEBUG.root'),
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

Outfile = "FinalAssociationAnalysis.root"	

print " Outfile set to " + Outfile

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START53_V7E::All'
		
### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### gen-reco association specific includes
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(Outfile),
    closeFileFast = cms.untracked.bool(True),
)	
	
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
   	src = cms.InputTag('offlinePrimaryVertices'),
	cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
	filter = cms.bool(True),
)
		
### IVF-specific includes
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
				
process.demo = cms.EDAnalyzer('FinalAssocTest',
    TrackCollection = cms.InputTag('generalTracks'),
    VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
	TrackingParticles = cms.InputTag("mergedtruth","MergedTrackTruth"),
	BeamSpot = cms.InputTag('offlineBeamSpot'),
    GetCleanedCollections = cms.bool(False),	
	ConversionsCollection = cms.InputTag('allConversions'),
	KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
	IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),
)


process.p = cms.Path(process.selectedPrimaryVertexQuality * process.inclusiveVertexing * process.demo)
