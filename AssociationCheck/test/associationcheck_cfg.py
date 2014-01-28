import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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
		
### IVF-specific includes
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")

process.TFileService = cms.Service('TFileService',
    fileName = cms.string("twcheck.root"),
    closeFileFast = cms.untracked.bool(True),
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/user/geisler/QCD_Pt-0to5_TuneZ2star_Summer12_PU_S10_START53_AODSIM.root'
    )
)

process.demo = cms.EDAnalyzer('AssociationCheck',		 
	  TrackCollection = cms.InputTag('generalTracks'),
      VertexCollection = cms.InputTag('offlinePrimaryVertices'),
      doReassociation = cms.bool(True),
      GetCleanedCollections = cms.bool(False),
      ConversionsCollection = cms.InputTag('allConversions'),
      V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
      V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
      NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
      IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),
	  FinalAssociation = cms.untracked.int32(1),			    
      ignoreMissingCollection = cms.bool(True),			    
      nTrackWeight = cms.double(0.001),	   	    
      MaxNumberOfAssociations = cms.int32(1),	
	  BeamSpot = cms.InputTag('offlineBeamSpot'),
)


process.p = cms.Path(process.inclusiveVertexing * process.demo)
