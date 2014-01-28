import FWCore.ParameterSet.Config as cms
import sys, os

OutfileAnlzr = "AssociationQualityAnlzr"

process = cms.Process("ASSOCQUALANLZR")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-RECO.root'),
    secondaryFileNames =  cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-DIGI-RAW-HLTDEBUG.root'),
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(190)
)

Outfile= "QualityAnalzr.root"	
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
		
### IVF-specific includes
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")

### validation-specific includes
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(Outfile),
    closeFileFast = cms.untracked.bool(True),
)

process.qualanlzrFV = cms.EDAnalyzer('QualityAnlzr',			 
	  TrackCollection = cms.InputTag('generalTracks'),
      VertexCollection = cms.InputTag('offlinePrimaryVertices'),
      doReassociation = cms.bool(True),
      GetCleanedCollections = cms.bool(False),
      ConversionsCollection = cms.InputTag('allConversions'),
      V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
      V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
      NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
      IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),
	  FinalAssociation = cms.untracked.int32(0),			    
      ignoreMissingCollection = cms.bool(True),			    
      nTrackWeight_z = cms.double(0.001),	  		    
      nTrackWeight_3D = cms.double(0.01),	   	    
      MaxNumberOfAssociations = cms.int32(3),	
	  BeamSpot = cms.InputTag('offlineBeamSpot'),
)

process.qualanlzrZ = cms.EDAnalyzer('QualityAnlzr',			 
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
      nTrackWeight_z = cms.double(0.001),	  		    
      nTrackWeight_3D = cms.double(0.01),	   	    
      MaxNumberOfAssociations = cms.int32(3),	
	  BeamSpot = cms.InputTag('offlineBeamSpot'),
)

process.qualanlzr3D = cms.EDAnalyzer('QualityAnlzr',			 
	  TrackCollection = cms.InputTag('generalTracks'),
      VertexCollection = cms.InputTag('offlinePrimaryVertices'),
      doReassociation = cms.bool(True),
      GetCleanedCollections = cms.bool(False),
      ConversionsCollection = cms.InputTag('allConversions'),
      V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
      V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
      NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
      IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),
	  FinalAssociation = cms.untracked.int32(2),			    
      ignoreMissingCollection = cms.bool(True),			    
      nTrackWeight_z = cms.double(0.001),	  		    
      nTrackWeight_3D = cms.double(0.01),	   	    
      MaxNumberOfAssociations = cms.int32(3),	
	  BeamSpot = cms.InputTag('offlineBeamSpot'),
)


process.p = cms.Path(
		process.inclusiveVertexing * 
		process.qualanlzrFV * 
		process.qualanlzrZ * 
		process.qualanlzr3D )
