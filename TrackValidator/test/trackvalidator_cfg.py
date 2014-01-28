import FWCore.ParameterSet.Config as cms

process = cms.Process("TRACKVALIDATOR")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-RECO.root'),
    secondaryFileNames =  cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-DIGI-RAW-HLTDEBUG.root'),
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

Outfile= "TrackValidator.root"	

print " Outfile set to " + Outfile

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START53_V26::All'
		
### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### validation-specific includes
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
process.load("MGeisler.TrackValidator.cutsRecoTracks_cfi")
from MGeisler.TrackValidator.TrackingParticleSelection_cfi import *

########### track selection configuration ########  
process.cutsRecoTracks.ptMin = cms.double(1.0)

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
		
process.out = cms.OutputModule("PoolOutputModule",
	outputCommands = cms.untracked.vstring('drop *'),
	fileName = cms.untracked.string('EmptyFile.root')
)
		

### Part of the AssociationMap
from CommonTools.RecoUtils.pf_pu_assomap_cfi import AssociationMaps	
		
process.assMap = AssociationMaps.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),   
)		
		
from CommonTools.RecoUtils.pf_pu_firstvertextracks_cfi import *	

process.amTracks = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMap'),	
        VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(3),
)
		
		  
### Part of the JetMet technique	

process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix="JetMet"
usePF2PAT(process,runPF2PAT=True, jetAlgo="AK5", runOnMC=True, postfix=postfix)
process.pfPileUpJetMet.Vertices = cms.InputTag('selectedPrimaryVertexQuality')
process.pfPileUpJetMet.checkClosestZVertex = cms.bool(False)
		
process.pfPUMuonJetMet = cms.EDProducer("TPPFCandidatesOnPFCandidates",
		bottomCollection = cms.InputTag("pfAllMuonsJetMet"),
		topCollection = cms.InputTag("pfMuonsFromVertexJetMet"),
		enable = cms.bool(True),
		name = cms.untracked.string("PUMuonJetMet"),
		verbose = cms.untracked.bool(False),		 
)
		
process.patPF2PATSequenceJetMet.replace(
	getattr(process,"pfMuonsFromVertexJetMet"),
	process.pfMuonsFromVertexJetMet * process.pfPUMuonJetMet 
)
		
process.pfNoMuonJetMet.topCollection = cms.InputTag("pfPUMuonJetMet", "", "")
		
		  
		
process.pfPUElectronJetMet = cms.EDProducer("TPPFCandidatesOnPFCandidates",
		bottomCollection = cms.InputTag("pfAllElectronsJetMet"),
		topCollection = cms.InputTag("pfElectronsFromVertexJetMet"),
		enable = cms.bool(True),
		name = cms.untracked.string("PUElectronJetMet"),
		verbose = cms.untracked.bool(False),		 
)
		
process.patPF2PATSequenceJetMet.replace(
	getattr(process,"pfElectronsFromVertexJetMet"),
	process.pfElectronsFromVertexJetMet * process.pfPUElectronJetMet 
)
		
process.pfNoElectronJetMet.topCollection = cms.InputTag("pfPUElectronJetMet", "", "")
		
process.jmTracks = cms.EDProducer("FromParticlesToTracks",
		ParticleCollection = cms.InputTag('pfNoElectronJetMet'),
		RefTrackCollection = cms.InputTag('cutsRecoTracks'),	
)
		
		  
### Part of the MuEg technique	
		
postfix="MuEg"
usePF2PAT(process,runPF2PAT=True, jetAlgo="AK5", runOnMC=True, postfix=postfix)
process.pfPileUpMuEg.Vertices = cms.InputTag('selectedPrimaryVertexQuality')
process.pfPileUpMuEg.checkClosestZVertex = cms.bool(True)
		
process.pfPUMuonMuEg = cms.EDProducer("TPPFCandidatesOnPFCandidates",
		bottomCollection = cms.InputTag("pfAllMuonsMuEg"),
		topCollection = cms.InputTag("pfMuonsFromVertexMuEg"),
		enable = cms.bool(True),
		name = cms.untracked.string("PUMuonMuEg"),
		verbose = cms.untracked.bool(False),		 
)
		
process.patPF2PATSequenceMuEg.replace(
	getattr(process,"pfMuonsFromVertexMuEg"),
	process.pfMuonsFromVertexMuEg * process.pfPUMuonMuEg 
)
		
process.pfNoMuonMuEg.topCollection = cms.InputTag("pfPUMuonMuEg", "", "")
		
		  
		
process.pfPUElectronMuEg = cms.EDProducer("TPPFCandidatesOnPFCandidates",
		bottomCollection = cms.InputTag("pfAllElectronsMuEg"),
		topCollection = cms.InputTag("pfElectronsFromVertexMuEg"),
		enable = cms.bool(True),
		name = cms.untracked.string("PUElectronMuEg"),
		verbose = cms.untracked.bool(False),		 
)
		
process.patPF2PATSequenceMuEg.replace(
	getattr(process,"pfElectronsFromVertexMuEg"),
	process.pfElectronsFromVertexMuEg * process.pfPUElectronMuEg 
)
		
process.pfNoElectronMuEg.topCollection = cms.InputTag("pfPUElectronMuEg", "", "")
		
process.meTracks = cms.EDProducer("FromParticlesToTracks",
		ParticleCollection = cms.InputTag('pfNoElectronMuEg'),
		RefTrackCollection = cms.InputTag('cutsRecoTracks'),	
)


### Part of the Validation
				
process.trackvalidator = cms.EDAnalyzer('TrackValidator',
		tcLabel = cms.VInputTag(cms.InputTag("cutsRecoTracks"), cms.InputTag("amTracks","T2V"), cms.InputTag("jmTracks"), cms.InputTag("meTracks")),
    	tcRefLabel = cms.InputTag("generalTracks"),
    	PULabel = cms.InputTag("addPileupInfo"),
    	TPLabel = cms.InputTag("mergedtruth","MergedTrackTruth"),
		ignoremissingtrackcollection=cms.bool(False),
		UseLogPt=cms.bool(False),
		generalTpSelector = TrackingParticleSelectionGeneral,
		useJetWeighting = cms.untracked.bool(False),
		genJetCollLabel = cms.untracked.string("ak5GenJets"),
		jetCollLabel = cms.untracked.string("ak5PFJets"),
)			

### sequences and paths

process.startingstuff = cms.Sequence(
		process.cutsRecoTracks *
		process.inclusiveVertexing *
		process.selectedPrimaryVertexQuality
)

process.assomaps = cms.Sequence(
    	process.assMap *
    	process.amTracks
)

process.jetmet = cms.Sequence(
    	process.patPF2PATSequenceJetMet *
		process.jmTracks
)

process.mueg = cms.Sequence(
    	process.patPF2PATSequenceMuEg *
		process.meTracks
)

process.validation = cms.Sequence(
      	process.trackvalidator 
)

process.p = cms.Path(
		process.startingstuff *
		process.assomaps *
		process.jetmet *
		process.mueg *
		process.validation
)
