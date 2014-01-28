import FWCore.ParameterSet.Config as cms

import sys


#_____________________
#
# READ INPUT PARAMETERS
#_____________________


# maximum number of events
MaxNumberOfEvents = 180 # reduce for testing

#if len(sys.argv) > 2: 
    #MaxNumberOfEvents = int(sys.argv[2])
	
	

process = cms.Process("DAUGHTERFINDER")

from PhysicsTools.PatAlgos.patSequences_cff import *
from TopQuarkAnalysis.Configuration.patRefSel_common_cfi import *
		
### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### Particle flow

postfix = 'PFRef'

# subtract charged hadronic pile-up particles (from wrong PVs)
# effects also JECs
usePFnoPU       = True # before any top projection
usePfIsoLessCHS = True # switch to new PF isolation with L1Fastjet CHS

# other switches for PF top projections (default: all 'True')
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = True # before MET top projection

# cuts used in top projections
# vertices
pfVertices  = 'goodOfflinePrimaryVertices'
pfD0Cut     = 0.2
pfDzCut     = 0.5
# muons
pfMuonSelectionCut = 'pt > 5.'
useMuonCutBasePF = False # use minimal (veto) muon selection cut on top of 'pfMuonSelectionCut'
pfMuonIsoConeR03 = False
pfMuonCombIsoCut = 0.2
# electrons
pfElectronSelectionCut  = 'pt > 5. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits < 2'
useElectronCutBasePF  = False # use minimal (veto) electron selection cut on top of 'pfElectronSelectionCut'
pfElectronIsoConeR03 = True
pfElectronCombIsoCut  = 0.2

### JEC levels

typeIMetCorrections = True

### Input

# list of input files
primInputFiles = cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-RECO.root')
secInputFiles = cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-DIGI-RAW-HLTDEBUG.root')

### Conditions

# maximum number of events
maxEvents = MaxNumberOfEvents # reduce for testing

# GlobalTags
globalTagMC   = 'START53_V18PR::All'

### Output

# output file
outputFile = 'Bdaughterfinder.root'

# event frequency of Fwk report
fwkReportEvery = 1

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = False


###
### Basic configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_basics_cff" )
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery

process.options.wantSummary = wantSummary
process.GlobalTag.globaltag = globalTagMC


###
### Input configuration
###
  
process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" )
process.source.fileNames = primInputFiles
process.source.secondaryFileNames = secInputFiles
#process.source.eventsToProcess = cms.untracked.VEventRange('1:3073')
process.maxEvents.input  = maxEvents


###
### Output configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_outputModule_cff" )
	
process.TFileService = cms.Service('TFileService',
    fileName = cms.string(outputFile),
    closeFileFast = cms.untracked.bool(True),
)

### Good vertex selection
process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
process.step0b = process.goodOfflinePrimaryVertices.clone( filter = True )


###
### AssociationMap configuration
###
		
### PFCandidate AssociationMap-specific includes
from CommonTools.RecoUtils.pfcand_assomap_cfi import PFCandAssoMap
		
process.PFCand2VertexAM = PFCandAssoMap.clone(
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
)
		
### PFCandidateCollection-specific includes
from CommonTools.RecoUtils.pfcand_nopu_witham_cfi import FirstVertexPFCandidates
		
process.PFCand = FirstVertexPFCandidates.clone(
	VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAM'),
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(outputFile),
    closeFileFast = cms.untracked.bool(True),
)
	
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
   	src = cms.InputTag('offlinePrimaryVertices'),
	cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
	filter = cms.bool(True),
)
		
### IVF-specific includes
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
		
### AssociationMap-specific includes
from CommonTools.RecoUtils.pf_pu_assomap_cfi import AssociationMaps, AssociationMapsJetMet
		
process.assMapWithOutIVF = AssociationMaps.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),   	    
		MaxNumberOfAssociations = cms.int32(3),
)
		
process.assMapWithIVF = AssociationMaps.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),	    
		MaxNumberOfAssociations = cms.int32(3),
		IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),   
)
						
process.assMapJetMet = AssociationMapsJetMet.clone(	
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),   	    
		MaxNumberOfAssociations = cms.int32(3),
)
	
process.cutsRecoTracks = cms.EDFilter("RecoTrackSelector",
    src = cms.InputTag("generalTracks"),
    maxChi2 = cms.double(5.0),
    tip = cms.double(0.2),
    minRapidity = cms.double(-2.4),
    lip = cms.double(17.0),
    ptMin = cms.double(1.),
    maxRapidity = cms.double(2.4),
    quality = cms.vstring('highPurity'),
    algorithm = cms.vstring(),
    minHit = cms.int32(8),
    minPixelHit = cms.int32(2), 
    min3DHit = cms.int32(0),
    beamSpot = cms.InputTag("offlineBeamSpot")
)

### GeneralTrack AssociationMap-specific includes		
from CommonTools.RecoUtils.pf_pu_firstvertextracks_cfi import *				  
				       
#without ivf
process.Q1FVTracksWithOutIVF = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapWithOutIVF'),	
        VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(1),
)				  
				       
process.Q2FVTracksWithOutIVF = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapWithOutIVF'),	
        VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(2),
)				  
				       
process.Q3FVTracksWithOutIVF = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapWithOutIVF'),	
        VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(3),
)	  
				       
#with ivf
process.Q1FVTracksWithIVF = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapWithIVF'),	
        VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(1),
)				  
				       
process.Q2FVTracksWithIVF = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapWithIVF'),	
        VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(2),
)				  
				       
process.Q3FVTracksWithIVF = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapWithIVF'),	
        VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(3),
)				  
			       
#jetmet				       
process.Q1FVTracksJetMet = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapJetMet'),	
        VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(1),
)				  
				       
process.Q2FVTracksJetMet = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapJetMet'),	
	    VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(2),
)				  
				       
process.Q3FVTracksJetMet = FirstVertexTracks.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		AssociationMap = cms.InputTag('assMapJetMet'),	
	    VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
		TrackCollection = cms.InputTag('cutsRecoTracks'),	   
		MinQuality = cms.int32(3),
)


###
### PAT/PF2PAT configuration
###

process.load( "PhysicsTools.PatAlgos.patSequences_cff" )

### Check JECs

# JEC set
jecSet = 'AK5PF'
if usePFnoPU:
  jecSet += 'chs'

# JEC levels
jecLevels = []
jecLevels.append( 'L1FastJet' )
jecLevels.append( 'L2Relative' )
jecLevels.append( 'L3Absolute' )

### Switch configuration

from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
usePF2PAT( process
         , runPF2PAT           = True
         , runOnMC             = True
         , jetAlgo             = 'AK5'
         , postfix             = postfix
         , jetCorrections      = ( jecSet
                                 , jecLevels
                                 )
         , typeIMetCorrections = typeIMetCorrections
         , pvCollection        = cms.InputTag( pfVertices )
         )

if useMuonCutBasePF:
  pfMuonSelectionCut += ' && %s'%( muonCut )
if useElectronCutBasePF:
  pfElectronSelectionCut += ' && %s'%( electronCut )

getattr( process, 'pfNoPileUp'   + postfix ).enable = usePFnoPU
getattr( process, 'pfNoMuon'     + postfix ).enable = useNoMuon
getattr( process, 'pfNoElectron' + postfix ).enable = useNoElectron
getattr( process, 'pfNoJet'      + postfix ).enable = useNoJet
getattr( process, 'pfNoTau'      + postfix ).enable = useNoTau

getattr( process, 'pfPileUpIso' + postfix ).checkClosestZVertex = usePfIsoLessCHS

getattr( process, 'pfMuonsFromVertex' + postfix ).d0Cut = pfD0Cut
getattr( process, 'pfMuonsFromVertex' + postfix ).dzCut = pfDzCut
getattr( process, 'pfSelectedMuons'   + postfix ).cut = pfMuonSelectionCut
getattr( process, 'pfIsolatedMuons'   + postfix ).isolationCut = pfMuonCombIsoCut

if pfMuonIsoConeR03:
  getattr( process, 'pfIsolatedMuons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix )
                                                                                            )
  getattr( process, 'pfIsolatedMuons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
  getattr( process, 'pfIsolatedMuons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                            , cms.InputTag( 'muPFIsoValueGamma03' + postfix )
                                                                                            )
  getattr( process, 'pfMuons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix )
                                                                                    )
  getattr( process, 'pfMuons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
  getattr( process, 'pfMuons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                    , cms.InputTag( 'muPFIsoValueGamma03' + postfix )
                                                                                    )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'muPFIsoValueChargedAll03' + postfix )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03' + postfix )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'muPFIsoValueGamma03' + postfix )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'muPFIsoValueCharged03' + postfix )

getattr( process, 'pfElectronsFromVertex' + postfix ).d0Cut = pfD0Cut
getattr( process, 'pfElectronsFromVertex' + postfix ).dzCut = pfDzCut
getattr( process, 'pfSelectedElectrons'   + postfix ).cut = pfElectronSelectionCut
getattr( process, 'pfIsolatedElectrons'   + postfix ).isolationCut = pfElectronCombIsoCut

if pfElectronIsoConeR03:
  getattr( process, 'pfIsolatedElectrons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
                                                                                                )
  getattr( process, 'pfIsolatedElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
  getattr( process, 'pfIsolatedElectrons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                                , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix )
                                                                                                )
  getattr( process, 'pfElectrons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
                                                                                        )
  getattr( process, 'pfElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
  getattr( process, 'pfElectrons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                        , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix )
                                                                                        )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )


### validation-specific includes
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
from MGeisler.TrackValidator.TrackingParticleSelection_cfi import *

TrackingParticleSelectionGeneral.ptMin = cms.double(0.8)
TrackingParticleSelectionGeneral.lip = cms.double(20.0)
TrackingParticleSelectionGeneral.tip = cms.double(0.5)
TrackingParticleSelectionGeneral.minHit = cms.int32(8)

process.bdf = cms.EDAnalyzer('ZDaughtersFinder',

    recoTrackColls = cms.VInputTag(cms.InputTag("generalTracks"), cms.InputTag("cutsRecoTracks"), cms.InputTag("Q1FVTracksWithOutIVF","T2V"), cms.InputTag("Q2FVTracksWithOutIVF","T2V"), cms.InputTag("Q3FVTracksWithOutIVF","T2V"), cms.InputTag("Q1FVTracksWithIVF","T2V"), cms.InputTag("Q2FVTracksWithIVF","T2V"), cms.InputTag("Q3FVTracksWithIVF","T2V"), cms.InputTag("Q1FVTracksJetMet","T2V"), cms.InputTag("Q2FVTracksJetMet","T2V"), cms.InputTag("Q3FVTracksJetMet","T2V")),
    
    recoJetColls = cms.VInputTag(cms.InputTag("patJetsPFRef"), cms.InputTag("patJetsPFRef"), cms.InputTag("patJetsPFNew"), cms.InputTag("patJetsPFNew"), cms.InputTag("patJetsPFNew"), cms.InputTag("patJetsPFRef"), cms.InputTag("patJetsPFRef"), cms.InputTag("patJetsPFRef"), cms.InputTag("patJetsPFRef"), cms.InputTag("patJetsPFRef"), cms.InputTag("patJetsPFRef")),
    
    pdgIds = cms.vint32(511,521,10511,10521,513,523,10513,10523,20513,20523,515,525,531,10531,533,10533,20533,535,541,10541,543,10543,20543,545, 5101,5103,5201,5203,5301,5303,5401,5403,5503,551,10551,100551,110551,200551,210551,553,10553,20553,30553,100553,110553,120553,13053,200553,210553, 220553,300553,9000553,9010553,555,10555,20555,100555,110555,120555,200555,557,100557,5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,5314, 5324,5332,5334,5142,5242,5412,5422,5414,5424,5342,5432,5434,5442,5444,5512,5522,5514,5524,5532,5534,5542,5544,5554),
    
    GenParticles = cms.InputTag("genParticles"),
    TrackingParticles = cms.InputTag("mergedtruth","MergedTrackTruth"),
    GenJets = cms.InputTag("ak5GenJets"),
    PUInfo = cms.InputTag("addPileupInfo"),
    maxDist = cms.double(0.2),
    maxDeltaR = cms.double(0.3),
    generalTpSelector = TrackingParticleSelectionGeneral,
)

# The paths

process.p = cms.Path()
process.p += process.goodOfflinePrimaryVertices
process.p += process.step0b

process.p += getattr( process, 'patPF2PATSequence' + postfix )

process.out.outputCommands =  cms.untracked.vstring('keep *')

delattr(process, 'outpath' )

newPostfix = 'PFNew'

sys.path.append('/user/geisler/CMSSW/Helpers/')
from cfg_tools import clonePath

OldNewdict = { "pfNoPileUpPFRef" : cms.InputTag("PFCand", "P2V") }
process.pNew = clonePath(process,process.p,"pfPileUpIsoPFRef",postfix,newPostfix,OldNewdict)

getattr(process,"pNew").replace(
    getattr(process,"pfPileUp"+newPostfix),
    process.PFCand2VertexAM 
)

getattr(process,"pNew").replace(
    getattr(process,"pfNoPileUp"+newPostfix),
    process.PFCand 
)

process.pNew += getattr( process, 'selectedPrimaryVertexQuality' )
process.pNew += getattr( process, 'inclusiveVertexing' )

process.pNew += getattr( process, 'assMapWithOutIVF' )
process.pNew += getattr( process, 'assMapWithIVF' )
process.pNew += getattr( process, 'assMapJetMet' )

process.pNew += getattr( process, 'cutsRecoTracks' )

process.pNew += getattr( process, 'Q1FVTracksWithOutIVF' )
process.pNew += getattr( process, 'Q2FVTracksWithOutIVF' )
process.pNew += getattr( process, 'Q3FVTracksWithOutIVF' )

process.pNew += getattr( process, 'Q1FVTracksWithIVF' )
process.pNew += getattr( process, 'Q2FVTracksWithIVF' )
process.pNew += getattr( process, 'Q3FVTracksWithIVF' )

process.pNew += getattr( process, 'Q1FVTracksJetMet' )
process.pNew += getattr( process, 'Q2FVTracksJetMet' )
process.pNew += getattr( process, 'Q3FVTracksJetMet' )

process.pNew += getattr( process, 'bdf' )
