import FWCore.ParameterSet.Config as cms

import sys


#_____________________
#
# READ INPUT PARAMETERS
#_____________________


# maximum number of events
MaxNumberOfEvents = 250 # reduce for testing

#if len(sys.argv) > 2: 
    #MaxNumberOfEvents = int(sys.argv[2])
	
	

process = cms.Process("JETVALIDATOR")

from PhysicsTools.PatAlgos.patSequences_cff import *
from TopQuarkAnalysis.Configuration.patRefSel_common_cfi import *

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
useRelVals = False # if 'False', "inputFiles" is used
#inputFiles = ["file:/user/geisler/QCD_Pt-0to5_TuneZ2star_Summer12_PU_S10_START53_AODSIM.root"] # overwritten, if "useRelVals" is 'True'
inputFiles = ["file:/user/geisler/RelValZmumuJets_Pt_20_300_CMSSW_5_3_4_cand1_PU_START53_GEN-SIM-RECO.root"] # overwritten, if "useRelVals" is 'True'

### Conditions

# maximum number of events
maxEvents = MaxNumberOfEvents # reduce for testing

# GlobalTags
globalTagMC   = 'START53_V18PR::All'

### Output

# output file
outputFile = 'jetmetvalidator.root'

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
process.source.fileNames = inputFiles
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
		
### IVF-specific includes
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
		
### PFCandidate AssociationMap-specific includes
from CommonTools.RecoUtils.pfcand_assomap_cfi import PFCandAssoMap,PFCandAssoMapJetMet
		
process.PFCand2VertexAMQ3 = PFCandAssoMap.clone(
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'), 
	MaxNumberOfAssociations = cms.int32(1),	
)
		
process.PFCand2VertexAMQ2 = PFCandAssoMap.clone(
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'), 
	MaxNumberOfAssociations = cms.int32(2),	
)
		
process.PFCand2VertexAMQ1 = PFCandAssoMap.clone(
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'), 
	MaxNumberOfAssociations = cms.int32(2),	
)

process.PFCand2VertexAMJetMET = PFCandAssoMapJetMet.clone(
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	IVFVertexCollection = cms.InputTag(''),
)
		
### PFCandidateCollection-specific includes
from CommonTools.RecoUtils.pfcand_nopu_witham_cfi import FirstVertexPFCandidates
		
process.PFCandQ3 = FirstVertexPFCandidates.clone(
	VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAMQ3'),
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	MinQuality = cms.int32(3),
)
		
process.PFCandQ2 = FirstVertexPFCandidates.clone(
	VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAMQ2'),
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	MinQuality = cms.int32(2),
)
		
process.PFCandQ1 = FirstVertexPFCandidates.clone(
	VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAMQ1'),
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	MinQuality = cms.int32(1),
)

process.PFCandJM = FirstVertexPFCandidates.clone(
	VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAMJetMET'),
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	MinQuality = cms.int32(3),
)

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

getattr( process, 'pfMETcorrType0' + postfix ).quality = cms.int32(3)
getattr( process, 'patPFMETtype0Corr' + postfix ).quality = cms.int32(3)

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

		
patType0CorrectedPFMet = cms.EDProducer("CorrectedPATMETProducer",
	src = cms.InputTag('patPFMet' + postfix ),
	applyType0Corrections = cms.bool(False),
	applyType1Corrections = cms.bool(True),
	srcType1Corrections = cms.VInputTag(
		cms.InputTag('patPFMETtype0Corr' + postfix ) 
	),
	applyType2Corrections = cms.bool(False),
)

setattr(process, "patType0CorrectedPFMet" + postfix, patType0CorrectedPFMet)

		
process.jetvalidator = cms.EDAnalyzer('JetValidator',
    genJetLabel = cms.InputTag("ak5GenJets"),
    recoJetLabels = cms.VInputTag(cms.InputTag("patJetsPFRef"), cms.InputTag("patJetsPFNewQ1"), cms.InputTag("patJetsPFNewQ2"), cms.InputTag("patJetsPFNewQ3")),
    PULabel = cms.InputTag("addPileupInfo"),
    UseLogPt = cms.bool(False),
    MaxNumberOfjetsPerEvent = cms.uint32(2),
)

		
process.metvalidator = cms.EDAnalyzer('METValidator',
	corrMETLabels = cms.VInputTag( cms.InputTag('patType0CorrectedPFMetPFRef'), cms.InputTag('patType0CorrectedPFMetPFNewQ3'), cms.InputTag('patType0CorrectedPFMetPFNewQ2'), cms.InputTag('patType0CorrectedPFMetPFNewQ1'), cms.InputTag('patType0CorrectedPFMetPFNewJM') ),
	rawMETLabels = cms.VInputTag( cms.InputTag('patPFMetPFRef'), cms.InputTag('patPFMetPFNewQ3'), cms.InputTag('patPFMetPFNewQ2'), cms.InputTag('patPFMetPFNewQ1'), cms.InputTag('patPFMetPFNewJM') ),
    PULabel = cms.InputTag("addPileupInfo"),
)

# The paths

process.p = cms.Path()
process.p += process.goodOfflinePrimaryVertices
process.p += process.inclusiveVertexing
process.p += process.step0b

process.p += getattr( process, 'patPF2PATSequence' + postfix )

process.p += getattr( process, 'patType0CorrectedPFMet' + postfix )

process.out.outputCommands =  cms.untracked.vstring('keep *')

delattr(process, 'outpath' )

sys.path.append('/user/geisler/CMSSW/Helpers/')
from cfg_tools import clonePath

newPostfixQ3 = 'PFNewQ3'

OldNewdictQ3 = { "pfNoPileUp" + postfix : cms.InputTag("PFCandQ3", "P2V"), "pfCandidateToVertexAssociation" + postfix : cms.InputTag("PFCand2VertexAMQ3") }
process.pNewQ3 = clonePath(process,process.p,"pfPileUpIsoPFRef",postfix,newPostfixQ3,OldNewdictQ3)

#delattr(process, 'pfCandidateToVertexAssociation' + newPostfixQ3 )

getattr(process,"pNewQ3").replace(
    getattr(process,"pfPileUp"+newPostfixQ3),
    process.PFCand2VertexAMQ3 
)

getattr(process,"pNewQ3").replace(
    getattr(process,"pfNoPileUp"+newPostfixQ3),
    process.PFCandQ3
)

getattr(process, "patPFMETtype0Corr" + newPostfixQ3).quality = cms.int32(3)

newPostfixQ2 = 'PFNewQ2'

OldNewdictQ2 = { "pfNoPileUp" + postfix : cms.InputTag("PFCandQ2", "P2V"), "pfCandidateToVertexAssociation" + postfix : cms.InputTag("PFCand2VertexAMQ2") }
process.pNewQ2 = clonePath(process,process.p,"pfPileUpIsoPFRef",postfix,newPostfixQ2,OldNewdictQ2)

#delattr(process, 'pfCandidateToVertexAssociation' + newPostfixQ2 )

getattr(process,"pNewQ2").replace(
    getattr(process,"pfPileUp"+newPostfixQ2),
    process.PFCand2VertexAMQ2 
)

getattr(process,"pNewQ2").replace(
    getattr(process,"pfNoPileUp"+newPostfixQ2),
    process.PFCandQ2
)

getattr(process, "patPFMETtype0Corr" + newPostfixQ2).quality = cms.int32(2)

newPostfixQ1 = 'PFNewQ1'

OldNewdictQ1 = { "pfNoPileUp" + postfix : cms.InputTag("PFCandQ1", "P2V"), "pfCandidateToVertexAssociation" + postfix : cms.InputTag("PFCand2VertexAMQ1") }
process.pNewQ1 = clonePath(process,process.p,"pfPileUpIsoPFRef",postfix,newPostfixQ1,OldNewdictQ1)

#delattr(process, 'pfCandidateToVertexAssociation' + newPostfixQ1 )

getattr(process,"pNewQ1").replace(
    getattr(process,"pfPileUp"+newPostfixQ1),
    process.PFCand2VertexAMQ1 
)

getattr(process,"pNewQ1").replace(
    getattr(process,"pfNoPileUp"+newPostfixQ1),
    process.PFCandQ1
)

getattr(process, "patPFMETtype0Corr" + newPostfixQ1).quality = cms.int32(1)

newPostfixJM = 'PFNewJM'

OldNewdictJM = { "pfNoPileUp" + postfix : cms.InputTag("PFCandJM", "P2V"), "pfCandidateToVertexAssociation" + postfix : cms.InputTag("PFCand2VertexAMJetMET") }
process.pNewJM = clonePath(process,process.p,"pfPileUpIsoPFRef",postfix,newPostfixJM,OldNewdictJM)

#delattr(process, 'pfCandidateToVertexAssociation' + newPostfixJM )

getattr(process,"pNewJM").replace(
    getattr(process,"pfPileUp"+newPostfixJM),
    process.PFCand2VertexAMJetMET 
)

getattr(process,"pNewJM").replace(
    getattr(process,"pfNoPileUp"+newPostfixJM),
    process.PFCandJM
)

getattr(process, "patPFMETtype0Corr" + newPostfixJM).quality = cms.int32(3)

process.pNewJM += getattr( process, 'jetvalidator' )
process.pNewJM += getattr( process, 'metvalidator' )
