import FWCore.ParameterSet.Config as cms

import sys

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')

options.register ('allinOne',
                  "",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,  
                  'Parsing all options in one')

options.parseArguments()
isData = False

if options.allinOne != "" :
	opt = options.allinOne.split(",")
	isData = bool(int(opt[0]))

print ""		
print " Is This Data?", isData


#_____________________
#
# READ INPUT PARAMETERS
#_____________________


# maximum number of events
MaxNumberOfEvents = 10 # reduce for testing

#if len(sys.argv) > 2: 
    #MaxNumberOfEvents = int(sys.argv[2])
	
	

process = cms.Process("METVALIDATOR")

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

### Conditions

# maximum number of events
maxEvents = MaxNumberOfEvents # reduce for testing

# GlobalTags

tag = 'START53_V18PR::All'
addOn = 'MC'
if isData:
	tag =  'GR_P_V42_AN4::All'
	addOn = 'Data'
	
print " Globaltag set to ", tag
print ""

### Output

# output file
outputFile = 'MetValidatorData_' + addOn + '.root'

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
process.GlobalTag.globaltag = tag


###
### Input configuration
###
  
process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" )
process.maxEvents.input  = maxEvents

try:
	import json
except:
	try:
		import simplejson as json
	except:
		print "Please use lxplus or set an environment (for example crab) with json lib available"
		sys.exit(1)
		
good_lumis = []
json_file = '/user/geisler/LumiCalc/goldenJSON_selected.txt'
if json_file is not None and json_file != '':
	jsonfile=file(json_file, 'r')
	jsondict = json.load(jsonfile)
	runs = jsondict.keys()
	runs.sort()
  	for run in runs:
		blocks = jsondict[run]
		blocks.sort()
		prevblock = [-2,-2]
		for lsrange in blocks:
			if lsrange[0] == prevblock[1]+1:
				prevblock[1] = lsrange[1]
				good_lumis[-1] = str("%s:%s-%s:%s" % (run, prevblock[0], run, prevblock[1]))
			else:
				good_lumis.append(str("%s:%s-%s:%s" % (run, lsrange[0], run, lsrange[1])))
			prevblock = lsrange

from MGeisler.TrackValidatorData.inputfiles import inputfilesDATA, inputfilesMC 

if isData:
	process.source.fileNames = inputfilesDATA
	process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*good_lumis)
else:
	process.source.fileNames = inputfilesMC


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
process.goodOfflinePrimaryVertices.filter = True
		

process.eventweight = cms.EDProducer("Reweight",
		recoVertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
)

	
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
	
process.doubleMuTrigger = hltHighLevel.clone(
	HLTPaths = cms.vstring( 'HLT_Mu17_Mu8_v*' ),
)
	
process.theTwoMuons = cms.EDProducer('MuonSelection',
	Muons = cms.InputTag("muons"),
	VertexCollection = cms.InputTag("goodOfflinePrimaryVertices"),
	minLayers = cms.int32(10),
	minPt = cms.double(20.),
	maxEta = cms.double(2.1),
)

process.muonFilter = cms.EDFilter('MuonFilter',
	Muons = cms.InputTag("theTwoMuons"),
	isolationFactor = cms.double(0.1),
	massCut = cms.double(30.),
)


###
### AssociationMap configuration
###
		
### IVF-specific includes
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
		
### PFCandidate AssociationMap-specific includes
from CommonTools.RecoUtils.pfcand_assomap_cfi import PFCandAssoMap,PFCandAssoMapJetMet
		
### PFCandidateCollection-specific includes
from CommonTools.RecoUtils.pfcand_nopu_witham_cfi import FirstVertexPFCandidates

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

getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patElectrons' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patMuons' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patPFCandidateIsoDepositSelection' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patPFTauIsolation' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patTaus' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patJetCorrections' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'jetTracksAssociatorAtVertex' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'impactParameterTagInfosAOD' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'secondaryVertexTagInfosAOD' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'softMuonTagInfosAOD' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'secondaryVertexNegativeTagInfosAOD' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'inclusiveVertexing' ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'inclusiveSecondaryVertexFinderTagInfosAOD' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'softElectronCands' ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'softElectronTagInfosAOD' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'btaggingJetTagsAOD' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patHPSPFTauDiscrimination' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'secondaryVertexNegativeTagInfosAOD' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patJetCharge' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'genParticlesForJetsNoNu' ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patJetFlavourId' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patJets' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patPFParticles' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patCandidateSummary' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'selectedPatElectrons' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'selectedPatMuons' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'selectedPatTaus' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'selectedPatJets' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'patMETs' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'selectedPatPFParticles' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'selectedPatCandidateSummary' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'countPatElectrons' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'countPatMuons' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'countPatTaus' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'countPatLeptons' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'countPatJets' + postfix ) )
getattr( process, 'patPF2PATSequence' + postfix ).remove( getattr( process, 'countPatPFParticles' + postfix ) )

getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'pfCandsNotInJet' + postfix ) )
getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'selectedPatJetsForMETtype1p2Corr' + postfix ) )
getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'selectedPatJetsForMETtype2Corr' + postfix ) )
getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'patPFJetMETtype1p2Corr' + postfix ) )
getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'patPFJetMETtype2Corr' + postfix ) )
#getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'type0PFMEtCorrection' + postfix ) )
#getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'patPFMETtype0Corr' + postfix ) )
getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'pfCandMETcorr' + postfix ) )
getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'patType1CorrectedPFMet' + postfix ) )
getattr( process, 'producePatPFMETCorrections' + postfix ).remove( getattr( process, 'patType1p2CorrectedPFMet' + postfix ) )

getattr( process, 'patPFMETtype0Corr' + postfix ).quality = cms.int32(3)
getattr( process, 'pfMETcorrType0' + postfix ).quality = cms.int32(3)
getattr( process, 'patPFMet' + postfix ).addGenMET = cms.bool(False)

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
                                                                                    , cms.InputTag( 'muPFIsoValueGamma03' + postfix ) )

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
                                                                                        , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix ) )
  
  

		
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

		
process.metvalidator = cms.EDAnalyzer('METValidator',
	corrMETLabels = cms.VInputTag( cms.InputTag('patType0CorrectedPFMetPFRef'), cms.InputTag('patType0CorrectedPFMetPFNewQ3'), cms.InputTag('patType0CorrectedPFMetPFNewQ2'), cms.InputTag('patType0CorrectedPFMetPFNewQ1'), cms.InputTag('patType0CorrectedPFMetPFJM') ),
	rawMETLabels = cms.VInputTag( cms.InputTag('patPFMetPFRef'), cms.InputTag('patPFMetPFNewQ3'), cms.InputTag('patPFMetPFNewQ2'), cms.InputTag('patPFMetPFNewQ1'), cms.InputTag('patPFMetPFJM') ),
    PULabel = cms.InputTag("addPileupInfo"),
	isData = cms.bool(isData),
	usePUInfo = cms.bool(False),
	Weight = cms.InputTag('eventweight'),
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
)

# The paths

process.p = cms.Path()
process.p += process.goodOfflinePrimaryVertices
process.p += process.doubleMuTrigger 
process.p += process.theTwoMuons 
process.p += process.muonFilter


if not isData:
	process.p += process.eventweight

process.p += process.inclusiveVertexing


for mod in getattr( process,  'patPF2PATSequence' + postfix ).moduleNames():
	if "atch" in mod or "Gen" in mod:
		getattr( process,  'patPF2PATSequence' + postfix ).remove( getattr( process, mod ) )

process.p += getattr( process, 'patPF2PATSequence' + postfix )

process.p += getattr( process, 'patType0CorrectedPFMet' + postfix )

process.out.outputCommands =  cms.untracked.vstring('keep *')

#delattr(process, 'outpath' )

sys.path.append('/user/geisler/CMSSW/Helpers/')
from cfg_tools import clonePath


## quality == 3

newPostfixQ3 = 'PFNewQ3'
		
process.PFCand2VertexAMPFNewQ3 = PFCandAssoMap.clone(
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'), 
	MaxNumberOfAssociations = cms.int32(1),	
	ignoreMissingCollection = cms.bool(False),	
)
		
process.pfNoPileUpAMPFNewQ3 = FirstVertexPFCandidates.clone(
	VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAMPFNewQ3'),
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	MinQuality = cms.int32(3),
)
		
process.pfPileUpAMPFNewQ3 = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfNoPileUpAMPFNewQ3", "P2V"),
    name = cms.untracked.string('pileUpOnPFCandidates'+newPostfixQ3),
    verbose = cms.untracked.bool(False)
)

process.myPUSequencePFNewQ3 = cms.Sequence(process.PFCand2VertexAMPFNewQ3*process.pfNoPileUpAMPFNewQ3*process.pfPileUpAMPFNewQ3)


OldNewdictQ3 = { "pfNoPileUp" + postfix : cms.InputTag("pfNoPileUpAM" + newPostfixQ3, "P2V"),
                 "pfNoPileUpIso" + postfix : cms.InputTag("pfNoPileUpAM" + newPostfixQ3, "P2V"),
				 "pfPileUp" + postfix : cms.InputTag("pfPileUpAM" + newPostfixQ3),
                 "pfPileUpIso" + postfix : cms.InputTag("pfPileUpAM" + newPostfixQ3),
				 "pfCandidateToVertexAssociation" + postfix : cms.InputTag('PFCand2VertexAM' + newPostfixQ3) }
				 
process.pNewQ3 = clonePath(process,process.p,"pfAllNeutralHadrons"+postfix,postfix,newPostfixQ3,OldNewdictQ3)

getattr( process, 'pNewQ3' ).remove( getattr( process, 'particleFlowDisplacedVertex' + newPostfixQ3 ) )
getattr( process, 'pNewQ3' ).remove( getattr( process, 'pfCandidateToVertexAssociation' + newPostfixQ3 ) )

getattr( process, 'pNewQ3' ).remove( getattr( process, 'pfPileUp' + newPostfixQ3 ) )
getattr( process, 'pNewQ3' ).remove( getattr( process, 'pfNoPileUp' + newPostfixQ3 ) )

getattr( process, 'pNewQ3' ).remove( getattr( process, 'pfNoPileUpIso' + newPostfixQ3 ) )

getattr(process,"pNewQ3").replace(
    getattr(process,"pfPileUpIso"+newPostfixQ3),
    process.myPUSequencePFNewQ3 
)

getattr(process, "patPFMETtype0Corr" + newPostfixQ3).quality = cms.int32(3)

## quality == 2

newPostfixQ2 = 'PFNewQ2'

OldNewdictQ2 = { "pfNoPileUpAM" + newPostfixQ3 : cms.InputTag("pfNoPileUpAM" + newPostfixQ2, "P2V"),
                 "pfNoPileUpIsoAM" + newPostfixQ3 : cms.InputTag("pfNoPileUpAM" + newPostfixQ2, "P2V"),
				 "pfPileUpAM" + newPostfixQ3 : cms.InputTag("pfPileUpAM" + newPostfixQ2),
                 "pfPileUpIsoAM" + newPostfixQ3 : cms.InputTag("pfPileUpAM" + newPostfixQ2),
				 "pfCandidateToVertexAssociation" + newPostfixQ3 : cms.InputTag('PFCand2VertexAM' + newPostfixQ2) }
				 
process.pNewQ2 = clonePath(process,process.pNewQ3,"PFCand2VertexAM"+newPostfixQ3,newPostfixQ3,newPostfixQ2,OldNewdictQ2)

getattr(process, "PFCand2VertexAM" + newPostfixQ2).MaxNumberOfAssociations = cms.int32(2)

getattr(process, "pfNoPileUpAM" + newPostfixQ2).VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAM' + newPostfixQ2)
getattr(process, "pfNoPileUpAM" + newPostfixQ2).MinQuality = cms.int32(2)

getattr(process, "pfPileUpAM" + newPostfixQ2).topCollection = cms.InputTag("pfNoPileUpAM" + newPostfixQ2, "P2V")
getattr(process, "pfPileUpAM" + newPostfixQ2).name = cms.untracked.string('pileUpOnPFCandidates'+newPostfixQ2)

getattr(process, "patPFMETtype0Corr" + newPostfixQ2).quality = cms.int32(2)


## quality == 1

newPostfixQ1 = 'PFNewQ1'

OldNewdictQ1 = { "pfNoPileUpAM" + newPostfixQ3 : cms.InputTag("pfNoPileUpAM" + newPostfixQ1, "P2V"),
                 "pfNoPileUpIsoAM" + newPostfixQ3 : cms.InputTag("pfNoPileUpAM" + newPostfixQ1, "P2V"),
				 "pfPileUpAM" + newPostfixQ3 : cms.InputTag("pfPileUpAM" + newPostfixQ1),
                 "pfPileUpIsoAM" + newPostfixQ3 : cms.InputTag("pfPileUpAM" + newPostfixQ1),
				 "pfCandidateToVertexAssociation" + newPostfixQ3 : cms.InputTag('PFCand2VertexAM' + newPostfixQ1) }
				 
process.pNewQ1 = clonePath(process,process.pNewQ3,"PFCand2VertexAM"+newPostfixQ3,newPostfixQ3,newPostfixQ1,OldNewdictQ1)

getattr(process, "PFCand2VertexAM" + newPostfixQ1).MaxNumberOfAssociations = cms.int32(2)

getattr(process, "pfNoPileUpAM" + newPostfixQ1).VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAM' + newPostfixQ1)
getattr(process, "pfNoPileUpAM" + newPostfixQ1).MinQuality = cms.int32(1)

getattr(process, "pfPileUpAM" + newPostfixQ1).topCollection = cms.InputTag("pfNoPileUpAM" + newPostfixQ1, "P2V")
getattr(process, "pfPileUpAM" + newPostfixQ1).name = cms.untracked.string('pileUpOnPFCandidates'+newPostfixQ1)

getattr(process, "patPFMETtype0Corr" + newPostfixQ1).quality = cms.int32(1)


## jet-met

newPostfixJM = 'PFJM'
		
process.PFCand2VertexAMJetMET = PFCandAssoMapJetMet.clone(
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	IVFVertexCollection = cms.InputTag(''),
)

process.pfNoPileUpAMPFJM = FirstVertexPFCandidates.clone(
	VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAMJetMET'),
	VertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
	AssociationType = cms.InputTag('PFCandsToVertex'),
	MinQuality = cms.int32(3),
)
		
process.pfPileUpAMPFJM = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfNoPileUpAMPFJM", "P2V"),
    name = cms.untracked.string('pileUpOnPFCandidates'+newPostfixJM),
    verbose = cms.untracked.bool(False)
)

process.myPUSequencePFJM = cms.Sequence(process.PFCand2VertexAMJetMET*process.pfNoPileUpAMPFJM*process.pfPileUpAMPFJM)


OldNewdictJM = { "pfNoPileUp" + postfix : cms.InputTag("pfNoPileUpAM" + newPostfixJM, "P2V"),
                 "pfNoPileUpIso" + postfix : cms.InputTag("pfNoPileUpAM" + newPostfixJM, "P2V"),
				 "pfPileUp" + postfix : cms.InputTag("pfPileUpAM" + newPostfixJM),
                 "pfPileUpIso" + postfix : cms.InputTag("pfPileUpAM" + newPostfixJM),
				 "pfCandidateToVertexAssociation" + postfix : cms.InputTag("PFCand2VertexAMJetMET") }
				 
process.pJM = clonePath(process,process.p,"pfAllNeutralHadrons"+postfix,postfix,newPostfixJM,OldNewdictJM)

getattr( process, 'pJM' ).remove( getattr( process, 'particleFlowDisplacedVertex' + newPostfixJM ) )
getattr( process, 'pJM' ).remove( getattr( process, 'pfCandidateToVertexAssociation' + newPostfixJM ) )

getattr( process, 'pJM' ).remove( getattr( process, 'pfPileUp' + newPostfixJM ) )
getattr( process, 'pJM' ).remove( getattr( process, 'pfNoPileUp' + newPostfixJM ) )

getattr( process, 'pJM' ).remove( getattr( process, 'pfNoPileUpIso' + newPostfixJM ) )

getattr(process,"pJM").replace(
    getattr(process,"pfPileUpIso"+newPostfixJM),
    process.myPUSequencePFJM 
)

getattr(process, "patPFMETtype0Corr" + newPostfixJM).quality = cms.int32(3)

process.pJM += getattr( process, 'metvalidator' )