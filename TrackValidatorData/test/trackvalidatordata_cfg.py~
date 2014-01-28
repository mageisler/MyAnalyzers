import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')

options.register ('allinOne',
                  "",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,  
                  'Parsing all options in one')

options.parseArguments()
isData = False
isZMuMu = False

if options.allinOne != "" :
	opt = options.allinOne.split(",")
	isData = bool(int(opt[0]))
	if len(opt)>1:
		isZMuMu = bool(int(opt[1]))

print ""		
print " Is This Data?", isData
print " Is This Z to MuMu?", isZMuMu
print ""
					

process = cms.Process("TRACKVALIDATORDATA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

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

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(),	
)

if isData:
	process.source.fileNames = inputfilesDATA
	process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(*good_lumis)
else:
	process.source.fileNames = inputfilesMC
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

tag = 'START53_V18PR::All'

addOn = "MC"
if isData:
	tag =  'GR_P_V42_AN4::All'
	addOn = "Data"
	
Outfile= "TrackValidatorData_" + addOn +  ".root"	

print " Outfile set to ", Outfile
print " Globaltag set to ", tag
print ""

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = tag
		
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
		

process.eventweight = cms.EDProducer("Reweight",
		recoVertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
)

ZMuMuInputTag = cms.InputTag("")
		
if isZMuMu:
	
	from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
	
	process.doubleMuTrigger = hltHighLevel.clone(
		HLTPaths = cms.vstring( 'HLT_Mu17_Mu8_v*' ),
	)
	
	process.theTwoMuons = cms.EDProducer('MuonSelection',
		Muons = cms.InputTag("muons"),
		VertexCollection = cms.InputTag("selectedPrimaryVertexQuality"),
  		minLayers = cms.int32(5),
 		minPt = cms.double(0.),
  		maxEta = cms.double(2.4),
	)
	
	process.muonFilter = cms.EDFilter('MuonFilter',
		Muons = cms.InputTag("theTwoMuons"),
  		isolationFactor = cms.double(0.15),
  		massCut = cms.double(40.),
	)
	
	process.ZMuMusequence = cms.Sequence(process.doubleMuTrigger * process.theTwoMuons * process.muonFilter)
	ZMuMuInputTag = cms.InputTag("theTwoMuons") 
		
### IVF-specific includes
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
		
### AssociationMap-specific includes
from CommonTools.RecoUtils.pf_pu_assomap_cfi import AssociationMaps, AssociationMapsJetMet
		
process.assMapA1 = AssociationMaps.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),	    
		MaxNumberOfAssociations = cms.int32(1),
		IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),   
)
		
process.assMapA2 = AssociationMaps.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),	    
		MaxNumberOfAssociations = cms.int32(2),
		IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),   
)
		
process.assMapA3 = AssociationMaps.clone(
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),	    
		MaxNumberOfAssociations = cms.int32(3),
		IVFVertexCollection = cms.InputTag('inclusiveMergedVertices'),   
)
						
process.assMapJetMetA1 = AssociationMapsJetMet.clone(	
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),   	    
		MaxNumberOfAssociations = cms.int32(1),
)
						
process.assMapJetMetA2 = AssociationMapsJetMet.clone(	
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),   	    
		MaxNumberOfAssociations = cms.int32(2),
)
						
process.assMapJetMetA3 = AssociationMapsJetMet.clone(	
		AssociationType = cms.InputTag('TracksToVertex'),
		VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),   	    
		MaxNumberOfAssociations = cms.int32(3),
)


process.trackvalidatordata = cms.EDAnalyzer('TrackValidatorData',
		amLabels = cms.VInputTag(cms.InputTag("assMapA1"), cms.InputTag("assMapA2"), cms.InputTag("assMapA3"), cms.InputTag("assMapJetMetA1"), cms.InputTag("assMapJetMetA2"), cms.InputTag("assMapJetMetA3")),
    	RefMuonCollection = ZMuMuInputTag,
		ignoremissingtrackcollection=cms.bool(False),
    	vcLabel = cms.InputTag("selectedPrimaryVertexQuality"),
		Weight=cms.InputTag("eventweight"),
		isData=cms.bool(isData),
		isZMuMu=cms.bool(isZMuMu),
)	  
		

### sequences and paths

process.startingstuff = cms.Sequence(
		process.inclusiveVertexing *
		process.selectedPrimaryVertexQuality
)

if not isData:
	process.startingstuff+= process.eventweight
		
if isZMuMu:
	process.startingstuff+= process.ZMuMusequence

process.assomaps = cms.Sequence(
    	process.assMapA1 *
    	process.assMapA2 *
    	process.assMapA3 *
    	process.assMapJetMetA1 *
    	process.assMapJetMetA2 *
    	process.assMapJetMetA3
)
		
		  

process.validation = cms.Sequence(
      	process.trackvalidatordata 
)

process.p = cms.Path(
		process.startingstuff *
		process.assomaps *
		process.validation
)

