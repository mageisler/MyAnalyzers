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

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50

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

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( 
	input = cms.untracked.int32(100000) 
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
	lumisToProcess = cms.untracked.VLuminosityBlockRange(*good_lumis),
)

from MGeisler.TrackValidatorData.inputfiles import inputfilesDATA, inputfilesMC 

if isData:
	process.source.fileNames = inputfilesDATA
else:
	process.source.fileNames = inputfilesMC

	
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
   		src = cms.InputTag('offlinePrimaryVertices'),
		cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
		filter = cms.bool(True),
)		

process.demo = cms.EDAnalyzer('FindReweight',
  	vertices = cms.InputTag('selectedPrimaryVertexQuality'),
  	puinfo  = cms.InputTag("addPileupInfo"),
  	isData = cms.bool(isData),
)


process.p = cms.Path(process.selectedPrimaryVertexQuality * process.demo)