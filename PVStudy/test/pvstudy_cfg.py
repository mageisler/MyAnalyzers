import FWCore.ParameterSet.Config as cms

process = cms.Process("PVStudy")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')

options.register ('runOnMC',
                  True, # default value, allowed : True, False
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,  
                  "run on Monte Carlo or real data")

options.parseArguments()
	
tag =  'GR_P_V42_AN4::All'
Outfile = "PVStudy_Run2012D.root"
Infile = ["/store/data/Run2012D/MinimumBias/AOD/PromptReco-v1/000/207/905/FC64E288-E938-E211-B71A-001D09F2437B.root"]
if options.runOnMC:
    Outfile = "PVStudy_MC-QCD.root"
    tag = 'START53_V18PR::All'
    Infile = ["file:/user/geisler/QCD_Pt-0to5_TuneZ2star_Summer12_PU_S10_START53_AODSIM.root"]

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = tag
    
print "Oufile is set to", Outfile, ", with tag", tag

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("MGeisler.PVStudy.pvstudy_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet(
    	input = cms.untracked.int32(-1)
)


process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(),
	secondaryFileNames = cms.untracked.vstring(), 
)

process.ana.verbose = False
process.ana.realData = options.runOnMC
process.ana.analyzeOnTheFly = True

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(Outfile),
    closeFileFast = cms.untracked.bool(True),
)

process.p = cms.Path(process.ana)
		
process.source.fileNames = Infile
