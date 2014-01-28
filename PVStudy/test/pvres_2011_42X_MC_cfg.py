import FWCore.ParameterSet.Config as cms

process = cms.Process("pvtxana")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.Geometry.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')

options.register ('runOnMC',
                  True, # default value, allowed : True, False
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,  
                  "run on Monte Carlo or real data")

options.parseArguments()
	
tag =  'GR_P_V42_AN4::All'
Outfile = "VertexAnalyzer_Resolution_Run2012D"
Infile = ["/store/data/Run2012D/MinimumBias/AOD/PromptReco-v1/000/207/905/FC64E288-E938-E211-B71A-001D09F2437B.root"]
if options.runOnMC:
    Outfile = "VertexAnalyzer_Resolution_MC-QCD"
    tag = 'START53_V18PR::All'
    Infile = ["file:/user/geisler/QCD_Pt-0to5_TuneZ2star_Summer12_PU_S10_START53_AODSIM.root"]
    
print "Oufile is set to", Outfile, ", with tag", tag

process.MessageLogger.cerr.FwkReport.reportEvery = 500

#======================================
# Global Tag
#======================================
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = tag

#======================================
# Input
#======================================

readFiles = cms.untracked.vstring(Infile)
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

process.maxEvents = cms.untracked.PSet(
    	input = cms.untracked.int32(-1)
)
 
#======================================
# Trigger Filter
#======================================
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlinePrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 0 && abs(z) <= 24 && position.Rho <= 2"),
                                           filter = cms.bool(True)
                                           )

#====================================
# VtxTrackSplitterProducer
#====================================
process.load("MGeisler.PVStudy.VtxTrackSplitterProducer_cff")

#====================================
# PVProducer based on SplittedTracks1
#====================================
process.load("RecoVertex.Configuration.RecoVertex_cff")
process.PVProducer1 = process.offlinePrimaryVertices.clone()
process.PVProducer1.TrackLabel = cms.InputTag("VtxTrackSplitterProducer","SplittedTracks1")

#====================================
# PVProducer based on SplittedTracks2
#====================================
process.PVProducer2 = process.offlinePrimaryVertices.clone()
process.PVProducer2.TrackLabel = cms.InputTag("VtxTrackSplitterProducer","SplittedTracks2")

#====================================
# PVT Analyzer 
#==================================== 
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("MGeisler.PVStudy.pvstudy_cfi")
process.ana.realData = options.runOnMC
process.ana.histoFileName = Outfile + '_histos.root'
process.ana.OutputFileName = Outfile + '_ntuples.root'

process.ana.nTrkMin = 5
process.ana.nTrkMax = 100
process.ana.ntrkdiffcut = cms.untracked.double(0.1)
process.ana.avgTrkPtInPVMin = cms.untracked.double(0)
process.ana.avgTrkPtInPVMax = cms.untracked.double(10000)

#======================================
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.allPath=cms.Path(process.primaryVertexFilter*process.VtxTrackSplitterProducer*(process.PVProducer1+process.PVProducer2)*process.ana)
