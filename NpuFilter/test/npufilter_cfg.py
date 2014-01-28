import FWCore.ParameterSet.Config as cms
import sys, os


#_____________________
#
# READ INPUT PARAMETERS
#_____________________


# maximum number of events
MaxNumberOfEvents = cms.untracked.int32(-1) # reduce for testing

process = cms.Process("METVALIDATOR")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:/user/geisler/RelValTTbar_PU_START53_GEN-SIM-DIGI-RAW-HLTDEBUG.root"),
)
		
process.maxEvents = cms.untracked.PSet(
    input = MaxNumberOfEvents
)	

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

process.npufilter = cms.EDFilter("NpuFilter",
    PULabel = cms.InputTag("addPileupInfo"),
    minPU = cms.int32(20),
    maxPU = cms.int32(21),
)

### sequences and paths

process.p = cms.Path(process.npufilter)