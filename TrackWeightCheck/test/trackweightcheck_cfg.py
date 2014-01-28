import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service('TFileService',
    fileName = cms.string("twcheck.root"),
    closeFileFast = cms.untracked.bool(True),
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/user/geisler/QCD_Pt-0to5_TuneZ2star_Summer12_PU_S10_START53_AODSIM.root'
    )
)

process.demo = cms.EDAnalyzer('TrackWeightCheck'
)


process.p = cms.Path(process.demo)
