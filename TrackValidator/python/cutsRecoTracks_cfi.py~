import FWCore.ParameterSet.Config as cms

import PhysicsTools.RecoAlgos.recoTrackSelector_cfi
cutsRecoTracks = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone(
	minRapidity = cms.double(-2.4),
	maxRapidity = cms.double(2.4),   
	minPixelHit = cms.int32(0),
    tip = cms.double(120.0),
    lip = cms.double(280.0),
	quality = cms.vstring('highPurity','tight','loose'),
)
