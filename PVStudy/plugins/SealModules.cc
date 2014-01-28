#include "FWCore/Framework/interface/MakerMacros.h"

#include "MGeisler/PVStudy/interface/LhcTrackAnalyzer.h"
#include "MGeisler/PVStudy/interface/PVStudy.h"
#include "MGeisler/PVStudy/interface/PVEffAnalyzer.h"
#include "MGeisler/PVStudy/interface/TrackSplitterProducer.h"
#include "MGeisler/PVStudy/interface/VtxTrackSplitterProducer.h" 

// DEFINE_SEAL_MODULE();
DEFINE_FWK_MODULE(LhcTrackAnalyzer);
DEFINE_FWK_MODULE(PVStudy);
DEFINE_FWK_MODULE(PVEffAnalyzer);
DEFINE_FWK_MODULE(TrackSplitterProducer);
DEFINE_FWK_MODULE(VtxTrackSplitterProducer);
