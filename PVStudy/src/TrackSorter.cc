#include "MGeisler/PVStudy/interface/TrackSorter.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "MGeisler/PVStudy/interface/TrackHigherPt.h"

reco::TrackCollection TrackSorter::sortedList(const reco::TrackCollection & unsortedTrackColl) const
{
  reco::TrackCollection tracksorted = unsortedTrackColl;
  std::sort(tracksorted.begin(), tracksorted.end(), TrackHigherPt());
  return tracksorted;
}

