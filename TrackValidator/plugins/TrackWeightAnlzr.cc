// -*- C++ -*-
//
// Package:    TrackWeightAnlzr
// Class:      TrackWeightAnlzr
// 
/**\class TrackWeightAnlzr TrackWeightAnlzr.cc MGeisler/TrackWeightAnlzr/src/TrackWeightAnlzr.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Tue Sep  4 15:43:59 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// root includes
#include "TH3F.h"
//
// class declaration
//

class TrackWeightAnlzr : public edm::EDAnalyzer {
   public:
      explicit TrackWeightAnlzr(const edm::ParameterSet&);
      ~TrackWeightAnlzr();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static reco::VertexRef TrackWeightAssociation(const reco::TrackBaseRef&, edm::Handle<reco::VertexCollection>);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      // ----------member data ---------------------------

      edm::InputTag tcLabel_;
      edm::InputTag vcLabel_;

      edm::InputTag tpLabel_;

      TH3F *h3_pur_tw;
};

//
// constants, enums and typedefs
//

using namespace std;
using namespace edm;
using namespace reco;

//
// static data member definitions
//

//
// constructors and destructor
//
TrackWeightAnlzr::TrackWeightAnlzr(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  tcLabel_ = iConfig.getParameter<InputTag>("TrackCollection");
  vcLabel_ = iConfig.getParameter<InputTag>("VertexCollection");

  tpLabel_ = iConfig.getParameter<InputTag>("TrackingParticles");

  //--------------


  Service<TFileService> tfs;

  h3_pur_tw = tfs->make<TH3F>("h3_pur_tw", "", 3, -0.5, 2.5, 2, -0.5, 1.5, 20, 0., 1.);

}


TrackWeightAnlzr::~TrackWeightAnlzr()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

VertexRef 
TrackWeightAnlzr::TrackWeightAssociation(const TrackBaseRef&  trackbaseRef, Handle<VertexCollection> vtxcollH) 
{

	VertexRef bestvertexref(vtxcollH,0);		
 	float bestweight = 0.;

	//loop over all vertices in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollH->size(); ++index_vtx){

          VertexRef vertexref(vtxcollH,index_vtx);

     	  //get the most probable vertex for the track
	  float weight = vertexref->trackWeight(trackbaseRef);
	  if(weight>bestweight){
  	    bestweight = weight;
	    bestvertexref = vertexref;
 	  } 

	}

  	return bestvertexref;

}

// ------------ method called for each event  ------------
void
TrackWeightAnlzr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //associate reco tracks to tracking particles
  ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get("quickTrackAssociatorByHits",theAssociator);
  TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product();

  //get the reco tracks   
  Handle<TrackCollection>  theTracksH;
  Handle<View<Track> >  theTracksV;
  iEvent.getByLabel(tcLabel_,theTracksH);
  iEvent.getByLabel(tcLabel_,theTracksV);

  //get the reco tracks   
  Handle<VertexCollection>  theVerticesH;
  iEvent.getByLabel(vcLabel_,theVerticesH);

  VertexRef firstVertex(theVerticesH, 0);

  //get the tracking particles   
  Handle<TrackingParticleCollection>  theTPsH;
  iEvent.getByLabel(tpLabel_,theTPsH);

  RecoToSimCollection recSimColl;
  recSimColl=theTrackAssociator_->associateRecoToSim(theTracksV,theTPsH,&iEvent,&iSetup);
	    
  //loop over all tracks of the track collection	
  for ( size_t idxTrack = 0; idxTrack < theTracksH->size(); ++idxTrack ) {

    TrackRef trackref(theTracksH, idxTrack);
    TrackBaseRef trackbaseref = TrackBaseRef(trackref);

    vector<pair<TrackingParticleRef, double> > tp;
    if(recSimColl.find(trackbaseref) != recSimColl.end()) tp = recSimColl[trackbaseref];
  
    int x_bin = 1;

    if (tp.size()==0){ 
      x_bin=2;
    }else{     
      TrackingParticleRef tpr = tp[0].first;
      if ((tpr->eventId().event() == 0) && (tpr->eventId().bunchCrossing() == 0)) x_bin=0;
    }

    VertexRef bestvertex = TrackWeightAnlzr::TrackWeightAssociation(trackbaseref, theVerticesH);
    double tweight = bestvertex->trackWeight(trackbaseref);

    if ( tweight > 1e-5 ) {

      int y_bin = 1;
      if ( bestvertex==firstVertex ) y_bin=0;

      h3_pur_tw->Fill( x_bin, y_bin, tweight);

    }

  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackWeightAnlzr::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackWeightAnlzr);
