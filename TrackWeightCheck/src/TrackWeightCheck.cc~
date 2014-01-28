// -*- C++ -*-
//
// Package:    TrackWeightCheck
// Class:      TrackWeightCheck
// 
/**\class TrackWeightCheck TrackWeightCheck.cc MGeisler/TrackWeightCheck/src/TrackWeightCheck.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Thu Aug  1 16:00:31 CEST 2013
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
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// root includes
#include "TH1F.h"

//
// class declaration
//

class TrackWeightCheck : public edm::EDAnalyzer {
   public:
      explicit TrackWeightCheck(const edm::ParameterSet&);
      ~TrackWeightCheck();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);


      // ----------member data ---------------------------

      std::vector<TH1F*> h_tw_vs_eta;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


using namespace edm;
using namespace std;
using namespace reco;
  
//
// constructors and destructor
//
TrackWeightCheck::TrackWeightCheck(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  Service<TFileService> tfs;

  for ( unsigned ite = 0; ite<3; ite++) {             
      
    char h_name[48];    
    
    sprintf(h_name,"h_tw_vs_eta_%i",ite);
     
    h_tw_vs_eta.push_back( tfs->make<TH1F>(h_name, " ; #eta; # entries", 50, -2.5, 2.5) );
    
  }

}


TrackWeightCheck::~TrackWeightCheck()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackWeightCheck::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //get track collection from the event
  Handle<TrackCollection>  tCH;
  iEvent.getByLabel("generalTracks", tCH);

  //get vertex collection from the event
  Handle<VertexCollection>  vCH;
  iEvent.getByLabel("offlinePrimaryVertices", vCH); 
  
  for ( unsigned trk_idx=0; trk_idx<tCH->size(); trk_idx++ ) {
  
    TrackRef trk_ref(tCH, trk_idx);
    
    double trk_eta = trk_ref->eta();
    
    bool hasTW = false;
    bool isSig = false;
    
    for ( unsigned vtx_idx=0; vtx_idx<vCH->size(); vtx_idx++ ) {
    
      VertexRef vtx_ref(vCH, vtx_idx);
    
      if ( vtx_ref->trackWeight( trk_ref )>0. ) {
        hasTW = true;
        
        if ( vtx_idx==0 ) {
          isSig = true;
        } 
        
        break;
        
      }
    
    }
    
    if ( hasTW ) {
      h_tw_vs_eta[1]->Fill(trk_eta);
      if ( isSig ) {
        h_tw_vs_eta[0]->Fill(trk_eta);
      }        
    } 
    
    h_tw_vs_eta[2]->Fill(trk_eta);    
  
  }
   
   
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackWeightCheck::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackWeightCheck);
