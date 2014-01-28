// -*- C++ -*-
//
// Package:    FromParticlesToTracks
// Class:      FromParticlesToTracks
// 
/**\class FromParticlesToTracks FromParticlesToTracks.cc MGeisler/FromParticlesToTracks/src/FromParticlesToTracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Thu Jan 23 10:43:52 CET 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"


//
// class declaration
//

class FromParticlesToTracks : public edm::EDProducer {
   public:
      explicit FromParticlesToTracks(const edm::ParameterSet&);
      ~FromParticlesToTracks();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
     
      virtual bool TrackMatch(reco::Track, reco::Track);
      

      // ----------member data ---------------------------
      
      edm::InputTag input_ParticleCollection_;
      edm::InputTag input_generalTracksCollection_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
FromParticlesToTracks::FromParticlesToTracks(const edm::ParameterSet& iConfig)
{

  //now do what ever other initialization is needed
  
  input_ParticleCollection_ = iConfig.getParameter<edm::InputTag>("ParticleCollection");
  input_generalTracksCollection_ = iConfig.getParameter<edm::InputTag>("RefTrackCollection");

  //register your products
  
  produces<reco::TrackCollection>();
  
}


FromParticlesToTracks::~FromParticlesToTracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
FromParticlesToTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace reco;
   
  std::auto_ptr<TrackCollection> output(new TrackCollection());
  
  //get the input pfCandidateCollection
  Handle<PFCandidateCollection> pfCandH;
  iEvent.getByLabel(input_ParticleCollection_, pfCandH);
 
  //get the input track collection
  Handle<TrackCollection> input_trckcollH;
  iEvent.getByLabel(input_generalTracksCollection_, input_trckcollH);
   
  for ( unsigned i=0; i<pfCandH->size(); i++ ) {

    PFCandidateRef candref(pfCandH, i);
    TrackRef trackref = candref->trackRef();

    if ( !trackref.isNull() && trackref.isAvailable() ) {

  	  for( unsigned int index_input_trck=0; index_input_trck<input_trckcollH->size(); index_input_trck++ ){

	    TrackRef input_trackref = TrackRef(input_trckcollH,index_input_trck);

   	    if( TrackMatch(*trackref, *input_trackref) ){
	  
	      output->push_back( *trackref );
	    
	    }
	    
	  }
	  
	}
  
  }
  
  iEvent.put( output );
 
}

bool 
FromParticlesToTracks::TrackMatch(reco::Track track1, reco::Track track2)
{

	return (
	  (track1).eta()  == (track2).eta() &&
	  (track1).phi()  == (track2).phi() &&
	  (track1).chi2() == (track2).chi2() &&
	  (track1).ndof() == (track2).ndof() &&
	  (track1).p()    == (track2).p()
	);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FromParticlesToTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FromParticlesToTracks);