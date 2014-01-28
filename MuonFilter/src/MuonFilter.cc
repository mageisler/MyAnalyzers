// -*- C++ -*-
//
// Package:    MuonFilter
// Class:      MuonFilter
// 
/**\class MuonFilter MuonFilter.cc MGeisler/MuonFilter/src/MuonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Tue Sep  3 13:53:29 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "TLorentzVector.h"

//
// class declaration
//

class MuonFilter : public edm::EDFilter {
   public:
      explicit MuonFilter(const edm::ParameterSet&);
      ~MuonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------
      
      edm::InputTag muonLabel_;
      
      double isolationFactor_;
      
      double massCut_;
};

using namespace edm;
using namespace std;
using namespace reco;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonFilter::MuonFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muonLabel_ = iConfig.getParameter<InputTag>("Muons");
  
  isolationFactor_ = iConfig.getParameter<double>("isolationFactor");
  
  massCut_ = iConfig.getParameter<double>("massCut");

}


MuonFilter::~MuonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MuonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //get the input muon collection
  Handle< MuonCollection > mCH;
  iEvent.getByLabel( muonLabel_, mCH );
  
  bool correctSize = (mCH->size()==2);
  
  if ( !correctSize ) return false;
   
  bool oneIsolated = false;
  TLorentzVector sum_p4;
  int totalCharge = 0;
  
  for ( unsigned mu_idx=0; mu_idx<mCH->size(); mu_idx++ ) {
  
    MuonRef mu_ref(mCH, mu_idx);
    
    totalCharge+= mu_ref->charge();
    
    sum_p4+= TLorentzVector( mu_ref->px(), mu_ref->py(), mu_ref->pz(), mu_ref->energy() );
    
    if ( oneIsolated ) continue;
    
    const MuonPFIsolation mI04 = mu_ref->pfIsolationR04();
    
    double iso_correction = mI04.sumNeutralHadronEt + mI04.sumPhotonEt - 0.5*mI04.sumPUPt;
         
    if ( iso_correction < 0. ) {
      iso_correction = 0.;
    }
    
    double isolation = mI04.sumChargedParticlePt + iso_correction;
    
    if ( isolation < isolationFactor_*mu_ref->pt() )  {
         
      oneIsolated = true;
       
    }
    
  }
  
  double invMass = sum_p4.M();
  
  return ( oneIsolated && ( totalCharge==0 ) && ( fabs( invMass-90. )<=massCut_ ) );
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonFilter);
