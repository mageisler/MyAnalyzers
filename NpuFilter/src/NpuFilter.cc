// -*- C++ -*-
//
// Package:    NpuFilter
// Class:      NpuFilter
// 
/**\class NpuFilter NpuFilter.cc MGeisler/NpuFilter/src/NpuFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Tue May 21 11:53:25 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class NpuFilter : public edm::EDFilter {
   public:
      explicit NpuFilter(const edm::ParameterSet&);
      ~NpuFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

  //input parameters

  edm::InputTag PULabel_;

  int minPU_;
  int maxPU_;

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
NpuFilter::NpuFilter(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed

  //read input parameter

  PULabel_ = iConfig.getParameter<edm::InputTag>("PULabel");

  minPU_ = iConfig.getParameter<int>("minPU");
  maxPU_ = iConfig.getParameter<int>("maxPU");

}


NpuFilter::~NpuFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
NpuFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //get the pileup information  
  edm::Handle< std::vector<PileupSummaryInfo> > puinfoH;
  iEvent.getByLabel("addPileupInfo",puinfoH);
  PileupSummaryInfo puinfo;    
  
  for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){ 
    if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
      puinfo=(*puinfoH)[puinfo_ite];
      break;
    }
  }
 
  int simpv = puinfo.getPU_NumInteractions();

  if ( ( simpv >= minPU_ ) &&
       ( simpv < maxPU_ ) ) {
       
    std::cout << "YES" << std::endl;    
    return true;

  }

  return false;


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NpuFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(NpuFilter);