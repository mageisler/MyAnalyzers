// -*- C++ -*-
//
// Package:    FindReweight
// Class:      FindReweight
// 
/**\class FindReweight FindReweight.cc MGeisler/FindReweight/src/FindReweight.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Fri Sep  6 14:57:46 CEST 2013
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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <TH1F.h>

//
// class declaration
//

class FindReweight : public edm::EDAnalyzer {
   public:
      explicit FindReweight(const edm::ParameterSet&);
      ~FindReweight();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      
      edm::InputTag vertices_;
      edm::InputTag puinfo_;
      
      bool isData_;
      
      TH1F* h_numVtx;
      TH1F* h_numVtx_z;
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
FindReweight::FindReweight(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  vertices_ = iConfig.getParameter<InputTag>("vertices");
  puinfo_ = iConfig.getParameter<InputTag>("puinfo");
  isData_ = iConfig.getParameter<bool>("isData");
  
  h_numVtx = new TH1F("h_numVtx","", 60, -0.5, 59.5);
  h_numVtx_z = new TH1F("h_numVtx_z","", 200, -50., 50.);
  
}


FindReweight::~FindReweight()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FindReweight::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  //get the vertex collection
  Handle< VertexCollection >  vcH ;
  iEvent.getByLabel(vertices_, vcH); 
  
  for ( unsigned vtx_idx=0; vtx_idx<vcH->size(); vtx_idx++) {
    VertexRef vtx_ref(vcH, vtx_idx);
    h_numVtx_z->Fill( vtx_ref->z() );
  }
  
  if ( !isData_ ) {

    //get the pileup information  
    Handle< vector<PileupSummaryInfo> > puinfoH;
    iEvent.getByLabel(puinfo_,puinfoH);
    PileupSummaryInfo puinfo; 
   
    for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){ 
      if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
        puinfo=(*puinfoH)[puinfo_ite];
        break;
      }
    }
    
    int ntrue = puinfo.getTrueNumInteractions();
    
    h_numVtx->Fill( ntrue );
      
  }
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FindReweight::endJob() 
{

  double nVinteg = h_numVtx->Integral();
  if ( nVinteg>0. ) {
    h_numVtx->Scale( 1./nVinteg );
  }

  double nVzinteg = h_numVtx_z->Integral();
  if ( nVzinteg>0. ) {
    h_numVtx_z->Scale( 1./nVzinteg );
  }
  
  cout << " Num vertices: " << endl;
  
  for (unsigned bin_ite=1; bin_ite<=60; bin_ite++ ) {
    cout << h_numVtx->GetBinContent(bin_ite) << " ," << endl;
  }
  
  cout << " " << endl;
  cout << " z dist.: " << endl;
  
  for (unsigned bin_ite=1; bin_ite<=200; bin_ite++ ) {
    cout << h_numVtx_z->GetBinContent(bin_ite) << " ," << endl;
  }


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FindReweight::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FindReweight);
