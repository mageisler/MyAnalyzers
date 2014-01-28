// -*- C++ -*-
//
// Package:    FirstVertexCheck
// Class:      FirstVertexCheck
// 
/**\class FirstVertexCheck FirstVertexCheck.cc MGeisler/FirstVertexCheck/src/FirstVertexCheck.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Wed Nov 27 14:58:10 CET 2013
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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include <TProfile.h>

//
// class declaration
//

class FirstVertexCheck : public edm::EDAnalyzer {
   public:
      explicit FirstVertexCheck(const edm::ParameterSet&);
      ~FirstVertexCheck();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------
      
      TProfile* p_vertexEff5_npu;
      TProfile* p_vertexEff4_npu;
      TProfile* p_vertexEff3_npu;
      TProfile* p_vertexEff2_npu;
      TProfile* p_vertexEff1_npu; 
      
      TProfile* p_vertexEffC_npu; 
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
FirstVertexCheck::FirstVertexCheck(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  //parameters for npu axis
  double minPUcount  = -0.5;  
  double maxPUcount = 59.5;
  int nintPUcount = 60;

  Service<TFileService> tfs;
  vector<TFileDirectory>* subDir(new vector<TFileDirectory>());
  subDir->push_back(tfs->mkdir("offlinePrimaryVertices"));

  p_vertexEff5_npu = subDir->at(0).make<TProfile>("p_vertexEff5_npu", "; number of simulated vertices; efficiency", nintPUcount, minPUcount, maxPUcount);
  p_vertexEff4_npu = subDir->at(0).make<TProfile>("p_vertexEff4_npu", "; number of simulated vertices; efficiency", nintPUcount, minPUcount, maxPUcount);
  p_vertexEff3_npu = subDir->at(0).make<TProfile>("p_vertexEff3_npu", "; number of simulated vertices; efficiency", nintPUcount, minPUcount, maxPUcount);
  p_vertexEff2_npu = subDir->at(0).make<TProfile>("p_vertexEff2_npu", "; number of simulated vertices; efficiency", nintPUcount, minPUcount, maxPUcount);
  p_vertexEff1_npu = subDir->at(0).make<TProfile>("p_vertexEff1_npu", "; number of simulated vertices; efficiency", nintPUcount, minPUcount, maxPUcount);
  
  p_vertexEff5_npu->Sumw2();
  p_vertexEff4_npu->Sumw2();
  p_vertexEff3_npu->Sumw2();
  p_vertexEff2_npu->Sumw2();
  p_vertexEff1_npu->Sumw2();

  p_vertexEffC_npu = subDir->at(0).make<TProfile>("p_vertexEffC_npu", "; number of simulated vertices; efficiency", nintPUcount, minPUcount, maxPUcount);
  
  p_vertexEffC_npu->Sumw2();

}


FirstVertexCheck::~FirstVertexCheck()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FirstVertexCheck::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //get the pileup information  
  Handle< vector<PileupSummaryInfo> > puinfoH;
  iEvent.getByLabel("addPileupInfo",puinfoH);
  PileupSummaryInfo puinfo; 
  
  for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){ 
    if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
      puinfo=(*puinfoH)[puinfo_ite];
      break;
    }
  }
  
  int npu = puinfo.getPU_NumInteractions();
  
  //get the tracking vertices 
  Handle< SimVertexContainer > svH;
  iEvent.getByLabel("g4SimHits",svH);
   
  SimVertexRef fvr_sim;
   
  for (unsigned int svi=0; svi<svH->size(); svi++) {
    SimVertexRef svr(svH, svi);
    if ( ( svr->noParent() ) && (svr->eventId().bunchCrossing()==0) && (svr->eventId().event()==0) && ( svr->vertexId()==0 ) ) {
      fvr_sim = svr;
      break; 
    }
  } 
  
  double z_sim = fvr_sim->position().z();
  
  //get vertex collection from the event
  Handle<VertexCollection>  vertexCollectionH;
  iEvent.getByLabel("offlinePrimaryVertices", vertexCollectionH);
  
  double min_dist = 10.;
     
  for (unsigned int vi=1; vi<vertexCollectionH->size(); vi++) {
  
    VertexRef vtx_ref(vertexCollectionH, vi); 
    
    double z_tmp = vtx_ref->position().z();
    
    double dist_tmp = fabs(z_tmp-z_sim);
    
    if ( dist_tmp<min_dist  ) {
    
      min_dist = dist_tmp;
    
    }
    
  
  }  
  
  VertexRef fvr_rec(vertexCollectionH, 0); 
    
  double z_rec = fvr_rec->position().z();
  
  double dist = fabs(z_rec-z_sim);
  
  double eff5 = 0.;
  double eff4 = 0.;
  double eff3 = 0.;
  double eff2 = 0.;
  double eff1 = 0.;
  
  double effC = 0.;
  
  if ( dist<=0.5 ) {
    eff5 = 1.;
  }
  
  if ( dist<=0.4 ) {
    eff4 = 1.;
  }
  
  if ( dist<=0.3 ) {
    eff3 = 1.;
  }
  
  if ( dist<=0.2 ) {
    eff2 = 1.;
  }
  
  if ( dist<=0.1 ) {
    eff1 = 1.;
  }
  
  if ( dist<min_dist ) {
    effC = 1.;  
  }

  p_vertexEff5_npu->Fill(npu, eff5);
  p_vertexEff4_npu->Fill(npu, eff4);
  p_vertexEff3_npu->Fill(npu, eff3);
  p_vertexEff2_npu->Fill(npu, eff2);
  p_vertexEff1_npu->Fill(npu, eff1);

  p_vertexEffC_npu->Fill(npu, effC);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FirstVertexCheck::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FirstVertexCheck);
