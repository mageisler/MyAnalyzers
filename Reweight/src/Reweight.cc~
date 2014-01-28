// -*- C++ -*-
//
// Package:    Reweight
// Class:      Reweight
// 
/**\class Reweight Reweight.cc MGeisler/Reweight/src/Reweight.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Tue Aug 27 10:51:03 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

// ROOT include files
#include <TH1F.h>


//
// class declaration
//

class Reweight : public edm::EDProducer {
   public:
      explicit Reweight(const edm::ParameterSet&);
      ~Reweight();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

  edm::InputTag vcLabel_;
      
  edm::LumiReWeighting LumiWeights_;
  
  std::vector< float > MC_distrV;
  std::vector< float > Data_distrV;
  
  std::vector< float > zWeightsV;
  
  TH1F* h_vertexNum_z;
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
Reweight::Reweight(const edm::ParameterSet& iConfig)
{
  //read input parameter

  vcLabel_ = iConfig.getParameter<edm::InputTag>("recoVertexCollection");

  //PileUp reweighting stuff
         
  float MC_distr[60] = {
    		4.00002e-06 ,
			1.2e-05 ,
			2.00001e-05 ,
			2.00001e-05 ,
			0.0001 ,
			0.000268001 ,
			0.00197601 ,
			0.00611202 ,
			0.010516 ,
			0.0139601 ,
			0.0171081 ,
			0.0203481 ,
			0.0259081 ,
			0.0330041 ,
			0.0408882 ,
			0.0491362 ,
			0.0556682 ,
			0.0566562 ,
			0.0557562 ,
			0.0532522 ,
			0.0491082 ,
			0.0486562 ,
			0.0456842 ,
			0.0433282 ,
			0.0416882 ,
			0.0400202 ,
			0.0366121 ,
			0.0343761 ,
			0.0309561 ,
			0.0281321 ,
			0.0255041 ,
			0.0224801 ,
			0.0193681 ,
			0.0166641 ,
			0.0142361 ,
			0.012012 ,
			0.010324 ,
			0.00838803 ,
			0.00682803 ,
			0.00573202 ,
			0.00476002 ,
			0.00328801 ,
			0.00268001 ,
			0.00218401 ,
			0.00151201 ,
			0.00129201 ,
			0.000888004 ,
			0.000672003 ,
			0.000532002 ,
			0.000404002 ,
			0.000304001 ,
			0.000212001 ,
			0.000180001 ,
			0.000112 ,
			7.20003e-05 ,
			4.00002e-05 ,
			2.80001e-05 ,
			1.60001e-05 ,
			8.00003e-06 ,
			4.00002e-06};
    			
  float Data_distr[60] = {
    		0.0 ,
			0.0 ,
			0.0 ,
			0.0 ,
			0.0 ,
			2.23444268035e-17 ,
			3.34675233933e-16 ,
			4.53722907965e-15 ,
			5.60983119342e-14 ,
			6.33131856256e-13 ,
			6.78111639186e-12 ,
			8.38302558008e-11 ,
			1.64133400682e-09 ,
			4.08270353508e-08 ,
			8.52644004591e-07 ,
			1.28283398605e-05 ,
			0.00013419169517 ,
			0.000969906504982 ,
			0.00485371940501 ,
			0.0169194995024 ,
			0.0415032949737 ,
			0.0729299790902 ,
			0.0950468931229 ,
			0.0981933182667 ,
			0.0900176256444 ,
			0.0829396029151 ,
			0.0800627768658 ,
			0.0776151244774 ,
			0.0728194922419 ,
			0.0656383631145 ,
			0.0567438532367 ,
			0.0465870764517 ,
			0.0357948623257 ,
			0.0253938222681 ,
			0.016480987393 ,
			0.00973762228851 ,
			0.00522703757513 ,
			0.00254826047949 ,
			0.00112884578536 ,
			0.000454736000275 ,
			0.00016668857559 ,
			5.56205103924e-05 ,
			1.6894847764e-05 ,
			4.67015014855e-06 ,
			1.17415552614e-06 ,
			2.68308684947e-07 ,
			5.56836132752e-08 ,
			1.04875898367e-08 ,
			1.79132446023e-09 ,
			2.77301019966e-10 ,
			3.8884250875e-11 ,
			4.93679630242e-12 ,
			5.67289246233e-13 ,
			5.89818615387e-14 ,
			5.54768331673e-15 ,
			4.71937660198e-16 ,
			3.68532141479e-17 ,
			2.21591572015e-18 ,
			0.0 ,
			0.0 };
    
  for (int i = 0; i<60; i++) {
    MC_distrV.push_back( MC_distr[i] );
    Data_distrV.push_back( Data_distr[i] );    
  }
    
  LumiWeights_ = edm::LumiReWeighting(MC_distrV, Data_distrV); 
  
  float MC_z[200] = {
    		0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			8.41369e-06 ,
			1.2748e-05 ,
			1.73373e-05 ,
			2.80456e-05 ,
			3.95188e-05 ,
			4.89524e-05 ,
			6.32301e-05 ,
			8.49017e-05 ,
			0.000152466 ,
			0.000164449 ,
			0.000186631 ,
			0.000264139 ,
			0.000340882 ,
			0.000459183 ,
			0.000555813 ,
			0.0007042 ,
			0.000901794 ,
			0.00114069 ,
			0.00145021 ,
			0.00176968 ,
			0.002146 ,
			0.00261997 ,
			0.00317527 ,
			0.00377826 ,
			0.00452478 ,
			0.00528838 ,
			0.00619834 ,
			0.00718299 ,
			0.0083502 ,
			0.00956152 ,
			0.0108371 ,
			0.0122039 ,
			0.0136605 ,
			0.0151643 ,
			0.0167381 ,
			0.0182921 ,
			0.0199282 ,
			0.0215459 ,
			0.0229844 ,
			0.0246001 ,
			0.0258397 ,
			0.0273185 ,
			0.0284347 ,
			0.029595 ,
			0.0304066 ,
			0.031211 ,
			0.0318336 ,
			0.032163 ,
			0.0322387 ,
			0.0321594 ,
			0.0320352 ,
			0.0314947 ,
			0.030946 ,
			0.0301753 ,
			0.0291269 ,
			0.0278478 ,
			0.0267076 ,
			0.0253792 ,
			0.0241282 ,
			0.0224648 ,
			0.0209409 ,
			0.0194588 ,
			0.0176947 ,
			0.0161517 ,
			0.0144777 ,
			0.0130473 ,
			0.0116945 ,
			0.010357 ,
			0.00908117 ,
			0.00794634 ,
			0.00678831 ,
			0.00587148 ,
			0.00493705 ,
			0.00416732 ,
			0.00356307 ,
			0.00296723 ,
			0.00249122 ,
			0.00199583 ,
			0.00166616 ,
			0.00133676 ,
			0.00105554 ,
			0.000842133 ,
			0.000658307 ,
			0.000500232 ,
			0.000399778 ,
			0.000303658 ,
			0.000239153 ,
			0.000155271 ,
			0.000155526 ,
			0.000122891 ,
			7.39385e-05 ,
			5.91508e-05 ,
			5.02272e-05 ,
			3.2125e-05 ,
			2.14167e-05 ,
			1.65724e-05 ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. };
  
  float Data_z[200] = {
    		0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			5.41714e-07 ,
			5.41714e-07 ,
			1.62514e-06 ,
			1.62514e-06 ,
			5.95885e-06 ,
			9.20914e-06 ,
			1.35428e-05 ,
			2.60023e-05 ,
			4.06285e-05 ,
			7.42148e-05 ,
			0.000105634 ,
			0.000149513 ,
			0.000225353 ,
			0.000318528 ,
			0.000431746 ,
			0.000638681 ,
			0.000921455 ,
			0.00120206 ,
			0.00155201 ,
			0.00207639 ,
			0.00264519 ,
			0.00341009 ,
			0.00441768 ,
			0.00540847 ,
			0.00664575 ,
			0.00801303 ,
			0.00954608 ,
			0.0113998 ,
			0.0133213 ,
			0.0153722 ,
			0.0176436 ,
			0.019707 ,
			0.0223322 ,
			0.0245721 ,
			0.0268993 ,
			0.029155 ,
			0.0314118 ,
			0.0331133 ,
			0.0350321 ,
			0.0365061 ,
			0.0377195 ,
			0.0385169 ,
			0.0391925 ,
			0.0393837 ,
			0.0390088 ,
			0.0387537 ,
			0.0378246 ,
			0.0365391 ,
			0.0351691 ,
			0.0334449 ,
			0.0315532 ,
			0.0292899 ,
			0.0270624 ,
			0.0246843 ,
			0.0222195 ,
			0.0201241 ,
			0.0176918 ,
			0.015545 ,
			0.0133392 ,
			0.0112433 ,
			0.00967068 ,
			0.00803199 ,
			0.00660295 ,
			0.00541985 ,
			0.00426329 ,
			0.00331529 ,
			0.00266361 ,
			0.00198213 ,
			0.00153305 ,
			0.00108451 ,
			0.000832073 ,
			0.000629471 ,
			0.000456665 ,
			0.000281691 ,
			0.000213977 ,
			0.000126219 ,
			9.31748e-05 ,
			5.74217e-05 ,
			3.73783e-05 ,
			2.00434e-05 ,
			1.57097e-05 ,
			7.58399e-06 ,
			3.25028e-06 ,
			2.70857e-06 ,
			1.08343e-06 ,
			1.08343e-06 ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. ,
			0. };
    
  for (int i = 0; i<200; i++) { 
    if (MC_z[i] != 0. ) {
      zWeightsV.push_back( Data_z[i] *1./ MC_z[i] );
    } else {
      zWeightsV.push_back( 0. );
    }   
  }  
  
  h_vertexNum_z = new TH1F("h_vertexNum_z", "z weights; z / cm; weights", 200, -50., 50.);
  
  produces<float>();
  
}


Reweight::~Reweight()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
Reweight::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace edm;
  using namespace std;
  using namespace reco;

  //get the reco vertices 
  Handle< VertexCollection > vcH;
  iEvent.getByLabel(vcLabel_, vcH);

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
     
  int nTrue = puinfo.getTrueNumInteractions();
  
  float pu_weight = LumiWeights_.weight( nTrue );
  float z_weight = 1.;
      
  for( unsigned vtx_idx=0; vtx_idx<vcH->size(); ++vtx_idx ) {
    VertexRef w_vtxref(vcH, vtx_idx);
    int z_bin_number = h_vertexNum_z->FindBin( w_vtxref->z() );
    if ( ( z_bin_number>0 ) && ( z_bin_number<201 ) ) {
      z_weight*= zWeightsV[ z_bin_number - 1 ];
    } else {
      z_weight*= 0.;
    }
  } 
  
  auto_ptr<float> weight_output(new float( pu_weight * z_weight ));  
  
  iEvent.put( weight_output ); 
 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Reweight::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Reweight);
