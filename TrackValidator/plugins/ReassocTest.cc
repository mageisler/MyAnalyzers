// -*- C++ -*-
//
// Package:    ReassocTest
// Class:      ReassocTest
// 
/**\class ReassocTest ReassocTest.cc MGeisler/ReassocTest/src/ReassocTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Mon Sep  3 15:36:13 CEST 2012
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

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"

// root includes
#include "TH1F.h"
#include "TH3F.h"

//
// class declaration
//

class ReassocTest : public edm::EDAnalyzer {
   public:
      explicit ReassocTest(const edm::ParameterSet&);
      ~ReassocTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static bool TrackWeightAssociation(const reco::TrackRef&, edm::Handle<reco::VertexCollection>);

      static reco::VertexRef FindVertexZ(reco::TrackRef, edm::Handle<reco::VertexCollection>, double);
      static reco::VertexRef FindVertex3(reco::TransientTrack, edm::Handle<reco::VertexCollection>, double);

      static std::auto_ptr<reco::ConversionCollection> GetCleanedConversions(edm::Handle<reco::ConversionCollection>, 
                                                                             edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromConversion(const reco::TrackRef, reco::ConversionCollection, reco::Conversion*);        
      static reco::VertexRef FindConversionVertex3(const reco::TrackRef, reco::Conversion, 	
                                                  edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                  edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, 
	                                          double);
      static reco::VertexRef FindConversionVertexZ(const reco::TrackRef, reco::Conversion, 
                                                   edm::Handle<reco::VertexCollection>, double);

      static std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedKshort(edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);
      static std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedLambda(edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);
      
      static bool ComesFromDecay(const reco::TrackRef, reco::VertexCompositeCandidateCollection, reco::VertexCompositeCandidate*);     
      static reco::VertexRef FindV0DecayVertex3(const reco::TrackRef, reco::VertexCompositeCandidate, 	
                                                  edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                  edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, 
	                                          double);
      static reco::VertexRef FindV0DecayVertexZ(const reco::TrackRef, reco::VertexCompositeCandidate, 
                                                   edm::Handle<reco::VertexCollection>, double);

      static std::auto_ptr<reco::PFDisplacedVertexCollection> GetCleanedNI(edm::Handle<reco::PFDisplacedVertexCollection>, edm::Handle<reco::BeamSpot>, bool); 
      static bool ComesFromNI(const reco::TrackRef, reco::PFDisplacedVertexCollection, reco::PFDisplacedVertex*);           
      static reco::VertexRef FindNIVertex3(const reco::TrackRef, reco::PFDisplacedVertex, edm::ESHandle<MagneticField>, const edm::EventSetup&, edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, double, reco::TransientTrack);
      static reco::VertexRef FindNIVertexZ(const reco::TrackRef, reco::PFDisplacedVertex, 
                                                   edm::Handle<reco::VertexCollection>, double);
                                                   
      static std::auto_ptr<reco::VertexCollection> GetCleanedIVF(edm::Handle<reco::VertexCollection>, edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromIVF(const reco::TrackRef, reco::VertexCollection, reco::Vertex*);     
      static reco::VertexRef FindIVFVertex3(const reco::TrackRef, reco::Vertex,edm::ESHandle<MagneticField>, const edm::EventSetup&, edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, double);
      static reco::VertexRef FindIVFVertexZ(const reco::TrackRef, reco::Vertex, edm::Handle<reco::VertexCollection>, double);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      edm::InputTag tcLabel_;
      edm::InputTag vcLabel_;

      edm::InputTag puLabel_;

      edm::InputTag tpLabel_;

      edm::InputTag beamSpotLabel_;

      bool cleanedColls_;

      edm::InputTag ConversionsCollection_;

      edm::InputTag KshortCollection_;
      edm::InputTag LambdaCollection_;

      edm::InputTag NIVertexCollection_;

      edm::InputTag IVFVertexCollection_;
      
      
      TH1F *h_numTrks_rhoVal, *h_numTrks_rhoSig;
      
      //conversion
      
      TH1F *h_conv_numTrks_rhoVal,    *h_conv_numTrks_rhoSig;

      TH3F *h3_conv_PrimZ_pur_rhoVal, *h3_conv_Prim3_pur_rhoVal;
      TH3F *h3_conv_SecZ_pur_rhoVal,  *h3_conv_Sec3_pur_rhoVal;

      TH3F *h3_conv_PrimZ_pur_rhoSig, *h3_conv_Prim3_pur_rhoSig;
      TH3F *h3_conv_SecZ_pur_rhoSig,  *h3_conv_Sec3_pur_rhoSig;
      
      //K decay
      
      TH1F *h_kdec_numTrks_rhoVal,    *h_kdec_numTrks_rhoSig;

      TH3F *h3_kdec_PrimZ_pur_rhoVal, *h3_kdec_Prim3_pur_rhoVal;
      TH3F *h3_kdec_SecZ_pur_rhoVal,  *h3_kdec_Sec3_pur_rhoVal;

      TH3F *h3_kdec_PrimZ_pur_rhoSig, *h3_kdec_Prim3_pur_rhoSig;
      TH3F *h3_kdec_SecZ_pur_rhoSig,  *h3_kdec_Sec3_pur_rhoSig;
      
      //L decay
      
      TH1F *h_ldec_numTrks_rhoVal,    *h_ldec_numTrks_rhoSig;

      TH3F *h3_ldec_PrimZ_pur_rhoVal, *h3_ldec_Prim3_pur_rhoVal;
      TH3F *h3_ldec_SecZ_pur_rhoVal,  *h3_ldec_Sec3_pur_rhoVal;

      TH3F *h3_ldec_PrimZ_pur_rhoSig, *h3_ldec_Prim3_pur_rhoSig;
      TH3F *h3_ldec_SecZ_pur_rhoSig,  *h3_ldec_Sec3_pur_rhoSig;
      
      //nuclear interaction
      
      TH1F *h_nuci_numTrks_rhoVal,    *h_nuci_numTrks_rhoSig;

      TH3F *h3_nuci_PrimZ_pur_rhoVal, *h3_nuci_Prim3_pur_rhoVal;
      TH3F *h3_nuci_SecZ_pur_rhoVal,  *h3_nuci_Sec3_pur_rhoVal;

      TH3F *h3_nuci_PrimZ_pur_rhoSig, *h3_nuci_Prim3_pur_rhoSig;
      TH3F *h3_nuci_SecZ_pur_rhoSig,  *h3_nuci_Sec3_pur_rhoSig;
      
      //inclusive vertex finder
      
      TH1F *h_invf_numTrks_rhoVal,    *h_invf_numTrks_rhoSig;

      TH3F *h3_invf_PrimZ_pur_rhoVal, *h3_invf_Prim3_pur_rhoVal;
      TH3F *h3_invf_SecZ_pur_rhoVal,  *h3_invf_Sec3_pur_rhoVal;

      TH3F *h3_invf_PrimZ_pur_rhoSig, *h3_invf_Prim3_pur_rhoSig;
      TH3F *h3_invf_SecZ_pur_rhoSig,  *h3_invf_Sec3_pur_rhoSig;
};

//
// constants, enums and typedefs
//

using namespace std;
using namespace edm;
using namespace reco;

const double eMass = 0.000511;
const double kMass = 0.49765;
const double lamMass = 1.11568;
const double piMass = 0.1396;

//
// static data member definitions
//

//
// constructors and destructor
//

ReassocTest::ReassocTest(const edm::ParameterSet& iConfig)
{

  double tpMin = -0.5;
  double tpMax = 2.5;
  int tpNbin = 3;
   
  double stepTp = (tpMax-tpMin)/tpNbin; 
  
  Double_t tpIntervalsa[31];  
  
  for (int k=0;k<tpNbin+1;k++) {
    double d = tpMin+k*stepTp;
    tpIntervalsa[k] = d;
  } 

  double rtMin = -0.5;
  double rtMax = 1.5;
  int rtNbin = 2;
   
  double stepRt = (rtMax-rtMin)/rtNbin; 
  
  Double_t rtIntervalsa[31];  
  
  for (int k=0;k<rtNbin+1;k++) {
    double d = rtMin+k*stepRt;
    rtIntervalsa[k] = d;
  } 

  double rhoMin = 0.01;
  double rhoMax = 150.;
  int rhoNbin = 30;

  double rhoMinLog = log10(rhoMin);
  double rhoMaxLog = log10(rhoMax);
   
  double stepRho = (rhoMaxLog-rhoMinLog)/rhoNbin; 
  
  Double_t rhoIntervalsa[31];      
  rhoIntervalsa[0] = rhoMin;
  
  for (int k=1;k<rhoNbin+1;k++) {
    double d=pow(10,rhoMinLog+k*stepRho);
    rhoIntervalsa[k] = d;
  }    


   //now do what ever initialization is needed

  tcLabel_ = iConfig.getParameter<InputTag>("TrackCollection");
  vcLabel_ = iConfig.getParameter<InputTag>("VertexCollection");

  puLabel_ = iConfig.getParameter<InputTag>("PileUp");

  tpLabel_ = iConfig.getParameter<InputTag>("TrackingParticles");

  beamSpotLabel_ = iConfig.getParameter<InputTag>("BeamSpot");

  cleanedColls_ = iConfig.getParameter<bool>("GetCleanedCollections");

  ConversionsCollection_= iConfig.getParameter<InputTag>("ConversionsCollection");

  KshortCollection_= iConfig.getParameter<InputTag>("KshortCollection");
  LambdaCollection_= iConfig.getParameter<InputTag>("LambdaCollection");

  NIVertexCollection_= iConfig.getParameter<InputTag>("NIVertexCollection");

  IVFVertexCollection_= iConfig.getParameter<InputTag>("IVFVertexCollection");

  //--------------


  Service<TFileService> tfs;
  
  h_numTrks_rhoVal = tfs->make<TH1F>("h_numTrks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_numTrks_rhoSig = tfs->make<TH1F>("h_numTrks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_numTrks_rhoVal->Sumw2();
  h_numTrks_rhoSig->Sumw2();  
      
  //conversion
  
  h_conv_numTrks_rhoVal = tfs->make<TH1F>("h_conv_numTrks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_conv_numTrks_rhoSig = tfs->make<TH1F>("h_conv_numTrks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_conv_numTrks_rhoVal->Sumw2();
  h_conv_numTrks_rhoSig->Sumw2();  

  h3_conv_PrimZ_pur_rhoVal = tfs->make<TH3F>("h3_conv_PrimZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Prim3_pur_rhoVal = tfs->make<TH3F>("h3_conv_Prim3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_SecZ_pur_rhoVal = tfs->make<TH3F>("h3_conv_SecZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Sec3_pur_rhoVal = tfs->make<TH3F>("h3_conv_Sec3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_conv_PrimZ_pur_rhoVal->Sumw2();
  h3_conv_Prim3_pur_rhoVal->Sumw2(); 
  h3_conv_SecZ_pur_rhoVal->Sumw2();
  h3_conv_Sec3_pur_rhoVal->Sumw2(); 

  h3_conv_PrimZ_pur_rhoSig = tfs->make<TH3F>("h3_conv_PrimZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Prim3_pur_rhoSig = tfs->make<TH3F>("h3_conv_Prim3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_SecZ_pur_rhoSig = tfs->make<TH3F>("h3_conv_SecZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Sec3_pur_rhoSig = tfs->make<TH3F>("h3_conv_Sec3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_conv_PrimZ_pur_rhoSig->Sumw2();
  h3_conv_Prim3_pur_rhoSig->Sumw2(); 
  h3_conv_SecZ_pur_rhoSig->Sumw2();
  h3_conv_Sec3_pur_rhoSig->Sumw2(); 
      
  //K decay
  
  h_kdec_numTrks_rhoVal = tfs->make<TH1F>("h_kdec_numTrks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_kdec_numTrks_rhoSig = tfs->make<TH1F>("h_kdec_numTrks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_kdec_numTrks_rhoVal->Sumw2();
  h_kdec_numTrks_rhoSig->Sumw2();  

  h3_kdec_PrimZ_pur_rhoVal = tfs->make<TH3F>("h3_kdec_PrimZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Prim3_pur_rhoVal = tfs->make<TH3F>("h3_kdec_Prim3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_SecZ_pur_rhoVal = tfs->make<TH3F>("h3_kdec_SecZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Sec3_pur_rhoVal = tfs->make<TH3F>("h3_kdec_Sec3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_kdec_PrimZ_pur_rhoVal->Sumw2();
  h3_kdec_Prim3_pur_rhoVal->Sumw2(); 
  h3_kdec_SecZ_pur_rhoVal->Sumw2();
  h3_kdec_Sec3_pur_rhoVal->Sumw2(); 

  h3_kdec_PrimZ_pur_rhoSig = tfs->make<TH3F>("h3_kdec_PrimZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Prim3_pur_rhoSig = tfs->make<TH3F>("h3_kdec_Prim3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_SecZ_pur_rhoSig = tfs->make<TH3F>("h3_kdec_SecZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Sec3_pur_rhoSig = tfs->make<TH3F>("h3_kdec_Sec3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_kdec_PrimZ_pur_rhoSig->Sumw2();
  h3_kdec_Prim3_pur_rhoSig->Sumw2(); 
  h3_kdec_SecZ_pur_rhoSig->Sumw2();
  h3_kdec_Sec3_pur_rhoSig->Sumw2(); 
      
  //L decay
  
  h_ldec_numTrks_rhoVal = tfs->make<TH1F>("h_ldec_numTrks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_ldec_numTrks_rhoSig = tfs->make<TH1F>("h_ldec_numTrks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_ldec_numTrks_rhoVal->Sumw2();
  h_ldec_numTrks_rhoSig->Sumw2();  

  h3_ldec_PrimZ_pur_rhoVal = tfs->make<TH3F>("h3_ldec_PrimZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Prim3_pur_rhoVal = tfs->make<TH3F>("h3_ldec_Prim3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_SecZ_pur_rhoVal = tfs->make<TH3F>("h3_ldec_SecZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Sec3_pur_rhoVal = tfs->make<TH3F>("h3_ldec_Sec3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_ldec_PrimZ_pur_rhoVal->Sumw2();
  h3_ldec_Prim3_pur_rhoVal->Sumw2(); 
  h3_ldec_SecZ_pur_rhoVal->Sumw2();
  h3_ldec_Sec3_pur_rhoVal->Sumw2(); 

  h3_ldec_PrimZ_pur_rhoSig = tfs->make<TH3F>("h3_ldec_PrimZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Prim3_pur_rhoSig = tfs->make<TH3F>("h3_ldec_Prim3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_SecZ_pur_rhoSig = tfs->make<TH3F>("h3_ldec_SecZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Sec3_pur_rhoSig = tfs->make<TH3F>("h3_ldec_Sec3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_ldec_PrimZ_pur_rhoSig->Sumw2();
  h3_ldec_Prim3_pur_rhoSig->Sumw2(); 
  h3_ldec_SecZ_pur_rhoSig->Sumw2();
  h3_ldec_Sec3_pur_rhoSig->Sumw2(); 
      
  //nuclear interaction
  
  h_nuci_numTrks_rhoVal = tfs->make<TH1F>("h_nuci_numTrks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_nuci_numTrks_rhoSig = tfs->make<TH1F>("h_nuci_numTrks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_nuci_numTrks_rhoVal->Sumw2();
  h_nuci_numTrks_rhoSig->Sumw2();  

  h3_nuci_PrimZ_pur_rhoVal = tfs->make<TH3F>("h3_nuci_PrimZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Prim3_pur_rhoVal = tfs->make<TH3F>("h3_nuci_Prim3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_SecZ_pur_rhoVal = tfs->make<TH3F>("h3_nuci_SecZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Sec3_pur_rhoVal = tfs->make<TH3F>("h3_nuci_Sec3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_nuci_PrimZ_pur_rhoVal->Sumw2();
  h3_nuci_Prim3_pur_rhoVal->Sumw2(); 
  h3_nuci_SecZ_pur_rhoVal->Sumw2();
  h3_nuci_Sec3_pur_rhoVal->Sumw2(); 

  h3_nuci_PrimZ_pur_rhoSig = tfs->make<TH3F>("h3_nuci_PrimZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Prim3_pur_rhoSig = tfs->make<TH3F>("h3_nuci_Prim3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_SecZ_pur_rhoSig = tfs->make<TH3F>("h3_nuci_SecZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Sec3_pur_rhoSig = tfs->make<TH3F>("h3_nuci_Sec3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_nuci_PrimZ_pur_rhoSig->Sumw2();
  h3_nuci_Prim3_pur_rhoSig->Sumw2(); 
  h3_nuci_SecZ_pur_rhoSig->Sumw2();
  h3_nuci_Sec3_pur_rhoSig->Sumw2(); 
      
  //inclusive vertex finder
  
  h_invf_numTrks_rhoVal = tfs->make<TH1F>("h_invf_numTrks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_invf_numTrks_rhoSig = tfs->make<TH1F>("h_invf_numTrks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_invf_numTrks_rhoVal->Sumw2();
  h_invf_numTrks_rhoSig->Sumw2();  

  h3_invf_PrimZ_pur_rhoVal = tfs->make<TH3F>("h3_invf_PrimZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Prim3_pur_rhoVal = tfs->make<TH3F>("h3_invf_Prim3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_SecZ_pur_rhoVal = tfs->make<TH3F>("h3_invf_SecZ_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Sec3_pur_rhoVal = tfs->make<TH3F>("h3_invf_Sec3_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_invf_PrimZ_pur_rhoVal->Sumw2();
  h3_invf_Prim3_pur_rhoVal->Sumw2(); 
  h3_invf_SecZ_pur_rhoVal->Sumw2();
  h3_invf_Sec3_pur_rhoVal->Sumw2(); 

  h3_invf_PrimZ_pur_rhoSig = tfs->make<TH3F>("h3_invf_PrimZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Prim3_pur_rhoSig = tfs->make<TH3F>("h3_invf_Prim3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_SecZ_pur_rhoSig = tfs->make<TH3F>("h3_invf_SecZ_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Sec3_pur_rhoSig = tfs->make<TH3F>("h3_invf_Sec3_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_invf_PrimZ_pur_rhoSig->Sumw2();
  h3_invf_Prim3_pur_rhoSig->Sumw2(); 
  h3_invf_SecZ_pur_rhoSig->Sumw2();
  h3_invf_Sec3_pur_rhoSig->Sumw2(); 
  
}


ReassocTest::~ReassocTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool 
ReassocTest::TrackWeightAssociation(const TrackRef&  trackRef, Handle<VertexCollection> vtxcollH) 
{

	VertexRef bestvertexref(vtxcollH,0);		
 	float bestweight = 0.;

	const TrackBaseRef& trackbaseRef = TrackBaseRef(trackRef);

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

  	return (bestweight > 1e-5);

}

VertexRef
ReassocTest::FindVertexZ(TrackRef trkref, Handle<VertexCollection> vtxCollH, double trackWeight)
{

	double ztrack = trkref->vertex().z();

	VertexRef foundVertexRef(vtxCollH,0);

	double dzmin = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

      VertexRef vertexref(vtxCollH,index_vtx);

	  double nTracks = vertexref->tracksSize();
 
	  //find and store the closest vertex in z
      double distance = fabs(ztrack - vertexref->z());

	  double weightedDistance = distance-trackWeight*nTracks;	

          if(weightedDistance<dzmin) {
            dzmin = weightedDistance; 
            foundVertexRef = vertexref;
          }
	
	}

	return foundVertexRef;

}

VertexRef
ReassocTest::FindVertex3(TransientTrack transtrk, Handle<VertexCollection> vtxCollH, double trackWeight)
{

	VertexRef foundVertexRef(vtxCollH,0);

	double d3min = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

      VertexRef vertexref(vtxCollH,index_vtx);

	  double nTracks = sqrt(vertexref->tracksSize());

      double distance = 1e5;	        
      pair<bool,Measurement1D> IpPair = IPTools::absoluteImpactParameter3D(transtrk, *vertexref);
 
	  if ( IpPair.first ) distance = IpPair.second.value();

	  double weightedDistance = distance-trackWeight*nTracks;	

      if(weightedDistance<d3min) {
        d3min = weightedDistance; 
        foundVertexRef = vertexref;
      }
	
	}

	return foundVertexRef;

}

auto_ptr<ConversionCollection> 
ReassocTest::GetCleanedConversions(Handle<ConversionCollection> convCollH, Handle<BeamSpot> bsH, bool cleanedColl)
{

    auto_ptr<ConversionCollection> cleanedConvColl(new ConversionCollection() );

	for (unsigned int convcoll_idx=0; convcoll_idx<convCollH->size(); convcoll_idx++){

	  ConversionRef convref(convCollH,convcoll_idx);

 	  if(!cleanedColl){   
        cleanedConvColl->push_back(*convref);
	    continue;
      }

	  if( (convref->nTracks()==2) &&
          (fabs(convref->pairInvariantMass())<=0.1) ){
    
        cleanedConvColl->push_back(*convref);

	  }

	}

  	return cleanedConvColl;

}

bool
ReassocTest::ComesFromConversion(const TrackRef trackref, ConversionCollection cleanedConvColl, Conversion* gamma)
{

	for(unsigned int convcoll_ite=0; convcoll_ite<cleanedConvColl.size(); convcoll_ite++){

	  if(ConversionTools::matchesConversion(trackref,cleanedConvColl.at(convcoll_ite))){
	
	    *gamma = cleanedConvColl.at(convcoll_ite);
	    return true;

  	  }

  	}

	return false;
}

VertexRef
ReassocTest::FindConversionVertex3(const reco::TrackRef trackref, reco::Conversion gamma, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

	math::XYZPoint conv_pos = gamma.conversionVertex().position();

	math::XYZVector conv_mom(gamma.refittedPair4Momentum().x(),
	                         gamma.refittedPair4Momentum().y(),
	                         gamma.refittedPair4Momentum().z());

	Track photon(trackref->chi2(), trackref->ndof(), conv_pos, conv_mom, 0, trackref->covariance());

    TransientTrack transpho(photon, &(*bFieldH) );
    transpho.setBeamSpot(*bsH);
   	transpho.setES(iSetup);

	return ReassocTest::FindVertex3(transpho, vtxcollH, tWeight); 

}

VertexRef
ReassocTest::FindConversionVertexZ(const reco::TrackRef trackref, reco::Conversion gamma, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

	math::XYZPoint conv_pos = gamma.conversionVertex().position();

	math::XYZVector conv_mom(gamma.refittedPair4Momentum().x(),
	                         gamma.refittedPair4Momentum().y(),
	                         gamma.refittedPair4Momentum().z());

    double ztrackfirst = conv_pos.z();
	double radius = conv_pos.rho();
	double tracktheta = conv_mom.theta();

	double ztrack = ztrackfirst - (radius/tan(tracktheta));

	VertexRef bestvertexref(vtxcollH, 0);

	double dzmin = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollH->size(); ++index_vtx){

      VertexRef vertexref(vtxcollH,index_vtx);

	  int nTracks = vertexref->tracksSize();
 
	  //find and store the closest vertex in z
      double distance = fabs(ztrack - vertexref->z());

	  double weightedDistance = distance-tWeight*nTracks;	

      if(weightedDistance<dzmin) {
        dzmin = weightedDistance; 
        bestvertexref = vertexref;
      }
	
	}

	return bestvertexref;

}

auto_ptr<VertexCompositeCandidateCollection>
ReassocTest::GetCleanedKshort(Handle<VertexCompositeCandidateCollection> KshortsH, Handle<BeamSpot> bsH, bool cleanedColl)
{

    auto_ptr<VertexCompositeCandidateCollection> cleanedKaonColl(new VertexCompositeCandidateCollection() );

	for (unsigned int kscoll_idx=0; kscoll_idx<KshortsH->size(); kscoll_idx++){

	  VertexCompositeCandidateRef ksref(KshortsH,kscoll_idx);

 	  if(!cleanedColl){   
        cleanedKaonColl->push_back(*ksref);
	    continue;
	  }

	}

	return cleanedKaonColl;

}

auto_ptr<VertexCompositeCandidateCollection>
ReassocTest::GetCleanedLambda(Handle<VertexCompositeCandidateCollection> LambdasH, Handle<BeamSpot> bsH, bool cleanedColl)
{

    auto_ptr<VertexCompositeCandidateCollection> cleanedLambdaColl(new VertexCompositeCandidateCollection() );

	for (unsigned int lambdacoll_idx=0; lambdacoll_idx<LambdasH->size(); lambdacoll_idx++){

	  VertexCompositeCandidateRef lambdaref(LambdasH,lambdacoll_idx);

 	  if(!cleanedColl){   
        cleanedLambdaColl->push_back(*lambdaref);
	    continue;
      }

	}

	return cleanedLambdaColl;
}

bool
ReassocTest::ComesFromDecay(const TrackRef trackref, VertexCompositeCandidateCollection cleanedVCCC, VertexCompositeCandidate* V0)
{

	//the part for the reassociation of particles from Kshort decays
	for(VertexCompositeCandidateCollection::const_iterator iKS=cleanedVCCC.begin(); iKS!=cleanedVCCC.end(); iKS++){

	  const RecoChargedCandidate *dauCand1 = dynamic_cast<const RecoChargedCandidate*>(iKS->daughter(0));
 	  TrackRef dauTk1 = dauCand1->track();
	  const RecoChargedCandidate *dauCand2 = dynamic_cast<const RecoChargedCandidate*>(iKS->daughter(1));
 	  TrackRef dauTk2 = dauCand2->track();

	  if((trackref==dauTk1) || (trackref==dauTk2)){
	  
	    *V0 = *iKS; 
	    return true;

	  }

	}

	return false;
}

VertexRef
ReassocTest::FindV0DecayVertex3(const reco::TrackRef trackref, reco::VertexCompositeCandidate V0_vtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

	math::XYZPoint dec_pos = V0_vtx.vertex();

	math::XYZVector dec_mom(V0_vtx.momentum().x(),
	                        V0_vtx.momentum().y(),
	                        V0_vtx.momentum().z());

	Track V0(trackref->chi2(), trackref->ndof(), dec_pos, dec_mom, 0, trackref->covariance());

    TransientTrack transV0(V0, &(*bFieldH) );
    transV0.setBeamSpot(*bsH);
    transV0.setES(iSetup);

	return ReassocTest::FindVertex3(transV0, vtxcollH, tWeight); 

}

VertexRef
ReassocTest::FindV0DecayVertexZ(const reco::TrackRef trackref, reco::VertexCompositeCandidate V0_vtx, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

	math::XYZPoint dec_pos = V0_vtx.vertex();

	math::XYZVector dec_mom(V0_vtx.momentum().x(),
	                        V0_vtx.momentum().y(),
	                        V0_vtx.momentum().z());

    double ztrackfirst = dec_pos.z();
	double radius = dec_pos.rho();
	double tracktheta = dec_mom.theta();

	double ztrack = ztrackfirst - (radius/tan(tracktheta));

	VertexRef bestvertexref(vtxcollH, 0);

	double dzmin = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollH->size(); ++index_vtx){

      VertexRef vertexref(vtxcollH,index_vtx);

	  int nTracks = vertexref->tracksSize();
 
	  //find and store the closest vertex in z
      double distance = fabs(ztrack - vertexref->z());

	  double weightedDistance = distance-tWeight*nTracks;	

      if(weightedDistance<dzmin) {
        dzmin = weightedDistance; 
        bestvertexref = vertexref;
      }
	
	}

	return bestvertexref;

}

auto_ptr<PFDisplacedVertexCollection>
ReassocTest::GetCleanedNI(Handle<PFDisplacedVertexCollection> NuclIntH, Handle<BeamSpot> bsH, bool cleanedColl)
{

    auto_ptr<PFDisplacedVertexCollection> cleanedNIColl(new PFDisplacedVertexCollection() );

	for (PFDisplacedVertexCollection::const_iterator niref=NuclIntH->begin(); niref!=NuclIntH->end(); niref++){

	  if(!cleanedColl){
	    cleanedNIColl->push_back(*niref);
	    continue;
      }

	}

	return cleanedNIColl;
}

bool
ReassocTest::ComesFromNI(const TrackRef trackref, PFDisplacedVertexCollection cleanedNI, PFDisplacedVertex* displVtx)
{

	//the part for the reassociation of particles from nuclear interactions
	for(PFDisplacedVertexCollection::const_iterator iDisplV=cleanedNI.begin(); iDisplV!=cleanedNI.end(); iDisplV++){

	  if(iDisplV->trackWeight(trackref)>1.e-2){
	  
	    *displVtx = *iDisplV; 
	    return true;

	  }

	}

	return false;
}

VertexRef
ReassocTest::FindNIVertex3(const reco::TrackRef trackref, reco::PFDisplacedVertex displVtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight, TransientTrack transhelp)
{ 

	TrackCollection refittedTracks = displVtx.refittedTracks();

	if((displVtx.isTherePrimaryTracks()) || (displVtx.isThereMergedTracks())){

	  for(TrackCollection::const_iterator trkcoll_ite=refittedTracks.begin(); trkcoll_ite!=refittedTracks.end(); trkcoll_ite++){
	
	    const TrackBaseRef retrackbaseref = displVtx.originalTrack(*trkcoll_ite); 

	    if(displVtx.isIncomingTrack(retrackbaseref)){

	      VertexRef bestvertexref(vtxcollH, 0);
	      float bestweight = 0.;

  	      for(unsigned int index_vtx=0;  index_vtx<vtxcollH->size(); ++index_vtx){

            VertexRef vertexref(vtxcollH,index_vtx);

     	    //get the most probable vertex for the track
	        float weight = vertexref->trackWeight(retrackbaseref);
	        if(weight>bestweight){
    	      bestweight = weight;
	          bestvertexref = vertexref;
 	        } 
 
	      }

	      if(bestweight>1.e-5) return bestvertexref;

    	  TransientTrack transIncom(*retrackbaseref, &(*bFieldH) );
    	  transIncom.setBeamSpot(*bsH);
    	  transIncom.setES(iSetup);

	      return FindVertex3(transIncom, vtxcollH, tWeight); 

	    }

	  }

	}

 	return FindVertex3(transhelp, vtxcollH, tWeight); 

}

VertexRef
ReassocTest::FindNIVertexZ(const reco::TrackRef trackref, reco::PFDisplacedVertex displVtx, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

	math::XYZPoint ni_pos = displVtx.position();

	math::XYZVector ni_mom(displVtx.primaryMomentum().x(),
	                       displVtx.primaryMomentum().y(),
	                       displVtx.primaryMomentum().z());

    double ztrackfirst = ni_pos.z();
	double radius = ni_pos.rho();
	double tracktheta = ni_mom.theta();

	double ztrack = ztrackfirst - (radius/tan(tracktheta));

	VertexRef bestvertexref(vtxcollH, 0);

	double dzmin = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollH->size(); ++index_vtx){

      VertexRef vertexref(vtxcollH,index_vtx);

	  int nTracks = vertexref->tracksSize();
 
	  //find and store the closest vertex in z
      double distance = fabs(ztrack - vertexref->z());

	  double weightedDistance = distance-tWeight*nTracks;	

      if(weightedDistance<dzmin) {
        dzmin = weightedDistance; 
        bestvertexref = vertexref;
      }
	
	}

	return bestvertexref;

}

auto_ptr<VertexCollection>
ReassocTest::GetCleanedIVF(Handle<VertexCollection> ifvH, Handle<BeamSpot> bsH, bool cleanedColl)
{

    auto_ptr<VertexCollection> cleanedIVFColl(new VertexCollection() );

	for (VertexCollection::const_iterator ivfref=ifvH->begin(); ivfref!=ifvH->end(); ivfref++){

	  if ( !cleanedColl ) {
	    cleanedIVFColl->push_back(*ivfref);
	    continue;
      }            

	}

	return cleanedIVFColl;
}

bool
ReassocTest::ComesFromIVF(const TrackRef trackref, VertexCollection cleanedIVF, Vertex* ivfVtx)
{

	for(VertexCollection::const_iterator iInclV=cleanedIVF.begin(); iInclV!=cleanedIVF.end(); iInclV++){

	  if(iInclV->trackWeight(trackref)>1.e-5){
	  
	    *ivfVtx = *iInclV; 
	    return true;

	  }

	}

	return false;
}

VertexRef
ReassocTest::FindIVFVertex3(const TrackRef trackref, Vertex ivfVtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, Handle<BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{

	math::XYZPoint iv_pos = ivfVtx.position();

	math::XYZVector iv_mom( ivfVtx.p4(0.1, 0.).x(),
	                        ivfVtx.p4(0.1, 0.).y(),
	                        ivfVtx.p4(0.1, 0.).z() );

	Track incom(trackref->chi2(), trackref->ndof(), iv_pos, iv_mom, 0, trackref->covariance());

	TransientTrack transIncom(incom, &(*bFieldH) );
    transIncom.setBeamSpot(*bsH);
    transIncom.setES(iSetup);

	return ReassocTest::FindVertex3(transIncom, vtxcollH, tWeight); 
}

VertexRef
ReassocTest::FindIVFVertexZ(const reco::TrackRef trackref, Vertex ivfVtx, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{

	math::XYZPoint iv_pos = ivfVtx.position();


	math::XYZVector iv_mom( ivfVtx.p4(0.1, 0.).x(),
	                        ivfVtx.p4(0.1, 0.).y(),
	                        ivfVtx.p4(0.1, 0.).z() );

    double ztrackfirst = iv_pos.z();
	double radius = iv_pos.rho();
	double tracktheta = iv_mom.theta();

	double ztrack = ztrackfirst - (radius/tan(tracktheta));

	VertexRef bestvertexref(vtxcollH, 0);

	double dzmin = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollH->size(); ++index_vtx){

      VertexRef vertexref(vtxcollH,index_vtx);

	  int nTracks = vertexref->tracksSize();
 
	  //find and store the closest vertex in z
      double distance = fabs(ztrack - vertexref->z());

	  double weightedDistance = distance-tWeight*nTracks;	

      if(weightedDistance<dzmin) {
        dzmin = weightedDistance; 
        bestvertexref = vertexref;
      }
	
	}

	return bestvertexref;
}

// ------------ method called for each event  ------------
void
ReassocTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //associate reco tracks to tracking particles
  ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get("quickTrackAssociatorByHits",theAssociator);
  TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product();

  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

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

  //get the offline beam spot
  Handle<BeamSpot>  bsH;
  iEvent.getByLabel(beamSpotLabel_,bsH);

  //get the conversion collection for the gamma conversions
  Handle<ConversionCollection> convCollH;
  iEvent.getByLabel(ConversionsCollection_, convCollH);
  auto_ptr<ConversionCollection> cleanedConvCollP = ReassocTest::GetCleanedConversions(convCollH,bsH,cleanedColls_);

  //get the vertex composite candidate collection for the Kshort's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH;
  iEvent.getByLabel(KshortCollection_, vertCompCandCollKshortH);
  auto_ptr<VertexCompositeCandidateCollection> cleanedKshortCollP = ReassocTest::GetCleanedKshort(vertCompCandCollKshortH,bsH,cleanedColls_);
  
  //get the vertex composite candidate collection for the Lambda's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
  iEvent.getByLabel(LambdaCollection_, vertCompCandCollLambdaH);
  auto_ptr<VertexCompositeCandidateCollection> cleanedLambdaCollP = ReassocTest::GetCleanedLambda(vertCompCandCollLambdaH,bsH,cleanedColls_);

  //get the displaced vertex collection for nuclear interactions
  Handle<PFDisplacedVertexCollection> displVertexCollH;
  iEvent.getByLabel(NIVertexCollection_,displVertexCollH);
  auto_ptr<PFDisplacedVertexCollection> cleanedNICollP = ReassocTest::GetCleanedNI(displVertexCollH,bsH,cleanedColls_);


  //get the vertex collection for inclusive vertex finder   
  Handle<VertexCollection> ivfVertexCollH;
  iEvent.getByLabel(IVFVertexCollection_, ivfVertexCollH);
  auto_ptr<VertexCollection> cleanedIVFCollP = ReassocTest::GetCleanedIVF(ivfVertexCollH, bsH, cleanedColls_);

	    
  //loop over all tracks of the track collection	
  for ( size_t idxTrack = 0; idxTrack < theTracksH->size(); ++idxTrack ) {

    TrackRef trackref(theTracksH, idxTrack);
    TrackBaseRef trackbaseref = TrackBaseRef(trackref);

    if ( ReassocTest::TrackWeightAssociation(trackref, theVerticesH) ) continue;

    TransientTrack transtrk(trackref, &(*bFieldHandle) );
    transtrk.setBeamSpot(*bsH);
    transtrk.setES(iSetup);

    vector<pair<TrackingParticleRef, double> > tp;
    if(recSimColl.find(trackbaseref) != recSimColl.end()) tp = recSimColl[trackbaseref];
  
    int x_bin = 1;

    if (tp.size()==0){ 
      x_bin=2;
    }else{     
      TrackingParticleRef tpr = tp[0].first;
      if ((tpr->eventId().event() == 0) && (tpr->eventId().bunchCrossing() == 0)) x_bin=0;
    }

    double rhoVal = transtrk.stateAtBeamLine().transverseImpactParameter().value();
    double rhoSig = transtrk.stateAtBeamLine().transverseImpactParameter().significance();
    
    h_numTrks_rhoVal->Fill( rhoVal );
    h_numTrks_rhoSig->Fill( rhoSig );

    //conversion
    Conversion gamma;
    if ( ReassocTest::ComesFromConversion(trackref, *cleanedConvCollP, &gamma) ){
    
      h_conv_numTrks_rhoVal->Fill( rhoVal );
      h_conv_numTrks_rhoSig->Fill( rhoSig );

      int y_bin_PrimZ = 1;
      if(ReassocTest::FindConversionVertexZ(trackref, gamma, theVerticesH, 0.001)==firstVertex) y_bin_PrimZ = 0;

      int y_bin_Prim3 = 1;
      if(ReassocTest::FindConversionVertex3(trackref, gamma, bFieldHandle, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_Prim3 = 0;

      int y_bin_SecZ = 1;
      if(ReassocTest::FindVertexZ(trackref,theVerticesH,0.001)==firstVertex) y_bin_SecZ=0;

      int y_bin_Sec3 = 1;
      if(ReassocTest::FindVertex3(transtrk,theVerticesH,0.01)==firstVertex) y_bin_Sec3=0;

      h3_conv_PrimZ_pur_rhoVal->Fill( x_bin, y_bin_PrimZ, rhoVal );
      h3_conv_Prim3_pur_rhoVal->Fill( x_bin, y_bin_Prim3, rhoVal );
      h3_conv_SecZ_pur_rhoVal->Fill( x_bin, y_bin_SecZ, rhoVal );
      h3_conv_Sec3_pur_rhoVal->Fill( x_bin, y_bin_Sec3, rhoVal );

      h3_conv_PrimZ_pur_rhoSig->Fill( x_bin, y_bin_PrimZ, rhoSig );
      h3_conv_Prim3_pur_rhoSig->Fill( x_bin, y_bin_Prim3, rhoSig );
      h3_conv_SecZ_pur_rhoSig->Fill( x_bin, y_bin_SecZ, rhoSig );
      h3_conv_Sec3_pur_rhoSig->Fill( x_bin, y_bin_Sec3, rhoSig );

    }

    //K decay
    VertexCompositeCandidate Ks;
    if ( ReassocTest::ComesFromDecay(trackref, *cleanedKshortCollP, &Ks) ) {
    
      h_kdec_numTrks_rhoVal->Fill( rhoVal );
      h_kdec_numTrks_rhoSig->Fill( rhoSig );

      int y_bin_PrimZ = 1;
      if(ReassocTest::FindV0DecayVertexZ(trackref, Ks, theVerticesH, 0.001)==firstVertex) y_bin_PrimZ = 0;

      int y_bin_Prim3 = 1;
      if(ReassocTest::FindV0DecayVertex3(trackref, Ks, bFieldHandle, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_Prim3 = 0;

      int y_bin_SecZ = 1;
      if(ReassocTest::FindVertexZ(trackref,theVerticesH,0.001)==firstVertex) y_bin_SecZ=0;

      int y_bin_Sec3 = 1;
      if(ReassocTest::FindVertex3(transtrk,theVerticesH,0.01)==firstVertex) y_bin_Sec3=0;

      h3_kdec_PrimZ_pur_rhoVal->Fill( x_bin, y_bin_PrimZ, rhoVal );
      h3_kdec_Prim3_pur_rhoVal->Fill( x_bin, y_bin_Prim3, rhoVal );
      h3_kdec_SecZ_pur_rhoVal->Fill( x_bin, y_bin_SecZ, rhoVal );
      h3_kdec_Sec3_pur_rhoVal->Fill( x_bin, y_bin_Sec3, rhoVal );

      h3_kdec_PrimZ_pur_rhoSig->Fill( x_bin, y_bin_PrimZ, rhoSig );
      h3_kdec_Prim3_pur_rhoSig->Fill( x_bin, y_bin_Prim3, rhoSig );
      h3_kdec_SecZ_pur_rhoSig->Fill( x_bin, y_bin_SecZ, rhoSig );
      h3_kdec_Sec3_pur_rhoSig->Fill( x_bin, y_bin_Sec3, rhoSig );

    }

    //L decay
    VertexCompositeCandidate La;
    if ( ReassocTest::ComesFromDecay(trackref, *cleanedLambdaCollP, &La) ) {
    
      h_ldec_numTrks_rhoVal->Fill( rhoVal );
      h_ldec_numTrks_rhoSig->Fill( rhoSig );

      int y_bin_PrimZ = 1;
      if(ReassocTest::FindV0DecayVertexZ(trackref, La, theVerticesH, 0.001)==firstVertex) y_bin_PrimZ = 0;

      int y_bin_Prim3 = 1;
      if(ReassocTest::FindV0DecayVertex3(trackref, La, bFieldHandle, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_Prim3 = 0;

      int y_bin_SecZ = 1;
      if(ReassocTest::FindVertexZ(trackref,theVerticesH,0.001)==firstVertex) y_bin_SecZ=0;

      int y_bin_Sec3 = 1;
      if(ReassocTest::FindVertex3(transtrk,theVerticesH,0.01)==firstVertex) y_bin_Sec3=0;

      h3_ldec_PrimZ_pur_rhoVal->Fill( x_bin, y_bin_PrimZ, rhoVal );
      h3_ldec_Prim3_pur_rhoVal->Fill( x_bin, y_bin_Prim3, rhoVal );
      h3_ldec_SecZ_pur_rhoVal->Fill( x_bin, y_bin_SecZ, rhoVal );
      h3_ldec_Sec3_pur_rhoVal->Fill( x_bin, y_bin_Sec3, rhoVal );

      h3_ldec_PrimZ_pur_rhoSig->Fill( x_bin, y_bin_PrimZ, rhoSig );
      h3_ldec_Prim3_pur_rhoSig->Fill( x_bin, y_bin_Prim3, rhoSig );
      h3_ldec_SecZ_pur_rhoSig->Fill( x_bin, y_bin_SecZ, rhoSig );
      h3_ldec_Sec3_pur_rhoSig->Fill( x_bin, y_bin_Sec3, rhoSig );
    }

    //Nuclear interactions
    PFDisplacedVertex displVtx;
    if ( ReassocTest::ComesFromNI(trackref, *cleanedNICollP, &displVtx) ) {
    
      h_nuci_numTrks_rhoVal->Fill( rhoVal );
      h_nuci_numTrks_rhoSig->Fill( rhoSig );

      int y_bin_PrimZ = 1;
      if(ReassocTest::FindNIVertexZ(trackref, displVtx, theVerticesH, 0.001)==firstVertex) y_bin_PrimZ = 0;

      int y_bin_Prim3 = 1;
      if(ReassocTest::FindNIVertex3(trackref, displVtx, bFieldHandle, iSetup, bsH, theVerticesH, 0.01, transtrk)==firstVertex) y_bin_Prim3 = 0;

      int y_bin_SecZ = 1;
      if(ReassocTest::FindVertexZ(trackref,theVerticesH,0.001)==firstVertex) y_bin_SecZ=0;

      int y_bin_Sec3 = 1;
      if(ReassocTest::FindVertex3(transtrk,theVerticesH,0.01)==firstVertex) y_bin_Sec3=0;

      h3_nuci_PrimZ_pur_rhoVal->Fill( x_bin, y_bin_PrimZ, rhoVal );
      h3_nuci_Prim3_pur_rhoVal->Fill( x_bin, y_bin_Prim3, rhoVal );
      h3_nuci_SecZ_pur_rhoVal->Fill( x_bin, y_bin_SecZ, rhoVal );
      h3_nuci_Sec3_pur_rhoVal->Fill( x_bin, y_bin_Sec3, rhoVal );

      h3_nuci_PrimZ_pur_rhoSig->Fill( x_bin, y_bin_PrimZ, rhoSig );
      h3_nuci_Prim3_pur_rhoSig->Fill( x_bin, y_bin_Prim3, rhoSig );
      h3_nuci_SecZ_pur_rhoSig->Fill( x_bin, y_bin_SecZ, rhoSig );
      h3_nuci_Sec3_pur_rhoSig->Fill( x_bin, y_bin_Sec3, rhoSig );

    }

    //Inclusive vertex
	Vertex ivfVtx;
    if ( ReassocTest::ComesFromIVF(trackref, *cleanedIVFCollP, &ivfVtx) ) {
    
      h_invf_numTrks_rhoVal->Fill( rhoVal );
      h_invf_numTrks_rhoSig->Fill( rhoSig );

      int y_bin_PrimZ = 1;
      if(ReassocTest::FindIVFVertexZ(trackref, ivfVtx, theVerticesH, 0.001)==firstVertex) y_bin_PrimZ = 0;

      int y_bin_Prim3 = 1;
      if(ReassocTest::FindIVFVertex3(trackref, ivfVtx, bFieldHandle, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_Prim3 = 0;

      int y_bin_SecZ = 1;
      if(ReassocTest::FindVertexZ(trackref,theVerticesH,0.001)==firstVertex) y_bin_SecZ=0;

      int y_bin_Sec3 = 1;
      if(ReassocTest::FindVertex3(transtrk,theVerticesH,0.01)==firstVertex) y_bin_Sec3=0;

      h3_invf_PrimZ_pur_rhoVal->Fill( x_bin, y_bin_PrimZ, rhoVal );
      h3_invf_Prim3_pur_rhoVal->Fill( x_bin, y_bin_Prim3, rhoVal );
      h3_invf_SecZ_pur_rhoVal->Fill( x_bin, y_bin_SecZ, rhoVal );
      h3_invf_Sec3_pur_rhoVal->Fill( x_bin, y_bin_Sec3, rhoVal );

      h3_invf_PrimZ_pur_rhoSig->Fill( x_bin, y_bin_PrimZ, rhoSig );
      h3_invf_Prim3_pur_rhoSig->Fill( x_bin, y_bin_Prim3, rhoSig );
      h3_invf_SecZ_pur_rhoSig->Fill( x_bin, y_bin_SecZ, rhoSig );
      h3_invf_Sec3_pur_rhoSig->Fill( x_bin, y_bin_Sec3, rhoSig );

    }

  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ReassocTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ReassocTest);