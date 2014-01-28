// -*- C++ -*-
//
// Package:    UnCleanedComp
// Class:      UnCleanedComp
// 
/**\class UnCleanedComp UnCleanedComp.cc MGeisler/UnCleanedComp/src/UnCleanedComp.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Thu Sep  6 14:15:40 CEST 2012
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
#include "TH3F.h"
#include "TH1F.h"
#include "TProfile3D.h"

//
// class declaration
//

class UnCleanedComp : public edm::EDAnalyzer {
   public:
      explicit UnCleanedComp(const edm::ParameterSet&);
      ~UnCleanedComp();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static bool TrackWeightAssociation(const reco::TrackRef&, edm::Handle<reco::VertexCollection>);

      static reco::VertexRef FindVertexZ(reco::TrackRef, edm::Handle<reco::VertexCollection>, double);
      static reco::VertexRef FindVertex3(reco::TransientTrack, edm::Handle<reco::VertexCollection>, double);

      static double dR(math::XYZPoint, math::XYZVector, edm::Handle<reco::BeamSpot>);

      virtual std::auto_ptr<reco::ConversionCollection> GetCleanedConversions(edm::Handle<reco::ConversionCollection>, 
                                                                             edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromConversion(const reco::TrackRef, reco::ConversionCollection, reco::Conversion*);        
      static reco::VertexRef FindConversionVertex3(const reco::TrackRef, reco::Conversion, 	
                                                  edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                  edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, 
	                                          double);

      virtual std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedKshort(edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);
      virtual std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedLambda(edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromV0Decay(const reco::TrackRef, reco::VertexCompositeCandidateCollection, 
	 	 	  	   reco::VertexCompositeCandidate*);        
      static reco::VertexRef FindV0DecayVertex3(const reco::TrackRef, reco::VertexCompositeCandidate, 	
                                                  edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                  edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, 
	                                          double);

      virtual std::auto_ptr<reco::PFDisplacedVertexCollection> GetCleanedNI(edm::Handle<reco::PFDisplacedVertexCollection>, edm::Handle<reco::BeamSpot>, bool); 
      static bool ComesFromNI(const reco::TrackRef, reco::PFDisplacedVertexCollection, reco::PFDisplacedVertex*);      
      static reco::VertexRef FindNIVertex3(const reco::TrackRef, reco::PFDisplacedVertex, edm::ESHandle<MagneticField>, const edm::EventSetup&, edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, double, reco::TransientTrack);
                                                   
      virtual std::auto_ptr<reco::VertexCollection> GetCleanedIVF(edm::Handle<reco::VertexCollection>, edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromIVF(const reco::TrackRef, reco::VertexCollection, reco::Vertex*);     
      static reco::VertexRef FindIVFVertex3(const reco::TrackRef, reco::Vertex,edm::ESHandle<MagneticField>, const edm::EventSetup&, edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, double);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      edm::InputTag tcLabel_;
      edm::InputTag vcLabel_;

      edm::InputTag tpLabel_;

      edm::InputTag beamSpotLabel_;

      edm::InputTag ConversionsCollection_;

      edm::InputTag KshortCollection_;
      edm::InputTag LambdaCollection_;

      edm::InputTag NIVertexCollection_;

      edm::InputTag IVFVertexCollection_;

      TH1F *h_conv_Clean_NumTracks_rhoVal, *h_conv_Clean_NumTracks_rhoSig;
      TH1F *h_conv_Unclean_NumTracks_rhoVal, *h_conv_Unclean_NumTracks_rhoSig;
      
      TH3F *h3_conv_Clean_pur_rhoVal, *h3_conv_Unclean_pur_rhoVal;
      TH3F *h3_conv_Clean_pur_rhoSig, *h3_conv_Unclean_pur_rhoSig;
      

      TH1F *h_kdec_Clean_NumTracks_rhoVal, *h_kdec_Clean_NumTracks_rhoSig;
      TH1F *h_kdec_Unclean_NumTracks_rhoVal, *h_kdec_Unclean_NumTracks_rhoSig;

      TH3F *h3_kdec_Clean_pur_rhoVal, *h3_kdec_Unclean_pur_rhoVal;
      TH3F *h3_kdec_Clean_pur_rhoSig, *h3_kdec_Unclean_pur_rhoSig;
  
      TProfile3D *p3_kdec_pur, *p3_kdec_eff;
      

      TH1F *h_ldec_Clean_NumTracks_rhoVal, *h_ldec_Clean_NumTracks_rhoSig;
      TH1F *h_ldec_Unclean_NumTracks_rhoVal, *h_ldec_Unclean_NumTracks_rhoSig;

      TH3F *h3_ldec_Clean_pur_rhoVal, *h3_ldec_Unclean_pur_rhoVal;
      TH3F *h3_ldec_Clean_pur_rhoSig, *h3_ldec_Unclean_pur_rhoSig;
  
      TProfile3D *p3_ldec_pur, *p3_ldec_eff;
      

      TH1F *h_nuci_Clean_NumTracks_rhoVal, *h_nuci_Clean_NumTracks_rhoSig;
      TH1F *h_nuci_Unclean_NumTracks_rhoVal, *h_nuci_Unclean_NumTracks_rhoSig;

      TH3F *h3_nuci_Clean_pur_rhoVal, *h3_nuci_Unclean_pur_rhoVal;
      TH3F *h3_nuci_Clean_pur_rhoSig, *h3_nuci_Unclean_pur_rhoSig;
  
      TProfile3D *p3_nuci_pur, *p3_nuci_eff;
      

      TH1F *h_invf_Clean_NumTracks_rhoVal, *h_invf_Clean_NumTracks_rhoSig;
      TH1F *h_invf_Unclean_NumTracks_rhoVal, *h_invf_Unclean_NumTracks_rhoSig;

      TH3F *h3_invf_Clean_pur_rhoVal, *h3_invf_Unclean_pur_rhoVal;
      TH3F *h3_invf_Clean_pur_rhoSig, *h3_invf_Unclean_pur_rhoSig;
      
      TProfile3D *p3_invf_pur, *p3_invf_eff;
      

      VertexDistanceXY distanceComputerXY;
      VertexDistance3D distanceComputer3D;
      
      VertexState BSVertexState;
      
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
UnCleanedComp::UnCleanedComp(const edm::ParameterSet& iConfig)
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

  tpLabel_ = iConfig.getParameter<InputTag>("TrackingParticles");

  beamSpotLabel_ = iConfig.getParameter<InputTag>("BeamSpot");

  ConversionsCollection_= iConfig.getParameter<InputTag>("ConversionsCollection");

  KshortCollection_= iConfig.getParameter<InputTag>("KshortCollection");
  LambdaCollection_= iConfig.getParameter<InputTag>("LambdaCollection");

  NIVertexCollection_= iConfig.getParameter<InputTag>("NIVertexCollection");

  IVFVertexCollection_= iConfig.getParameter<InputTag>("IVFVertexCollection");

  //--------------


  Service<TFileService> tfs;

  //conversion
  
  h_conv_Clean_NumTracks_rhoVal = tfs->make<TH1F>("h_conv_Clean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_conv_Clean_NumTracks_rhoSig = tfs->make<TH1F>("h_conv_Clean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  h_conv_Unclean_NumTracks_rhoVal = tfs->make<TH1F>("h_conv_Unclean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_conv_Unclean_NumTracks_rhoSig = tfs->make<TH1F>("h_conv_Unclean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_conv_Clean_NumTracks_rhoVal->Sumw2(); 
  h_conv_Clean_NumTracks_rhoSig->Sumw2(); 
  h_conv_Unclean_NumTracks_rhoVal->Sumw2(); 
  h_conv_Unclean_NumTracks_rhoSig->Sumw2(); 
  
  h3_conv_Clean_pur_rhoVal = tfs->make<TH3F>("h3_conv_Clean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Unclean_pur_rhoVal = tfs->make<TH3F>("h3_conv_Unclean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Clean_pur_rhoSig = tfs->make<TH3F>("h3_conv_Clean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Unclean_pur_rhoSig = tfs->make<TH3F>("h3_conv_Unclean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_conv_Clean_pur_rhoVal->Sumw2();
  h3_conv_Unclean_pur_rhoVal->Sumw2(); 
  h3_conv_Clean_pur_rhoSig->Sumw2();
  h3_conv_Unclean_pur_rhoSig->Sumw2(); 

  //K decay
  
  h_kdec_Clean_NumTracks_rhoVal = tfs->make<TH1F>("h_kdec_Clean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_kdec_Clean_NumTracks_rhoSig = tfs->make<TH1F>("h_kdec_Clean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  h_kdec_Unclean_NumTracks_rhoVal = tfs->make<TH1F>("h_kdec_Unclean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_kdec_Unclean_NumTracks_rhoSig = tfs->make<TH1F>("h_kdec_Unclean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_kdec_Clean_NumTracks_rhoVal->Sumw2(); 
  h_kdec_Clean_NumTracks_rhoSig->Sumw2(); 
  h_kdec_Unclean_NumTracks_rhoVal->Sumw2(); 
  h_kdec_Unclean_NumTracks_rhoSig->Sumw2(); 
  
  h3_kdec_Clean_pur_rhoVal = tfs->make<TH3F>("h3_kdec_Clean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Unclean_pur_rhoVal = tfs->make<TH3F>("h3_kdec_Unclean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Clean_pur_rhoSig = tfs->make<TH3F>("h3_kdec_Clean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Unclean_pur_rhoSig = tfs->make<TH3F>("h3_kdec_Unclean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_kdec_Clean_pur_rhoVal->Sumw2();
  h3_kdec_Unclean_pur_rhoVal->Sumw2(); 
  h3_kdec_Clean_pur_rhoSig->Sumw2();
  h3_kdec_Unclean_pur_rhoSig->Sumw2(); 
  
  p3_kdec_pur = tfs->make<TProfile3D>("p3_kdec_pur", "", 50,0.,50., 28,-0.49,0.49, 45,0.,45. );
  p3_kdec_eff = tfs->make<TProfile3D>("p3_kdec_eff", "", 50,0.,50., 28,-0.49,0.49, 45,0.,45. );
  p3_kdec_pur->Sumw2();
  p3_kdec_eff->Sumw2();

  //L decay
  
  h_ldec_Clean_NumTracks_rhoVal = tfs->make<TH1F>("h_ldec_Clean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_ldec_Clean_NumTracks_rhoSig = tfs->make<TH1F>("h_ldec_Clean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  h_ldec_Unclean_NumTracks_rhoVal = tfs->make<TH1F>("h_ldec_Unclean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_ldec_Unclean_NumTracks_rhoSig = tfs->make<TH1F>("h_ldec_Unclean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_ldec_Clean_NumTracks_rhoVal->Sumw2(); 
  h_ldec_Clean_NumTracks_rhoSig->Sumw2(); 
  h_ldec_Unclean_NumTracks_rhoVal->Sumw2(); 
  h_ldec_Unclean_NumTracks_rhoSig->Sumw2(); 
  
  h3_ldec_Clean_pur_rhoVal = tfs->make<TH3F>("h3_ldec_Clean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Unclean_pur_rhoVal = tfs->make<TH3F>("h3_ldec_Unclean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Clean_pur_rhoSig = tfs->make<TH3F>("h3_ldec_Clean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Unclean_pur_rhoSig = tfs->make<TH3F>("h3_ldec_Unclean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_ldec_Clean_pur_rhoVal->Sumw2();
  h3_ldec_Unclean_pur_rhoVal->Sumw2(); 
  h3_ldec_Clean_pur_rhoSig->Sumw2();
  h3_ldec_Unclean_pur_rhoSig->Sumw2(); 
  
  p3_ldec_pur = tfs->make<TProfile3D>("p3_ldec_pur", "", 50,0.,50., 40,-0.5,0.5, 45,0.,45. );
  p3_ldec_eff = tfs->make<TProfile3D>("p3_ldec_eff", "", 50,0.,50., 40,-0.5,0.5, 45,0.,45. );
  p3_ldec_pur->Sumw2();
  p3_ldec_eff->Sumw2();

  //Nuclear interaction
  
  h_nuci_Clean_NumTracks_rhoVal = tfs->make<TH1F>("h_nuci_Clean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_nuci_Clean_NumTracks_rhoSig = tfs->make<TH1F>("h_nuci_Clean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  h_nuci_Unclean_NumTracks_rhoVal = tfs->make<TH1F>("h_nuci_Unclean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_nuci_Unclean_NumTracks_rhoSig = tfs->make<TH1F>("h_nuci_Unclean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_nuci_Clean_NumTracks_rhoVal->Sumw2(); 
  h_nuci_Clean_NumTracks_rhoSig->Sumw2(); 
  h_nuci_Unclean_NumTracks_rhoVal->Sumw2(); 
  h_nuci_Unclean_NumTracks_rhoSig->Sumw2(); 
  
  h3_nuci_Clean_pur_rhoVal = tfs->make<TH3F>("h3_nuci_Clean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Unclean_pur_rhoVal = tfs->make<TH3F>("h3_nuci_Unclean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Clean_pur_rhoSig = tfs->make<TH3F>("h3_nuci_Clean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Unclean_pur_rhoSig = tfs->make<TH3F>("h3_nuci_Unclean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_nuci_Clean_pur_rhoVal->Sumw2();
  h3_nuci_Unclean_pur_rhoVal->Sumw2(); 
  h3_nuci_Clean_pur_rhoSig->Sumw2();
  h3_nuci_Unclean_pur_rhoSig->Sumw2();
  
  p3_nuci_pur = tfs->make<TProfile3D>("p3_nuci_pur", "", 50,0.,50., 12,0.,12., 45,0.,45. );
  p3_nuci_eff = tfs->make<TProfile3D>("p3_nuci_eff", "", 50,0.,50., 12,0.,12., 45,0.,45. );
  p3_nuci_pur->Sumw2();
  p3_nuci_eff->Sumw2();

  //Inclusive vertex finder
  
  h_invf_Clean_NumTracks_rhoVal = tfs->make<TH1F>("h_invf_Clean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_invf_Clean_NumTracks_rhoSig = tfs->make<TH1F>("h_invf_Clean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  h_invf_Unclean_NumTracks_rhoVal = tfs->make<TH1F>("h_invf_Unclean_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);
  h_invf_Unclean_NumTracks_rhoSig = tfs->make<TH1F>("h_invf_Unclean_NumTracks_rhoSig", "", rhoNbin, rhoIntervalsa);
  
  h_invf_Clean_NumTracks_rhoVal->Sumw2(); 
  h_invf_Clean_NumTracks_rhoSig->Sumw2(); 
  h_invf_Unclean_NumTracks_rhoVal->Sumw2(); 
  h_invf_Unclean_NumTracks_rhoSig->Sumw2(); 
  
  h3_invf_Clean_pur_rhoVal = tfs->make<TH3F>("h3_invf_Clean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Unclean_pur_rhoVal = tfs->make<TH3F>("h3_invf_Unclean_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Clean_pur_rhoSig = tfs->make<TH3F>("h3_invf_Clean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Unclean_pur_rhoSig = tfs->make<TH3F>("h3_invf_Unclean_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_invf_Clean_pur_rhoVal->Sumw2();
  h3_invf_Unclean_pur_rhoVal->Sumw2(); 
  h3_invf_Clean_pur_rhoSig->Sumw2();
  h3_invf_Unclean_pur_rhoSig->Sumw2(); 
  
  p3_invf_pur = tfs->make<TProfile3D>("p3_invf_pur", "", 50,0.,50., 12,0.,12., 45,0.,45. );
  p3_invf_eff = tfs->make<TProfile3D>("p3_invf_eff", "", 50,0.,50., 12,0.,12., 45,0.,45. );
  p3_invf_pur->Sumw2();
  p3_invf_eff->Sumw2();

}


UnCleanedComp::~UnCleanedComp()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool 
UnCleanedComp::TrackWeightAssociation(const TrackRef&  trackRef, Handle<VertexCollection> vtxcollH) 
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
UnCleanedComp::FindVertexZ(TrackRef trkref, Handle<VertexCollection> vtxCollH, double trackWeight)
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
UnCleanedComp::FindVertex3(TransientTrack transtrk, Handle<VertexCollection> vtxCollH, double trackWeight)
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
UnCleanedComp::GetCleanedConversions(Handle<ConversionCollection> convCollH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<ConversionCollection> cleanedConvColl(new ConversionCollection() );

  for ( unsigned int convcoll_idx=0; convcoll_idx<convCollH->size(); convcoll_idx++ ){

    ConversionRef convref(convCollH,convcoll_idx);

    if ( !cleanedColl ) {   
      cleanedConvColl->push_back(*convref);
      continue;
    }

    if ( convref->quality( Conversion::highPurity ) ){

      cleanedConvColl->push_back(*convref);

    }

  }

  return cleanedConvColl;

}

bool
UnCleanedComp::ComesFromConversion(const TrackRef trackref, ConversionCollection cleanedConvColl, Conversion* gamma)
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
UnCleanedComp::FindConversionVertex3(const reco::TrackRef trackref, reco::Conversion gamma, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 
  
  math::XYZPoint conv_pos = gamma.conversionVertex().position();

  math::XYZVector conv_mom(gamma.refittedPair4Momentum().x(),
   						   gamma.refittedPair4Momentum().y(),
   						   gamma.refittedPair4Momentum().z());

  Track photon(trackref->chi2(), trackref->ndof(), conv_pos, conv_mom, 0, trackref->covariance());

  TransientTrack transpho(photon, &(*bFieldH) );
  transpho.setBeamSpot(*bsH);
  transpho.setES(iSetup);

  return UnCleanedComp::FindVertex3(transpho, vtxcollH, tWeight);  

}

auto_ptr<VertexCompositeCandidateCollection>
UnCleanedComp::GetCleanedKshort(Handle<VertexCompositeCandidateCollection> KshortsH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<VertexCompositeCandidateCollection> cleanedKaonColl(new VertexCompositeCandidateCollection() );

  for ( unsigned int kscoll_idx=0; kscoll_idx<KshortsH->size(); kscoll_idx++ ) {

    VertexCompositeCandidateRef ksref(KshortsH,kscoll_idx);

    if ( !cleanedColl ) {   
      cleanedKaonColl->push_back(*ksref);
      continue;
    }

    GlobalPoint dec_pos = RecoVertex::convertPos(ksref->vertex());    

    GlobalError decayVertexError = GlobalError(ksref->vertexCovariance(0,0), ksref->vertexCovariance(0,1), ksref->vertexCovariance(1,1), ksref->vertexCovariance(0,2), ksref->vertexCovariance(1,2), ksref->vertexCovariance(2,2));  

    double kaon_significance = ( distanceComputerXY.distance( BSVertexState, VertexState( dec_pos, decayVertexError ) ) ).significance();

    if ( ( ksref->vertexNormalizedChi2()<=7. ) &&
         ( fabs(ksref->mass() - kMass)<=0.06 ) &&
         ( kaon_significance>25. ) ) {

      cleanedKaonColl->push_back(*ksref);

    }

  }

  return cleanedKaonColl;

}

auto_ptr<VertexCompositeCandidateCollection>
UnCleanedComp::GetCleanedLambda(Handle<VertexCompositeCandidateCollection> LambdasH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<VertexCompositeCandidateCollection> cleanedLambdaColl(new VertexCompositeCandidateCollection() );

  for ( unsigned int lambdacoll_idx=0; lambdacoll_idx<LambdasH->size(); lambdacoll_idx++ ) {

    VertexCompositeCandidateRef lambdaref(LambdasH,lambdacoll_idx);

    if ( !cleanedColl ) {   
      cleanedLambdaColl->push_back(*lambdaref);
      continue;
    }

    GlobalPoint dec_pos = RecoVertex::convertPos(lambdaref->vertex());    

    GlobalError decayVertexError = GlobalError(lambdaref->vertexCovariance(0,0), lambdaref->vertexCovariance(0,1), lambdaref->vertexCovariance(1,1), lambdaref->vertexCovariance(0,2), lambdaref->vertexCovariance(1,2), lambdaref->vertexCovariance(2,2));  

    double lambda_significance = ( distanceComputerXY.distance( BSVertexState, VertexState( dec_pos, decayVertexError ) ) ).significance();

    if ( ( lambdaref->vertexNormalizedChi2()<=7. ) &&
         ( fabs(lambdaref->mass() - lamMass)<=0.04 ) &&
         ( lambda_significance>26. ) ){

      cleanedLambdaColl->push_back(*lambdaref);

    }

  }

  return cleanedLambdaColl;
	
}

bool
UnCleanedComp::ComesFromV0Decay(const TrackRef trackref, VertexCompositeCandidateCollection cleanedVCCC, VertexCompositeCandidate* V0)
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
UnCleanedComp::FindV0DecayVertex3(const reco::TrackRef trackref, reco::VertexCompositeCandidate V0_vtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

  math::XYZPoint dec_pos = V0_vtx.vertex();

  math::XYZVector dec_mom(V0_vtx.momentum().x(),
						  V0_vtx.momentum().y(),
						  V0_vtx.momentum().z());

  Track V0(trackref->chi2(), trackref->ndof(), dec_pos, dec_mom, 0, trackref->covariance());
 
  TransientTrack transV0(V0, &(*bFieldH) );
  transV0.setBeamSpot(*bsH);
  transV0.setES(iSetup);

  return UnCleanedComp::FindVertex3(transV0, vtxcollH, tWeight);

}

auto_ptr<PFDisplacedVertexCollection>
UnCleanedComp::GetCleanedNI(Handle<PFDisplacedVertexCollection> NuclIntH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<PFDisplacedVertexCollection> cleanedNIColl(new PFDisplacedVertexCollection() );

  for ( PFDisplacedVertexCollection::const_iterator niref=NuclIntH->begin(); niref!=NuclIntH->end(); niref++ ) {

    if ( !cleanedColl ) {
      cleanedNIColl->push_back(*niref);
      continue;
    }

    GlobalPoint ni_pos = RecoVertex::convertPos( niref->position() );    
    GlobalError interactionVertexError = RecoVertex::convertError( niref->error() );

    double nuclint_distance = ( distanceComputerXY.distance( BSVertexState, VertexState( ni_pos, interactionVertexError ) ) ).value();

    if ( ( !niref->isFake() ) &&
         ( niref->isNucl() ) &&
         ( niref->normalizedChi2()<=2. ) &&
         ( niref->tracksSize()>=2 ) &&
         ( nuclint_distance>3. ) ) {

      cleanedNIColl->push_back(*niref);
 
    }

  }

  return cleanedNIColl;
	
}

bool
UnCleanedComp::ComesFromNI(const TrackRef trackref, PFDisplacedVertexCollection cleanedNI, PFDisplacedVertex* displVtx)
{

  //the part for the reassociation of particles from nuclear interactions
  for ( PFDisplacedVertexCollection::const_iterator iDisplV=cleanedNI.begin(); iDisplV!=cleanedNI.end(); iDisplV++ ) {

    if ( iDisplV->trackWeight(trackref)>1.e-2 ) {

      *displVtx = *iDisplV; 
      return true;

    }

  }

  return false;
	
}

VertexRef
UnCleanedComp::FindNIVertex3(const reco::TrackRef trackref, reco::PFDisplacedVertex displVtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight, TransientTrack transhelp)
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

auto_ptr<VertexCollection>
UnCleanedComp::GetCleanedIVF(Handle<VertexCollection> ifvH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<VertexCollection> cleanedIVFColl(new VertexCollection() );
 
  for ( VertexCollection::const_iterator ivfref=ifvH->begin(); ivfref!=ifvH->end(); ivfref++ ) {

    if ( !cleanedColl ) {
      cleanedIVFColl->push_back(*ivfref);
      continue;
    } 

    GlobalPoint iv_pos = RecoVertex::convertPos( ivfref->position() );    
    GlobalError iv_err = RecoVertex::convertError( ivfref->error() );  
    
    double ivf_significance = ( distanceComputerXY.distance( BSVertexState, VertexState( iv_pos, iv_err ))).significance();
    
    if ( ( ivfref->isValid() ) && 
  	     ( !ivfref->isFake() ) && 
  	     ( ivfref->chi2()<=10. ) && 
  	     ( ivfref->nTracks(0.)>=2 )&& 
  	     ( ivf_significance>=5. )  ) {
  	     
      cleanedIVFColl->push_back(*ivfref);
 
    }            
 
  }

  return cleanedIVFColl;
	
}

bool
UnCleanedComp::ComesFromIVF(const TrackRef trackref, VertexCollection cleanedIVF, Vertex* ivfVtx)
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
UnCleanedComp::FindIVFVertex3(const TrackRef trackref, Vertex ivfVtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, Handle<BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{

  math::XYZPoint iv_pos = ivfVtx.position();

  math::XYZVector iv_mom( ivfVtx.p4(0.1, 0.).x(),
						  ivfVtx.p4(0.1, 0.).y(),
						  ivfVtx.p4(0.1, 0.).z() );
  
  Track incom(trackref->chi2(), trackref->ndof(), iv_pos, iv_mom, 0, trackref->covariance());

  TransientTrack transIncom(incom, &(*bFieldH) );
  transIncom.setBeamSpot(*bsH);
  transIncom.setES(iSetup);

  return UnCleanedComp::FindVertex3(transIncom, vtxcollH, tWeight);
   
}

// ------------ method called for each event  ------------
void
UnCleanedComp::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //associate reco tracks to tracking particles
  ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get("quickTrackAssociatorByHits",theAssociator);
  TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product();

  ESHandle<MagneticField> bfH;
  iSetup.get<IdealMagneticFieldRecord>().get(bfH);

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

  GlobalPoint bsPosition = RecoVertex::convertPos(bsH->position());
  GlobalError bsError = RecoVertex::convertError(bsH->covariance3D());
  
  BSVertexState = VertexState(bsPosition,bsError);

  //get the conversion collection for the gamma conversions
  Handle<ConversionCollection> convCollH;
  iEvent.getByLabel(ConversionsCollection_, convCollH);

  auto_ptr<ConversionCollection> cleanedConvCollP = UnCleanedComp::GetCleanedConversions(convCollH,bsH,true);
  auto_ptr<ConversionCollection> uncleanedConvCollP = UnCleanedComp::GetCleanedConversions(convCollH,bsH,false);

  //get the vertex composite candidate collection for the Kshort's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH;
  iEvent.getByLabel(KshortCollection_, vertCompCandCollKshortH);

  auto_ptr<VertexCompositeCandidateCollection> cleanedKshortCollP = UnCleanedComp::GetCleanedKshort(vertCompCandCollKshortH,bsH,true);
  auto_ptr<VertexCompositeCandidateCollection> uncleanedKshortCollP = UnCleanedComp::GetCleanedKshort(vertCompCandCollKshortH,bsH,false);
  
  //get the vertex composite candidate collection for the Lambda's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
  iEvent.getByLabel(LambdaCollection_, vertCompCandCollLambdaH);

  auto_ptr<VertexCompositeCandidateCollection> cleanedLambdaCollP = UnCleanedComp::GetCleanedLambda(vertCompCandCollLambdaH,bsH,true);
  auto_ptr<VertexCompositeCandidateCollection> uncleanedLambdaCollP = UnCleanedComp::GetCleanedLambda(vertCompCandCollLambdaH,bsH,false);

  //get the displaced vertex collection for nuclear interactions
  Handle<PFDisplacedVertexCollection> displVertexCollH;
  iEvent.getByLabel(NIVertexCollection_,displVertexCollH);

  auto_ptr<PFDisplacedVertexCollection> cleanedNICollP = UnCleanedComp::GetCleanedNI(displVertexCollH,bsH,true);
  auto_ptr<PFDisplacedVertexCollection> uncleanedNICollP = UnCleanedComp::GetCleanedNI(displVertexCollH,bsH,false);


  //get the vertex collection for inclusive vertex finder   
  Handle<VertexCollection> ivfVertexCollH;
  iEvent.getByLabel(IVFVertexCollection_, ivfVertexCollH);
  
  auto_ptr<VertexCollection> cleanedIVFCollP = UnCleanedComp::GetCleanedIVF(ivfVertexCollH, bsH, true);
  auto_ptr<VertexCollection> uncleanedIVFCollP = UnCleanedComp::GetCleanedIVF(ivfVertexCollH, bsH, false);

	    
  //loop over all tracks of the track collection	
  for ( size_t idxTrack = 0; idxTrack < theTracksH->size(); ++idxTrack ) {

    TrackRef trackref(theTracksH, idxTrack);
    TrackBaseRef trackbaseref = TrackBaseRef(trackref);

    if(UnCleanedComp::TrackWeightAssociation(trackref, theVerticesH)) continue;

    TransientTrack transtrk(trackref, &(*bfH) );
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
    
    //conversion
    Conversion gamma;
    if ( UnCleanedComp::ComesFromConversion(trackref, *uncleanedConvCollP, &gamma) ){
      
      h_conv_Unclean_NumTracks_rhoVal->Fill( rhoVal );
      h_conv_Unclean_NumTracks_rhoSig->Fill( rhoSig );

      int y_bin_uncleaned = 1;
      if(UnCleanedComp::FindConversionVertex3(trackref, gamma, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_cleaned = 1;
      if ( UnCleanedComp::ComesFromConversion(trackref, *cleanedConvCollP, &gamma) ){
    
        h_conv_Clean_NumTracks_rhoVal->Fill( rhoVal );
        h_conv_Clean_NumTracks_rhoSig->Fill( rhoSig );
        
        if(UnCleanedComp::FindConversionVertex3(trackref, gamma, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_cleaned = 0;
        
      }else{
        if(UnCleanedComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_cleaned = 0;
      }

      h3_conv_Unclean_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_conv_Unclean_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_conv_Clean_pur_rhoVal->Fill( x_bin, y_bin_cleaned, rhoVal );
      h3_conv_Clean_pur_rhoSig->Fill( x_bin, y_bin_cleaned, rhoSig );

    }

    //K decay
    VertexCompositeCandidate Ks;
    if( UnCleanedComp::ComesFromV0Decay(trackref, *uncleanedKshortCollP, &Ks) ){
      
      h_kdec_Unclean_NumTracks_rhoVal->Fill( rhoVal );
      h_kdec_Unclean_NumTracks_rhoSig->Fill( rhoSig );   
    
      double chi2 = Ks.vertexNormalizedChi2();
      double massD = fabs(Ks.mass() - kMass);

      GlobalPoint ks_pos = RecoVertex::convertPos( Ks.vertex() );    
      GlobalError ks_err = RecoVertex::convertError( Ks.vertexCovariance() ); 
    
      double s = (distanceComputerXY.distance(BSVertexState, VertexState(ks_pos, ks_err))).significance();

      int y_bin_uncleaned = 1;
      if(UnCleanedComp::FindV0DecayVertex3(trackref, Ks, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_cleaned = 1;
      if ( UnCleanedComp::ComesFromV0Decay(trackref, *cleanedKshortCollP, &Ks) ){
    
        h_kdec_Clean_NumTracks_rhoVal->Fill( rhoVal );
        h_kdec_Clean_NumTracks_rhoSig->Fill( rhoSig );
                 
        if(UnCleanedComp::FindV0DecayVertex3(trackref, Ks, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_cleaned = 0;
        
      }else{
        if(UnCleanedComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_cleaned = 0;
      }

      h3_kdec_Unclean_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_kdec_Unclean_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_kdec_Clean_pur_rhoVal->Fill( x_bin, y_bin_cleaned, rhoVal );
      h3_kdec_Clean_pur_rhoSig->Fill( x_bin, y_bin_cleaned, rhoSig );
      
      if ( y_bin_uncleaned==0 ) {
        if ( x_bin==0 ) {
          p3_kdec_pur->Fill( chi2, massD, s, 1. );
        } else {
          p3_kdec_pur->Fill( chi2, massD, s, 0. );
        }
      }
      
      if ( x_bin==0 ) {
        if ( y_bin_uncleaned==0 ) {
          p3_kdec_eff->Fill( chi2, massD, s, 1. );
        } else {
          p3_kdec_eff->Fill( chi2, massD, s, 0. );
        }
      } 

    }

    //Lambda
    VertexCompositeCandidate La;
    if( UnCleanedComp::ComesFromV0Decay(trackref, *uncleanedLambdaCollP, &La) ){
      
      h_ldec_Unclean_NumTracks_rhoVal->Fill( rhoVal );
      h_ldec_Unclean_NumTracks_rhoSig->Fill( rhoSig );
    
      double chi2 = La.vertexNormalizedChi2();
      double massD = fabs(La.mass() - lamMass);

      GlobalPoint la_pos = RecoVertex::convertPos( La.vertex() );    
      GlobalError la_err = RecoVertex::convertError( La.vertexCovariance() );  
    
      double s = (distanceComputerXY.distance(BSVertexState, VertexState(la_pos, la_err))).significance();

      int y_bin_uncleaned = 1;
      if(UnCleanedComp::FindV0DecayVertex3(trackref, La, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_cleaned = 1;
      if ( UnCleanedComp::ComesFromV0Decay(trackref, *cleanedLambdaCollP, &La) ){
    
        h_ldec_Clean_NumTracks_rhoVal->Fill( rhoVal );
        h_ldec_Clean_NumTracks_rhoSig->Fill( rhoSig );
        
        if(UnCleanedComp::FindV0DecayVertex3(trackref, La, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_cleaned = 0;
        
      }else{
        if(UnCleanedComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_cleaned = 0;
      }

      h3_ldec_Unclean_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_ldec_Unclean_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_ldec_Clean_pur_rhoVal->Fill( x_bin, y_bin_cleaned, rhoVal );
      h3_ldec_Clean_pur_rhoSig->Fill( x_bin, y_bin_cleaned, rhoSig );
      
      if ( y_bin_uncleaned==0 ) {
        if ( x_bin==0 ) {
          p3_ldec_pur->Fill( chi2, massD, s, 1. );
        } else {
          p3_ldec_pur->Fill( chi2, massD, s, 0. );
        }
      }
      
      if ( x_bin==0 ) {
        if ( y_bin_uncleaned==0 ) {
          p3_ldec_eff->Fill( chi2, massD, s, 1. );
        } else {
          p3_ldec_eff->Fill( chi2, massD, s, 0. );
        }
      }

    }
    
    //Nuclear interaction
    PFDisplacedVertex displVtx;
    if ( UnCleanedComp::ComesFromNI(trackref, *uncleanedNICollP, &displVtx) ){
      
      h_nuci_Unclean_NumTracks_rhoVal->Fill( rhoVal );
      h_nuci_Unclean_NumTracks_rhoSig->Fill( rhoSig );
    
      double chi2 = displVtx.normalizedChi2();
      int ntrks = displVtx.tracksSize();

      GlobalPoint ni_pos = RecoVertex::convertPos( displVtx.position() );    
      GlobalError ni_err = RecoVertex::convertError( displVtx.error() );
    
      double s = (distanceComputerXY.distance(BSVertexState, VertexState(ni_pos, ni_err))).value();

      int y_bin_uncleaned = 1;
      if(UnCleanedComp::FindNIVertex3(trackref, displVtx, bfH, iSetup, bsH, theVerticesH, 0.01, transtrk)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_cleaned = 1;
      if ( UnCleanedComp::ComesFromNI(trackref, *cleanedNICollP, &displVtx) ){
      
        h_nuci_Clean_NumTracks_rhoVal->Fill( rhoVal );
        h_nuci_Clean_NumTracks_rhoSig->Fill( rhoSig );
        
        if(UnCleanedComp::FindNIVertex3(trackref, displVtx, bfH, iSetup, bsH, theVerticesH, 0.01, transtrk)==firstVertex) y_bin_cleaned = 0;
        
      }else{
        if(UnCleanedComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_cleaned = 0;
      }

      h3_nuci_Unclean_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_nuci_Unclean_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_nuci_Clean_pur_rhoVal->Fill( x_bin, y_bin_cleaned, rhoVal );
      h3_nuci_Clean_pur_rhoSig->Fill( x_bin, y_bin_cleaned, rhoSig );
      
      if ( y_bin_uncleaned==0 ) {
        if ( x_bin==0 ) {
          p3_nuci_pur->Fill( chi2, ntrks, s, 1. );
        } else {
          p3_nuci_pur->Fill( chi2, ntrks, s, 0. );
        }
      }
      
      if ( x_bin==0 ) {
        if ( y_bin_uncleaned==0 ) {
          p3_nuci_eff->Fill( chi2, ntrks, s, 1. );
        } else {
          p3_nuci_eff->Fill( chi2, ntrks, s, 0. );
        }
      }

    }
    
    //Inclusive vertex finder
	Vertex ivfVtx;
    if ( UnCleanedComp::ComesFromIVF(trackref, *uncleanedIVFCollP, &ivfVtx) ) {
      
      h_invf_Unclean_NumTracks_rhoVal->Fill( rhoVal );
      h_invf_Unclean_NumTracks_rhoSig->Fill( rhoSig );
    
      double chi2 = ivfVtx.normalizedChi2();
      int ntrks = ivfVtx.nTracks(0.);

      GlobalPoint iv_pos = RecoVertex::convertPos( ivfVtx.position() );    
      GlobalError iv_err = RecoVertex::convertError( ivfVtx.error() );  
    
      double s = (distanceComputerXY.distance(BSVertexState, VertexState(iv_pos, iv_err))).significance();

      int y_bin_uncleaned = 1;
      if(UnCleanedComp::FindIVFVertex3(trackref, ivfVtx, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_cleaned = 1;
      if ( UnCleanedComp::ComesFromIVF(trackref, *cleanedIVFCollP, &ivfVtx) ){
      
        h_invf_Clean_NumTracks_rhoVal->Fill( rhoVal );
        h_invf_Clean_NumTracks_rhoSig->Fill( rhoSig );

        if(UnCleanedComp::FindIVFVertex3(trackref, ivfVtx, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_cleaned = 0;
        
      }else{
        if(UnCleanedComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_cleaned = 0;
      }

      h3_invf_Unclean_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_invf_Unclean_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_invf_Clean_pur_rhoVal->Fill( x_bin, y_bin_cleaned, rhoVal );
      h3_invf_Clean_pur_rhoSig->Fill( x_bin, y_bin_cleaned, rhoSig );
      
      if ( y_bin_uncleaned==0 ) {
        if ( x_bin==0 ) {
          p3_invf_pur->Fill( chi2, ntrks, s, 1. );
        } else {
          p3_invf_pur->Fill( chi2, ntrks, s, 0. );
        }
      }
      
      if ( x_bin==0 ) {
        if ( y_bin_uncleaned==0 ) {
          p3_invf_eff->Fill( chi2, ntrks, s, 1. );
        } else {
          p3_invf_eff->Fill( chi2, ntrks, s, 0. );
        }
      }

    }

  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
UnCleanedComp::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(UnCleanedComp);
