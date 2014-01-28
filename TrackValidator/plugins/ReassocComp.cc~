// -*- C++ -*-
//
// Package:    ReassocComp
// Class:      ReassocComp
// 
/**\class ReassocComp ReassocComp.cc MGeisler/ReassocComp/src/ReassocComp.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Tue Sep 25 12:12:54 CEST 2012
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

//
// class declaration
//

class ReassocComp : public edm::EDAnalyzer {
   public:
      explicit ReassocComp(const edm::ParameterSet&);
      ~ReassocComp();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static bool TrackWeightAssociation(reco::TrackBaseRef&, edm::Handle<reco::VertexCollection>);

      static reco::VertexRef FindVertexZ(reco::TrackRef, edm::Handle<reco::VertexCollection>, double);
      static reco::VertexRef FindVertex3D(reco::TransientTrack, edm::Handle<reco::VertexCollection>, double);

      static double dR(math::XYZPoint, math::XYZVector, edm::Handle<reco::BeamSpot>);

      virtual std::auto_ptr<reco::ConversionCollection> GetCleanedConversions(edm::Handle<reco::ConversionCollection>, 
                                                                             edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromConversion(const reco::TrackRef, reco::ConversionCollection, reco::Conversion*);        
      static reco::VertexRef FindConversionVertex3D(const reco::TrackRef, reco::Conversion, edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                                    edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, double);

      virtual std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedKshort(edm::Handle<reco::VertexCompositeCandidateCollection>,
                                                                                       edm::Handle<reco::BeamSpot>, bool);
      virtual std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedLambda(edm::Handle<reco::VertexCompositeCandidateCollection>, 
                                                                                       edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromV0Decay(const reco::TrackRef, reco::VertexCompositeCandidateCollection,                
                                   reco::VertexCompositeCandidate*);         
      static reco::VertexRef FindV0DecayVertex3D(const reco::TrackRef, reco::VertexCompositeCandidate, edm::ESHandle<MagneticField>, 
                                                 const edm::EventSetup&, edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, 
	                                             double);

      virtual std::auto_ptr<reco::PFDisplacedVertexCollection> GetCleanedNI(edm::Handle<reco::PFDisplacedVertexCollection>, 
                                                                            edm::Handle<reco::BeamSpot>, bool); 
      static bool ComesFromNI(const reco::TrackRef, reco::PFDisplacedVertexCollection, reco::PFDisplacedVertex*);      
      static reco::VertexRef FindNIVertex3D(const reco::TrackRef, reco::PFDisplacedVertex, edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                            edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, double, reco::TransientTrack);
                                                   
      virtual std::auto_ptr<reco::VertexCollection> GetCleanedIVF(edm::Handle<reco::VertexCollection>, edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromIVF(const reco::TrackRef, reco::VertexCollection, reco::Vertex*);     
      static reco::VertexRef FindIVFVertex3D(const reco::TrackRef, reco::Vertex,edm::ESHandle<MagneticField>, const edm::EventSetup&, 
                                            edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, double);


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
      
      TH1F *h_conv_numTracks_nMissHits; 
      TH3F *h3_conv_Sec_pur_rhoVal, *h3_conv_Sec_pur_rhoSig, *h3_conv_Sec_pur_nMissHits;
      TH3F *h3_conv_Fin_pur_rhoVal, *h3_conv_Fin_pur_rhoSig, *h3_conv_Fin_pur_nMissHits;

      TH1F *h_kdec_numTracks_nMissHits;
      TH3F *h3_kdec_Sec_pur_rhoVal, *h3_kdec_Sec_pur_rhoSig, *h3_kdec_Sec_pur_nMissHits;
      TH3F *h3_kdec_Fin_pur_rhoVal, *h3_kdec_Fin_pur_rhoSig, *h3_kdec_Fin_pur_nMissHits;

      TH1F *h_ldec_numTracks_nMissHits;
      TH3F *h3_ldec_Sec_pur_rhoVal, *h3_ldec_Sec_pur_rhoSig, *h3_ldec_Sec_pur_nMissHits;
      TH3F *h3_ldec_Fin_pur_rhoVal, *h3_ldec_Fin_pur_rhoSig, *h3_ldec_Fin_pur_nMissHits;

      TH1F *h_nuci_numTracks_nMissHits;
      TH3F *h3_nuci_Sec_pur_rhoVal, *h3_nuci_Sec_pur_rhoSig, *h3_nuci_Sec_pur_nMissHits;
      TH3F *h3_nuci_Fin_pur_rhoVal, *h3_nuci_Fin_pur_rhoSig, *h3_nuci_Fin_pur_nMissHits;

      TH1F *h_invf_numTracks_nMissHits;
      TH3F *h3_invf_Sec_pur_rhoVal, *h3_invf_Sec_pur_rhoSig, *h3_invf_Sec_pur_nMissHits;
      TH3F *h3_invf_Fin_pur_rhoVal, *h3_invf_Fin_pur_rhoSig, *h3_invf_Fin_pur_nMissHits;
      

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
ReassocComp::ReassocComp(const edm::ParameterSet& iConfig)
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

  double missMin = -0.5;
  double missMax = 4.5;
  int missNbin = 5;  
   
  double stepNmh = (missMax-missMin)/missNbin; 
  
  Double_t nmhIntervalsa[6];  
  
  for (int k=0;k<missNbin+1;k++) {
    double d = missMin+k*stepNmh;
    nmhIntervalsa[k] = d;
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
  
  h_conv_numTracks_nMissHits = tfs->make<TH1F>("h_conv_numTracks_nMissHits", "", missNbin, nmhIntervalsa);
  h_conv_numTracks_nMissHits->Sumw2();
  
  h3_conv_Sec_pur_rhoVal = tfs->make<TH3F>("h3_conv_Sec_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Sec_pur_rhoSig = tfs->make<TH3F>("h3_conv_Sec_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Sec_pur_nMissHits = tfs->make<TH3F>("h3_conv_Sec_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  h3_conv_Fin_pur_rhoVal = tfs->make<TH3F>("h3_conv_Fin_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Fin_pur_rhoSig = tfs->make<TH3F>("h3_conv_Fin_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_conv_Fin_pur_nMissHits = tfs->make<TH3F>("h3_conv_Fin_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  
  h3_conv_Sec_pur_rhoVal->Sumw2();
  h3_conv_Sec_pur_rhoSig->Sumw2(); 
  h3_conv_Sec_pur_nMissHits->Sumw2(); 
  h3_conv_Fin_pur_rhoVal->Sumw2();
  h3_conv_Fin_pur_rhoSig->Sumw2(); 
  h3_conv_Fin_pur_nMissHits->Sumw2(); 

  //K decay
  
  h_kdec_numTracks_nMissHits = tfs->make<TH1F>("h_kdec_numTracks_nMissHits", "", missNbin, nmhIntervalsa);
  h_kdec_numTracks_nMissHits->Sumw2();
  
  h3_kdec_Sec_pur_rhoVal = tfs->make<TH3F>("h3_kdec_Sec_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Sec_pur_rhoSig = tfs->make<TH3F>("h3_kdec_Sec_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Sec_pur_nMissHits = tfs->make<TH3F>("h3_kdec_Sec_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  h3_kdec_Fin_pur_rhoVal = tfs->make<TH3F>("h3_kdec_Fin_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Fin_pur_rhoSig = tfs->make<TH3F>("h3_kdec_Fin_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_kdec_Fin_pur_nMissHits = tfs->make<TH3F>("h3_kdec_Fin_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  
  h3_kdec_Sec_pur_rhoVal->Sumw2();
  h3_kdec_Sec_pur_rhoSig->Sumw2(); 
  h3_kdec_Sec_pur_nMissHits->Sumw2(); 
  h3_kdec_Fin_pur_rhoVal->Sumw2();
  h3_kdec_Fin_pur_rhoSig->Sumw2();
  h3_kdec_Fin_pur_nMissHits->Sumw2(); 

  //L decay
  
  h_ldec_numTracks_nMissHits = tfs->make<TH1F>("h_ldec_numTracks_nMissHits", "", missNbin, nmhIntervalsa);
  h_ldec_numTracks_nMissHits->Sumw2();
  
  h3_ldec_Sec_pur_rhoVal = tfs->make<TH3F>("h3_ldec_Sec_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Sec_pur_rhoSig = tfs->make<TH3F>("h3_ldec_Sec_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Sec_pur_nMissHits = tfs->make<TH3F>("h3_ldec_Sec_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  h3_ldec_Fin_pur_rhoVal = tfs->make<TH3F>("h3_ldec_Fin_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Fin_pur_rhoSig = tfs->make<TH3F>("h3_ldec_Fin_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_ldec_Fin_pur_nMissHits = tfs->make<TH3F>("h3_ldec_Fin_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  
  h3_ldec_Sec_pur_rhoVal->Sumw2();
  h3_ldec_Sec_pur_rhoSig->Sumw2(); 
  h3_ldec_Sec_pur_nMissHits->Sumw2(); 
  h3_ldec_Fin_pur_rhoVal->Sumw2();
  h3_ldec_Fin_pur_rhoSig->Sumw2(); 
  h3_ldec_Fin_pur_nMissHits->Sumw2(); 

  //Nuclear interaction
  
  h_nuci_numTracks_nMissHits = tfs->make<TH1F>("h_nuci_numTracks_nMissHits", "", missNbin, nmhIntervalsa);
  h_nuci_numTracks_nMissHits->Sumw2();
  
  h3_nuci_Sec_pur_rhoVal = tfs->make<TH3F>("h3_nuci_Sec_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Sec_pur_rhoSig = tfs->make<TH3F>("h3_nuci_Sec_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Sec_pur_nMissHits = tfs->make<TH3F>("h3_nuci_Sec_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  h3_nuci_Fin_pur_rhoVal = tfs->make<TH3F>("h3_nuci_Fin_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Fin_pur_rhoSig = tfs->make<TH3F>("h3_nuci_Fin_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_nuci_Fin_pur_nMissHits = tfs->make<TH3F>("h3_nuci_Fin_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  
  h3_nuci_Sec_pur_rhoVal->Sumw2();
  h3_nuci_Sec_pur_rhoSig->Sumw2();
  h3_nuci_Sec_pur_nMissHits->Sumw2();  
  h3_nuci_Fin_pur_rhoVal->Sumw2();
  h3_nuci_Fin_pur_rhoSig->Sumw2(); 
  h3_nuci_Fin_pur_nMissHits->Sumw2(); 

  //Inclusive vertex
  
  h_invf_numTracks_nMissHits = tfs->make<TH1F>("h_invf_numTracks_nMissHits", "", missNbin, nmhIntervalsa);
  h_invf_numTracks_nMissHits->Sumw2();
  
  h3_invf_Sec_pur_rhoVal = tfs->make<TH3F>("h3_invf_Sec_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Sec_pur_rhoSig = tfs->make<TH3F>("h3_invf_Sec_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Sec_pur_nMissHits = tfs->make<TH3F>("h3_invf_Sec_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  h3_invf_Fin_pur_rhoVal = tfs->make<TH3F>("h3_invf_Fin_pur_rhoVal", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Fin_pur_rhoSig = tfs->make<TH3F>("h3_invf_Fin_pur_rhoSig", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_invf_Fin_pur_nMissHits = tfs->make<TH3F>("h3_invf_Fin_pur_nMissHits", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, missNbin, nmhIntervalsa);
  
  h3_invf_Sec_pur_rhoVal->Sumw2();
  h3_invf_Sec_pur_rhoSig->Sumw2(); 
  h3_invf_Sec_pur_nMissHits->Sumw2(); 
  h3_invf_Fin_pur_rhoVal->Sumw2();
  h3_invf_Fin_pur_rhoSig->Sumw2(); 
  h3_invf_Fin_pur_nMissHits->Sumw2(); 

}


ReassocComp::~ReassocComp()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool 
ReassocComp::TrackWeightAssociation(TrackBaseRef&  trackbaseRef, Handle<VertexCollection> vtxcollH) 
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

  return (bestweight > 1e-5);

}

VertexRef
ReassocComp::FindVertexZ(TrackRef trkref, Handle<VertexCollection> vtxCollH, double trackWeight)
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
ReassocComp::FindVertex3D(TransientTrack transtrk, Handle<VertexCollection> vtxCollH, double trackWeight)
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
ReassocComp::GetCleanedConversions(Handle<ConversionCollection> convCollH, Handle<BeamSpot> bsH, bool cleanedColl)
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
ReassocComp::ComesFromConversion(const TrackRef trackref, ConversionCollection cleanedConvColl, Conversion* gamma)
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
ReassocComp::FindConversionVertex3D(const reco::TrackRef trackref, reco::Conversion gamma, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 
  
  math::XYZPoint conv_pos = gamma.conversionVertex().position();

  math::XYZVector conv_mom(gamma.refittedPair4Momentum().x(),
   						   gamma.refittedPair4Momentum().y(),
   						   gamma.refittedPair4Momentum().z());

  Track photon(trackref->chi2(), trackref->ndof(), conv_pos, conv_mom, 0, trackref->covariance());

  TransientTrack transpho(photon, &(*bFieldH) );
  transpho.setBeamSpot(*bsH);
  transpho.setES(iSetup);

  return ReassocComp::FindVertex3D(transpho, vtxcollH, tWeight);

}

auto_ptr<VertexCompositeCandidateCollection>
ReassocComp::GetCleanedKshort(Handle<VertexCompositeCandidateCollection> KshortsH, Handle<BeamSpot> bsH, bool cleanedColl)
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
ReassocComp::GetCleanedLambda(Handle<VertexCompositeCandidateCollection> LambdasH, Handle<BeamSpot> bsH, bool cleanedColl)
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
ReassocComp::ComesFromV0Decay(const TrackRef trackref, VertexCompositeCandidateCollection cleanedVCCC, VertexCompositeCandidate* V0)
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
ReassocComp::FindV0DecayVertex3D(const reco::TrackRef trackref, reco::VertexCompositeCandidate V0_vtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

  math::XYZPoint dec_pos = V0_vtx.vertex();

  math::XYZVector dec_mom(V0_vtx.momentum().x(),
						  V0_vtx.momentum().y(),
						  V0_vtx.momentum().z());

  Track V0(trackref->chi2(), trackref->ndof(), dec_pos, dec_mom, 0, trackref->covariance());
 
  TransientTrack transV0(V0, &(*bFieldH) );
  transV0.setBeamSpot(*bsH);
  transV0.setES(iSetup);

  return ReassocComp::FindVertex3D(transV0, vtxcollH, tWeight);

}

auto_ptr<PFDisplacedVertexCollection>
ReassocComp::GetCleanedNI(Handle<PFDisplacedVertexCollection> NuclIntH, Handle<BeamSpot> bsH, bool cleanedColl)
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
         ( nuclint_distance>10. ) ) {

      cleanedNIColl->push_back(*niref);
 
    }

  }

  return cleanedNIColl;
  
}

bool
ReassocComp::ComesFromNI(const TrackRef trackref, PFDisplacedVertexCollection cleanedNI, PFDisplacedVertex* displVtx)
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
ReassocComp::FindNIVertex3D(const reco::TrackRef trackref, reco::PFDisplacedVertex displVtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight, TransientTrack transhelp)
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

        return FindVertex3D(transIncom, vtxcollH, tWeight); 

      }

    }

  }

  return ReassocComp::FindVertex3D(transhelp, vtxcollH, tWeight);

}

auto_ptr<VertexCollection>
ReassocComp::GetCleanedIVF(Handle<VertexCollection> ifvH, Handle<BeamSpot> bsH, bool cleanedColl)
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
  	     ( ivf_significance>=5. ) && 
  	     ( ivfref->nTracks(0.)>=2 ) ) {
  	     
      cleanedIVFColl->push_back(*ivfref);
 
    }            
 
  }

  return cleanedIVFColl;
  
}

bool
ReassocComp::ComesFromIVF(const TrackRef trackref, VertexCollection cleanedIVF, Vertex* ivfVtx)
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
ReassocComp::FindIVFVertex3D(const TrackRef trackref, Vertex ivfVtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, Handle<BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{

  math::XYZPoint iv_pos = ivfVtx.position();

  math::XYZVector iv_mom( ivfVtx.p4(0.1, 0.).x(),
						  ivfVtx.p4(0.1, 0.).y(),
						  ivfVtx.p4(0.1, 0.).z() );
  
  Track incom(trackref->chi2(), trackref->ndof(), iv_pos, iv_mom, 0, trackref->covariance());

  TransientTrack transIncom(incom, &(*bFieldH) );
  transIncom.setBeamSpot(*bsH);
  transIncom.setES(iSetup);

  return ReassocComp::FindVertex3D(transIncom, vtxcollH, tWeight); 
	
}

// ------------ method called for each event  ------------
void
ReassocComp::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  auto_ptr<ConversionCollection> uncleanedConvCollP = ReassocComp::GetCleanedConversions(convCollH,bsH,false);

  //get the vertex composite candidate collection for the Kshort's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH;
  iEvent.getByLabel(KshortCollection_, vertCompCandCollKshortH);
  auto_ptr<VertexCompositeCandidateCollection> uncleanedKshortCollP = ReassocComp::GetCleanedKshort(vertCompCandCollKshortH,bsH,false);
  
  //get the vertex composite candidate collection for the Lambda's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
  iEvent.getByLabel(LambdaCollection_, vertCompCandCollLambdaH);
  auto_ptr<VertexCompositeCandidateCollection> uncleanedLambdaCollP = ReassocComp::GetCleanedLambda(vertCompCandCollLambdaH,bsH,false);

  //get the displaced vertex collection for nuclear interactions
  Handle<PFDisplacedVertexCollection> displVertexCollH;
  iEvent.getByLabel(NIVertexCollection_,displVertexCollH);
  auto_ptr<PFDisplacedVertexCollection> uncleanedNICollP = ReassocComp::GetCleanedNI(displVertexCollH,bsH,false);

  //get the vertex collection for inclusive vertex finder   
  Handle<VertexCollection> ivfVertexCollH;
  iEvent.getByLabel(IVFVertexCollection_, ivfVertexCollH);  
  auto_ptr<VertexCollection> uncleanedIVFCollP = ReassocComp::GetCleanedIVF(ivfVertexCollH, bsH, false);

	    
  //loop over all tracks of the track collection	
  for ( size_t idxTrack = 0; idxTrack < theTracksH->size(); ++idxTrack ) {

    TrackRef trackref(theTracksH, idxTrack);
    TrackBaseRef trackbaseref = TrackBaseRef(trackref);

    if(ReassocComp::TrackWeightAssociation(trackbaseref, theVerticesH)) continue;

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
    
    int missHits = trackref->trackerExpectedHitsInner().numberOfLostHits();

    //Conversion
    Conversion gamma;
    if ( ReassocComp::ComesFromConversion(trackref, *uncleanedConvCollP, &gamma) ){

      h_conv_numTracks_nMissHits->Fill( missHits );

      int y_bin_uncleaned = 1;
      if(ReassocComp::FindConversionVertex3D(trackref, gamma, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_final = 1;
      if(ReassocComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_final = 0;

      h3_conv_Sec_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_conv_Sec_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_conv_Sec_pur_nMissHits->Fill( x_bin, y_bin_uncleaned, missHits );
      h3_conv_Fin_pur_rhoVal->Fill( x_bin, y_bin_final, rhoVal );
      h3_conv_Fin_pur_rhoSig->Fill( x_bin, y_bin_final, rhoSig );
      h3_conv_Fin_pur_nMissHits->Fill( x_bin, y_bin_final, missHits );

    }

    //Kaon Decay
    VertexCompositeCandidate V0;
    if ( ReassocComp::ComesFromV0Decay(trackref, *uncleanedKshortCollP, &V0) ){

      h_kdec_numTracks_nMissHits->Fill( missHits );

      int y_bin_uncleaned = 1;
      if(ReassocComp::FindV0DecayVertex3D(trackref, V0, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_final = 1;
      if(ReassocComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_final = 0;

      h3_kdec_Sec_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_kdec_Sec_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_kdec_Sec_pur_nMissHits->Fill( x_bin, y_bin_uncleaned, missHits );
      h3_kdec_Fin_pur_rhoVal->Fill( x_bin, y_bin_final, rhoVal );
      h3_kdec_Fin_pur_rhoSig->Fill( x_bin, y_bin_final, rhoSig );
      h3_kdec_Fin_pur_nMissHits->Fill( x_bin, y_bin_final, missHits );

    }

    //Lambda Decay
    if ( ReassocComp::ComesFromV0Decay(trackref, *uncleanedLambdaCollP, &V0) ){

      h_ldec_numTracks_nMissHits->Fill( missHits );

      int y_bin_uncleaned = 1;
      if(ReassocComp::FindV0DecayVertex3D(trackref, V0, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_final = 1;
      if(ReassocComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_final = 0;

      h3_ldec_Sec_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_ldec_Sec_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_ldec_Sec_pur_nMissHits->Fill( x_bin, y_bin_uncleaned, missHits );
      h3_ldec_Fin_pur_rhoVal->Fill( x_bin, y_bin_final, rhoVal );
      h3_ldec_Fin_pur_rhoSig->Fill( x_bin, y_bin_final, rhoSig );
      h3_ldec_Fin_pur_nMissHits->Fill( x_bin, y_bin_final, missHits );

    }

    //Nuclear Interaction
    PFDisplacedVertex displVtx;
    if ( ReassocComp::ComesFromNI(trackref, *uncleanedNICollP, &displVtx) ) {

      h_nuci_numTracks_nMissHits->Fill( missHits );

      int y_bin_uncleaned = 1;
      if(ReassocComp::FindNIVertex3D(trackref, displVtx, bfH, iSetup, bsH, theVerticesH, 0.01, transtrk)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_final = 1;
      if(ReassocComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_final = 0;

      h3_nuci_Sec_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_nuci_Sec_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_nuci_Sec_pur_nMissHits->Fill( x_bin, y_bin_uncleaned, missHits );
      h3_nuci_Fin_pur_rhoVal->Fill( x_bin, y_bin_final, rhoVal );
      h3_nuci_Fin_pur_rhoSig->Fill( x_bin, y_bin_final, rhoSig );
      h3_nuci_Fin_pur_nMissHits->Fill( x_bin, y_bin_final, missHits );

    }

    //Inclusive vertex
	Vertex ivfVtx;
    if ( ReassocComp::ComesFromIVF(trackref, *uncleanedIVFCollP, &ivfVtx) ) {

      h_invf_numTracks_nMissHits->Fill( missHits );

      int y_bin_uncleaned = 1;
      if(ReassocComp::FindIVFVertex3D(trackref, ivfVtx, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex) y_bin_uncleaned = 0;

      int y_bin_final = 1;
      if(ReassocComp::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex) y_bin_final = 0;

      h3_invf_Sec_pur_rhoVal->Fill( x_bin, y_bin_uncleaned, rhoVal );
      h3_invf_Sec_pur_rhoSig->Fill( x_bin, y_bin_uncleaned, rhoSig );
      h3_invf_Sec_pur_nMissHits->Fill( x_bin, y_bin_uncleaned, missHits );
      h3_invf_Fin_pur_rhoVal->Fill( x_bin, y_bin_final, rhoVal );
      h3_invf_Fin_pur_rhoSig->Fill( x_bin, y_bin_final, rhoSig );
      h3_invf_Fin_pur_nMissHits->Fill( x_bin, y_bin_final, missHits );

    }

  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ReassocComp::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ReassocComp);
