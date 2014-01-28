// -*- C++ -*-
//
// Package:    FinalAssocTest
// Class:      FinalAssocTest
// 
/**\class FinalAssocTest FinalAssocTest.cc MGeisler/FinalAssocTest/src/FinalAssocTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Wed Sep  5 09:15:18 CEST 2012
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

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

// root includes
#include "TH1F.h"
#include "TH3F.h"

//
// class declaration
//

class FinalAssocTest : public edm::EDAnalyzer {
   public:
      explicit FinalAssocTest(const edm::ParameterSet&);
      ~FinalAssocTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static bool TrackWeightAssociation(const reco::TrackRef&, edm::Handle<reco::VertexCollection>);

      static std::auto_ptr<reco::ConversionCollection> GetCleanedConversions(edm::Handle<reco::ConversionCollection>, 
                                                                             edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromConversion(const reco::TrackRef, reco::ConversionCollection, reco::Conversion*); 

      static std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedKshort(edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);
      static std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedLambda(edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromV0Decay(const reco::TrackRef, reco::VertexCompositeCandidateCollection, 
	 	 	  	   reco::VertexCompositeCandidateCollection, reco::VertexCompositeCandidate*); 

      static std::auto_ptr<reco::PFDisplacedVertexCollection> GetCleanedNI(edm::Handle<reco::PFDisplacedVertexCollection>, edm::Handle<reco::BeamSpot>, bool); 
      static bool ComesFromNI(const reco::TrackRef, reco::PFDisplacedVertexCollection, reco::PFDisplacedVertex*);
                                                   
      static std::auto_ptr<reco::VertexCollection> GetCleanedIVF(edm::Handle<reco::VertexCollection>, edm::Handle<reco::BeamSpot>, bool);
      static bool ComesFromIVF(const reco::TrackRef, reco::VertexCollection, reco::Vertex*);  

      static reco::VertexRef FindVertexZS(reco::TrackRef, edm::Handle<reco::VertexCollection>, double);
      static reco::VertexRef FindVertex3S(reco::TransientTrack, edm::Handle<reco::VertexCollection>, double);

      static reco::VertexRef FindVertexZ(reco::TrackRef, edm::Handle<reco::VertexCollection>, double);
      static reco::VertexRef FindVertex3(reco::TransientTrack, edm::Handle<reco::VertexCollection>, double);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      edm::InputTag tcLabel_;
      edm::InputTag vcLabel_;

      edm::InputTag tpLabel_;

      edm::InputTag beamSpotLabel_;

      bool cleanedColls_;

      edm::InputTag ConversionsCollection_;

      edm::InputTag KshortCollection_;
      edm::InputTag LambdaCollection_;

      edm::InputTag NIVertexCollection_;

      edm::InputTag IVFVertexCollection_;

      TH1F *h_NumTracks_rhoVal;
      
      TH3F *h3_Fin_pur_z_0, *h3_Fin_pur_3_0, *h3_Fin_pur_1_0;
      
      TH3F *h3_Fin_pur_z_01, *h3_Fin_pur_z_001;
      TH3F *h3_Fin_pur_z_0001, *h3_Fin_pur_z_00001;
             
      TH3F *h3_Fin_pur_3_01, *h3_Fin_pur_3_001;
      TH3F *h3_Fin_pur_3_0001, *h3_Fin_pur_3_00001;
      
      TH3F *h3_Fin_pur_z_01S, *h3_Fin_pur_z_001S;
      TH3F *h3_Fin_pur_z_0001S, *h3_Fin_pur_z_00001S;
             
      TH3F *h3_Fin_pur_3_01S, *h3_Fin_pur_3_001S;
      TH3F *h3_Fin_pur_3_0001S, *h3_Fin_pur_3_00001S; 
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
FinalAssocTest::FinalAssocTest(const edm::ParameterSet& iConfig)
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

  cleanedColls_ = iConfig.getParameter<bool>("GetCleanedCollections");

  ConversionsCollection_= iConfig.getParameter<InputTag>("ConversionsCollection");

  KshortCollection_= iConfig.getParameter<InputTag>("KshortCollection");
  LambdaCollection_= iConfig.getParameter<InputTag>("LambdaCollection");

  NIVertexCollection_= iConfig.getParameter<InputTag>("NIVertexCollection");

  IVFVertexCollection_= iConfig.getParameter<InputTag>("IVFVertexCollection");
  
  //--------------


  Service<TFileService> tfs;
  
  h_NumTracks_rhoVal = tfs->make<TH1F>("h_NumTracks_rhoVal", "", rhoNbin, rhoIntervalsa);  
  h_NumTracks_rhoVal->Sumw2();
  
  h3_Fin_pur_z_0 = tfs->make<TH3F>("h3_Fin_pur_z_0", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_3_0 = tfs->make<TH3F>("h3_Fin_pur_3_0", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_1_0 = tfs->make<TH3F>("h3_Fin_pur_1_0", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
    
  h3_Fin_pur_z_0->Sumw2();  
  h3_Fin_pur_3_0->Sumw2(); 
  h3_Fin_pur_1_0->Sumw2(); 
  
  
  h3_Fin_pur_z_01 = tfs->make<TH3F>("h3_Fin_pur_z_01", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_z_001 = tfs->make<TH3F>("h3_Fin_pur_z_001", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_z_0001 = tfs->make<TH3F>("h3_Fin_pur_z_0001", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_z_00001 = tfs->make<TH3F>("h3_Fin_pur_z_00001", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
    
  h3_Fin_pur_z_01->Sumw2();  
  h3_Fin_pur_z_001->Sumw2();  
  h3_Fin_pur_z_0001->Sumw2(); 
  h3_Fin_pur_z_00001->Sumw2();  
  
  
  h3_Fin_pur_3_01 = tfs->make<TH3F>("h3_Fin_pur_3_01", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_3_001 = tfs->make<TH3F>("h3_Fin_pur_3_001", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_3_0001 = tfs->make<TH3F>("h3_Fin_pur_3_0001", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_3_00001 = tfs->make<TH3F>("h3_Fin_pur_3_00001", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_Fin_pur_3_01->Sumw2();  
  h3_Fin_pur_3_001->Sumw2();  
  h3_Fin_pur_3_0001->Sumw2();
  h3_Fin_pur_3_00001->Sumw2();
  
  
  h3_Fin_pur_z_01S = tfs->make<TH3F>("h3_Fin_pur_z_01S", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_z_001S = tfs->make<TH3F>("h3_Fin_pur_z_001S", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_z_0001S = tfs->make<TH3F>("h3_Fin_pur_z_0001S", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_z_00001S = tfs->make<TH3F>("h3_Fin_pur_z_00001S", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
    
  h3_Fin_pur_z_01S->Sumw2();  
  h3_Fin_pur_z_001S->Sumw2();  
  h3_Fin_pur_z_0001S->Sumw2(); 
  h3_Fin_pur_z_00001S->Sumw2();  
  
  
  h3_Fin_pur_3_01S = tfs->make<TH3F>("h3_Fin_pur_3_01S", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_3_001S = tfs->make<TH3F>("h3_Fin_pur_3_001S", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_3_0001S = tfs->make<TH3F>("h3_Fin_pur_3_0001S", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  h3_Fin_pur_3_00001S = tfs->make<TH3F>("h3_Fin_pur_3_00001S", "", tpNbin, tpIntervalsa, rtNbin, rtIntervalsa, rhoNbin, rhoIntervalsa);
  
  h3_Fin_pur_3_01S->Sumw2();  
  h3_Fin_pur_3_001S->Sumw2();  
  h3_Fin_pur_3_0001S->Sumw2();
  h3_Fin_pur_3_00001S->Sumw2();

}


FinalAssocTest::~FinalAssocTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool 
FinalAssocTest::TrackWeightAssociation(const TrackRef&  trackRef, Handle<VertexCollection> vtxcollH) 
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

auto_ptr<ConversionCollection> 
FinalAssocTest::GetCleanedConversions(Handle<ConversionCollection> convCollH, Handle<BeamSpot> bsH, bool cleanedColl)
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
FinalAssocTest::ComesFromConversion(const TrackRef trackref, ConversionCollection cleanedConvColl, Conversion* gamma)
{

	for(unsigned int convcoll_ite=0; convcoll_ite<cleanedConvColl.size(); convcoll_ite++){

	  if(ConversionTools::matchesConversion(trackref,cleanedConvColl.at(convcoll_ite))){
	
	    *gamma = cleanedConvColl.at(convcoll_ite);
	    return true;

  	  }

  	}

	return false;
}

auto_ptr<VertexCompositeCandidateCollection>
FinalAssocTest::GetCleanedKshort(Handle<VertexCompositeCandidateCollection> KshortsH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<VertexCompositeCandidateCollection> cleanedKaonColl(new VertexCompositeCandidateCollection() );

  for (unsigned int kscoll_idx=0; kscoll_idx<KshortsH->size(); kscoll_idx++){

	VertexCompositeCandidateRef ksref(KshortsH,kscoll_idx);

	if(!cleanedColl){   
	  cleanedKaonColl->push_back(*ksref);
	  continue;
	}

	VertexDistanceXY distanceComputer;

	GlobalPoint dec_pos = RecoVertex::convertPos(ksref->vertex());    

	GlobalError decayVertexError = GlobalError(ksref->vertexCovariance(0,0), ksref->vertexCovariance(0,1), ksref->vertexCovariance(1,1), ksref->vertexCovariance(0,2), ksref->vertexCovariance(1,2), ksref->vertexCovariance(2,2));

	math::XYZVector dec_mom(ksref->momentum().x(),
							ksref->momentum().y(),
							ksref->momentum().z());    

	GlobalPoint bsPosition = RecoVertex::convertPos(bsH->position());
	GlobalError bsError = RecoVertex::convertError(bsH->covariance3D());
	
	double kaon_significance = (distanceComputer.distance(VertexState(bsPosition,bsError), VertexState(dec_pos, decayVertexError))).significance();

	if ( (ksref->vertex().rho()>=3.) &&
         (ksref->vertexNormalizedChi2()<=7.) &&
		 (fabs(ksref->mass() - kMass)<=0.07) &&
		 (kaon_significance>15.) ){

	  cleanedKaonColl->push_back(*ksref);

	}

  }
 
  return cleanedKaonColl;

}


auto_ptr<VertexCompositeCandidateCollection>
FinalAssocTest::GetCleanedLambda(Handle<VertexCompositeCandidateCollection> LambdasH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<VertexCompositeCandidateCollection> cleanedLambdaColl(new VertexCompositeCandidateCollection() );

  for (unsigned int lambdacoll_idx=0; lambdacoll_idx<LambdasH->size(); lambdacoll_idx++){

    VertexCompositeCandidateRef lambdaref(LambdasH,lambdacoll_idx);

    if(!cleanedColl){   
      cleanedLambdaColl->push_back(*lambdaref);
      continue;
    }

    VertexDistanceXY distanceComputer;
  
    GlobalPoint dec_pos = RecoVertex::convertPos(lambdaref->vertex());    

    GlobalError decayVertexError = GlobalError(lambdaref->vertexCovariance(0,0), lambdaref->vertexCovariance(0,1), lambdaref->vertexCovariance(1,1), lambdaref->vertexCovariance(0,2), lambdaref->vertexCovariance(1,2), lambdaref->vertexCovariance(2,2));

    math::XYZVector dec_mom(lambdaref->momentum().x(),
						    lambdaref->momentum().y(),
						    lambdaref->momentum().z());    

    GlobalPoint bsPosition = RecoVertex::convertPos(bsH->position());
    GlobalError bsError = RecoVertex::convertError(bsH->covariance3D());

    double lambda_significance = (distanceComputer.distance(VertexState(bsPosition,bsError), VertexState(dec_pos, decayVertexError))).significance();

    if ( (lambdaref->vertex().rho()>=3.) &&
	     (lambdaref->vertexNormalizedChi2()<=7.) &&
	     (fabs(lambdaref->mass() - lamMass)<=0.05) &&
	     (lambda_significance>15.) ){

      cleanedLambdaColl->push_back(*lambdaref);

    }

  }

  return cleanedLambdaColl;
	
}

bool
FinalAssocTest::ComesFromV0Decay(const TrackRef trackref, VertexCompositeCandidateCollection cleanedKshort, VertexCompositeCandidateCollection cleanedLambda, VertexCompositeCandidate* V0)
{

	//the part for the reassociation of particles from Kshort decays
	for(VertexCompositeCandidateCollection::const_iterator iKS=cleanedKshort.begin(); iKS!=cleanedKshort.end(); iKS++){

	  const RecoChargedCandidate *dauCand1 = dynamic_cast<const RecoChargedCandidate*>(iKS->daughter(0));
 	  TrackRef dauTk1 = dauCand1->track();
	  const RecoChargedCandidate *dauCand2 = dynamic_cast<const RecoChargedCandidate*>(iKS->daughter(1));
 	  TrackRef dauTk2 = dauCand2->track();

	  if((trackref==dauTk1) || (trackref==dauTk2)){
	  
	    *V0 = *iKS; 
	    return true;

	  }

	}

	//the part for the reassociation of particles from Lambda decays
	for(VertexCompositeCandidateCollection::const_iterator iLambda=cleanedLambda.begin(); iLambda!=cleanedLambda.end(); iLambda++){

	  const RecoChargedCandidate *dauCand1 = dynamic_cast<const RecoChargedCandidate*>(iLambda->daughter(0));
 	  TrackRef dauTk1 = dauCand1->track();
	  const RecoChargedCandidate *dauCand2 = dynamic_cast<const RecoChargedCandidate*>(iLambda->daughter(1));
 	  TrackRef dauTk2 = dauCand2->track();

   	  if((trackref==dauTk1) || (trackref==dauTk2)){
	  
	    *V0 = *iLambda; 
	    return true;

	  }

	}

	return false;
}

auto_ptr<PFDisplacedVertexCollection>
FinalAssocTest::GetCleanedNI(Handle<PFDisplacedVertexCollection> NuclIntH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<PFDisplacedVertexCollection> cleanedNIColl(new PFDisplacedVertexCollection() );

  for (PFDisplacedVertexCollection::const_iterator niref=NuclIntH->begin(); niref!=NuclIntH->end(); niref++){

    if(!cleanedColl){
      cleanedNIColl->push_back(*niref);
      continue;
    }

    VertexDistanceXY distanceComputer;
   
    GlobalPoint ni_pos = RecoVertex::convertPos(niref->position());    
    GlobalError interactionVertexError = RecoVertex::convertError(niref->error());

    math::XYZVector ni_mom(niref->primaryMomentum().x(),
						   niref->primaryMomentum().y(),
						   niref->primaryMomentum().z());

    GlobalPoint bsPosition = RecoVertex::convertPos(bsH->position());
    GlobalError bsError = RecoVertex::convertError(bsH->covariance3D());

    double nuclint_significance = (distanceComputer.distance(VertexState(bsPosition,bsError), VertexState(ni_pos, interactionVertexError))).significance();

    if ( (niref->position().rho()>=3.) &&
    	!(niref->isFake()) &&
    	 (niref->isNucl()) &&
	     (nuclint_significance>15.) ){

      cleanedNIColl->push_back(*niref);

    }

  }

  return cleanedNIColl;
  
}

bool
FinalAssocTest::ComesFromNI(const TrackRef trackref, PFDisplacedVertexCollection cleanedNI, PFDisplacedVertex* displVtx)
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

auto_ptr<VertexCollection>
FinalAssocTest::GetCleanedIVF(Handle<VertexCollection> ifvH, Handle<BeamSpot> bsH, bool cleanedColl)
{

  auto_ptr<VertexCollection> cleanedIVFColl(new VertexCollection() );

  for (VertexCollection::const_iterator ivfref=ifvH->begin(); ivfref!=ifvH->end(); ivfref++){

    if ( !cleanedColl ) {
      cleanedIVFColl->push_back(*ivfref);
      continue;
    } 
 
    if ( ( ivfref->isValid() ) && 
         ( !ivfref->isFake() ) && 
         ( ivfref->ndof()>4. ) && 
         ( ivfref->chi2()>0. ) && 
         ( ivfref->nTracks(0.)>0. ) ) {

      cleanedIVFColl->push_back(*ivfref);

    }            

  }

  return cleanedIVFColl;
	
}

bool
FinalAssocTest::ComesFromIVF(const TrackRef trackref, VertexCollection cleanedIVF, Vertex* ivfVtx)
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
FinalAssocTest::FindVertexZS(TrackRef trkref, Handle<VertexCollection> vtxCollH, double trackWeight)
{

	double ztrack = trkref->vertex().z();

	VertexRef foundVertexRef(vtxCollH,0);

	double dzmin = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

      VertexRef vertexref(vtxCollH,index_vtx);

	  double nTracks = sqrt(vertexref->tracksSize());
 
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
FinalAssocTest::FindVertex3S(TransientTrack transtrk, Handle<VertexCollection> vtxCollH, double trackWeight)
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

VertexRef
FinalAssocTest::FindVertexZ(TrackRef trkref, Handle<VertexCollection> vtxCollH, double trackWeight)
{

	double ztrack = trkref->vertex().z();

	VertexRef foundVertexRef(vtxCollH,0);

	double dzmin = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

      VertexRef vertexref(vtxCollH,index_vtx);

	  int nTracks = vertexref->tracksSize();
 
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
FinalAssocTest::FindVertex3(TransientTrack transtrk, Handle<VertexCollection> vtxCollH, double trackWeight)
{

	VertexRef foundVertexRef(vtxCollH,0);

	double d3min = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

      VertexRef vertexref(vtxCollH,index_vtx);

	  int nTracks = vertexref->tracksSize();

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

// ------------ method called for each event  ------------
void
FinalAssocTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  //get the conversion collection for the gamma conversions
  Handle<ConversionCollection> convCollH;
  iEvent.getByLabel(ConversionsCollection_, convCollH);
  auto_ptr<ConversionCollection> cleanedConvCollP = FinalAssocTest::GetCleanedConversions(convCollH,bsH,cleanedColls_);

  //get the vertex composite candidate collection for the Kshort's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH;
  iEvent.getByLabel(KshortCollection_, vertCompCandCollKshortH);
  auto_ptr<VertexCompositeCandidateCollection> cleanedKshortCollP = FinalAssocTest::GetCleanedKshort(vertCompCandCollKshortH,bsH,cleanedColls_);
  
  //get the vertex composite candidate collection for the Lambda's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
  iEvent.getByLabel(LambdaCollection_, vertCompCandCollLambdaH);
  auto_ptr<VertexCompositeCandidateCollection> cleanedLambdaCollP = FinalAssocTest::GetCleanedLambda(vertCompCandCollLambdaH,bsH,cleanedColls_);

  //get the displaced vertex collection for nuclear interactions
  Handle<PFDisplacedVertexCollection> displVertexCollH;
  iEvent.getByLabel(NIVertexCollection_,displVertexCollH);
  auto_ptr<PFDisplacedVertexCollection> cleanedNICollP = FinalAssocTest::GetCleanedNI(displVertexCollH,bsH,cleanedColls_);


  //get the vertex collection for inclusive vertex finder   
  Handle<VertexCollection> ivfVertexCollH;
  iEvent.getByLabel(IVFVertexCollection_, ivfVertexCollH);
  auto_ptr<VertexCollection> uncleanedIVFCollP = FinalAssocTest::GetCleanedIVF(ivfVertexCollH, bsH, cleanedColls_);

	    
  //loop over all tracks of the track collection	
  for ( size_t idxTrack = 0; idxTrack < theTracksH->size(); ++idxTrack ) {

    TrackRef trackref(theTracksH, idxTrack);
    TrackBaseRef trackbaseref = TrackBaseRef(trackref);

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

    Conversion gamma;
    VertexCompositeCandidate V0;
    PFDisplacedVertex displVtx;
	Vertex ivfVtx;

    if( ( FinalAssocTest::TrackWeightAssociation(trackref, theVerticesH) ) || 
        ( FinalAssocTest::ComesFromConversion(trackref, *cleanedConvCollP, &gamma) ) || 
        ( FinalAssocTest::ComesFromV0Decay(trackref, *cleanedKshortCollP, *cleanedLambdaCollP, &V0) ) || 
        ( FinalAssocTest::ComesFromNI(trackref, *cleanedNICollP, &displVtx) ) || 
        ( FinalAssocTest::ComesFromIVF(trackref, *uncleanedIVFCollP, &ivfVtx) ) ) continue;

    double rhoVal = transtrk.stateAtBeamLine().transverseImpactParameter().value();
    
    h_NumTracks_rhoVal->Fill( rhoVal );
  
    int y_bin_z_0 = 1;
    if( FinalAssocTest::FindVertexZ(trackref, theVerticesH, 0.)==firstVertex ) y_bin_z_0 = 0;
    int y_bin_z_01 = 1;
    if( FinalAssocTest::FindVertexZ(trackref, theVerticesH, 0.1)==firstVertex ) y_bin_z_01 = 0;
    int y_bin_z_001 = 1;
    if( FinalAssocTest::FindVertexZ(trackref, theVerticesH, 0.01)==firstVertex ) y_bin_z_001 = 0;
    int y_bin_z_0001 = 1;
    if( FinalAssocTest::FindVertexZ(trackref, theVerticesH, 0.001)==firstVertex ) y_bin_z_0001 = 0;
    int y_bin_z_00001 = 1;
    if( FinalAssocTest::FindVertexZ(trackref, theVerticesH, 0.0001)==firstVertex ) y_bin_z_00001 = 0;
  
    int y_bin_3_0 = 1;
    if( FinalAssocTest::FindVertex3(transtrk, theVerticesH, 0.)==firstVertex ) y_bin_3_0 = 0;
    int y_bin_3_01 = 1;
    if( FinalAssocTest::FindVertex3(transtrk, theVerticesH, 0.1)==firstVertex ) y_bin_3_01 = 0;
    int y_bin_3_001 = 1;
    if( FinalAssocTest::FindVertex3(transtrk, theVerticesH, 0.01)==firstVertex ) y_bin_3_001 = 0;
    int y_bin_3_0001 = 1;
    if( FinalAssocTest::FindVertex3(transtrk, theVerticesH, 0.001)==firstVertex ) y_bin_3_0001 = 0;
    int y_bin_3_00001 = 1;
    if( FinalAssocTest::FindVertex3(transtrk, theVerticesH, 0.0001)==firstVertex ) y_bin_3_00001 = 0;

    h3_Fin_pur_z_0->Fill( x_bin, y_bin_z_0, rhoVal );
    h3_Fin_pur_z_01->Fill( x_bin, y_bin_z_01, rhoVal );
    h3_Fin_pur_z_001->Fill( x_bin, y_bin_z_001, rhoVal );
    h3_Fin_pur_z_0001->Fill( x_bin, y_bin_z_0001, rhoVal );
    h3_Fin_pur_z_00001->Fill( x_bin, y_bin_z_00001, rhoVal );

    h3_Fin_pur_3_0->Fill( x_bin, y_bin_3_0, rhoVal );
    h3_Fin_pur_3_01->Fill( x_bin, y_bin_3_01, rhoVal );
    h3_Fin_pur_3_001->Fill( x_bin, y_bin_3_001, rhoVal );
    h3_Fin_pur_3_0001->Fill( x_bin, y_bin_3_0001, rhoVal );
    h3_Fin_pur_3_00001->Fill( x_bin, y_bin_3_00001, rhoVal ); 
  
    int y_bin_z_01S = 1;
    if( FinalAssocTest::FindVertexZS(trackref, theVerticesH, 0.1)==firstVertex ) y_bin_z_01S = 0;
    int y_bin_z_001S = 1;
    if( FinalAssocTest::FindVertexZS(trackref, theVerticesH, 0.01)==firstVertex ) y_bin_z_001S = 0;
    int y_bin_z_0001S = 1;
    if( FinalAssocTest::FindVertexZS(trackref, theVerticesH, 0.001)==firstVertex ) y_bin_z_0001S = 0;
    int y_bin_z_00001S = 1;
    if( FinalAssocTest::FindVertexZS(trackref, theVerticesH, 0.0001)==firstVertex ) y_bin_z_00001S = 0;
  
    int y_bin_3_01S = 1;
    if( FinalAssocTest::FindVertex3S(transtrk, theVerticesH, 0.1)==firstVertex ) y_bin_3_01S = 0;
    int y_bin_3_001S = 1;
    if( FinalAssocTest::FindVertex3S(transtrk, theVerticesH, 0.01)==firstVertex ) y_bin_3_001S = 0;
    int y_bin_3_0001S = 1;
    if( FinalAssocTest::FindVertex3S(transtrk, theVerticesH, 0.001)==firstVertex ) y_bin_3_0001S = 0;
    int y_bin_3_00001S = 1;
    if( FinalAssocTest::FindVertex3S(transtrk, theVerticesH, 0.0001)==firstVertex ) y_bin_3_00001S = 0;

    h3_Fin_pur_z_01S->Fill( x_bin, y_bin_z_01S, rhoVal );
    h3_Fin_pur_z_001S->Fill( x_bin, y_bin_z_001S, rhoVal );
    h3_Fin_pur_z_0001S->Fill( x_bin, y_bin_z_0001S, rhoVal );
    h3_Fin_pur_z_00001S->Fill( x_bin, y_bin_z_00001S, rhoVal );

    h3_Fin_pur_3_01S->Fill( x_bin, y_bin_3_01S, rhoVal );
    h3_Fin_pur_3_001S->Fill( x_bin, y_bin_3_001S, rhoVal );
    h3_Fin_pur_3_0001S->Fill( x_bin, y_bin_3_0001S, rhoVal );
    h3_Fin_pur_3_00001S->Fill( x_bin, y_bin_3_00001S, rhoVal );      

    int y_bin_1_0 = 0;
    h3_Fin_pur_1_0->Fill( x_bin, y_bin_1_0, rhoVal );

  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FinalAssocTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FinalAssocTest);
