// -*- C++ -*-
//
// Package:    FilterFinder
// Class:      FilterFinder
// 
/**\class FilterFinder FilterFinder.cc MGeisler/FilterFinder/src/FilterFinder.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Wed Sep 19 09:27:08 CEST 2012
// $Id$
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

//
// class declaration
//

class FilterFinder : public edm::EDAnalyzer {
   public:
      explicit FilterFinder(const edm::ParameterSet&);
      ~FilterFinder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static bool TrackWeightAssociation(reco::TrackBaseRef&, edm::Handle<reco::VertexCollection>);

      static std::auto_ptr<reco::ConversionCollection> GetCleanedConversions(edm::Handle<reco::ConversionCollection>);
      static bool ComesFromConversion(const reco::TrackRef, reco::ConversionCollection, reco::Conversion*);         
      static reco::VertexRef FindConversionVertex3D(const reco::TrackRef, reco::Conversion, 	
                                                    edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                    edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, 
	                                            double);     

      static std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedV0(edm::Handle<reco::VertexCompositeCandidateCollection>);
      static bool ComesFromV0Decay(const reco::TrackRef, reco::VertexCompositeCandidateCollection,                
                                   reco::VertexCompositeCandidate*);      
      static reco::VertexRef FindV0DecayVertex3D(const reco::TrackRef, reco::VertexCompositeCandidate, 	
                                                 edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                 edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, 
	                                         double);

      static std::auto_ptr<reco::PFDisplacedVertexCollection> GetCleanedNI(edm::Handle<reco::PFDisplacedVertexCollection>);
      static bool ComesFromNI(const reco::TrackRef, reco::PFDisplacedVertexCollection, reco::PFDisplacedVertex*);    
      static reco::VertexRef FindNIVertex3D(const reco::TrackRef, reco::PFDisplacedVertex, 	
                                            edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				            edm::Handle<reco::BeamSpot>, edm::Handle<reco::VertexCollection>, 
	                                    double);

      static reco::VertexRef FindVertex3D(reco::TransientTrack, edm::Handle<reco::VertexCollection>, double);

      static double dR(math::XYZPoint, math::XYZVector, edm::Handle<reco::BeamSpot>);


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

  //--------------

      TH3F* conversion_nTracks, *conversion_pairInvariantMass, *conversion_radius, *conversion_deltaR;

      TH3F* kaonDecay_chi2, *kaonDecay_deltaMass, *kaonDecay_significance, *kaonDecay_deltaR;
      TH3F* lambdaDecay_chi2, *lambdaDecay_deltaMass, *lambdaDecay_significance, *lambdaDecay_deltaR;

      TH3F* nuclearInteraction_significance, *nuclearInteraction_deltaR;
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
FilterFinder::FilterFinder(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  tcLabel_ = iConfig.getParameter<InputTag>("TrackCollection");
  vcLabel_ = iConfig.getParameter<InputTag>("VertexCollection");

  tpLabel_ = iConfig.getParameter<InputTag>("TrackingParticles");

  beamSpotLabel_ = iConfig.getParameter<InputTag>("BeamSpot");

  ConversionsCollection_= iConfig.getParameter<InputTag>("ConversionsCollection");

  KshortCollection_= iConfig.getParameter<InputTag>("KshortCollection");
  LambdaCollection_= iConfig.getParameter<InputTag>("LambdaCollection");

  NIVertexCollection_= iConfig.getParameter<InputTag>("NIVertexCollection");

  //--------------


  Service<TFileService> tfs;

  conversion_nTracks = tfs->make<TH3F>("conversion_nTracks", "conversion_nTracks; sim association; reco association; Number of Tracks", 3, -0.5, 2.5, 2, -0.5, 1.5, 4, -0.5, 3.5);
  conversion_pairInvariantMass = tfs->make<TH3F>("conversion_pairInvariantMass", "conversion_pairInvariantMass; sim association; reco association; Pair Invariant Mass / GeV", 3, -0.5, 2.5, 2, -0.5, 1.5, 100, 0., 1.);
  conversion_radius = tfs->make<TH3F>("conversion_radius", "conversion_radius; sim association; reco association; d_{xy} / cm", 3, -0.5, 2.5, 2, -0.5, 1.5, 100, 0., 10.);
  conversion_deltaR = tfs->make<TH3F>("conversion_deltaR", "conversion_nTracks; sim association; reco association; deltaR(#gamma,Vtx-BS)", 3, -0.5, 2.5, 2, -0.5, 1.5, 200, 0., 0.8);


  kaonDecay_chi2 = tfs->make<TH3F>("kaonDecay_chi2", "kaonDecay_chi2; sim association; reco association; #chi^{2}", 3, -0.5, 2.5, 2, -0.5, 1.5, 80, 0., 8.);
  kaonDecay_deltaMass = tfs->make<TH3F>("kaonDecay_deltaMass", "kaonDecay_deltaMass; sim association; reco association; #Delta_{Mass} / GeV", 3, -0.5, 2.5, 2, -0.5, 1.5, 150, 0., 0.15);
  kaonDecay_significance = tfs->make<TH3F>("kaonDecay_significance", "kaonDecay_significance; sim association; reco association; #sigma_{Distance}", 3, -0.5, 2.5, 2, -0.5, 1.5, 25, 0., 50.);
  kaonDecay_deltaR = tfs->make<TH3F>("kaonDecay_deltaR", "kaonDecay_deltaR; sim association; reco association; deltaR(K,Vtx-BS)", 3, -0.5, 2.5, 2, -0.5, 1.5, 100, 0., 1.);

  lambdaDecay_chi2 = tfs->make<TH3F>("lambdaDecay_chi2", "lambdaDecay_chi2; sim association; reco association; #chi^{2}", 3, -0.5, 2.5, 2, -0.5, 1.5, 80, 0., 8.);
  lambdaDecay_deltaMass = tfs->make<TH3F>("lambdaDecay_deltaMass", "lambdaDecay_deltaMass; sim association; reco association; #Delta_{Mass} / GeV", 3, -0.5, 2.5, 2, -0.5, 1.5, 150, 0., 0.15);
  lambdaDecay_significance = tfs->make<TH3F>("lambdaDecay_significance", "lambdaDecay_significance; sim association; reco association; #sigma_{Distance}", 3, -0.5, 2.5, 2, -0.5, 1.5, 50, 0., 50.);
  lambdaDecay_deltaR = tfs->make<TH3F>("lambdaDecay_deltaR", "lambdaDecay_deltaR; sim association; reco association; deltaR(#lambda,Vtx-BS)", 3, -0.5, 2.5, 2, -0.5, 1.5, 100, 0., 1.);


  nuclearInteraction_significance = tfs->make<TH3F>("nuclearInteraction_significance", "nuclearInteraction_significance; sim association; reco association; #sigma_{Distance}", 3, -0.5, 2.5, 2, -0.5, 1.5, 50, 0., 50.);
  nuclearInteraction_deltaR = tfs->make<TH3F>("nuclearInteraction_deltaR", "nuclearInteraction_deltaR; sim association; reco association; deltaR(Incom,Vtx-BS)", 3, -0.5, 2.5, 2, -0.5, 1.5, 100, 0., 1.);

}


FilterFinder::~FilterFinder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool 
FilterFinder::TrackWeightAssociation(TrackBaseRef&  trackbaseRef, Handle<VertexCollection> vtxcollH) 
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
FilterFinder::FindVertex3D(TransientTrack transtrk, Handle<VertexCollection> vtxCollH, double trackWeight)
{

	VertexRef foundVertexRef(vtxCollH,0);

	double d3min = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

          VertexRef vertexref(vtxCollH,index_vtx);

	  double nTracks = sqrt(vertexref->tracksSize());

          double distance = 10001.;	        
          pair<bool,Measurement1D> IpPair = IPTools::absoluteImpactParameter3D(transtrk, *vertexref);
 
	  if(IpPair.first)
            distance = IpPair.second.value();

	  double weightedDistance = distance-trackWeight*nTracks;	

          if(weightedDistance<d3min) {
            d3min = weightedDistance; 
            foundVertexRef = vertexref;
          }
	
	}

	return foundVertexRef;

}

auto_ptr<ConversionCollection> 
FilterFinder::GetCleanedConversions(Handle<ConversionCollection> convCollH)
{
     	auto_ptr<ConversionCollection> cleanedConvColl(new ConversionCollection() );

	for (unsigned int convcoll_idx=0; convcoll_idx<convCollH->size(); convcoll_idx++){

	  ConversionRef convref(convCollH,convcoll_idx);
  
          cleanedConvColl->push_back(*convref);

      	}

  	return cleanedConvColl;

}

bool
FilterFinder::ComesFromConversion(const TrackRef trackref, ConversionCollection cleanedConvColl, Conversion* gamma)
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
FilterFinder::FindConversionVertex3D(const reco::TrackRef trackref, reco::Conversion gamma, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

	math::XYZPoint conv_pos = gamma.conversionVertex().position();

	math::XYZVector conv_mom(gamma.refittedPair4Momentum().x(),
	                         gamma.refittedPair4Momentum().y(),
	                         gamma.refittedPair4Momentum().z());

	Track photon(trackref->chi2(), trackref->ndof(), conv_pos, conv_mom, 0, trackref->covariance());

    	TransientTrack transpho(photon, &(*bFieldH) );
    	transpho.setBeamSpot(*bsH);
    	transpho.setES(iSetup);

	return FilterFinder::FindVertex3D(transpho, vtxcollH, tWeight); 

}

auto_ptr<VertexCompositeCandidateCollection>
FilterFinder::GetCleanedV0(Handle<VertexCompositeCandidateCollection> V0sH)
{

     	auto_ptr<VertexCompositeCandidateCollection> cleanedColl(new VertexCompositeCandidateCollection() );

	for (unsigned int V0coll_idx=0; V0coll_idx<V0sH->size(); V0coll_idx++){

	  VertexCompositeCandidateRef V0ref(V0sH,V0coll_idx);

	  if (V0ref->vertex().rho()<=3.) continue;
  
          cleanedColl->push_back(*V0ref);

	}

	return cleanedColl;

}

bool
FilterFinder::ComesFromV0Decay(const TrackRef trackref, VertexCompositeCandidateCollection cleanedV0, VertexCompositeCandidate* V0)
{

	//the part for the reassociation of particles from V0 decays
	for(VertexCompositeCandidateCollection::const_iterator iV0=cleanedV0.begin(); iV0!=cleanedV0.end(); iV0++){

	  const RecoChargedCandidate *dauCand1 = dynamic_cast<const RecoChargedCandidate*>(iV0->daughter(0));
 	  TrackRef dauTk1 = dauCand1->track();
	  const RecoChargedCandidate *dauCand2 = dynamic_cast<const RecoChargedCandidate*>(iV0->daughter(1));
 	  TrackRef dauTk2 = dauCand2->track();

	  if((trackref==dauTk1) || (trackref==dauTk2)){
	  
	    *V0 = *iV0; 
	    return true;

	  }

	}

	return false;
}

VertexRef
FilterFinder::FindV0DecayVertex3D(const reco::TrackRef trackref, reco::VertexCompositeCandidate V0_vtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
{ 

	math::XYZPoint dec_pos = V0_vtx.vertex();

	math::XYZVector dec_mom(V0_vtx.momentum().x(),
	                        V0_vtx.momentum().y(),
	                        V0_vtx.momentum().z());

	Track V0(V0_vtx.vertexChi2(), V0_vtx.vertexNdof(), dec_pos, dec_mom, 0, trackref->covariance());

    	TransientTrack transV0(V0, &(*bFieldH) );
    	transV0.setBeamSpot(*bsH);
    	transV0.setES(iSetup);

	return FilterFinder::FindVertex3D(transV0, vtxcollH, tWeight); 

}

auto_ptr<PFDisplacedVertexCollection>
FilterFinder::GetCleanedNI(Handle<PFDisplacedVertexCollection> NuclIntH)
{

     	auto_ptr<PFDisplacedVertexCollection> cleanedNIColl(new PFDisplacedVertexCollection() );

	for (PFDisplacedVertexCollection::const_iterator niref=NuclIntH->begin(); niref!=NuclIntH->end(); niref++){

	  if ( (niref->position().rho()<=3.) || (niref->isFake()) || !(niref->isNucl()) ) continue;

	  cleanedNIColl->push_back(*niref);

	}

	return cleanedNIColl;
}

bool
FilterFinder::ComesFromNI(const TrackRef trackref, PFDisplacedVertexCollection cleanedNI, PFDisplacedVertex* displVtx)
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
FilterFinder::FindNIVertex3D(const reco::TrackRef trackref, reco::PFDisplacedVertex displVtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, edm::Handle<reco::VertexCollection> vtxcollH, double tWeight)
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

	      if(bestweight>1.e-5) 
                return bestvertexref;

    	      TransientTrack transIncom(*retrackbaseref, &(*bFieldH) );
    	      transIncom.setBeamSpot(*bsH);
    	      transIncom.setES(iSetup);

	      return FindVertex3D(transIncom, vtxcollH, tWeight); 

	    }

	  }

	}

	math::XYZPoint ni_pos = displVtx.position();

	math::XYZVector ni_mom(displVtx.primaryMomentum().x(),
	                       displVtx.primaryMomentum().y(),
	                       displVtx.primaryMomentum().z());

	Track incom(trackref->chi2(), trackref->ndof(), ni_pos, ni_mom, 0, trackref->covariance());

    	TransientTrack transIncom(incom, &(*bFieldH) );
    	transIncom.setBeamSpot(*bsH);
    	transIncom.setES(iSetup);

 	return FindVertex3D(transIncom, vtxcollH, tWeight); 

}

double
FilterFinder::dR(math::XYZPoint vtx_pos, math::XYZVector vtx_mom, edm::Handle<reco::BeamSpot> bsH)
{

	double bs_x = bsH->x0();
	double bs_y = bsH->y0();
	double bs_z = bsH->z0();

     	double connVec_x = vtx_pos.x() - bs_x;
	double connVec_y = vtx_pos.y() - bs_y;
	double connVec_z = vtx_pos.z() - bs_z;

     	double connVec_r = sqrt(connVec_x*connVec_x + connVec_y*connVec_y + connVec_z*connVec_z);
	double connVec_theta = acos(connVec_z*1./connVec_r);

	double connVec_eta = -1.*log(tan(connVec_theta*1./2.));
	double connVec_phi = atan2(connVec_y,connVec_x);

	return deltaR(vtx_mom.eta(),vtx_mom.phi(),connVec_eta,connVec_phi);
    
}

// ------------ method called for each event  ------------
void
FilterFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  //get the reco vertices   
  Handle<VertexCollection>  theVerticesH;
  iEvent.getByLabel(vcLabel_,theVerticesH);

  VertexRef firstVertex(theVerticesH, 0);

  VertexDistance3D distanceComputer;

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
  auto_ptr<ConversionCollection> cleanedConvCollP = FilterFinder::GetCleanedConversions(convCollH);

  //get the vertex composite candidate collection for the Kshort's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH;
  iEvent.getByLabel(KshortCollection_, vertCompCandCollKshortH);
  auto_ptr<VertexCompositeCandidateCollection> cleanedKshortCollP = FilterFinder::GetCleanedV0(vertCompCandCollKshortH);
  
  //get the vertex composite candidate collection for the Lambda's
  Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
  iEvent.getByLabel(LambdaCollection_, vertCompCandCollLambdaH);
  auto_ptr<VertexCompositeCandidateCollection> cleanedLambdaCollP = FilterFinder::GetCleanedV0(vertCompCandCollLambdaH);

  //get the displaced vertex collection for nuclear interactions
  Handle<PFDisplacedVertexCollection> displVertexCollH;
  iEvent.getByLabel(NIVertexCollection_,displVertexCollH);
  auto_ptr<PFDisplacedVertexCollection> cleanedNICollP = FilterFinder::GetCleanedNI(displVertexCollH);

	    
  //loop over all tracks of the track collection	
  for ( size_t idxTrack = 0; idxTrack < theTracksH->size(); ++idxTrack ) {

    TrackRef trackref(theTracksH, idxTrack);
    TrackBaseRef trackbaseref = TrackBaseRef(trackref);

    if(FilterFinder::TrackWeightAssociation(trackbaseref, theVerticesH)) continue;

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

    //Conversion
    Conversion gamma;
    if ( FilterFinder::ComesFromConversion(trackref, *cleanedConvCollP, &gamma) ){

      int y_bin_Conversion = 1;
      if( FilterFinder::FindConversionVertex3D(trackref, gamma, bfH, iSetup, bsH, theVerticesH, 0.01)==firstVertex ) 
        y_bin_Conversion = 0;


      math::XYZPoint conv_pos = gamma.conversionVertex().position();

      math::XYZVector conv_mom(gamma.refittedPair4Momentum().x(),
	                       gamma.refittedPair4Momentum().y(),
	                       gamma.refittedPair4Momentum().z());
  
      double conversion_dR = FilterFinder::dR(conv_pos,conv_mom,bsH);

      conversion_nTracks->Fill( x_bin,y_bin_Conversion,gamma.nTracks() );
      conversion_pairInvariantMass->Fill( x_bin,y_bin_Conversion,gamma.pairInvariantMass() ); 
      conversion_radius->Fill( x_bin,y_bin_Conversion,gamma.dxy(bsH->position()) );
      conversion_deltaR->Fill( x_bin,y_bin_Conversion,conversion_dR );

    }

    //Kaon
    VertexCompositeCandidate V0;
    if ( FilterFinder::ComesFromV0Decay(trackref, *cleanedKshortCollP, &V0) ) {

      VertexRef foundPrimaryVertex = FilterFinder::FindV0DecayVertex3D(trackref, V0, bfH, iSetup, bsH, theVerticesH, 0.01);

      int y_bin_Kaon = 1;
      if( foundPrimaryVertex==firstVertex ) y_bin_Kaon = 0;

      GlobalPoint dec_pos = RecoVertex::convertPos(V0.vertex());    

      GlobalError decayVertexError = GlobalError(V0.vertexCovariance(0,0), V0.vertexCovariance(0,1), V0.vertexCovariance(1,1), V0.vertexCovariance(0,2), V0.vertexCovariance(1,2), V0.vertexCovariance(2,2));

      math::XYZVector dec_mom(V0.momentum().x(),
	                      V0.momentum().y(),
	                      V0.momentum().z());    

      GlobalPoint primaryVertexPosition = RecoVertex::convertPos(foundPrimaryVertex->position());
      GlobalError primaryVertexError = RecoVertex::convertError(foundPrimaryVertex->error());

      double kaon_chi2 = V0.vertexNormalizedChi2();
      double kaon_dM = fabs(V0.mass() - kMass);
      double kaon_significance = (distanceComputer.distance(VertexState(primaryVertexPosition,primaryVertexError), VertexState(dec_pos, decayVertexError))).significance();
      double kaon_dR = FilterFinder::dR(V0.vertex(),dec_mom,bsH);

      kaonDecay_chi2->Fill( x_bin,y_bin_Kaon,kaon_chi2 );
      kaonDecay_deltaMass->Fill( x_bin,y_bin_Kaon,kaon_dM ); 
      kaonDecay_significance->Fill( x_bin,y_bin_Kaon,kaon_significance );
      kaonDecay_deltaR->Fill( x_bin,y_bin_Kaon,kaon_dR );

    }

    //Lambda
    if ( FilterFinder::ComesFromV0Decay(trackref, *cleanedLambdaCollP, &V0) ) {

      VertexRef foundPrimaryVertex = FilterFinder::FindV0DecayVertex3D(trackref, V0, bfH, iSetup, bsH, theVerticesH, 0.01);

      int y_bin_Lambda = 1;
      if( foundPrimaryVertex==firstVertex ) y_bin_Lambda = 0;

      GlobalPoint dec_pos = RecoVertex::convertPos(V0.vertex());    

      GlobalError decayVertexError = GlobalError(V0.vertexCovariance(0,0), V0.vertexCovariance(0,1), V0.vertexCovariance(1,1), V0.vertexCovariance(0,2), V0.vertexCovariance(1,2), V0.vertexCovariance(2,2));

      math::XYZVector dec_mom(V0.momentum().x(),
	                      V0.momentum().y(),
	                      V0.momentum().z()); 

      GlobalPoint primaryVertexPosition = RecoVertex::convertPos(foundPrimaryVertex->position());
      GlobalError primaryVertexError = RecoVertex::convertError(foundPrimaryVertex->error());
  
      double lambda_chi2 = V0.vertexNormalizedChi2();
      double lambda_dM = fabs(V0.mass() - lamMass);
      double lambda_significance = (distanceComputer.distance(VertexState(primaryVertexPosition,primaryVertexError), VertexState(dec_pos, decayVertexError))).significance();
      double lambda_dR = FilterFinder::dR(V0.vertex(),dec_mom,bsH);

      lambdaDecay_chi2->Fill( x_bin,y_bin_Lambda,lambda_chi2 );
      lambdaDecay_deltaMass->Fill( x_bin,y_bin_Lambda,lambda_dM ); 
      lambdaDecay_significance->Fill( x_bin,y_bin_Lambda,lambda_significance );
      lambdaDecay_deltaR->Fill( x_bin,y_bin_Lambda,lambda_dR );

    }

    //Nuclear Interaction
    PFDisplacedVertex displVtx;
    if ( FilterFinder::ComesFromNI(trackref, *cleanedNICollP, &displVtx) ) {
  
      VertexRef foundPrimaryVertex = FilterFinder::FindNIVertex3D(trackref, displVtx, bfH, iSetup, bsH, theVerticesH, 0.01);

      int y_bin_NuclearInteraction = 1;
      if( foundPrimaryVertex==firstVertex ) 
        y_bin_NuclearInteraction = 0;

      GlobalPoint ni_pos = RecoVertex::convertPos(displVtx.position());    
      GlobalError interactionVertexError = RecoVertex::convertError(displVtx.error());

      math::XYZVector ni_mom(displVtx.primaryMomentum().x(),
	                     displVtx.primaryMomentum().y(),
	                     displVtx.primaryMomentum().z());

      GlobalPoint primaryVertexPosition = RecoVertex::convertPos(foundPrimaryVertex->position());
      GlobalError primaryVertexError = RecoVertex::convertError(foundPrimaryVertex->error());
  
      double nuclInteraction_significance = (distanceComputer.distance(VertexState(primaryVertexPosition,primaryVertexError), VertexState(ni_pos, interactionVertexError))).significance();
      double nuclearInteraction_dR = FilterFinder::dR(displVtx.position(),ni_mom,bsH);

      nuclearInteraction_significance->Fill( x_bin,y_bin_NuclearInteraction,nuclInteraction_significance );
      nuclearInteraction_deltaR->Fill( x_bin,y_bin_NuclearInteraction,nuclearInteraction_dR );

    }

  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FilterFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FilterFinder);
