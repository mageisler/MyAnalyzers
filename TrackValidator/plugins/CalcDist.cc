// -*- C++ -*-
//
// Package:    CalcDist
// Class:      CalcDist
// 
/**\class CalcDist CalcDist.cc MGeisler/CalcDist/src/CalcDist.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Tue Aug 21 16:20:13 CEST 2012
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

#include "DataFormats/Math/interface/Point3D.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// ROOT include files
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
//
// class declaration
//

class CalcDist : public edm::EDAnalyzer {
   public:
      explicit CalcDist(const edm::ParameterSet&);
      ~CalcDist();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      static bool TrackWeightAssociation(const reco::TrackRef&, edm::Handle<reco::VertexCollection>);
      static reco::VertexRef FindVertexZ(reco::TrackRef, edm::Handle<reco::VertexCollection>, double);
      static reco::VertexRef FindVertex3(reco::TransientTrack, edm::Handle<reco::VertexCollection>, double);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      edm::InputTag tcLabel_;
      edm::InputTag vcLabel_;

      edm::InputTag puLabel_;

      edm::InputTag tpLabel_;

      edm::InputTag beamSpotLabel_;

      TH2F* dist_Vtxdiff_hist;
      TH2F* dist_absVtxdiff_hist;

      TH2F* dist3_Find3_hist;
      TH2F* dist3_FindZ_hist;

      TH2F* distZ_Find3_hist;
      TH2F* distZ_FindZ_hist;

      TH3F *assoc_FindZ_hist, *assoc_Find3_hist;

};

using namespace std;
using namespace edm;
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
CalcDist::CalcDist(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  tcLabel_ = iConfig.getParameter<InputTag>("TrackCollection");
  vcLabel_ = iConfig.getParameter<InputTag>("VertexCollection");

  puLabel_ = iConfig.getParameter<InputTag>("PileUp");

  tpLabel_ = iConfig.getParameter<InputTag>("TrackingParticles");

  beamSpotLabel_ = iConfig.getParameter<InputTag>("BeamSpot");

  //--------------


  Service<TFileService> tfs;

  dist_Vtxdiff_hist = tfs->make<TH2F>("dist_Vtxdiff_hist", "dist_Vtxdiff_hist; (wrtBS - noBS) / cm; #rho__{BS}^{Reco} / cm", 101, -5.05, 5.05, 100, 0., 10.);
  dist_absVtxdiff_hist = tfs->make<TH2F>("dist_absVtxdiff_hist", "dist_absVtxdiff_hist; (wrtBS - noBS) / cm; #rho__{BS}^{Reco} / cm", 101, -5.05, 5.05, 100, 0., 10.);

  dist3_Find3_hist = tfs->make<TH2F>("dist3_Find3_hist", "dist3_Find3_hist; 3D distance / cm; #rho__{BS}^{Reco} / cm", 51, -0.05, 5.05, 100, 0., 10.);
  dist3_FindZ_hist = tfs->make<TH2F>("dist3_FindZ_hist", "dist3_FindZ_hist; 3D distance / cm; #rho__{BS}^{Reco} / cm", 51, -0.05, 5.05, 100, 0., 10.);

  distZ_Find3_hist = tfs->make<TH2F>("distZ_Find3_hist", "distZ_Find3_hist; z distance / cm; #rho__{BS}^{Reco} / cm", 51, -0.05, 5.05, 100, 0., 10.);
  distZ_FindZ_hist = tfs->make<TH2F>("distZ_FindZ_hist", "distZ_FindZ_hist; z distance / cm; #rho__{BS}^{Reco} / cm", 51, -0.05, 5.05, 100, 0., 10.);

  assoc_FindZ_hist = tfs->make<TH3F>("assoc_FindZ_hist", "assoc_FindZ_hist; sim association; reco association; #rho_{BS}^{Reco} / cm", 2, -0.5, 1.5, 2, -0.5, 1.5, 100, 0., 10.);
  assoc_Find3_hist = tfs->make<TH3F>("assoc_Find3_hist", "assoc_Find3_hist; sim association; reco association; #rho_{BS}^{Reco} / cm", 2, -0.5, 1.5, 2, -0.5, 1.5, 100, 0., 10.);

}


CalcDist::~CalcDist()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

bool 
CalcDist::TrackWeightAssociation(const TrackRef&  trackRef, Handle<VertexCollection> vtxcollH) 
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
CalcDist::FindVertexZ(TrackRef trkref, Handle<VertexCollection> vtxCollH, double trackWeight)
{

	double ztrack = trkref->vertex().z();

	VertexRef foundVertexRef(vtxCollH,0);

	double dzmin = 5.;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

          VertexRef vertexref(vtxCollH,index_vtx);

	  int nTracks = sqrt(vertexref->tracksSize());
 
	  //find and store the closest vertex in z
          double distance = fabs(ztrack - vertexref->z());

	  double weightedDistance = max(0.0,distance-trackWeight*nTracks);	

          if(weightedDistance<dzmin) {
            dzmin = weightedDistance; 
            foundVertexRef = vertexref;
          }
	
	}

	return foundVertexRef;

}

VertexRef
CalcDist::FindVertex3(TransientTrack transtrk, Handle<VertexCollection> vtxCollH, double trackWeight)
{

	VertexRef foundVertexRef(vtxCollH,0);

	double d3min = 5.;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

          VertexRef vertexref(vtxCollH,index_vtx);

	  GlobalPoint vtxpos(vertexref->x(),vertexref->y(),vertexref->z());
	  GlobalPoint closestPoint = transtrk.trajectoryStateClosestToPoint(vtxpos).position();

	  int nTracks = sqrt(vertexref->tracksSize());
 
	  //find and store the closest vertex in z
	  double x_dist = vtxpos.x() - closestPoint.x();
	  double y_dist = vtxpos.y() - closestPoint.y();
	  double z_dist = vtxpos.z() - closestPoint.z();

          double distance = sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist);

	  double weightedDistance = max(0.0,distance-trackWeight*nTracks);	

          if(weightedDistance<d3min) {
            d3min = weightedDistance; 
            foundVertexRef = vertexref;
          }
	
	}

	return foundVertexRef;

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CalcDist::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  //get the pileup information  
  Handle< vector<PileupSummaryInfo> > puinfoH;
  iEvent.getByLabel(puLabel_,puinfoH);
  PileupSummaryInfo puinfo;      
  
  for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){ 
    if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
      puinfo=(*puinfoH)[puinfo_ite];
      break;
    }
  }

  //get the tracking vertices   
  Handle<TrackingVertexCollection> theTVsH ;
  iEvent.getByLabel(tpLabel_,theTVsH);

  vector< vector<TrackingVertexRef> > realTrackingV;
  vector< vector<bool> > realTrackingVAvail;
  
  for (unsigned int puinfo_ite=0;puinfo_ite<puinfoH->size();++puinfo_ite){ 
  
    vector<TrackingVertexRef> help;
	
    for(int coll_ite=0; coll_ite<=puinfoH->at(puinfo_ite).getPU_NumInteractions(); coll_ite++){

      TrackingVertexRef tv_help;
       help.push_back(tv_help);

    }

    realTrackingV.push_back(help);

  }

  int minBC = puinfoH->at(0).getBunchCrossing();
  
  for(TrackingVertexCollection::size_type tv_ite=0; tv_ite<theTVsH->size(); tv_ite++){

    TrackingVertexRef tvr(theTVsH, tv_ite);

    if( (tvr->eventId().event()<(signed int)realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).size()) &&
        (realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).at(tvr->eventId().event()).isNull()) &&
        (tvr->nSourceTracks()==0) ){

      realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).at(tvr->eventId().event()) = tvr;

    }

  }
	    
  //loop over all tracks of the track collection	
  for ( size_t idxTrack = 0; idxTrack < theTracksH->size(); ++idxTrack ) {

    TrackRef trackref = TrackRef(theTracksH, idxTrack);
    TrackBaseRef trackbaseref = TrackBaseRef(trackref);

    TransientTrack transtrk(trackref, &(*bFieldHandle) );
    transtrk.setBeamSpot(*bsH);

    vector<pair<TrackingParticleRef, double> > tp;
    if(recSimColl.find(trackbaseref) != recSimColl.end()) tp = recSimColl[trackbaseref];

    if (tp.size()==0) continue;

    TrackingParticleRef tpr = tp[0].first;

    if (realTrackingV.at(tpr->eventId().bunchCrossing()-minBC).at(tpr->eventId().event()).isNull()) continue;

    double z_gen = (realTrackingV.at(tpr->eventId().bunchCrossing()-minBC).at(tpr->eventId().event()))->position().z();
    double x_gen = (realTrackingV.at(tpr->eventId().bunchCrossing()-minBC).at(tpr->eventId().event()))->position().x();
    double y_gen = (realTrackingV.at(tpr->eventId().bunchCrossing()-minBC).at(tpr->eventId().event()))->position().y();

    if ( !(CalcDist::TrackWeightAssociation(trackref, theVerticesH)) ){

      int x_bin = 1;
      if((tpr->eventId().event() == 0) && (tpr->eventId().bunchCrossing() == 0)) x_bin=0;

      int yZ_bin = 1;
      if(CalcDist::FindVertexZ(trackref,theVerticesH,0.01)==firstVertex) yZ_bin=0;

      int y3_bin = 1;
      if(CalcDist::FindVertex3(transtrk,theVerticesH,0.01)==firstVertex) y3_bin=0;

      double z_recoZ = (CalcDist::FindVertexZ(trackref,theVerticesH,0.01))->position().z();
      double x_recoZ = (CalcDist::FindVertexZ(trackref,theVerticesH,0.01))->position().x();
      double y_recoZ = (CalcDist::FindVertexZ(trackref,theVerticesH,0.01))->position().y();

      double z_reco3 = (CalcDist::FindVertex3(transtrk,theVerticesH,0.01))->position().z();
      double x_reco3 = (CalcDist::FindVertex3(transtrk,theVerticesH,0.01))->position().x();
      double y_reco3 = (CalcDist::FindVertex3(transtrk,theVerticesH,0.01))->position().y();

      double rho = transtrk.stateAtBeamLine().transverseImpactParameter().value();

      double ZdistZ = fabs(z_gen - z_recoZ);
      double XdistZ = fabs(x_gen - x_recoZ);
      double YdistZ = fabs(y_gen - y_recoZ);

      double Zdist3 = fabs(z_gen - z_reco3);
      double Xdist3 = fabs(x_gen - x_reco3);
      double Ydist3 = fabs(y_gen - y_reco3);

      double AbsdistZ = sqrt(ZdistZ*ZdistZ + XdistZ*XdistZ + YdistZ*YdistZ);
      double Absdist3 = sqrt(Zdist3*Zdist3 + Xdist3*Xdist3 + Ydist3*Ydist3);

      dist_Vtxdiff_hist->Fill(ZdistZ - Zdist3,rho);
      dist_absVtxdiff_hist->Fill(AbsdistZ - Absdist3,rho);

      dist3_Find3_hist->Fill(Absdist3,rho);
      dist3_FindZ_hist->Fill(AbsdistZ,rho);

      distZ_Find3_hist->Fill(Zdist3,rho);
      distZ_FindZ_hist->Fill(ZdistZ,rho);

      assoc_Find3_hist->Fill(x_bin,y3_bin,rho);
      assoc_FindZ_hist->Fill(x_bin,yZ_bin,rho);

    }

  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CalcDist::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CalcDist);