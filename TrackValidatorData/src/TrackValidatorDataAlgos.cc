#include "MGeisler/TrackValidatorData/interface/TrackValidatorDataAlgos.h"

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "TLorentzVector.h"
   
using namespace edm;
using namespace std;
using namespace reco;

TrackValidatorDataAlgos::TrackValidatorDataAlgos(const edm::ParameterSet& iConfig)
{ 
  
  //parameters for vs_eta plots
  minEta  = -2.5;  
  maxEta  = 2.5;
  nintEta = 50;

  //parameters for vs_pt plots
  minPt  = 0.1;
  maxPt  = 100.;
  nintPt = 40;

  //parameters for vs_minDist plots
  minD  = 0.;
  maxD  = 5.;
  maxDZ  = 0.5;
  nintD = 500;

  //parameters for vs_quality plots
  minQ  = -0.5;
  maxQ  = 5.5;
  nintQ = 6;

  //parameters for vs_tweight plots
  minTW  = 0.;
  maxTW  = 1.;
  nintTW = 100;
  
  //parameters for track number plots
  minTrackcount  = -0.5;
  maxTrackcount  = 999.5;
  nintTrackcount = 250;
  
  //parameters for Pileup plots
  minVertcount  = -0.5;
  maxVertcount  = 59.5;
  nintVertcount = 60;
  
  //parameters for invariant mass plots
  minM  = 50.;
  maxM  = 130.;
  nintM = 40;

}

void 
TrackValidatorDataAlgos::BookHistos(TFileDirectory subDir) 
{ 

  //histograms

  h_numVtx.push_back(subDir.make<TH1F>("h_numVtx", "number of reconstructed vertices per event; number of reco vertices; number of events", nintVertcount, minVertcount, maxVertcount));
  
  h_numVtx.back()->Sumw2();

  h_numTrks.push_back(subDir.make<TH1F>("h_numTrks", "number of considered signal tracks per event; number of signal tracks; number of events", nintTrackcount, minTrackcount, maxTrackcount));
  
  h_numTrks.back()->Sumw2();

  h_minDist.push_back(subDir.make<TH1F>("h_minDist", "minimum distance to signal vertex; d_{min} / cm; number of tracks", nintD, minD, maxD));
  h_quality.push_back(subDir.make<TH1F>("h_quality", "quality; quality class; number of associations", nintQ, minQ, maxQ));
  h_tweight.push_back(subDir.make<TH1F>("h_tweight", "track weight; signal track weight; number of associations", nintTW, minTW, maxTW));
  
  h_minDist.back()->Sumw2();
  h_quality.back()->Sumw2();
  h_tweight.back()->Sumw2();

  h_numTrks_eta.push_back(subDir.make<TH1F>("h_numTrks_eta", "number of considered signal tracks vs eta; #eta; number of tracks", nintEta, minEta, maxEta));
  h_numTrks_pt.push_back(subDir.make<TH1F>("h_numTrks_pt", "number of considered signal tracks vs pt; p_{t} / GeV; number of tracks", nintPt, minPt, maxPt));
  
  h_numTrks_eta.back()->Sumw2();
  h_numTrks_pt.back()->Sumw2();

  h_zmumunum_eta.push_back(subDir.make<TH1F>("h_zmumunum_eta", "found muons from Z decays vs eta; #eta; number of muons", nintEta, minEta, maxEta));
  h_zmumunum_pt.push_back(subDir.make<TH1F>("h_zmumunum_pt", "found muons from Z decays vs pt; p_{t} / GeV; number of muons", nintPt, minPt, maxPt));
  h_zmumunum_nrv.push_back(subDir.make<TH1F>("h_zmumunum_nrv", "found muons from Z decays vs nrv; number of reco  vertices; number of muons", nintVertcount, minVertcount, maxVertcount));
  h_zmumunum_minDist.push_back(subDir.make<TH1F>("h_zmumunum_minDist", "found muons from Z decays vs minimal distance; d_{min} / cm; number of muons", nintD, minD, maxDZ));
  h_zmumunum_quality.push_back(subDir.make<TH1F>("h_zmumunum_quality", "found muons from Z decays vs quality; quality class; number of muons", nintQ, minQ, maxQ));
  
  h_zmumunum_eta.back()->Sumw2();
  h_zmumunum_pt.back()->Sumw2();
  h_zmumunum_nrv.back()->Sumw2();
  h_zmumunum_minDist.back()->Sumw2();
  h_zmumunum_quality.back()->Sumw2();

  h_zmumuden_eta.push_back(subDir.make<TH1F>("h_zmumuden_eta", "triggered muons from Z decays vs eta; #eta; number of muons", nintEta, minEta, maxEta));
  h_zmumuden_pt.push_back(subDir.make<TH1F>("h_zmumuden_pt", "triggered muons from Z decays vs pt; p_{t} / GeV; number of muons", nintPt, minPt, maxPt));
  h_zmumuden_nrv.push_back(subDir.make<TH1F>("h_zmumuden_nrv", "triggered muons from Z decays vs nrv; number of reco  vertices; number of muons", nintVertcount, minVertcount, maxVertcount));
  h_zmumuden_minDist.push_back(subDir.make<TH1F>("h_zmumuden_minDist", "triggered muons from Z decays vs minimal distance; d_{min} / cm; number of muons", nintD, minD, maxDZ));
  h_zmumuden_quality.push_back(subDir.make<TH1F>("h_zmumuden_quality", "triggered muons from Z decays vs quality; quality class; number of muons", nintQ, minQ, maxQ));
  
  h_zmumuden_eta.back()->Sumw2();
  h_zmumuden_pt.back()->Sumw2();
  h_zmumuden_nrv.back()->Sumw2();
  h_zmumuden_minDist.back()->Sumw2();
  h_zmumuden_quality.back()->Sumw2();

  h_invariantMass_passed.push_back(subDir.make<TH1F>("h_invariantMass_passed", "invariant mass from the passed muon pairs; invariant mass / GeV; number of events", nintM, minM, maxM));
  h_invariantMass_failed.push_back(subDir.make<TH1F>("h_invariantMass_failed", "invariant mass from the failed muon pairs; invariant mass / GeV; number of events", nintM, minM, maxM));
  
  h_invariantMass_passed.back()->Sumw2();
  h_invariantMass_failed.back()->Sumw2();

  //profiles

  p_numTrks_nrv.push_back(subDir.make<TProfile>("p_numTrks_nrv", "number of considered signal tracks vs nrv; number of reco vertices; number of tracks", nintVertcount, minVertcount, maxVertcount));
  
  p_numTrks_nrv.back()->Sumw2();

  p_minDist_eta.push_back(subDir.make<TProfile>("p_minDist_eta", "minimum distance to signal vertex vs eta; eta; d_{min} / cm", nintEta, minEta, maxEta));
  p_minDist_pt.push_back(subDir.make<TProfile>("p_minDist_pt", "minimum distance to signal vertex vs pt; p_{t} / GeV; d_{min} / cm", nintPt, minPt, maxPt));
  p_minDist_nrv.push_back(subDir.make<TProfile>("p_minDist_nrv", "minimum distance to signal vertex vs nrv; number of reco vertices; d_{min} / cm", nintVertcount, minVertcount, maxVertcount));
  
  p_minDist_eta.back()->Sumw2();
  p_minDist_pt.back()->Sumw2();
  p_minDist_nrv.back()->Sumw2(); 

  p_quality_eta.push_back(subDir.make<TProfile>("p_quality_eta", "quality class vs eta; eta; quality class", nintEta, minEta, maxEta));
  p_quality_pt.push_back(subDir.make<TProfile>("p_quality_pt", "quality class vs pt; p_{t} / GeV; quality class", nintPt, minPt, maxPt));
  p_quality_nrv.push_back(subDir.make<TProfile>("p_quality_nrv", "quality class vs nrv; number of reco vertices; quality class", nintVertcount, minVertcount, maxVertcount));
  
  p_quality_eta.back()->Sumw2();
  p_quality_pt.back()->Sumw2();
  p_quality_nrv.back()->Sumw2(); 

  p_tweight_eta.push_back(subDir.make<TProfile>("p_tweight_eta", "track weight vs eta; eta; track weight", nintEta, minEta, maxEta));
  p_tweight_pt.push_back(subDir.make<TProfile>("p_tweight_pt", "track weight vs pt; p_{t} / GeV; track weight", nintPt, minPt, maxPt));
  p_tweight_nrv.push_back(subDir.make<TProfile>("p_tweight_nrv", "track weight vs nrv; number of reco vertices; track weight", nintVertcount, minVertcount, maxVertcount));
  p_tweight_quality.push_back(subDir.make<TProfile>("p_tweight_quality", "track weight vs quality; quality class; track weight", nintQ, minQ, maxQ));
  p_tweight_minDist.push_back(subDir.make<TProfile>("p_tweight_minDist", "track weight vs minimal distance; d_{min} / cm; track weight", nintD, minD, maxD));
  
  p_tweight_eta.back()->Sumw2();
  p_tweight_pt.back()->Sumw2();
  p_tweight_nrv.back()->Sumw2();
  p_tweight_quality.back()->Sumw2();
  p_tweight_minDist.back()->Sumw2(); 

  p_zmumueff_eta.push_back(subDir.make<TProfile>("p_zmumueff_eta", "efficiency of muons from Z Decays vs eta; eta; efficiency", nintEta, minEta, maxEta));
  p_zmumueff_pt.push_back(subDir.make<TProfile>("p_zmumueff_pt", "efficiency of muons from Z Decays vs pt; p_{t} / GeV; efficiency", nintPt, minPt, maxPt));
  p_zmumueff_nrv.push_back(subDir.make<TProfile>("p_zmumueff_nrv", "efficiency of muons from Z Decays vs nrv; number of reco vertices; efficiency", nintVertcount, minVertcount, maxVertcount));
  p_zmumueff_minDist.push_back(subDir.make<TProfile>("p_zmumueff_minDist", "efficiency of muons from Z Decays vs minimal distance; d_{min} / cm; efficiency", nintD, minD, maxDZ));
  p_zmumueff_quality.push_back(subDir.make<TProfile>("p_zmumueff_quality", "efficiency of muons from Z Decays vs quality; quality class; efficiency", nintQ, minQ, maxQ));
  
  p_zmumueff_eta.back()->Sumw2();
  p_zmumueff_pt.back()->Sumw2();
  p_zmumueff_nrv.back()->Sumw2();
  p_zmumueff_minDist.back()->Sumw2();
  p_zmumueff_quality.back()->Sumw2();

}

void 
TrackValidatorDataAlgos::get_input_collections(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //get the beamspot
  iEvent.getByLabel("offlineBeamSpot", bsH);

  //get the magnetic field
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldH);

}

void 
TrackValidatorDataAlgos::fill_independent_histos(int counter, int nrv, int rt, float w)
{

  h_numVtx.at(counter)->Fill(nrv, w);

  h_numTrks.at(counter)->Fill(rt, w);
  
  p_numTrks_nrv.at(counter)->Fill(nrv, rt, w);

}

void 
TrackValidatorDataAlgos::fill_track_dependend_histos(int counter, const TrackRef trk_ref, int q, const edm::EventSetup& iSetup, const VertexRef vtx_ref, int nrv, float w)
{

  double trk_eta = trk_ref->eta();
  double trk_pt = trk_ref->pt();

  TransientTrack transtrk(trk_ref, &(*bFieldH) );
  transtrk.setBeamSpot(*bsH);
  transtrk.setES(iSetup);
	    	    
  double min_distance = IPTools::absoluteImpactParameter3D( transtrk, *(vtx_ref ) ).second.value();
  double track_weight = vtx_ref->trackWeight(trk_ref);
  
  h_minDist.at(counter)->Fill(min_distance, w);
  h_quality.at(counter)->Fill(q, w);
  h_tweight.at(counter)->Fill(track_weight, w);
  
  h_numTrks_eta.at(counter)->Fill(trk_eta, w);
  h_numTrks_pt.at(counter)->Fill(trk_pt, w);  
  
  p_minDist_eta.at(counter)->Fill(trk_eta, min_distance, w);
  p_minDist_pt.at(counter)->Fill(trk_pt, min_distance, w);
  p_minDist_nrv.at(counter)->Fill(nrv, min_distance, w);  
  
  p_quality_eta.at(counter)->Fill(trk_eta, q, w);
  p_quality_pt.at(counter)->Fill(trk_pt, q, w);
  p_quality_nrv.at(counter)->Fill(nrv, q, w);
  
  p_tweight_eta.at(counter)->Fill(trk_eta, track_weight, w);
  p_tweight_pt.at(counter)->Fill(trk_pt, track_weight, w);
  p_tweight_nrv.at(counter)->Fill(nrv, track_weight, w);
  p_tweight_quality.at(counter)->Fill(q, track_weight, w);
  p_tweight_minDist.at(counter)->Fill(min_distance, track_weight, w);      

}

void 
TrackValidatorDataAlgos::fill_zmumu_histos(int counter, Handle< MuonCollection > refMuonsH, const TrackQualityPairVector tqpV, const edm::EventSetup& iSetup, const VertexRef vtx_ref, int nrv, float w)
{   

  if ( refMuonsH->size()!=2 ) {
    cout << "Not exactly two muons in the reference collection, which is bad!" << endl;
    return;
  }
  
  bool bothRecoed = true;  
  TLorentzVector sum_p4;
  
  for ( unsigned rmu_idx=0; rmu_idx<refMuonsH->size(); rmu_idx++ ) {
  
    MuonRef mu_ref(refMuonsH, rmu_idx);
    TrackRef mu_trkref = mu_ref->innerTrack();

    TransientTrack transtrk(mu_trkref, &(*bFieldH) );
    transtrk.setBeamSpot(*bsH);
    transtrk.setES(iSetup);
	    	    
    double min_distance = IPTools::absoluteImpactParameter3D( transtrk, *(vtx_ref ) ).second.value();
    
    sum_p4+= TLorentzVector( mu_ref->px(), mu_ref->py(), mu_ref->pz(), mu_ref->energy() );
  
    double rm_eta = mu_ref->eta();
    double rm_pt = mu_ref->pt();        
    
    int quality = find_muon_in_collection(mu_trkref, tqpV);
    bool recoed = (quality<10); 
  
    h_zmumuden_eta[counter]->Fill(rm_eta, w); 
    h_zmumuden_pt[counter]->Fill(rm_pt, w); 
    h_zmumuden_nrv[counter]->Fill(nrv, w);  
    h_zmumunum_minDist[counter]->Fill(min_distance, w);  
    h_zmumunum_quality[counter]->Fill(quality, w);
    
    if ( recoed ) {
  
      h_zmumunum_eta[counter]->Fill(rm_eta, w); 
      h_zmumunum_pt[counter]->Fill(rm_pt, w); 
      h_zmumunum_nrv[counter]->Fill(nrv, w);
      h_zmumuden_minDist[counter]->Fill(min_distance, w);
      h_zmumuden_quality[counter]->Fill(quality, w);
  
      p_zmumueff_eta[counter]->Fill(rm_eta, 1., w); 
      p_zmumueff_pt[counter]->Fill(rm_pt, 1., w); 
      p_zmumueff_nrv[counter]->Fill(nrv, 1., w); 
      p_zmumueff_minDist[counter]->Fill(min_distance, 1., w);
      p_zmumueff_quality[counter]->Fill(quality, 1., w);         
    
    } else {
    
      bothRecoed = false;
      
      p_zmumueff_eta[counter]->Fill(rm_eta, 0., w); 
      p_zmumueff_pt[counter]->Fill(rm_pt, 0., w); 
      p_zmumueff_nrv[counter]->Fill(nrv, 0., w); 
      p_zmumueff_minDist[counter]->Fill(min_distance, 0., w); 
      p_zmumueff_quality[counter]->Fill(quality, 0., w); 
      
    }
  
  }
  
  double invariantMass = sum_p4.M();
   
  if ( bothRecoed ) {
    h_invariantMass_passed[counter]->Fill(invariantMass, w);
  } else { 
    h_invariantMass_failed[counter]->Fill(invariantMass, w);
  }
  
}

int
TrackValidatorDataAlgos::find_muon_in_collection(const TrackRef muon, const TrackQualityPairVector tqpV)
{

  for (unsigned trk_idx=0; trk_idx<tqpV.size(); trk_idx++ ) {
  
    TrackRef trk_ref = tqpV[trk_idx].first;
    
    if ( muon->eta()  == trk_ref->eta() &&
	     muon->phi()  == trk_ref->phi() &&
	     muon->chi2() == trk_ref->chi2() &&
	     muon->ndof() == trk_ref->ndof() &&
	     muon->p()    == trk_ref->p() ) return tqpV[trk_idx].second;
  
  }

  return 15;
  
}