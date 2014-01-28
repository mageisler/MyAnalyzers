// -*- C++ -*-
//
// Package:    TrackResolutionCheck
// Class:      TrackResolutionCheck
// 
/**\class TrackResolutionCheck TrackResolutionCheck.cc MGeisler/TrackResolutionCheck/src/TrackResolutionCheck.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Fri Nov 29 13:51:22 CET 2013
// $Id$
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"

#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
//
// class declaration
//

class TrackResolutionCheck : public edm::EDAnalyzer {
   public:
      explicit TrackResolutionCheck(const edm::ParameterSet&);
      ~TrackResolutionCheck();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      double get_gauss_error(TH1D*);
      std::vector< double > get_RMS_error(TH1D*);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();

      // ----------member data ---------------------------
  	   
      std::vector< double > pTintervalsv;
      
      TH1F* h_TracksNum_pt;
      
      TH1F* p_dZres_barrel;
      TH1F* p_dZres_transition;
      TH1F* p_dZres_endcap;
      
      TH2F* h2_dZres_barrel;
      TH2F* h2_dZres_transition;
      TH2F* h2_dZres_endcap;
      
      TProfile* p_Z0error_barrel;
      TProfile* p_Z0error_transition;
      TProfile* p_Z0error_endcap; 
      
      double eta_bar;
      double eta_tra;
      double eta_end;
      
      int nintPt;
      
  	  double minRes;  
      double maxRes;      
      int nintRes;
      
};

//
// constants, enums and typedefs
//

const double cm2um = 10000;

//
// static data member definitions
//

using namespace edm;
using namespace std;
using namespace reco;

//
// constructors and destructor
//
TrackResolutionCheck::TrackResolutionCheck(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
      
  minRes = -0.05;  
  maxRes = +0.05;      
  nintRes = 200;

  //parameters for pt axis
//   double minPt  = 0.;  
  double maxPt = 150.;
  double minPtLog = log10(0.1);  
  double maxPtLog = log10(maxPt); 
  nintPt = 30;
  
  eta_bar = 0.9;
  eta_tra = 1.4;
  eta_end = 2.5;
   
  Double_t pTintervalsa[31]; 
    
  double stepPt = (maxPtLog-minPtLog)/nintPt;
  pTintervalsa[0] = minPtLog;
  for (int k=1;k<nintPt+1;k++) {
    double d=pow(10,minPtLog+k*stepPt);
    pTintervalsa[k] = d;
  }    

  Service<TFileService> tfs;
  vector<TFileDirectory>* subDir(new vector<TFileDirectory>());
  subDir->push_back(tfs->mkdir("generalTracks"));

  h_TracksNum_pt = subDir->at(0).make<TH1F>( "h_TracksNum_pt", "; p_{T} / GeV; number of reco tracks", nintPt, pTintervalsa );
  h_TracksNum_pt->Sumw2();

  p_dZres_barrel = subDir->at(0).make<TH1F>( "p_dZres_barrel", "; p_{T} / GeV; resolution d_{Z} / #mum", nintPt, pTintervalsa );
  p_dZres_transition = subDir->at(0).make<TH1F>("p_dZres_transition", "; p_{T} / GeV; resolution d_{Z} / #mum", nintPt, pTintervalsa);
  p_dZres_endcap = subDir->at(0).make<TH1F>("p_dZres_endcap", "; p_{T} / GeV; resolution d_{Z} / #mum", nintPt, pTintervalsa);
  
  p_dZres_barrel->Sumw2();
  p_dZres_transition->Sumw2();
  p_dZres_endcap->Sumw2();

  h2_dZres_barrel = new TH2F( "h2_dZres_barrel", "", nintPt, pTintervalsa, nintRes,minRes,maxRes );
  h2_dZres_transition =  new TH2F("h2_dZres_transition", "", nintPt, pTintervalsa, nintRes,minRes,maxRes );
  h2_dZres_endcap =  new TH2F("h2_dZres_endcap", "", nintPt, pTintervalsa, nintRes,minRes,maxRes );
  
  h2_dZres_barrel->Sumw2();
  h2_dZres_transition->Sumw2();
  h2_dZres_endcap->Sumw2();

  p_Z0error_barrel = subDir->at(0).make<TProfile>("p_Z0error_barrel", "; p_{T} / GeV; #sigma_{Z_{0}} / #mum", nintPt, pTintervalsa);
  p_Z0error_transition = subDir->at(0).make<TProfile>("p_Z0error_transition", "; p_{T} / GeV; #sigma_{Z_{0}} / #mum", nintPt, pTintervalsa);
  p_Z0error_endcap = subDir->at(0).make<TProfile>("p_Z0error_endcap", "; p_{T} / GeV; #sigma_{Z_{0}} / #mum", nintPt, pTintervalsa);
  
  p_Z0error_barrel->Sumw2();
  p_Z0error_transition->Sumw2();
  p_Z0error_endcap->Sumw2(); 
  
}


TrackResolutionCheck::~TrackResolutionCheck()
{
 
   // do anyyhing here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackResolutionCheck::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get the beamspot
  Handle<BeamSpot>  bsH ;
  iEvent.getByLabel("offlineBeamSpot", bsH);

  //get the tracking particles   
  Handle<TrackingParticleCollection>  TPCollectionH ;
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionH);

  //get the reco tracks   
  Handle< TrackCollection>  tcH ;
  Handle< View<Track> >  tvH ;
  iEvent.getByLabel("generalTracks",tcH);
  iEvent.getByLabel("generalTracks",tvH);
      

  //associate reco tracks to tracking particles
  ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get("quickTrackAssociatorByHits",theAssociator);
  TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product();
    
  RecoToSimCollection rsC = theTrackAssociator_->associateRecoToSim(tvH,TPCollectionH,&iEvent,&iSetup);
  
  ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
  iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP",parametersDefinerTP);


  for ( unsigned int tc_ite=0; tc_ite<tcH->size(); ++tc_ite ) {

    TrackRef trk_ref(tcH, tc_ite);
    TrackBaseRef track_bref(trk_ref);
    
    double trk_eta = fabs(trk_ref->eta());
    double trk_pt = trk_ref->pt();
    
	h_TracksNum_pt->Fill( trk_pt );
    
    double trk_z0Err = trk_ref->dzError();
    
	if ( trk_eta<= eta_bar ) {
	  p_Z0error_barrel->Fill(trk_pt, trk_z0Err );
	} else {
	  if ( trk_eta<= eta_tra ) {
		p_Z0error_transition->Fill(trk_pt, trk_z0Err );
	  } else {
		p_Z0error_endcap->Fill(trk_pt, trk_z0Err );
	  }
	}
          
    vector<pair<TrackingParticleRef, double> > tpd;
    if(rsC.find(track_bref) != rsC.end()) tpd = rsC[track_bref];

    if (tpd.size()!=0) {
      
      TrackingParticleRef tpr = tpd.begin()->first;      
      
	  TrackingParticle::Vector tpm = parametersDefinerTP->momentum(iEvent,iSetup, *(tpr.get()) );
	  TrackingParticle::Point  tpv = parametersDefinerTP->vertex(iEvent,iSetup, *(tpr.get()) );
      
      double dz_rec = trk_ref->dz( bsH->position() );
      double dz_sim = tpv.z() - (tpv.x()*tpm.x()+tpv.y()*tpm.y())/sqrt(tpm.perp2()) * tpm.z()/sqrt(tpm.perp2());
      
      double dist = dz_rec-dz_sim;
      
  	  cout << dz_rec << " - " << dz_sim << endl;
  	  cout << dist << endl;
  	  cout << "" << endl;
      
	  if ( trk_eta<= eta_bar ) {
 		h2_dZres_barrel->Fill( trk_pt, dist );
	  } else {
		if ( trk_eta<= eta_tra ) {
		  h2_dZres_transition->Fill( trk_pt, dist );
		} else {
		  h2_dZres_endcap->Fill( trk_pt, dist );
		}
	  }
      
        
    }      
  
  } 

}

vector< double >
TrackResolutionCheck::get_RMS_error(TH1D* in_histo)
{

  vector< double > output;
  
  output.push_back( in_histo->GetRMS() );
  output.push_back( in_histo->GetRMSError() );
  
  return output;

}

double
TrackResolutionCheck::get_gauss_error(TH1D* in_histo)
{

  double result = 0.;

  TAxis *axis0 = in_histo->GetXaxis();
  int nbin = axis0->GetLast();
  double nOF = in_histo->GetBinContent(nbin+1);
  double nUF = in_histo->GetBinContent(0);
  
  double fitRange = 2.*in_histo->GetRMS();
  double sigMax[2] = {0.,0.};
  
//   cout << "Jetzt: " << endl;
//   cout << nOF<< endl;
//   cout << nUF << endl;
//   cout << in_histo->GetEntries()<< endl;
//   cout << "" << endl;

  if ( in_histo->GetEntries() - nOF - nUF >= 10) {
    for (int bi = 0; bi < nbin; bi++) {
      if ( (axis0->GetBinLowEdge(bi) < 0) && (in_histo->GetBinContent(bi) > 0) ) {
	    sigMax[0] = axis0->GetBinLowEdge(bi);
	    if ( abs(sigMax[0]) > abs(sigMax[1]) ) sigMax[1] = abs(sigMax[0]);
      }
      if ( (axis0->GetBinLowEdge(bi) >= 0) && (in_histo->GetBinContent(bi) > 0) ) {
	    sigMax[0] = axis0->GetBinUpEdge(bi);
	    if ( abs(sigMax[0]) > abs(sigMax[1]) ) sigMax[1] = abs(sigMax[0]);
      }
    }

    TF1 *fgaus = new TF1("fgaus","gaus",-fitRange, fitRange);
    fgaus->SetParameter(1, 0.);
    fgaus->SetParLimits(1, -fitRange/10., fitRange/10.);
    fgaus->SetParLimits(2, 0., sqrt(2.)*sigMax[1]);
    in_histo->Fit(fgaus,"QLRM");
    result = ((fgaus->GetParameter(2))?fgaus->GetParameter(2):0.);
  }
  
//   cout << result << endl;
  
  return result / sqrt(2.);

}

// ------------ method called at the end of the job  ------------
void
TrackResolutionCheck::endJob()
{
  
  cout << h2_dZres_barrel->Integral() << endl;
  
  for ( int bn=1; bn<=nintPt ; bn++ ) {
      p_dZres_barrel->Fill(bn, get_gauss_error(h2_dZres_barrel->ProjectionY("pSrn",bn,bn,"e") ) ); 
      p_dZres_transition->Fill(bn, get_gauss_error(h2_dZres_transition->ProjectionY("pSrn",bn,bn,"e") ) ); 
      p_dZres_endcap->Fill(bn, get_gauss_error(h2_dZres_endcap->ProjectionY("pSrn",bn,bn,"e") ) ); 
  } 
  
  cout << "" << endl;


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackResolutionCheck::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackResolutionCheck);
