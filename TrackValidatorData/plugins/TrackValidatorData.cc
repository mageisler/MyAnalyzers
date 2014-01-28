#include "MGeisler/TrackValidatorData/interface/TrackValidatorData.h"

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToManyWithQuality.h"
#include "DataFormats/Common/interface/OneToManyWithQualityGeneric.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// constants, enums and typedefs
//

typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, int> > TrackToVertexAssMap;

//
// static data member definitions
//

//
// constructors and destructor
//
TrackValidatorData::TrackValidatorData(const edm::ParameterSet& iConfig):TrackValidatorDataAlgos(iConfig)
{
   //now do what ever initialization is needed
  refMCLabel_ = iConfig.getParameter<InputTag>("RefMuonCollection");
  
  amLabels_ = iConfig.getParameter<vector<InputTag> >("amLabels");

  vcLabel_ = iConfig.getParameter<InputTag>("vcLabel");

  wLabel_ = iConfig.getParameter<InputTag>("Weight");

  Service<TFileService> tfs;
  vector<TFileDirectory>* subDir(new vector<TFileDirectory>());

  for(unsigned aml=0; aml<amLabels_.size(); aml++){

    InputTag aMap = amLabels_[aml];
    string dirName = "";
    dirName += aMap.label();

    subDir->push_back(tfs->mkdir(dirName));
    TrackValidatorDataAlgos::BookHistos(subDir->at(aml));

  }


  isData_ = iConfig.getParameter<bool>("isData");
  isZMuMu_ = iConfig.getParameter<bool>("isZMuMu");
  ignoremissingtkcollection_ = iConfig.getParameter<bool>("ignoremissingtrackcollection");

}

TrackValidatorData::~TrackValidatorData()
{
}

//
// member functions
//

// ------------ method called for each event  ------------
void
TrackValidatorData::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  TrackValidatorDataAlgos::get_input_collections(iEvent, iSetup); 
  
  //get the vertex collection
  Handle< VertexCollection >  vcH ;
  iEvent.getByLabel(vcLabel_, vcH);
  
  int numRecoVtx = vcH->size();
  VertexRef fVR( vcH, 0 );
  
  float weight = 1.;
  
  if ( !isData_ ) {  

    //get reweighting weight from the event
    Handle<float>  wH;
    iEvent.getByLabel(wLabel_, wH);
    
    weight = *(wH.product());
  
  }

  //get reference track collection from the event
  Handle< MuonCollection >  rmcH;
  iEvent.getByLabel(refMCLabel_, rmcH);

  // ########################################################
  // part of the charged particle analysis
  // ########################################################

  //loop over input collections
  for(unsigned aml=0; aml<amLabels_.size(); aml++){ 

    //get track collection from the event
    Handle< TrackToVertexAssMap >  amH;
    iEvent.getByLabel(amLabels_[aml], amH);

	const TrackQualityPairVector trckcoll = amH->begin()->val;   
    

    // ###########################
    // fill independent histograms
    // ###########################

    TrackValidatorDataAlgos::fill_independent_histos(aml, numRecoVtx, trckcoll.size(), weight); 


    // #####################################################
    // fill reconstruction histograms (LOOP OVER RECOTRACKS)
    // #####################################################

    for ( unsigned int trckcoll_ite = 0; trckcoll_ite < trckcoll.size(); trckcoll_ite++ ){

      TrackRef trk_ref = trckcoll[trckcoll_ite].first;
      int quality = trckcoll[trckcoll_ite].second;

      TrackValidatorDataAlgos::fill_track_dependend_histos(aml, trk_ref, quality, iSetup, fVR, numRecoVtx, weight);

    }
    
    if ( !isZMuMu_ ) continue;


    // ##################################################
    // check if both muons have been correctly associated
    // ##################################################

    TrackValidatorDataAlgos::fill_zmumu_histos(aml, rmcH, trckcoll, iSetup, fVR, numRecoVtx, weight);
  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackValidatorData::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackValidatorData);
