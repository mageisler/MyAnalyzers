#include "MGeisler/JetValidator/interface/JetValidator.h"

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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JetValidator::JetValidator(const edm::ParameterSet& iConfig):JetValidatorAlgos(iConfig)
{
   //now do what ever initialization is needed

  genJetLabel_ = iConfig.getParameter<InputTag>("genJetLabel");
  recoJetLabels_ = iConfig.getParameter<vector<InputTag> >("recoJetLabels");

  puLabel_ = iConfig.getParameter<InputTag>("PULabel");

  Service<TFileService> tfs;
  vector<TFileDirectory>* subDir(new vector<TFileDirectory>());

  JetValidatorAlgos::CreateIntervalVectors();

  for(unsigned jcl=0; jcl<recoJetLabels_.size(); jcl++){

    JetValidatorAlgos::initialize();

    InputTag jColl = recoJetLabels_[jcl];
    string dirName = "";
    dirName += jColl.label();

    subDir->push_back(tfs->mkdir(dirName));
    JetValidatorAlgos::BookHistos(subDir->at(jcl));

  }

}


JetValidator::~JetValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JetValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  //get gen jet collection from the event
  Handle<GenJetCollection> genJetCollH;
  iEvent.getByLabel(genJetLabel_,genJetCollH);

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

  int npu = puinfo.getPU_NumInteractions();

  //loop over input reco jet collections
  for (unsigned jcl=0; jcl<recoJetLabels_.size(); jcl++){

    //get reco jet collection from the event
    Handle<JetCollection> recoJetCollH;
    iEvent.getByLabel(recoJetLabels_[jcl],recoJetCollH);

    // ############################
    // fill independent histograms
    // ############################

    JetValidatorAlgos::fill_independent_histos(jcl,npu,recoJetCollH->size(),genJetCollH->size()); 

    // ###############################################
    // fill simulation histograms (LOOP OVER GENJETS)
    // ###############################################

    for (unsigned gj_ite=0; gj_ite<genJetCollH->size(); gj_ite++){

      GenJetRef gjr(genJetCollH, gj_ite);

      JetRef matchedRecoJet = JetValidatorAlgos::get_best_matching_recoJet(gjr,recoJetCollH);

      JetValidatorAlgos::fill_recoAssociated_genJet_histos(jcl,gjr,matchedRecoJet,npu,gj_ite);

    }

    // ####################################################
    // fill reconstruction histograms (LOOP OVER RECOJETS)
    // ####################################################

    for (unsigned rj_ite=0; rj_ite<recoJetCollH->size(); rj_ite++){

      JetRef rjr(recoJetCollH, rj_ite);

      GenJetRef matchedGenJet = JetValidatorAlgos::get_best_matching_genJet(rjr,genJetCollH);

      JetValidatorAlgos::fill_genAssociated_recoJet_histos(jcl,rjr,matchedGenJet,npu);

    }

  }

}



// ------------ method called when ending the processing of a run  ------------
void 
JetValidator::endRun(edm::Run const&, edm::EventSetup const&)
{

  //loop over input jet collections
  for(unsigned jcl=0; jcl<recoJetLabels_.size(); jcl++){

    JetValidatorAlgos::fillHistosFromVectors(jcl); 
 
  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetValidator);