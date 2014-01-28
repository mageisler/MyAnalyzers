#include "MGeisler/METValidator/interface/METValidator.h"

// user include files

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/CorrMETData.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "JetMETCorrections/Type1MET/interface/CorrectedMETProducerT.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// root includes
#include "TMath.h"


using namespace std;
using namespace edm;
using namespace pat;
using namespace CorrectedMETProducer_namespace;
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
METValidator::METValidator(const edm::ParameterSet& iConfig)
{
  
  //parameters for Pileup plots
  minVertcount  = -0.5;
  maxVertcount  = 59.5;
  nintVertcount = 60;
  
  //parameters for MET plots
  minMET  = -0.5;
  maxMET  = 99.5;
  nintMET = 200;

  //now do what ever initialization is needed

  corrMETLabels_ = iConfig.getParameter<vector<InputTag> >("corrMETLabels");
  rawMETLabels_ = iConfig.getParameter<vector<InputTag> >("rawMETLabels");

  isData_ = iConfig.getParameter<bool>("isData");

  usePUInfo_ = iConfig.getParameter<bool>("usePUInfo");

  puLabel_ = iConfig.getParameter<InputTag>("PULabel");

  wLabel_ = iConfig.getParameter<InputTag>("Weight");

  VertexCollectionLabel_ = iConfig.getParameter<InputTag>("VertexCollection");

  Service<TFileService> tfs;
  vector<TFileDirectory>* subDir(new vector<TFileDirectory>());

  for(unsigned mtl=0; mtl<corrMETLabels_.size(); mtl++){

    InputTag mtag = corrMETLabels_[mtl];
    string dirName = "";
    dirName += mtag.label();

    subDir->push_back(tfs->mkdir(dirName));
     
    //Book met histograms

    h_numVtx.push_back(subDir->at(mtl).make<TH1F>("h_numVtx", "Number of reconstructed vertices per event", nintVertcount, minVertcount, maxVertcount));
  
    h_numVtx.back()->Sumw2();

    h_corrMET_Type0.push_back(subDir->at(mtl).make<TH2F>("h_corrMET_Type0", "MET distribution vs number of generated vertices", nintMET, minMET, maxMET, nintVertcount, minVertcount, maxVertcount));
  
    h_corrMET_Type0.back()->Sumw2();
    
    h_rawMET_Type0.push_back(subDir->at(mtl).make<TH2F>("h_rawMET_Type0", "MET distribution vs number of generated vertices", nintMET, minMET, maxMET, nintVertcount, minVertcount, maxVertcount));
  
    h_rawMET_Type0.back()->Sumw2();

  }

}


METValidator::~METValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
METValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if ( corrMETLabels_.size()!=rawMETLabels_.size() ) return;
  
  //get the vertex collection
  Handle< VertexCollection >  vcH ;
  iEvent.getByLabel(VertexCollectionLabel_, vcH);
  
  int npu = vcH->size();
      
  float weight = 1.;
  
  if ( !isData_ ) {  

    //get reweighting weight from the event
    Handle<float>  wH;
    iEvent.getByLabel(wLabel_, wH);
    
    weight = *(wH.product());
  
    if ( usePUInfo_ ) {

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

      npu = puinfo.getPU_NumInteractions();
      
    }
  
  }

  //loop over input reco met
  for (unsigned mtl=0; mtl<corrMETLabels_.size(); mtl++){

    h_numVtx.at(mtl)->Fill(npu, weight);

    //get raw met from the event
    edm::Handle< std::vector<pat::MET> > rawMETH;
    iEvent.getByLabel(rawMETLabels_[mtl], rawMETH);
    pat::MET rawMET = rawMETH->at(0);
    
    double met_raw = rawMET.et();

    //get corrected met from the event
    edm::Handle< std::vector<pat::MET> > corrMETH;
    iEvent.getByLabel(corrMETLabels_[mtl], corrMETH);
    pat::MET corrMET = corrMETH->at(0);
    
    double met_corr = corrMET.et();

    h_rawMET_Type0[mtl]->Fill(met_raw, npu, weight);
    h_corrMET_Type0[mtl]->Fill(met_corr, npu, weight);

  }


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
METValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(METValidator);