#include "MGeisler/IsolationValidator/interface/IsolationValidator.h"

// system include files
#include <string>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// constants, enums and typedefs
//

typedef ValueMap<double> IsoMap;

//
// static data member definitions
//

//
// constructors and destructor
//
IsolationValidator::IsolationValidator(const edm::ParameterSet& iConfig)
{
  
  //parameters for Pileup plots
  minVertcount  = -0.5;
  maxVertcount  = 59.5;
  nintVertcount = 60;
  
  //parameters for Iso plots
  minIso  = 0.;
  maxIso  = 3.;
  nintIso = 60;
  
  
  //now do what ever initialization is needed
  sufLabels_ = iConfig.getParameter<vector<InputTag> >("Suffixes");
  
  wLabel_ = iConfig.getParameter<InputTag>("Weight");
  
  isData_ = iConfig.getParameter<bool>("isData");

  VertexCollectionLabel_ = iConfig.getParameter<InputTag>("VertexCollection");

  Service<TFileService> tfs;
  vector<TFileDirectory>* subDir(new vector<TFileDirectory>());

  for ( unsigned il_idx=0; il_idx<sufLabels_.size(); il_idx++ ){

    InputTag itag = sufLabels_[il_idx];
    string dirName = "Isolation";
    dirName += itag.label();

    subDir->push_back(tfs->mkdir(dirName));
     
    //Book histograms

    h_numVtx.push_back(subDir->at(il_idx).make<TH1F>("h_numVtx", "Number of reconstructed vertices per event", nintVertcount, minVertcount, maxVertcount));
  
    h_numVtx.back()->Sumw2();


    h_phIso_numVtx.push_back(subDir->at(il_idx).make<TH2F>("h_phIso_numVtx", "Photon isolation vs. number of reconstructed vertices per event", nintIso, minIso, maxIso, nintVertcount, minVertcount, maxVertcount));
    h_elIso_numVtx.push_back(subDir->at(il_idx).make<TH2F>("h_elIso_numVtx", "Electron isolation vs. number of reconstructed vertices per event", nintIso, minIso, maxIso, nintVertcount, minVertcount, maxVertcount));
    h_muIso_numVtx.push_back(subDir->at(il_idx).make<TH2F>("h_muIso_numVtx", "Muon isolation vs. number of reconstructed vertices per event", nintIso, minIso, maxIso, nintVertcount, minVertcount, maxVertcount));
  
    h_phIso_numVtx.back()->Sumw2();
    h_elIso_numVtx.back()->Sumw2();
    h_muIso_numVtx.back()->Sumw2();

  }

}


IsolationValidator::~IsolationValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
IsolationValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
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
  
  }

  //loop over input reco met
  for ( unsigned il_idx=0; il_idx<sufLabels_.size(); il_idx++ ){    

    h_numVtx.at(il_idx)->Fill(npu, weight);
    
    string suffix = sufLabels_[il_idx].label();
    
    // photons
    
    string photonsL = "pfSelectedPhotons" + suffix;
    Handle< PFCandidateCollection > photonsH;
    iEvent.getByLabel( photonsL, photonsH);
    
    string phIVChargedL = "phPFIsoValueCharged04PFId" + suffix;
    Handle< IsoMap > phIVChargedH;
    iEvent.getByLabel( phIVChargedL, phIVChargedH);
    
    string phIVGammaL = "phPFIsoValueGamma04PFId" + suffix;
    string phIVNeutralL = "phPFIsoValueNeutral04PFId" + suffix;
    vector< Handle< IsoMap > > phIVNeutralH( 2 );
    iEvent.getByLabel( phIVGammaL, phIVNeutralH[0] );
    iEvent.getByLabel( phIVNeutralL, phIVNeutralH[1] );
    
    for ( unsigned ph_idx=0; ph_idx<photonsH->size(); ph_idx++ ){
     
      PFCandidateRef ph_ref(photonsH, ph_idx);
        
      const IsoMap & isoMapC = *(phIVChargedH);  
      double isoSumCharged = isoMapC[ph_ref];
      
      double isoSumNeutral=0.0;
      
      for(unsigned iMap = 0; iMap<phIVNeutralH.size(); ++iMap) {
        
        const IsoMap & isoMapN = *(phIVNeutralH[iMap]);
        double val = isoMapN[ph_ref];
        isoSumNeutral+=val;
        
      }
      
      double ph_iso= ( isoSumCharged+isoSumNeutral ) / ph_ref->pt();
      
      h_phIso_numVtx.at(il_idx)->Fill(ph_iso, npu, weight);
     
    } 
    
    // electrons    
    
    string electronsL = "pfSelectedElectrons" + suffix;  
    Handle< PFCandidateCollection > electronsH; 
    iEvent.getByLabel( electronsL, electronsH);
    
    string elIVChargedL = "elPFIsoValueCharged03PFId" + suffix;
    Handle< IsoMap > elIVChargedH;
    iEvent.getByLabel( elIVChargedL, elIVChargedH);
    
    string elIVGammaL = "elPFIsoValueGamma03PFId" + suffix;
    string elIVNeutralL = "elPFIsoValueNeutral03PFId" + suffix;
    vector< Handle< IsoMap > > elIVNeutralH( 2 );
    iEvent.getByLabel( elIVGammaL, elIVNeutralH[0] );
    iEvent.getByLabel( elIVNeutralL, elIVNeutralH[1] );
    
    for ( unsigned el_idx=0; el_idx<electronsH->size(); el_idx++ ){
     
      PFCandidateRef el_ref(electronsH, el_idx);
        
      const IsoMap & isoMapC = *(elIVChargedH);  
      double isoSumCharged = isoMapC[el_ref];
      
      double isoSumNeutral=0.0;
      
      for(unsigned iMap = 0; iMap<elIVNeutralH.size(); ++iMap) {
        
        const IsoMap & isoMapN = *(elIVNeutralH[iMap]);
        double val = isoMapN[el_ref];
        isoSumNeutral+=val;
        
      }
      
      double el_iso = ( isoSumCharged+isoSumNeutral ) / el_ref->pt();
      
      h_elIso_numVtx.at(il_idx)->Fill(el_iso, npu, weight);
     
    } 
    
    // muons    
    
    string muonsL = "pfSelectedMuons" + suffix;    
    Handle< PFCandidateCollection > muonsH;     
    iEvent.getByLabel( muonsL, muonsH);
    
    string muIVChargedL = "muPFIsoValueCharged04" + suffix;
    Handle< IsoMap > muIVChargedH;
    iEvent.getByLabel( muIVChargedL, muIVChargedH);
    
    string muIVGammaL = "muPFIsoValueGamma04" + suffix;
    string muIVNeutralL = "muPFIsoValueNeutral04" + suffix;
    vector< Handle< IsoMap > > muIVNeutralH( 2 );
    iEvent.getByLabel( muIVGammaL, muIVNeutralH[0] );
    iEvent.getByLabel( muIVNeutralL, muIVNeutralH[1] );
    
    for ( unsigned mu_idx=0; mu_idx<muonsH->size(); mu_idx++ ){
     
      PFCandidateRef mu_ref(muonsH, mu_idx);
        
      const IsoMap & isoMapC = *(muIVChargedH);  
      double isoSumCharged = isoMapC[mu_ref];
      
      double isoSumNeutral=0.0;
      
      for(unsigned iMap = 0; iMap<muIVNeutralH.size(); ++iMap) {
        
        const IsoMap & isoMapN = *(muIVNeutralH[iMap]);
        double val = isoMapN[mu_ref];
        isoSumNeutral+=val;
        
      }
      
      double mu_iso = ( isoSumCharged+isoSumNeutral ) / mu_ref->pt();
      
      h_muIso_numVtx.at(il_idx)->Fill(mu_iso, npu, weight);
     
    } 
    
  }
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
IsolationValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsolationValidator);
