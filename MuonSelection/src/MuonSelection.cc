
#include "MGeisler/MuonSelection/interface/MuonSelection.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace edm;
using namespace std;
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
MuonSelection::MuonSelection(const edm::ParameterSet& iConfig)
{

  //register your products
  produces< MuonCollection >();
   
  //now do what ever other initialization is needed
  muonLabel_ = iConfig.getParameter<InputTag>("Muons");
  
  vertexCollectionLabel_ = iConfig.getParameter<InputTag>("VertexCollection");
  
  minLayers_ = iConfig.getParameter<int>("minLayers");
  
  minPt_ = iConfig.getParameter<double>("minPt");
  
  maxEta_ = iConfig.getParameter<double>("maxEta");
  
}


MuonSelection::~MuonSelection()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonSelection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  auto_ptr<MuonCollection> muon_output(new MuonCollection);

  //get the input muon collection
  Handle< MuonCollection > mCH;
  iEvent.getByLabel( muonLabel_, mCH );

  //get the vertex collection
  Handle< VertexCollection > vCH;
  iEvent.getByLabel( vertexCollectionLabel_, vCH );
  
  VertexRef fV_ref(vCH, 0);
  
  for ( unsigned mu_idx=0; mu_idx<mCH->size(); mu_idx++ ) {
  
    MuonRef mu_ref(mCH, mu_idx);
    
    if ( ( muon::isTightMuon(*mu_ref, *fV_ref) ) &&
         ( mu_ref->innerTrack()->hitPattern().trackerLayersWithMeasurement() >= minLayers_ ) &&
         ( mu_ref->pt() >= minPt_ ) &&
         ( fabs( mu_ref->eta() ) <= maxEta_ ) ) {
         
      muon_output->push_back( *mu_ref );
      
    }        
  
  }
  
  iEvent.put( muon_output ); 
 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonSelection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonSelection);
