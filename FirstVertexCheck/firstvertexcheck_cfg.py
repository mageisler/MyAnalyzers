import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
	
tag =  'GR_P_V42_AN4::All'
Outfile = "FirstVertexCheck_MC-QCD.root"
tag = 'START53_V18PR::All'
#Infile = ["file:/user/geisler/RelValZmumuJets_Pt_20_300_CMSSW_5_3_4_cand1_PU_START53_GEN-SIM-RECO.root"]
Infile = ['/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/001288D3-2AE3-E111-A943-0030487D5E5F.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/002D1EF6-33E3-E111-96E1-0030487E54B7.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0043E78A-56E3-E111-8743-002481E76052.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/006D0D9B-91E3-E111-8F8B-0025904B1342.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0086EC5D-72E3-E111-98B3-003048F0E3D2.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/008AE9FF-92E3-E111-A402-002481E0EA70.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/009B6F46-60E3-E111-B90C-003048D43730.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/00A80552-8EE3-E111-9F11-0025904B12A8.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/00B5BF5E-29E3-E111-B714-0030487E4EC7.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/00B5F5BC-8FE3-E111-B862-00266CF32EAC.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/00C6DDDF-50E3-E111-8A00-0030487F92A7.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/00D78FF9-55E3-E111-952F-003048D4DCD8.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/00E7C790-5DE3-E111-91CD-0030487E4B8D.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/02071AEA-80E3-E111-9B14-00266CF2AACC.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/021E4722-91E3-E111-B64A-003048D437C4.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0257F3A9-38E3-E111-83C2-003048C69328.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/027606E5-8EE3-E111-874B-002481E101DC.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/028242AC-57E3-E111-8F00-003048D3C880.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/028461FD-49E3-E111-BA3B-003048C69310.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/02AE2566-61E3-E111-AFC1-0030487D5D5B.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/02B1ADAF-96E3-E111-94A3-00266CF32920.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/02BE8D1B-3EE3-E111-8336-0030487F1655.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/02F3ED75-35E3-E111-A47D-003048C693C8.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0409D946-93E3-E111-A6B2-00266CF33340.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0418ADB9-38E3-E111-9A1C-0030487D857D.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0420C508-3DE3-E111-BF2A-003048C692CA.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/04597A88-88E3-E111-86B5-00266CFFA2EC.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/045E8E57-3CE3-E111-BC55-003048C692DE.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/046563E3-97E3-E111-ACE1-00266CFFA124.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/046ED408-99E3-E111-B547-003048F0E186.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/046F5762-7FE3-E111-90A5-00266CF32E70.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0477032B-99E3-E111-AC44-002481E0D6EE.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0479907F-7EE3-E111-8C5A-0025901D42C0.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0486042A-26E3-E111-A803-0030487F910D.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/049B820F-99E3-E111-9A82-003048F0E81E.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/04B2F1DC-27E3-E111-A82B-003048C692A4.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/04C8EF53-29E3-E111-B10B-003048C68A9A.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/04CCFF09-8CE3-E111-98FB-00266CFFA678.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/04DCDD88-74E3-E111-A252-00266CFFA124.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/04FE7D8B-2DE3-E111-9280-0030487D5EA7.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/06032774-53E3-E111-8F63-0030487D8633.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/060783FE-8BE3-E111-8F2C-0030487D5DA9.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/0653D2B8-8BE3-E111-A577-0030487CF3F7.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/06567F36-5EE3-E111-BC27-003048D4DEDC.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/06A2F24D-91E3-E111-B58F-003048D4365C.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/06E9E82A-68E3-E111-B926-0030487E54B5.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/087BAEBF-38E3-E111-8403-0030487F91D9.root',
       '/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/088BCD34-3AE3-E111-9B6B-0030487F929D.root']
    
print "Oufile is set to", Outfile, ", with tag", tag

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = tag

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring()
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(Outfile),
    closeFileFast = cms.untracked.bool(True),
)

process.demo = cms.EDAnalyzer('FirstVertexCheck'
)


process.p = cms.Path(process.demo)
		
process.source.fileNames = Infile
