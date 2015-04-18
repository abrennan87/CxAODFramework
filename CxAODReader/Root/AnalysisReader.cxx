#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <CxAODReader/AnalysisReader.h>

#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODBTagging/BTagging.h" 

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
using namespace std;


// this is needed to distribute the algorithm to the workers
ClassImp(AnalysisReader)



AnalysisReader :: AnalysisReader () :
m_analysisType(-1),
    m_eventCounter(0),
    m_isMC(false),
    m_weight(1.),
    m_sumOfWeights(1.),
    m_isSherpaVJets(0),
    m_isSherpaPt0VJets(0),
    m_SherpaPt0VJetsCut(70000)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode AnalysisReader :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  job.useXAOD ();

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init( "MyxAODAnalysis" ).ignore(); // call before opening first file


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode AnalysisReader :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  // In absence of hist manager suggest m_hist_Xlep_histname for pointers
  // name on file can be the same

  Info("histInitialize()", "Analysis type %i", m_analysisType ); 

  if (m_analysisType==e0Lep) 
    histInitialize_0Lep();
  else if (m_analysisType==e1Lep) 
    histInitialize_1Lep();
  else if (m_analysisType==e2Lep) 
    histInitialize_2Lep();
  else if (m_analysisType==VBFa)
    histInitialize_VBFa();
  else if (m_analysisType==monoWZH)
	histInitialize_monoWZH();
  else {
    Warning("histInitialize()", "Invalid analysis type %i",m_analysisType);
    return EL::StatusCode::FAILURE;
  }

  // some histograms useful for all analyses
  // Truth pt
  m_hist_VPtTruth = new TH1F("VptTruth","V pt truth; V pt truth [GeV]; Events",100,0,1000);
  m_hist_VPtTruth->Sumw2();
  wk()->addOutput(m_hist_VPtTruth);
      
  // b-tagging flavour plots
  bTaggerHists(0); // 0 for booking histograms
  return EL::StatusCode::SUCCESS;

}

EL::StatusCode AnalysisReader :: bTaggerHists(int bookfill, const xAOD::Jet *jet)
{
  // bookfill 0/1 for booking/filling
  // fills histograms per jet and jet flavour for a list of b-taggers

  if(bookfill==0) {
    TString taggers[6]={"SV1_IP3D","MV1","MV2c00","MV2c10","MV2c20","MVb"};
    TString flavs[3]={" L"," C"," B"};
    
    for (auto &tagger : taggers) {
      TString hbasename(tagger);
      float xmin=-1.1; float xmax=1.1; 
      if (hbasename.Contains("SV1_IP3D")) {xmin=-30; xmax=80;}
      for (auto &flav : flavs) {
	TString hname(hbasename); 
	hname+=flav;
	m_hist_taggers[hname] = new TH1F(TString("Taggers/"+hname),TString(";"+hname+";Events"),110,xmin,xmax);
	m_hist_taggers[hname]->Sumw2();
	wk()->addOutput(m_hist_taggers[hname]);
      }
    }

  } else if (bookfill==1) {

    if (!jet) return EL::StatusCode::SUCCESS;
    
    int label(1000);
    jet->getAttribute<int>("TruthLabelID", label);
    //    std::cout << " label " << label << std::endl;
    float sv1_plus_IP3D=m_superDecorator.get(jet,JetFloatProps::SV1_IP3D);
    float mv1=m_superDecorator.get(jet,JetFloatProps::MV1);
    float mv2c00=m_superDecorator.get(jet,JetFloatProps::MV2c00);
    float mv2c10=m_superDecorator.get(jet,JetFloatProps::MV2c10);
    float mv2c20=m_superDecorator.get(jet,JetFloatProps::MV2c20);
    float mvb=m_superDecorator.get(jet,JetFloatProps::MVb);
    
    TString flav(" L");
    if (label==4) flav=" C";
    else if (label==5) flav=" B";
    m_hist_taggers["SV1_IP3D"+flav]->Fill(sv1_plus_IP3D,m_weight);
    m_hist_taggers["MV1"+flav]->Fill(mv1,m_weight);
    m_hist_taggers["MV2c00"+flav]->Fill(mv2c00,m_weight);
    m_hist_taggers["MV2c10"+flav]->Fill(mv2c10,m_weight);
    m_hist_taggers["MV2c20"+flav]->Fill(mv2c20,m_weight);
    m_hist_taggers["MVb"+flav]->Fill(mvb,m_weight);
    
  }

  return EL::StatusCode::SUCCESS;

}


EL::StatusCode AnalysisReader :: histInitialize_0Lep()
{
  //m_hist_0lep_
  //m_hist_0lep_MET = new TH1F("MET","MET",100,0,600);
  //m_hist_olep_MET->Sumw2();
  //wk()->addOutput(m_hist_0lep_MET);

  TH1::SetDefaultSumw2();
  TString buffer;
  int nbins;
  double min;
  double max;
  TString prefix = "2j_";

  buffer = prefix+"dRBB"; nbins = 100; min = 0.; max = 5.;
  m_hist2j_0lep_dRBB = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist2j_0lep_dRBB);
  //
  buffer = prefix+"dEtaBB"; nbins = 100; min = 0.; max = 5.;
  m_hist2j_0lep_dEtaBB = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist2j_0lep_dEtaBB);
  //
  buffer = prefix+"dPhiVBB"; nbins = 100; min = 0.; max = 3.15;
  m_hist2j_0lep_dPhiVBB = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist2j_0lep_dPhiVBB);
  //
  buffer = prefix+"mBB"; nbins = 100; min = 0.; max = 500.;
  m_hist2j_0lep_mBB = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist2j_0lep_mBB);
  //
  buffer = prefix+"MEff"; nbins = 100; min = 0.; max = 1000.;
  m_hist2j_0lep_MEff = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist2j_0lep_MEff);
  //
  buffer = prefix+"pTB1"; nbins = 100; min = 0.; max = 500.;
  m_hist2j_0lep_pTB1 = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist2j_0lep_pTB1);
  //
  buffer = prefix+"pTB2"; nbins = 100; min = 0.; max = 500.;
  m_hist2j_0lep_pTB2 = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist2j_0lep_pTB2);
  //
  buffer = prefix+"MET"; nbins = 100; min = 0.; max = 500.;
  m_hist2j_0lep_MET = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist2j_0lep_MET);

  prefix = "3j_";
  
  buffer = prefix+"dRBB"; nbins = 100; min = 0.; max = 5.;
  m_hist3j_0lep_dRBB = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_dRBB);
  //
  buffer = prefix+"dEtaBB"; nbins = 100; min = 0.; max = 5.;
  m_hist3j_0lep_dEtaBB = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_dEtaBB);
  //
  buffer = prefix+"dPhiVBB"; nbins = 100; min = 0.; max = 3.15;
  m_hist3j_0lep_dPhiVBB = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_dPhiVBB);
  //
  buffer = prefix+"mBB"; nbins = 100; min = 0.; max = 500.;
  m_hist3j_0lep_mBB = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_mBB);
  //
  buffer = prefix+"MEff3"; nbins = 100; min = 0.; max = 1000.;
  m_hist3j_0lep_MEff3 = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_MEff3);
  //
  buffer = prefix+"pTB1"; nbins = 100; min = 0.; max = 500.;
  m_hist3j_0lep_pTB1 = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_pTB1);
  //
  buffer = prefix+"pTB2"; nbins = 100; min = 0.; max = 500.;
  m_hist3j_0lep_pTB2 = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_pTB2);
  //
  buffer = prefix+"MET"; nbins = 100; min = 0.; max = 500.;
  m_hist3j_0lep_MET = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_MET);
  //
  buffer = prefix+"MBBJ"; nbins = 100; min = 0.; max = 750.;
  m_hist3j_0lep_MBBJ = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_MBBJ);
  //
  buffer = prefix+"pTJ3"; nbins = 100; min = 0.; max = 150.;
  m_hist3j_0lep_pTJ3 = new TH1F(TString(buffer), TString(buffer), nbins, min, max);
  wk()->addOutput(m_hist3j_0lep_pTJ3);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader :: histInitialize_1Lep()
{
  //m_hist_1lep_
  m_hist_1lep_NJets = new TH1F("NJets","Number of signal jets; N^{sig}_{Jets}; Jets",10,0,10);
  m_hist_1lep_NJets->Sumw2();
  wk()->addOutput(m_hist_1lep_NJets);
  m_hist_1lep_j1Pt = new TH1F("j1Pt","Pt of leading Jet",100,0,400);
  m_hist_1lep_j1Pt->Sumw2();
  wk()->addOutput(m_hist_1lep_j1Pt);
  m_hist_1lep_fatJetPt = new TH1F("FatJetPt","Pt of leading Fat Jet",100,200,600);
  m_hist_1lep_fatJetPt->Sumw2();
  wk()->addOutput(m_hist_1lep_fatJetPt);
  m_hist_1lep_fatJetm = new TH1F("FatJetm","Mass of leading Fat Jet",100,0,400);
  m_hist_1lep_fatJetm->Sumw2();
  wk()->addOutput(m_hist_1lep_fatJetm);
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode AnalysisReader :: histInitialize_2Lep()
{
  m_hist_2lep_ZPV = new TH1F("ZPV","z vertex; z [cm]; Events",100,-200,200);
  m_hist_2lep_ZPV->Sumw2();
  wk()->addOutput(m_hist_2lep_ZPV);
  m_hist_2lep_MET = new TH1F("MET","MET",100,0,100);
  m_hist_2lep_MET->Sumw2();
  wk()->addOutput(m_hist_2lep_MET);
  m_hist_2lep_NJets = new TH1F("NJets","Number of signal jets; N^{sig}_{Jets}; Jets",10,0,10);
  m_hist_2lep_NJets->Sumw2();
  wk()->addOutput(m_hist_2lep_NJets);
  m_hist_2lep_Mmumu = new TH1F("Mmumu","Z mass muon; M_{#mu#mu} [GeV]; Events",60,60,120);
  m_hist_2lep_Mmumu->Sumw2();
  wk()->addOutput(m_hist_2lep_Mmumu);
  m_hist_2lep_Mee = new TH1F("Mee","Z mass elec; M_{ee} [GeV]; Events",60,60,120);
  m_hist_2lep_Mee->Sumw2();
  wk()->addOutput(m_hist_2lep_Mee);
  m_hist_2lep_PtJets = new TH1F("PtJets","PtJets",100,0,100);
  m_hist_2lep_PtJets->Sumw2();
  wk()->addOutput(m_hist_2lep_PtJets);
  m_hist_2lep_EtaJets = new TH1F("EtaJets","EtaJets",60,-3,3);
  m_hist_2lep_EtaJets->Sumw2();
  wk()->addOutput(m_hist_2lep_EtaJets);
  m_hist_2lep_DPhiJets = new TH1F("DPhiJets","DPhiJets",50,0,M_PI);
  m_hist_2lep_DPhiJets->Sumw2();
  wk()->addOutput(m_hist_2lep_DPhiJets);
  m_hist_2lep_PtZ = new TH1F("PtZ","PtZ",100,0,500);
  m_hist_2lep_PtZ->Sumw2();
  wk()->addOutput(m_hist_2lep_PtZ);
  m_hist_2lep_PtZee = new TH1F("PtZee","PtZee",100,0,500);
  m_hist_2lep_PtZee->Sumw2();
  wk()->addOutput(m_hist_2lep_PtZee);
  m_hist_2lep_PtZmumu = new TH1F("PtZmumu","PtZmumu",100,0,500);
  m_hist_2lep_PtZmumu->Sumw2();
  wk()->addOutput(m_hist_2lep_PtZmumu);

  // Cut flow hist
  static std::string cuts [18] = {"All", "GRL", "ZPV", "2Lep", "2Jet", "2SigJet", "Mll", "MET", "DRmin", "Pt2045",  "DRjj",  "Bin0", "Bin1","Bin2","Bin3","Bin4","MjjUncor","Mjj"};

  m_hist_2lep_CutFlowMu = new TH1F("CutFlowMuon","CutFlowMuon", 18, 0.5, 18.5);
  m_hist_2lep_CutFlowMu->Sumw2();
  for(int i=0; i<18; i++) {
    m_hist_2lep_CutFlowMu->GetXaxis()->SetBinLabel(i+1,cuts[i].c_str());
  }
  m_hist_2lep_CutFlowMu_noWght = new TH1F("CutFlowMuon_noWeight","CutFlowMuon_noWeight", 18, 0.5, 18.5);
  m_hist_2lep_CutFlowMu_noWght->Sumw2();
  for(int i=0; i<18; i++) {
    m_hist_2lep_CutFlowMu_noWght->GetXaxis()->SetBinLabel(i+1,cuts[i].c_str());
  }
  wk()->addOutput(m_hist_2lep_CutFlowMu);
  wk()->addOutput(m_hist_2lep_CutFlowMu_noWght);
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode AnalysisReader :: histInitialize_VBFa()
{
  //m_hist_VBFa_
 
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader :: histInitialize_monoWZH()
{
  
  TH1::SetDefaultSumw2();
  TString buffer;
  int nbins;
  double min;
  double max;
  static std::string cuts [15] = {"All", "pre-selection", "met250", "met350", "one fatjet", "el veto", "mu veto", "addjetveto", "metjetOL", "500met" ,"mJ(WZ)", "mJ(H)", "1subjet", "2subjet", "2btag"};
  // before any pre-selection
  buffer = "Before pre-selection fatjet Mass"; nbins = 100; min = 0; max = 1000 ;
  m_hist_mono_prepre_fatjetm = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_prepre_fatjetm);

  buffer = "Before pre-selection fatjet p_{T}"; nbins = 100; min = 0; max = 1000 ;
  m_hist_mono_prepre_fatjetpt = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_prepre_fatjetpt);

  //pre-selection
  buffer = "Pre-selection E_{T}^{miss}"; nbins = 60; min = 200.; max = 2000. ;
  m_hist_mono_pre_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_MET);
  
  buffer = "Pre-selection #Phi(MET)"; nbins = 16; min = -4.; max = 4. ;
  m_hist_mono_pre_METPhi = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_METPhi);

  buffer = "Pre-selection Fatjet Multiplicity"; nbins = 7; min = -0.5; max = 6.5 ;
  m_hist_mono_pre_Nfatjet = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_Nfatjet);
  
  buffer = "Pre-selection fatjet Mass"; nbins = 100; min = 0; max = 1000 ;
  m_hist_mono_pre_fatjetm = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_fatjetm);
 
  buffer = "Pre-selection fatjet p_{T}"; nbins = 100; min = 0; max = 1000 ;
  m_hist_mono_pre_fatjetpt = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_fatjetpt);
  
  buffer = "Pre-selection fatjet Eta"; nbins = 10; min = 0; max = 10 ;
  m_hist_mono_pre_fatjeteta = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_fatjeteta);
  
   buffer = "Pre-selection jet Mass"; nbins = 50; min = 0; max = 500 ;
  m_hist_mono_pre_jetm = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_jetm);
  
  buffer = "Pre-selection jet p_{T}"; nbins = 100; min = 0; max = 1000 ;
  m_hist_mono_pre_jetpt = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_jetpt);
  
    buffer = "Pre-selection jet Eta"; nbins = 10; min = 0; max = 5 ;
  m_hist_mono_pre_jeteta = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_pre_jeteta);
  
  
  //cutflow
  buffer = "Cutflow"; nbins = 15; min = 0.5; max = 15.5 ;
  m_hist_mono_cutflow = new TH1D(TString(buffer),TString(buffer),nbins, min, max);
  for(unsigned int i=0; i<15; i++) {
    m_hist_mono_cutflow->GetXaxis()->SetBinLabel(i+1,cuts[i].c_str());
  }
  wk()->addOutput(m_hist_mono_cutflow);

  buffer = "Cutflow_noweight"; nbins = 15; min = 0.5; max = 15.5 ;
  m_hist_mono_cutflow_noweight = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  for(unsigned int i=0; i<15; i++) {
    m_hist_mono_cutflow_noweight->GetXaxis()->SetBinLabel(i+1,cuts[i].c_str());
  }
  wk()->addOutput(m_hist_mono_cutflow_noweight);
 
  buffer = "MET250 cutflow MET"; nbins = 180; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_met250_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_met250_MET);

  buffer = "MET350 cutflow MET"; nbins = 180; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_met350_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_met350_MET);
 
  buffer = "fatjet cutflow MET"; nbins = 180; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_fjet_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_fjet_MET);
  
  buffer = "elveto cutflow MET"; nbins = 180; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_elveto_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_elveto_MET);

  buffer = "muveto cutflow MET"; nbins = 180; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_muveto_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_muveto_MET);

  buffer = "jetveto cutflow MET"; nbins = 180; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_jetveto_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_jetveto_MET);
 
  buffer = "metjetOL cutflow MET"; nbins = 180; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_metjetOL_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_metjetOL_MET);
 
  buffer = "MET cutflow MET"; nbins = 180; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_met_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_met_MET);
  
  buffer = "m(fatjet)(WZ) cutflow MET"; nbins = 60; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_mjWZ_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_mjWZ_MET);

  buffer = "m(fatjet)(H) cutflow MET"; nbins = 60; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_mjH_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_mjH_MET);

  buffer = "1 subjet cutflow MET"; nbins = 60; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_1subjet_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_1subjet_MET);

  buffer = "2 subjet cutflow MET"; nbins = 60; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_2subjet_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_2subjet_MET);
 
  buffer = "2 b-tag cutflow MET"; nbins = 60; min = 200.; max = 2000. ;
  m_hist_mono_cutflow_2btag_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_2btag_MET);
 
  buffer = "Signal Fatjet Multiplicity"; nbins = 7; min = -0.5; max = 6.5 ;
  m_hist_mono_cutflow_Nsigfjet = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_Nsigfjet);
  
  buffer = "#Delta R(fatjet ,addjet)"; nbins = 10; min = 0.; max = 3. ;
  m_hist_mono_cutflow_fjetaddjetDeltaR = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_cutflow_fjetaddjetDeltaR);
  
  
  //n-1 histogram

  buffer = "met250 n-1"; nbins = 60; min = 0.; max = 1800. ;
  m_hist_nminusone_met250 = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_met250);

  buffer = "met350 n-1"; nbins = 60; min = 0.; max = 1800. ;
  m_hist_nminusone_met350 = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_met350);
 
  buffer = "fatjet n-1 signal fatjet Multiplicity"; nbins = 6; min = -0.5; max = 5.5 ;
  m_hist_nminusone_fjet = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_fjet);
  
  buffer = "lepveto n-1 sig lepton Multiplicity"; nbins = 10; min = -0.5; max = 9.5 ;
  m_hist_nminusone_lepveto = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_lepveto);

  buffer = "photveto n-1 sig photon Multiplicity"; nbins = 10; min = -0.5; max = 9.5 ;
  m_hist_nminusone_photveto = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_photveto);

  buffer = "addjetveto n-1 sig addjet Multiplicity"; nbins = 10; min = -0.5; max = 9.5 ;
  m_hist_nminusone_jetveto = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_jetveto);  

  buffer = "met-jet overlap n-1 (overlapping multiplicity)"; nbins = 7; min = -3.5; max = 3.5 ;
  m_hist_nminusone_metjetOL = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_metjetOL);

  buffer = "met500 n-1"; nbins = 60; min = 0.; max = 1800. ;
  m_hist_nminusone_met500 = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_met500);  
  
  buffer = "mJ n-1"; nbins = 50; min = 0.; max = 200. ;
  m_hist_nminusone_mj = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_nminusone_mj);    
  
  buffer = "Efficiency E_{T}^{miss}"; nbins = 60; min = 200.; max = 2000. ;
  m_hist_mono_eff_MET = new TH1F(TString(buffer),TString(buffer),nbins, min, max);
  wk()->addOutput(m_hist_mono_eff_MET);
  
	
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  //std::cout << "fileExecute Input file is " << wk()->inputFile()->GetName() << std::endl;

  // For everyfile we need to find the original sum of weights from the xAOD stored in the meta-data histogram
  TFile* inputfile=wk()->inputFile();

  TH1* metahist=(inputfile) ?  (TH1*)inputfile->Get("MetaData_EventCount") : 0;
  
  m_sumOfWeights=1.0;
  if (metahist) {
    float nevents=metahist->GetBinContent(7);
    m_sumOfWeights=metahist->GetBinContent(8);
    std::cout << " Number of Events/Sum of weights on xAOD " << nevents << " " << m_sumOfWeights << std::endl;
  }

  // Flag if is Pt0 slice, to kill overlapping events. PtV slices only in 13TeV xAOD samples
  TString filename(inputfile->GetName());
  m_isSherpaPt0VJets= (filename.Contains("Pt0") && m_comEnergy=="13TeV") ? true : false;
  if(filename.Contains("Sherpa_CT10_W") && filename.Contains("Pt0")  && m_comEnergy=="13TeV") m_SherpaPt0VJetsCut = 40000;

  // general Sherpa flag
  m_isSherpaVJets = (filename.Contains("Sherpa")) ? true : false;
  

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode AnalysisReader :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode AnalysisReader :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  std::cout<< "initial" <<std::endl; 
  m_event = wk()->xaodEvent();

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events on first file = %lli", m_event->getEntries() ); // print long long int
  
  m_eventCounter = 0;

  Info("initialize()", "Analysis type %i", m_analysisType ); 

  if (m_analysisType==e0Lep) 
    m_zerolepSel = new VHbb0lepEvtSelection(); 
  else if (m_analysisType==e1Lep) 
    m_onelepSel = new VHbb1lepEvtSelection(); 
  else if (m_analysisType==e2Lep) 
    m_twolepSel = new VHbb2lepEvtSelection(); 
  else if (m_analysisType==VBFa)
    m_VBFaSel = new VBFbbAEvtSelection();
  else if (m_analysisType==monoWZH)
	m_zerolepSel = new VHbb0lepEvtSelection();
  m_overlapRemoval = new OverlapRemoval();

  //
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo___Nominal").isSuccess() ){
    Error("initialize()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // get MC flag - different info on data/MC files
  m_isMC=m_superDecorator.get(eventInfo, EventInfoIntProps::isMC);
  Info("initialize()", "isMC = %i", m_isMC );

  // initilise cross section provider - steerable centre of mass?
  if (m_isMC) {
    std::string xSectionFile = gSystem->Getenv("ROOTCOREBIN");
    if(m_comEnergy=="13TeV") 
      xSectionFile += "/data/CxAODTools/XSections_13TeV.txt";
    else if (m_comEnergy=="8TeV") 
      xSectionFile += "/data/CxAODTools/XSections_8TeV.txt";
    else {
      Error("initialize()", "Unknown COM energy. Exiting. %s",m_comEnergy.Data());
      return EL::StatusCode::FAILURE;
    }
    m_xSection = new XSectionProvider(xSectionFile);
  }

  // initilise sumOfWeights provider
  if (m_isMC) {
    // which analysis
    std::string ana_read = "";
    if (m_analysisType==e0Lep) ana_read = "0lep";
    else if (m_analysisType==e1Lep) ana_read = "1lep";
    else if (m_analysisType==e2Lep) ana_read = "2lep";
	else if (m_analysisType==monoWZH) ana_read = "0lep";
    else {
      Warning("histInitialize()", "Invalid analysis type %i",m_analysisType);
      return EL::StatusCode::FAILURE;
    }
    //
    std::string sumOfWeightsFile = gSystem->Getenv("ROOTCOREBIN");
    sumOfWeightsFile += "/data/CxAODTools/yields.";
    sumOfWeightsFile += ana_read;
    if(m_comEnergy=="13TeV") sumOfWeightsFile += ".13TeV.txt";
    else if (m_comEnergy=="8TeV") sumOfWeightsFile += ".8TeV.txt";
    else {
      Error("initialize()", "Unknown COM energy. Exiting. %s",m_comEnergy.Data());
      return EL::StatusCode::FAILURE;
    }
    m_sumOfWeights_fix = new sumOfWeightsProvider(sumOfWeightsFile);
    cout << "sumOfWeightsFile = " << sumOfWeightsFile << endl;
  }
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  bool debug = true;
  //std::cout<< "execute" <<std::endl; 
  // print every 10000 events
  if( (m_eventCounter % 10000) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
  m_eventCounter++;

  if( m_eventCounter==1){
    double weighted_value = m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("All"));
    cout << "weighted_value: " << weighted_value << endl;
  }

  if( m_eventCounter==41){
    double weighted_value = m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("All"));
    cout << "weighted_value: " << weighted_value << endl;
  }

  //----------------------------
  // Event information
  //--------------------------- 

  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo___Nominal").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // get MC flag - different info on data/MC files - now set in job initialise
  // m_isMC=m_superDecorator.get(eventInfo, EventInfoIntProps::isMC);

  // detemine if it is a Sherpa Pt0 sample, and remove events overlapping with slice (flag set for 13TeV only)
  float vpt_truth(-999);
  if (m_isSherpaVJets) {
    const xAOD::TruthParticleContainer* truthParts = 0;
    if( ! m_event->retrieve( truthParts, "TruthParticle___Nominal").isSuccess() ){
      Error("execute()", "Failed to retrieve  TruthParticle collection. Exiting." );
      return EL::StatusCode::FAILURE;
    }
    const xAOD::TruthParticle* lep0(0);
    const xAOD::TruthParticle* lep1(0);
    for (unsigned int iPart(0) ; iPart < truthParts->size(); ++iPart) {
      const xAOD::TruthParticle* part = truthParts->at(iPart);
      // assume leptons are the only status 3 particles in SherpaV+jets
      if (part->status()==3) {
	if (!lep0) lep0=part;
	else lep1=part;
      }
    }  
    if (lep0 && lep1) {
      TLorentzVector V(lep0->p4());
      V+=lep1->p4();
      vpt_truth=V.Pt();
      if (m_isSherpaPt0VJets &&  vpt_truth > m_SherpaPt0VJetsCut)  return EL::StatusCode::SUCCESS;
    }
  }

 

  // reset event weight
  m_weight=1.0; 


  // Query - reading the cross section every event.
  // This is because on the first event fileExecute and changeInput are called before initialize
  // so there is no event available to read datasetid
  if (m_isMC) {
    int datasetid=eventInfo->mcChannelNumber();
    // fix m_sumOfWeights when production done on xAOD
    m_sumOfWeights = (m_sumOfWeights_fix) ? m_sumOfWeights_fix->getsumOfWeights(datasetid) : 1.0;
    //
    float sigmaEff = (m_xSection) ? m_xSection->getXSection(datasetid) : 1.0;
    //std::cout << "Cross section times eff. for dataset id " << datasetid << " = " << sigmaEff << std::endl;

    // we are normalising to MC lumi, sumOfWeights calculated per file.
    m_weight= (m_sumOfWeights) ? sigmaEff/m_sumOfWeights : 1.0;			// 
  }

  // multiply by MC generator weight 
  if (m_isMC) 
    m_weight *=m_superDecorator.get(eventInfo,EventInfoFloatProps::MCEventWeight);

  //std::cout << "m_weight " << m_weight << std::endl;
  
  // Sherpa Vpt  
  if (m_isSherpaVJets) 
    m_hist_VPtTruth->Fill(vpt_truth/1000, m_weight);
  
  if (m_analysisType==e0Lep) 
    fill_0Lep();
  else if (m_analysisType==e1Lep) 
    fill_1Lep();
  else if (m_analysisType==e2Lep) 
    fill_2Lep();
  else if (m_analysisType==VBFa)
    fill_VBFa();
  else if (m_analysisType==monoWZH)
	fill_monoWZH();
 
  
  
  return EL::StatusCode::SUCCESS;

}

EL::StatusCode AnalysisReader :: fill_0Lep()
{

  // data GRL
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo___Nominal").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  if (!m_isMC) {
    bool passGRL=m_superDecorator.get(eventInfo, EventInfoIntProps::passGRL);
    if(!passGRL) return EL::StatusCode::SUCCESS;
  }

  // lookout for events with no z-vertex
  float ZPV=m_superDecorator.get(eventInfo, EventInfoFloatProps::ZPV);
  if (ZPV<-9998) {
    std::cout << " Event without primary vertex " << ZPV  << " " << eventInfo->runNumber() <<  " " << eventInfo->eventNumber() << std::endl;
    return EL::StatusCode::SUCCESS;
  }

  //----------------------------
  // Jets
  //--------------------------- 
  
  const xAOD::JetContainer* jets = 0;
  if ( !m_event->retrieve( jets, "AntiKt4EMTopoJets___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve AntiKt4LCTopoJets container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //----------------------------
  // Muons
  //--------------------------- 

  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
  //if ( !m_event->retrieve( muons, "Muons___Original" ).isSuccess() ){ // retrieve arguments: container type, container key
  //if ( !m_event->retrieve( muons, "Muons___MUONS_SCALE__1down" ).isSuccess() ){ // retrieve arguments: container type, container key
  //if ( !m_event->retrieve( muons, "Muons___MUONS_SCALE__1up" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //----------------------------
  // Elecrons
  //--------------------------- 
  const xAOD::ElectronContainer* elecs = 0;
  if ( !m_event->retrieve( elecs, "ElectronCollection___Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //---------------------------
  //Photons
  //---------------------------
  const xAOD::PhotonContainer* phots = 0;
  if( !m_event->retrieve( phots, "PhotonCollection__Nominal" ).isSuccess() ){
  //  Error("execute()", "Failed to retrieve Photons container. Exiting." );
  //  return EL::StatusCode::FAILURE;
  }


  //----------------------------
  // Missing Et
  //--------------------------- 
  const xAOD::MissingETContainer* METNominal = 0;
  if ( !m_event->retrieve( METNominal, "MET_RefFinal___Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Missing Et container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  const xAOD::MissingET * met = 0;
  if(METNominal->size() > 0) 
  {
  met = METNominal->at(0);
  //m_hist_monoWZH_MET->Fill(met->met(),m_weight);
  }
  //----------------------------
  // Overlap Removal
  //--------------------------- 
  // Do here for Nominal, systematics
  m_overlapRemoval->removeOverlap(elecs, phots, muons, jets);

  //----------------------------
  // Event Selection
  //--------------------------- 
  
  
  bool pass0lep=m_zerolepSel->passSelection(eventInfo,
                                            met,
                                            elecs,
					    phots,
                                            muons,
                                            jets);

  //  if(debug)  std::cout << " Passed 0 lepton "  << pass0lep << std::endl;

  ResultVHbb0lep zerolepResult=m_zerolepSel->result();

  if (zerolepResult.pass != pass0lep)
    std::cout << " 0 lep ERROR "  << pass0lep << " " << zerolepResult.pass << std::endl;

  // select events
  if (!pass0lep) return EL::StatusCode::SUCCESS;

  //
  //const xAOD::Electron* el1=zerolepResult.el1;
  //const xAOD::Electron* el2=zerolepResult.el2;
  //const xAOD::Muon* mu1=zerolepResult.mu1;
  //const xAOD::Muon* mu2=zerolepResult.mu2;
  std::vector<const xAOD::Jet*> selJets=zerolepResult.jets;

  // the derivation/selection keeps high pt 1-jet events, reject here
  int njets=selJets.size();
  if(njets<2) return EL::StatusCode::SUCCESS;
  
  // jets are selected with at least 2 signal jets, but selection vector includes forward jets
  int nsigjets=0;
  std::vector<const xAOD::Jet*> sigJets;
  for (int i(0) ; i<njets ; ++i) {
    if (fabs(selJets.at(i)->eta())>2.5) continue;
    if (selJets.at(i)->pt()<20000) continue;
    //
    float jvf0;
    selJets.at(i)->getAttribute("jvf0",jvf0);
    if (selJets.at(i)->pt()<50000 && fabs(selJets.at(i)->eta())<2.4 && fabs(jvf0)<=0.5) continue;
    //
    nsigjets++;
    sigJets.push_back(selJets.at(i));
  }
  //
  int nforwardjets=0;
  for (int i(0) ; i<njets ; ++i) {
    if (fabs(selJets.at(i)->eta())<2.5) continue;
    if (selJets.at(i)->pt()<30000) continue;
    nforwardjets++;
  }

  // require 2 or 3 signal jets
  if(nsigjets<2 || nsigjets>3) return EL::StatusCode::SUCCESS;

  // no forward jets
  if(nforwardjets>0) return EL::StatusCode::SUCCESS;
                                      
  // require leading jet 45
  if (sigJets.at(0)->pt()<45000) return EL::StatusCode::SUCCESS;


  // and finally cut on MET
  double MET = met->met()/1000.;
  if(MET<120) return EL::StatusCode::SUCCESS;

  // LV
  TLorentzVector j1Vec, j2Vec, j3Vec, metVec;
  j1Vec.SetPtEtaPhiM(sigJets.at(0)->pt()/1000., sigJets.at(0)->eta(), sigJets.at(0)->phi(), sigJets.at(0)->m()/1000.);
  j2Vec.SetPtEtaPhiM(sigJets.at(1)->pt()/1000., sigJets.at(1)->eta(), sigJets.at(1)->phi(), sigJets.at(1)->m()/1000.);
  if (nsigjets==3) j3Vec.SetPtEtaPhiM(sigJets.at(2)->pt()/1000., sigJets.at(2)->eta(), sigJets.at(2)->phi(), sigJets.at(2)->m()/1000.);
  metVec.SetPxPyPzE(met->mpx()/1000., met->mpy()/1000., 0, met->met()/1000.);

  // BDT input vars (for later should be jet corrected quantities)
  double bdt_dRBB = j1Vec.DeltaR(j2Vec);
  double bdt_dEtaBB = fabs( j1Vec.Eta() - j2Vec.Eta() );
  double bdt_dPhiVBB = fabs( metVec.DeltaPhi(j1Vec+j2Vec) );
  double bdt_mBB = (j1Vec+j2Vec).M();
  double bdt_MEff = j1Vec.Pt() + j2Vec.Pt() + metVec.Pt();
  double bdt_MEff3 = 0;
  if (nsigjets==3) bdt_MEff3 = j1Vec.Pt() + j2Vec.Pt() + j3Vec.Pt() + metVec.Pt();
  double bdt_pTB1 = j1Vec.Pt();
  double bdt_pTB2 = j2Vec.Pt();
  double bdt_MET = metVec.Pt();
  double bdt_MBBJ = 0;
  if (nsigjets==3) bdt_MBBJ = (j1Vec+j2Vec+j3Vec).M();
  double bdt_pTJ3 = 0;
  if (nsigjets==3) bdt_pTJ3 = j3Vec.Pt();

  // MVA selection (for later should be on jet uncorrected quantities)
  //-------------------------------------------------------

  // dR cut
  if (MET<200 && bdt_dRBB<0.7) return EL::StatusCode::SUCCESS;
  
  // mindPhi cut
  double mindPhi1 = fabs( j1Vec.DeltaPhi(metVec) );  
  double mindPhi2 = fabs( j2Vec.DeltaPhi(metVec) );       
  double mindPhi3 = 1000;
  if (nsigjets==3) fabs( j3Vec.DeltaPhi(metVec) );       
  //
  double mindPhi = mindPhi1;
  if (mindPhi2<mindPhi) mindPhi = mindPhi2;
  if (mindPhi3<mindPhi) mindPhi = mindPhi3;
  //
  if (mindPhi<1.5) return EL::StatusCode::SUCCESS;

  // dPhiMETMPT cut -> not possible yet

  // fill histos
  if (nsigjets==2) {
    m_hist2j_0lep_dRBB->Fill(bdt_dRBB, m_weight);
    m_hist2j_0lep_dEtaBB->Fill(bdt_dEtaBB, m_weight);
    m_hist2j_0lep_dPhiVBB->Fill(bdt_dPhiVBB, m_weight);
    m_hist2j_0lep_mBB->Fill(bdt_mBB, m_weight);
    m_hist2j_0lep_MEff->Fill(bdt_MEff, m_weight);
    m_hist2j_0lep_pTB1->Fill(bdt_pTB1, m_weight);
    m_hist2j_0lep_pTB2->Fill(bdt_pTB2, m_weight);
    m_hist2j_0lep_MET->Fill(bdt_MET, m_weight);
  }
  if (nsigjets==3) {
    m_hist3j_0lep_dRBB->Fill(bdt_dRBB, m_weight);
    m_hist3j_0lep_dEtaBB->Fill(bdt_dEtaBB, m_weight);
    m_hist3j_0lep_dPhiVBB->Fill(bdt_dPhiVBB, m_weight);
    m_hist3j_0lep_mBB->Fill(bdt_mBB, m_weight);
    m_hist3j_0lep_MEff3->Fill(bdt_MEff3, m_weight);
    m_hist3j_0lep_pTB1->Fill(bdt_pTB1, m_weight);
    m_hist3j_0lep_pTB2->Fill(bdt_pTB2, m_weight);
    m_hist3j_0lep_MET->Fill(bdt_MET, m_weight);
    m_hist3j_0lep_MBBJ->Fill(bdt_MBBJ, m_weight);
    m_hist3j_0lep_pTJ3->Fill(bdt_pTJ3, m_weight);
  }

  return EL::StatusCode::SUCCESS;
} 

EL::StatusCode AnalysisReader :: fill_1Lep()
{

  bool debug = false;

//need to add cutflow hist paul
  //data GRL
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo___Nominal").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  if (!m_isMC) {
    bool passGRL=m_superDecorator.get(eventInfo, EventInfoIntProps::passGRL);
    if(!passGRL) return EL::StatusCode::SUCCESS;
  }

  // lookout for events with no z-vertex
  float ZPV=m_superDecorator.get(eventInfo, EventInfoFloatProps::ZPV);
  if (ZPV<-9998) {
    std::cout << " Events without primary vertex " << ZPV  << " " << eventInfo->runNumber() <<  " " << eventInfo->eventNumber() << std::endl;
    return EL::StatusCode::SUCCESS;
  }

  //----------------------
  //Load Objects
  //----------------------


  const xAOD::JetContainer* jets = 0;
  if ( !m_event->retrieve( jets, "AntiKt4EMTopoJets___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve AntiKt4LCTopoJets container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //-----------------------
  //FatJets
  //----------------------


  const xAOD::JetContainer* fatjets = 0;
  if ( !m_event->retrieve( fatjets, "CamKt12LCTopoJets___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve CamKt12LCTopoJets container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::ElectronContainer* elecs = 0;
  if ( !m_event->retrieve( elecs, "ElectronCollection___Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::MissingETContainer* METNominal = 0;
  if ( !m_event->retrieve( METNominal, "MET_RefFinal___Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Missing Et container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //---------------------------
  //Photons
  //---------------------------
  const xAOD::PhotonContainer* phots = 0;
  if( !m_event->retrieve( phots, "PhotonCollection__Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Photons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::MissingET * met = 0;
  //if(METNominal->size() > 0) met = METNominal->at(0);
  met = METNominal->at(0);//is this right paul?
  //nope, want met from  VHbbEvntSelectio

  //----------------------------
  // Overlap Removal
  //--------------------------- 
  // Do here for Nominal, systematics
  //m_overlapRemoval->removeOverlapConst(elecs, muons, jets);
  m_overlapRemoval->removeOverlap(elecs, phots, muons, jets);

  //------------------------
  //Event Selection
  //-------------------------

  //make mass plot fo fat jets before any selection
  bool isFatEvent = false;
  std::vector<const xAOD::Jet*> fatsigJets;
  for (unsigned int iJet(0) ; iJet < fatjets->size(); ++iJet) {
    const xAOD::Jet * fatjet = fatjets->at(iJet);
      if(m_superDecorator.get(fatjet,FatJetIntProps::isFatJet)){
          isFatEvent = true;
          m_hist_1lep_fatJetPt->Fill(fatjet->pt()/1000, m_weight);
          m_hist_1lep_fatJetm->Fill(fatjet->m()/1000, m_weight);
      }else{
       //   std::cout << "not a fatjet" << std::endl;
      }

  }

  bool pass1lep=m_onelepSel->passSelection(eventInfo,
                                           met,
                                           elecs,
                                           phots,
                                           muons,
                                           jets);

  if(debug)  std::cout << " Passed 1 lepton "  << pass1lep << std::endl;

  ResultVHbb1lep onelepResult=m_onelepSel->result();

  if (onelepResult.pass != pass1lep)
    std::cout << " 1 lep ERROR "  << pass1lep << " " << onelepResult.pass << std::endl;

  if (!pass1lep) return EL::StatusCode::SUCCESS;

  //do final event selection
  const xAOD::Electron* el1=onelepResult.el;
//  const xAOD::Electron* el2=onelepResult.el2;
  const xAOD::Muon* mu1=onelepResult.mu;
//  const xAOD::Muon* mu2=onelepResult.mu2;
//  el2 and mu2 are just junk in the 1 lep analysis.
//  access met like el, mu?
//  probs not, doesnt need ordered
  std::vector<const xAOD::Jet*> selJets=onelepResult.jets;


  if (pass1lep&&debug) {
    std::cout << " el1 "  << el1 << std::endl;
    std::cout << " mu1 "  << mu1 << std::endl;
    
    //add met here?

    int njet=selJets.size();
    std::cout << " found 1 lep jets " <<  njet << std::endl;
    for (int i(0) ; i<njet ; ++i) {
      std::cout << i << " jet pt " << selJets.at(i)->pt() << std::endl;
    }        
  }



  //only doing muon for now
  if (!mu1) return EL::StatusCode::SUCCESS;

  // the derivation/selection keeps high pt 1-jet events, reject here
  // change for fat jet stuff?
  int njets=selJets.size();
  if(njets<2) return EL::StatusCode::SUCCESS;

  // jets are selected with at least 2 signal jets, but selection vector includes forward jets
  int nsigjets=0;
  std::vector<const xAOD::Jet*> sigJets;
  for (int i(0) ; i<njets ; ++i) {
    if (fabs(selJets.at(i)->eta())>2.5) continue;
    if (selJets.at(i)->pt()<20000) continue;
    nsigjets++;
    sigJets.push_back(selJets.at(i));
  } 

  // require 2 signal jets
  if(nsigjets<2) return EL::StatusCode::SUCCESS;

  // check selection jets with xAOD jets
  int njets1(0);
  for (unsigned int iJet(0) ; iJet < jets->size(); ++iJet) {
    const xAOD::Jet * jet = jets->at(iJet);
    if (!m_superDecorator.get(jet,JetIntProps::passOR)) continue; // should be updated by overlapRemover
    if (m_superDecorator.get(jet,JetIntProps::isVetoJet)) { // applied to keep CxAOD jets already
      if (m_superDecorator.get(jet,JetIntProps::isSignalJet)) // pt>20 GeV, |eta|<2.5
        njets1++;
    }
  }

  if (njets1 != nsigjets) 
    std::cout << "njets1 " << njets1 << " nsigjets " << nsigjets << std::endl;

  TLorentzVector mu1Vec;
  mu1Vec.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), mu1->m());
  //make jets and fill?
  //
  if (selJets.size() > 1  ){
    if(selJets.at(0)->pt()> 45000){
      m_hist_1lep_j1Pt->Fill(sigJets.at(0)->pt()/1000, m_weight);
    }
  }
  m_hist_1lep_NJets->Fill(nsigjets,m_weight);

  return EL::StatusCode::SUCCESS;
} 

EL::StatusCode AnalysisReader :: fill_monoWZH()
{

  bool debug = true;
  if(debug) std::cout<< "Filling" <<std::endl; 

//need to add cutflow hist paul
  //data GRL
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo___Nominal").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //cout << "eventInfo->eventNumber(): " << eventInfo->eventNumber() << endl;

  if (!m_isMC) {
    bool passGRL=m_superDecorator.get(eventInfo, EventInfoIntProps::passGRL);
    if(!passGRL) return EL::StatusCode::SUCCESS;
  }

  // lookout for events with no z-vertex
  float ZPV=m_superDecorator.get(eventInfo, EventInfoFloatProps::ZPV);
  if (ZPV<-9998) {
    std::cout << " Events without primary vertex " << ZPV  << " " << eventInfo->runNumber() <<  " " << eventInfo->eventNumber() << std::endl;
    return EL::StatusCode::SUCCESS;
  }

  //----------------------
  //Load Objects
  //----------------------


  const xAOD::JetContainer* jets = 0;
  if ( !m_event->retrieve( jets, "AntiKt4LCTopoJets___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key		
  //if ( !m_event->retrieve( jets, "AntiKt4EMTopoJets___Nominal" ).isSuccess() ){								// AntiKt4EMTopoJets___Nominal for initial sensitivity testing, update to LC on Mar 31
    Error("execute()", "Failed to retrieve AntiKt4EMTopoJets container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //cout << "jets->size(): " << jets->size() << endl;
  //for(int i(0); i < jets->size() ;i++){
  //   cout << "jets->at(i)->pt(): " << jets->at(i)->pt() << ", jets->at(i)->eta(): " << jets->at(i)->eta() << endl;
  //}

  //-----------------------
  //FatJets
  //----------------------


  const xAOD::JetContainer* fatjets = 0;
  //if ( !m_event->retrieve( fatjets, "CamKt12LCTopoJets___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
  if ( !m_event->retrieve( fatjets, "AntiKt10LCTopoJets___TrimmedPtFrac5SmallR20" ).isSuccess() ){ 						// CamKt12LCTopoJets___Nominal for standard CxAOD checking, AntiKt10LCTopoJets___TrimmedPtFrac5SmallR20 for testing Arturo's
    Error("execute()", "Failed to retrieve AntiKt10LCTopoJets___TrimmedPtFrac5SmallR20 container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

//  const xAOD::JetContainer* fatjets_Cam = 0;			// if running on my (old) CxAODs, change this name to fatjets and comment out above
//  if ( !m_event->retrieve( fatjets, "CamKt12LCTopoJets___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
//    Error("execute()", "Failed to retrieve AntiKt10LCTopoJets___TrimmedPtFrac5SmallR20 container. Exiting." );
//    return EL::StatusCode::FAILURE;
//  } 

//  cout << "fatjets->size(): " << fatjets->size() << endl;
//  for(int i(0); i < fatjets->size() ;i++){
//     cout << "fatjets->at(i)->pt(): " << fatjets->at(i)->pt() << ", fatjets->at(i)->eta(): " << fatjets->at(i)->eta() << endl;
//  }

  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::ElectronContainer* elecs = 0;
  if ( !m_event->retrieve( elecs, "ElectronCollection___Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::MissingETContainer* METNominal = 0;
  if ( !m_event->retrieve( METNominal, "MET_RefFinal___Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Missing Et container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //cout << "METNominal->at(0)->met(): " << METNominal->at(0)->met() << endl;  

  //---------------------------
  //Photons
  //---------------------------
  const xAOD::PhotonContainer* phots = 0;
  //if( !m_event->retrieve( phots, "PhotonCollection__Nominal" ).isSuccess() ){
  //  Error("execute()", "Failed to retrieve Photons container. Exiting." );
  //  return EL::StatusCode::FAILURE;
  //}

  const xAOD::MissingET * met = 0;
  if(METNominal->size() > 0) met = METNominal->at(0);
  
  //----------------------------
  // Overlap Removal
  //--------------------------- 
  // Do here for Nominal, systematics
  //m_overlapRemoval->removeOverlap(elecs, phots, muons, jets);

  //------------------------
  //Event Selection
  //-------------------------

  //make mass plot fo fat jets before any selection
  //bool isFatEvent = false;
  //std::vector<const xAOD::Jet*> fatsigJets;
  // for (unsigned int iJet(0) ; iJet < fatjets->size(); ++iJet) {
   // const xAOD::Jet * fatjet = fatjets->at(iJet);
      // if(m_superDecorator.get(fatjet,FatJetIntProps::isFatJet)){
          // isFatEvent = true;
		  //Need to add new hists
          //m_hist_1lep_fatJetPt->Fill(fatjet->pt()/1000, m_weight);
          //m_hist_1lep_fatJetm->Fill(fatjet->m()/1000, m_weight);
      // }else{
         //std::cout << "not a fatjet" << std::endl;
      // }

  // }
  //Need to def new selection
  bool passmonoWZH=true;
//  bool passmonoWZH=m_zerolepSel->passSelection(eventInfo,
//                                            met,
//                                            elecs,
//											phots,
//                                            muons,
//                                            jets);

  //  if(debug)  std::cout << " Passed 0 lepton "  << pass0lep << std::endl;
  //Following code need to be check
  ResultVHbb0lep zerolepResult=m_zerolepSel->result();
    
  //if (zerolepResult.pass != pass0lep)
  //  std::cout << " 0 lep ERROR "  << pass0lep << " " << zerolepResult.pass << std::endl;

  // select events
  if (!passmonoWZH) return EL::StatusCode::SUCCESS;

  //std::vector<const xAOD::Jet*> selJets=zerolepResult.jets;

  std::vector<const xAOD::Jet*> selJets;
  for(unsigned int iJet(0); iJet < jets->size();iJet++){
    const xAOD::Jet* Jet = jets->at(iJet);
    selJets.push_back(Jet);
  }
  
  //pre-selection 
  //could be implement in CxAODMaker  
  bool pass_presele = true;
  if(METNominal->size() == 0) pass_presele = false;
  // MET > 200GeV
  //if((met->met()/1000.) < 200 ) pass_presele = false;
  // at least one fatjet 
  //if(fatjets->size() < 1) pass_presele = false;
  if(debug) std::cout << "Done pre_selection" << std::endl; 
  //Event selection
  bool pass_met250 = false;
  bool pass_fatjet = false;
  bool pass_elveto = false;
  bool pass_muveto = false;
  bool pass_lepveto = false;
  bool pass_photveto = false;
  bool pass_jetveto = false;
  bool pass_metcut = false;
  bool pass_jetmass_WZ = false;
  bool pass_jetmass_H = false;
  bool pass_jetOL = false;
  bool pass_1subjet = false;
  bool pass_2subjet = false;
  bool pass_2btag = false; 
 
  m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("All"), m_weight);
  m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("All"));

if(pass_presele)
{
  // MET cut at 250 GeV
  if (met->met()/1000. > 250.) pass_met250 = true;

  //fat jet Pt > 250GeV |eta| < 1.2
  int nfatjet = fatjets->size();
  std::vector<const xAOD::Jet*> fatsigJets;
  if(debug) std::cout << "Done fatjet cut start" << std::endl;
  for(int i(0); i < nfatjet ;i++){
          if(fatjets->at(i)->pt() < 250000.) continue;
          if(fabs(fatjets->at(i)->eta()) > 1.2) continue;
          fatsigJets.push_back(fatjets->at(i));
  }
  if(fatsigJets.size() > 0) pass_fatjet = true;

  // define leading fat jet
  float max_pt = 0.;
  int lead_i;
  TLorentzVector lead_fjet;             // define the leading fat jet
  if(fatsigJets.size() > 0){
          for (int i(0); i< fatsigJets.size(); i++){
                  if (fatsigJets.at(i)->pt() > max_pt){
                           max_pt = fatsigJets.at(i)->pt();
                           lead_i = i;
                  }
          }
          lead_fjet.SetPtEtaPhiM(fatsigJets.at(lead_i)->pt(), fatsigJets.at(lead_i)->eta(), fatsigJets.at(lead_i)->phi(), fatsigJets.at(lead_i)->m());
  }
  if(debug) std::cout << "Done fatjet cut" << std::endl;

  int n_el = elecs->size();
  int n_mu = muons->size();
  //int n_ph = phots->size();
  if(debug) cout << "n_el: " << n_el << ", n_mu: " << n_mu << endl; // ", n_ph: " << n_ph << endl;

  // veto electrons which pass: isVHLooseElectron, ptcone20/pt < 0.1
  std::vector<const xAOD::Electron*> vetoEls;
  for (int i(0); i < n_el; i++){
    const xAOD::Electron* elec = elecs->at(i);
    bool looseEle = m_superDecorator.get(elec, ElecIntProps::isVHLooseElectron);
    float ptEle = elec->pt();
    float ptconeEle = elec->isolationValue(xAOD::Iso::ptcone20);
    if(debug) cout << "looseEle: " << looseEle << ", ptEle: " << ptEle << ", ptconeEle: " << ptconeEle << endl;
    if (looseEle == 0) continue;
    if (ptconeEle/ptEle > 0.1) continue;
    vetoEls.push_back(elec);
  }

  // veto muons which pass: isVHLooseMuon, ptcone20/pt < 0.1
  std::vector<const xAOD::Muon*> vetoMus;
  for (int i(0); i < n_mu; i++){
    const xAOD::Muon* muon = muons->at(i);
    bool looseMu = m_superDecorator.get(muon, MuonIntProps::isVHLooseMuon);
    float ptMu = muon->pt();
    float ptconeMu = -999.;
    muon->isolation(ptconeMu,xAOD::Iso::ptcone20); 
    if(debug) cout << "looseMu: " << looseMu << ", ptMu: " << ptMu << ", ptconeMu: " << ptconeMu << endl;
    if (looseMu == 0) continue;
    if (ptconeMu/ptMu > 0.1) continue;
    vetoMus.push_back(muon);
  }

  // veto photos - ignore for now
  //std::vector<const xAOD::Photon*> vetoPhots;
  //for (int i(0); i < n_ph; i++){
  //  if (phots->at(i)->pt() < 10000.) continue;
  //  if (fabs(phots->at(i)->eta()) > 2.37) continue;
  //  vetoPhots.push_back(phots->at(i));
  //}

  if(debug) cout << "n_vetoEl: " << vetoEls.size() << ", n_vetoMu: " << vetoMus.size() << endl; // ", n_vetoPhots: " << vetoPhots.size() << endl;

  if (vetoEls.size() == 0) pass_elveto = true;
  if (vetoMus.size() == 0) pass_muveto = true;

  // create set of narrow signal jets: 
  // pT > 20, eta < 4.5 (loosest pt, eta cuts), 
  // deltaR(electron) > 0.4
  std::vector<const xAOD::Jet*> sigJets_base;
  int nbasejets=selJets.size();						// note that selJets come from 'std::vector<const xAOD::Jet*> selJets=zerolepResult.jets;', is this correct?
  for(int i(0);i < nbasejets;i++){                                          // loop over antiKt jets
    const xAOD::Jet* jet_1 = selJets.at(i);
    if(jet_1->pt() < 20000.) continue;
    if(fabs(jet_1->eta()) > 4.5) continue;
    TLorentzVector narjet;
    narjet.SetPtEtaPhiM(jet_1->pt(), jet_1->eta(), jet_1->phi(), jet_1->m());
    bool jetElOL = false;
    if (vetoEls.size() > 0){
      for (int i(0); i < vetoEls.size(); i++){				// loop over all veto-electrons
        const xAOD::Electron* vetelec = vetoEls.at(i);
        TLorentzVector vetoEl;
  	vetoEl.SetPtEtaPhiM(vetelec->pt(), vetelec->eta(), vetelec->phi(), vetelec->m());
        float deltaR = 20.;
        deltaR = vetoEl.DeltaR(narjet);			// check the overlap
        if (deltaR < 0.4) jetElOL = true;				// if overlaps with any electron, we drop the jet
      }	  
    }
    if (jetElOL == true) continue; 
    sigJets_base.push_back(jet_1);
  }

  // create one stream of jets with high pT (40 GeV)
  std::vector<const xAOD::Jet*> sigJets;
  int njets=sigJets_base.size();
  for(int i(0);i < njets;i++){
    const xAOD::Jet* jet_2 = sigJets_base.at(i);
    if(jet_2->pt() < 40000.) continue;
    // check overlap with lead fat jet
    TLorentzVector addjet;
    addjet.SetPtEtaPhiM(jet_2->pt(), jet_2->eta(), jet_2->phi(), jet_2->m());
    float deltaR_fj = 20.;
    if (fatsigJets.size() > 0) deltaR_fj = lead_fjet.DeltaR(addjet);
    if (deltaR_fj < 0.9) continue;
    sigJets.push_back(jet_2);
  }

  // create second stream of jets with:
  // close to fat jet
  // eta < 2.5
  std::vector<const xAOD::Jet*> subJets;
  for(int i(0);i < njets;i++){
    const xAOD::Jet* jet_3 = sigJets_base.at(i);
    if(sigJets_base.at(i)->eta() > 2.5) continue;
    // require overlap with fat jet
    TLorentzVector addjet;
    addjet.SetPtEtaPhiM(jet_3->pt(), jet_3->eta(), jet_3->phi(), jet_3->m());
    float deltaR_fj2 = 20.;								// some redundancy here with the code above :(
    if (fatsigJets.size() > 0) deltaR_fj2 = lead_fjet.DeltaR(addjet);
    if (deltaR_fj2 > 1.0) continue;
    subJets.push_back(jet_3);
  }

  // check for b-tag
  std::vector<const xAOD::Jet*> b_subJets;
  for(int i(0);i < subJets.size();i++){
    const xAOD::Jet* jet_5 = subJets.at(i);
    float btag = m_superDecorator.get(jet_5, JetFloatProps::MV1);
    if (btag < 0.971966) continue;
    b_subJets.push_back(jet_5);
  }

  // veto event if more than 1 jet passing first stream
  if (sigJets.size() < 2) pass_jetveto = true;

  // veto events if (base) narrow jet overlaps with MET
  int num_close_jets = 0;
  TLorentzVector metVec;
  metVec.SetPxPyPzE(met->mpx()/1000., met->mpy()/1000., 0, met->met()/1000.);
  for (int i(0); i < sigJets_base.size(); i++){						// loop over the 'base' jets, 
   const xAOD::Jet* jet_4 = sigJets_base.at(i);						// this should be sigJets_base, but change for comparison with YT
    TLorentzVector seljet;
    seljet.SetPtEtaPhiM(jet_4->pt(), jet_4->eta(), jet_4->phi(), jet_4->m());
    float delPhi = fabs(seljet.DeltaPhi(metVec));					
    if (delPhi <= 0.4) num_close_jets++;
  }
  if (num_close_jets == 0) pass_jetOL = true;

  // MET > 500 GeV
  if((met->met()/1000.) > 500 ) pass_metcut = true;

  // cut on mass of the fat jet:
  // 50 < m_fj < 120 for mono-WZ channel
  // 90 < m_fj < 145 for mono-H channel
  float lead_mass = 0.; 
  if(fatsigJets.size() > 0){
    lead_mass = fatsigJets.at(lead_i)->m();
    if (lead_mass > 50000. && lead_mass < 120000.) pass_jetmass_WZ = true;
    if (lead_mass > 90000. && lead_mass < 145000.) pass_jetmass_H = true;
  }

  if (subJets.size() > 0) pass_1subjet = true;
  if (subJets.size() > 1) pass_2subjet = true;
  if (b_subJets.size() == 2) pass_2btag = true;    

  if(debug) std::cout << "Done monoWZ selection" << std::endl;

  //fill cutflow hist 
  
  //m_hist_mono_pre_MET->Fill(met->met()/1000., m_weight);
  m_hist_mono_pre_MET->Fill(met->met()/1000.);
  m_hist_mono_pre_Nfatjet->Fill(nfatjet);
  for(int i(0); i < nfatjet ;i++){
	m_hist_mono_pre_fatjetm->Fill(fatjets->at(i)->m()/1000., m_weight); 
	m_hist_mono_pre_fatjetpt->Fill(fatjets->at(i)->pt()/1000., m_weight); 
	m_hist_mono_pre_fatjeteta->Fill(fatjets->at(i)->eta()/1000., m_weight); 
  }
  for(int i(0); i < njets ;i++){
	m_hist_mono_pre_jetm->Fill(selJets.at(i)->m()/1000., m_weight); 
	m_hist_mono_pre_jetpt->Fill(selJets.at(i)->pt()/1000., m_weight); 
  	m_hist_mono_pre_jeteta->Fill(selJets.at(i)->eta()/1000., m_weight); 
  }
  m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("pre-selection"), m_weight);
  m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("pre-selection"));
  m_hist_mono_cutflow_Nsigfjet->Fill(fatsigJets.size());
 
  if(pass_met250)
  {
        //m_hist_mono_cutflow_met250_MET->Fill(met->met()/1000., m_weight);
	m_hist_mono_cutflow_met250_MET->Fill(met->met()/1000.);
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("met250"), m_weight);				// fine
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("met350"), m_weight);
	m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("met250"));
 	m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("met350"));		// inserted these otherwise weird stuff happens with yields
  }
  if(pass_met250 && pass_fatjet)
  {
	//m_hist_mono_cutflow_fjet_MET->Fill(met->met()/1000., m_weight);
	m_hist_mono_cutflow_fjet_MET->Fill(met->met()/1000.);									// fine
	m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("one fatjet"), m_weight);
	m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("one fatjet"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto)
  {
        //m_hist_mono_cutflow_lepveto_MET->Fill(met->met()/1000., m_weight);
	m_hist_mono_cutflow_elveto_MET->Fill(met->met()/1000.);								// fine
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("el veto"), m_weight);
	m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("el veto"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto)
  {
        //m_hist_mono_cutflow_photveto_MET->Fill(met->met()/1000., m_weight);
        m_hist_mono_cutflow_muveto_MET->Fill(met->met()/1000.);								// fine
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("mu veto"), m_weight);
	m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("mu veto"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto)
  {
	//m_hist_mono_cutflow_jetveto_MET->Fill(met->met()/1000., m_weight);
	m_hist_mono_cutflow_jetveto_MET->Fill(met->met()/1000.);								// fine
	m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("addjetveto"), m_weight);
	m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("addjetveto"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL)
  {
        //m_hist_mono_cutflow_metjetOL_MET->Fill(met->met()/1000., m_weight);
        m_hist_mono_cutflow_metjetOL_MET->Fill(met->met()/1000.);								// fine
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("metjetOL"), m_weight);
        m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("metjetOL"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut)
  {
	//m_hist_mono_cutflow_met_MET->Fill(met->met()/1000., m_weight);
	m_hist_mono_cutflow_met_MET->Fill(met->met()/1000.);									// fine
	m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("500met"), m_weight);
	m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("500met"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut && pass_jetmass_WZ)	// fine
  { 
	//m_hist_mono_cutflow_mj_MET->Fill(met->met()/1000., m_weight);
	m_hist_mono_cutflow_mjWZ_MET->Fill(met->met()/1000.);
	m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("mJ(WZ)"), m_weight);
	m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("mJ(WZ)"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut && pass_jetmass_H)		// fine
  {
        //m_hist_mono_cutflow_mj_MET->Fill(met->met()/1000., m_weight);
        m_hist_mono_cutflow_mjH_MET->Fill(met->met()/1000.);
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("mJ(H)"), m_weight);
        m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("mJ(H)"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut && pass_jetmass_H && pass_1subjet)  
  {
        m_hist_mono_cutflow_1subjet_MET->Fill(met->met()/1000.);
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("1subjet"), m_weight);
        m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("1subjet"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut && pass_jetmass_H && pass_1subjet && pass_2subjet)  
  {
        m_hist_mono_cutflow_2subjet_MET->Fill(met->met()/1000.);
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("2subjet"), m_weight);
        m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("2subjet"));
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut && pass_jetmass_H && pass_1subjet && pass_2subjet && pass_2btag)  
  {
        m_hist_mono_cutflow_2btag_MET->Fill(met->met()/1000.);
        m_hist_mono_cutflow->Fill(m_hist_mono_cutflow->GetXaxis()->FindBin("2btag"), m_weight);
        m_hist_mono_cutflow_noweight->Fill(m_hist_mono_cutflow_noweight->GetXaxis()->FindBin("2btag"));
  }

  if(debug) std::cout << "Done fill_hist" << std::endl;
 
  //fill n-1 histogram
  if(pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut && pass_jetmass_WZ)
  {
	m_hist_nminusone_met250->Fill(met->met()/1000., m_weight);
  }								
  if(pass_met250 && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut && pass_jetmass_WZ)
  {
	m_hist_nminusone_fjet->Fill(fatsigJets.size(), m_weight);
  }
  if(pass_met250 && pass_fatjet && pass_jetveto && pass_jetOL && pass_metcut && pass_jetmass_WZ)				// split into two with el and mu
  {
	m_hist_nminusone_lepveto->Fill((vetoEls.size()+vetoMus.size()), m_weight);
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetOL && pass_metcut && pass_jetmass_WZ)
  {
	m_hist_nminusone_jetveto->Fill(sigJets.size(), m_weight);
  }				
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_metcut && pass_jetmass_WZ)		// This should be improved to show the overlap between all jets and MET
  {
	m_hist_nminusone_metjetOL->Fill(num_close_jets, m_weight); 
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_jetmass_WZ)
  {
        m_hist_nminusone_met500->Fill(met->met()/1000., m_weight);
  }
  if(pass_met250 && pass_fatjet && pass_elveto && pass_muveto && pass_jetveto && pass_jetOL && pass_metcut)
  {
        m_hist_nminusone_mj->Fill(lead_mass/1000., m_weight);
  }
}


  return EL::StatusCode::SUCCESS;
} 

EL::StatusCode AnalysisReader :: fill_2Lep()
{
  bool debug = false;
  
  m_hist_2lep_CutFlowMu->Fill(m_hist_2lep_CutFlowMu->GetXaxis()->FindBin("All"),m_weight);
  m_hist_2lep_CutFlowMu_noWght->Fill(m_hist_2lep_CutFlowMu_noWght->GetXaxis()->FindBin("All"));

  // data GRL
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo___Nominal").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  if (!m_isMC) {
    bool passGRL=m_superDecorator.get(eventInfo, EventInfoIntProps::passGRL);
    if(!passGRL) return EL::StatusCode::SUCCESS;
  }

  m_hist_2lep_CutFlowMu->Fill(m_hist_2lep_CutFlowMu->GetXaxis()->FindBin("GRL"),m_weight);
  m_hist_2lep_CutFlowMu_noWght->Fill(m_hist_2lep_CutFlowMu_noWght->GetXaxis()->FindBin("GRL"));


  // lookout for events with no z-vertex
  float ZPV=m_superDecorator.get(eventInfo, EventInfoFloatProps::ZPV);
  if (ZPV<-9998) {
    std::cout << " Event without primary vertex " << ZPV  << " " << eventInfo->runNumber() <<  " " << eventInfo->eventNumber() << std::endl;
    return EL::StatusCode::SUCCESS;
  }
  
  m_hist_2lep_CutFlowMu->Fill(m_hist_2lep_CutFlowMu->GetXaxis()->FindBin("ZPV"),m_weight);
  m_hist_2lep_CutFlowMu_noWght->Fill(m_hist_2lep_CutFlowMu_noWght->GetXaxis()->FindBin("ZPV"));

  //----------------------------
  // Jets
  //--------------------------- 

  const xAOD::JetContainer* jets = 0;
  if ( !m_event->retrieve( jets, "AntiKt4EMTopoJets___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve AntiKt4LCTopoJets container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  //  Info("execute()", "  number of jets = %lu", jets->size());

  // loop over the jets in the container
  /*
  xAOD::JetContainer::const_iterator jet_itr = jets->begin();
  xAOD::JetContainer::const_iterator jet_end = jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    Info("execute()", "  jet pt = %.2f GeV", ((*jet_itr)->pt() * 0.001));
    int jor=m_superDecorator.get(*jet_itr,JetIntProps::passOR);
  } // end for loop over jets
  */

  //----------------------------
  // Muons
  //--------------------------- 

  // always need to read nominal, since variations are shallow copies
//  const xAOD::MuonContainer* muonsNom = 0;
//  if ( !m_event->retrieve( muonsNom, "Muons___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
//    Error("execute()", "Failed to retrieve Muons container. Exiting." );
//    return EL::StatusCode::FAILURE;
//  }

  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons___Nominal" ).isSuccess() ){ // retrieve arguments: container type, container key
  //if ( !m_event->retrieve( muons, "Muons___Original" ).isSuccess() ){ // retrieve arguments: container type, container key
  //if ( !m_event->retrieve( muons, "Muons___MUONS_SCALE__1down" ).isSuccess() ){ // retrieve arguments: container type, container key
  //if ( !m_event->retrieve( muons, "Muons___MUONS_SCALE__1up" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }


  // loop over the muons in the container
  /*
  xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
  xAOD::MuonContainer::const_iterator muon_end = muons->end();
  for( ; muon_itr != muon_end; ++muon_itr ) {
    Info("execute()", "  original muon pt = %.2f GeV", ((*muon_itr)->pt() * 0.001));
    // Test a decorator that doesn't exist on CxAOD
    if (!m_superDecorator.get((*muon_itr),MuonFloatProps::z0)) continue;
    // should give Error in <ObjectDecorator::getAuxValue()>: Accessor 'z0' not available.
  } 
  */


  //----------------------------
  // Electrons
  //--------------------------- 
  const xAOD::ElectronContainer* elecs = 0;
  if ( !m_event->retrieve( elecs, "ElectronCollection___Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  

  // loop over the electrons in the container
/*
  xAOD::ElectronContainer::const_iterator elec_itr = elecs->begin();
  xAOD::ElectronContainer::const_iterator elec_end = elecs->end();
  for( ; elec_itr != elec_end; ++elec_itr ) {
    Info("execute()", "  original elec pt = %.2f GeV", ((*elec_itr)->pt() * 0.001)); // just to print out something
  } 
*/

  //---------------------------
  //Photons
  //--------------------------- 
  const xAOD::PhotonContainer* phots = 0;
  if( !m_event->retrieve( phots, "PhotonCollection__Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Photons container. Exiting.");
    return EL::StatusCode::FAILURE;
  }


  //----------------------------
  // Missing Et
  //--------------------------- 
  const xAOD::MissingETContainer* METNominal = 0;
  if ( !m_event->retrieve( METNominal, "MET_RefFinal___Nominal" ).isSuccess() ){
    Error("execute()", "Failed to retrieve Missing Et container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  

  const xAOD::MissingET * met = 0;
  if(METNominal->size() > 0) met = METNominal->at(0);
  /*
  double sumET = met->sumet();
  double MET = met->met();
  std::cout << " sumET, MET " << sumET << " " << MET << std::endl;
  */


  //----------------------------
  // Overlap Removal
  //--------------------------- 
  // Do here for Nominal, systematics
  if(debug)  std::cout << " Overlap Removal "  << std::endl;
  //std::cout << " Run Event " << eventInfo->runNumber() <<  " " << eventInfo->eventNumber() << std::endl;
  m_overlapRemoval->removeOverlap(elecs, phots, muons, jets);
  if(debug)  std::cout << " Overlap Removal (done)"  << std::endl;

  //----------------------------
  // Event Selection
  //--------------------------- 

  if(debug)  std::cout << " Event Slection "  << std::endl;
  bool pass2lep=m_twolepSel->passSelection(eventInfo,
					   met,
					   elecs,
					   phots,
					   muons,
					   jets);

  if(debug)  std::cout << " Passed 2 lepton "  << pass2lep << std::endl;

  ResultVHbb2lep twolepResult=m_twolepSel->result();
  
  if (twolepResult.pass != pass2lep)
    std::cout << " 2 lep ERROR "  << pass2lep << " " << twolepResult.pass << std::endl;

  // select events
  if (!pass2lep) return EL::StatusCode::SUCCESS;

  //
  const xAOD::Electron* el1=twolepResult.el1;
  const xAOD::Electron* el2=twolepResult.el2;
  const xAOD::Muon* mu1=twolepResult.mu1;
  const xAOD::Muon* mu2=twolepResult.mu2;
  std::vector<const xAOD::Jet*> selJets=twolepResult.jets;

  if (pass2lep&&debug) {
    std::cout << " el1 "  << el1 << " el2 " << el2 << std::endl;
    std::cout << " mu1 "  << mu1 << " mu2 " << mu2 << std::endl;
    
    if (el1 && el2) 
      std::cout << " el1 "  << el1->pt() << " el2 " << el2->pt() << std::endl;
    if (mu1 && mu2) 
      std::cout << " mu1 "  << mu1->pt() << " mu2 " << mu2->pt() << std::endl;

    int njet=selJets.size();
    std::cout << " found 2 lep jets " <<  njet << std::endl;
    for (int i(0) ; i<njet ; ++i) {
      std::cout << i << " jet pt " << selJets.at(i)->pt() << std::endl;
    }        
  }



  m_hist_2lep_CutFlowMu->Fill(m_hist_2lep_CutFlowMu->GetXaxis()->FindBin("2Lep"),m_weight);
  m_hist_2lep_CutFlowMu_noWght->Fill(m_hist_2lep_CutFlowMu_noWght->GetXaxis()->FindBin("2Lep"));

  //
  bool isMu=(mu1&&mu2);
  bool isE=(el1&&el2);

  // We only look at mu-mu
  //  if (!(mu1&&mu2)) return EL::StatusCode::SUCCESS;
  if (!isMu&&!isE) return EL::StatusCode::SUCCESS;

  // reject e-mu
  // if (!(mu1&&mu2)&&!(el1&&el2)) return EL::StatusCode::SUCCESS;

  // the derivation/selection keeps high pt 1-jet events, reject here
  int njets=selJets.size();
  if(njets<2) return EL::StatusCode::SUCCESS;

  m_hist_2lep_CutFlowMu->Fill(m_hist_2lep_CutFlowMu->GetXaxis()->FindBin("2Jet"),m_weight);
  m_hist_2lep_CutFlowMu_noWght->Fill(m_hist_2lep_CutFlowMu_noWght->GetXaxis()->FindBin("2Jet"));


  // jets are selected with at least 2 signal jets, but selection vector includes forward jets
  int nsigjets=0;
  std::vector<const xAOD::Jet*> sigJets;
  for (int i(0) ; i<njets ; ++i) {
    if (fabs(selJets.at(i)->eta())>2.5) continue;
    if (selJets.at(i)->pt()<20000) continue;
    nsigjets++;
    sigJets.push_back(selJets.at(i));
  } 

  // require 2 signal jets
  if(nsigjets<2) return EL::StatusCode::SUCCESS;

  // leading jet pt cut. sorting _should_ be done by event selection
  if (sigJets.at(0)->pt()/1000.<45.) return EL::StatusCode::SUCCESS;

  m_hist_2lep_CutFlowMu->Fill(m_hist_2lep_CutFlowMu->GetXaxis()->FindBin("2SigJet"),m_weight);
  m_hist_2lep_CutFlowMu_noWght->Fill(m_hist_2lep_CutFlowMu_noWght->GetXaxis()->FindBin("2SigJet"));

  // and finally cut on MET
  double MET = met->met()/1000.;
  bool metcut=(MET<60);

  // check selection jets with xAOD jets
  int njets2(0);
  for (unsigned int iJet(0) ; iJet < jets->size(); ++iJet) {
    const xAOD::Jet * jet = jets->at(iJet);
    if (!m_superDecorator.get(jet,JetIntProps::passOR)) continue; // should be updated by overlapRemover
    if (m_superDecorator.get(jet,JetIntProps::isVetoJet)) { // applied to keep CxAOD jets already
      if (m_superDecorator.get(jet,JetIntProps::isSignalJet)) // pt>20 GeV, |eta|<2.5
        njets2++;
    }
  }

  if (njets2 != nsigjets) 
    std::cout << "njets2 " << njets2 << " nsigjets " << nsigjets << std::endl;


  // fill mZ
  TLorentzVector l1Vec;
  if(isMu)
    l1Vec.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), mu1->m());
  else
    l1Vec.SetPtEtaPhiM(el1->pt(), el1->eta(), el1->phi(), el1->m());

  TLorentzVector l2Vec;
  if(isMu)
    l2Vec.SetPtEtaPhiM(mu2->pt(), mu2->eta(), mu2->phi(), mu2->m());
  else
    l2Vec.SetPtEtaPhiM(el2->pt(), el2->eta(), el2->phi(), el2->m());

  TLorentzVector ZVec(l1Vec);
  ZVec+=l2Vec;

  if(isMu && metcut)
    m_hist_2lep_Mmumu->Fill(ZVec.M()/1000.,m_weight);
  else if(isE && metcut)
    m_hist_2lep_Mee->Fill(ZVec.M()/1000.,m_weight);

  if (ZVec.M()/1000. < 71.) return EL::StatusCode::SUCCESS;
  if (ZVec.M()/1000. >= 121.) return EL::StatusCode::SUCCESS;

  m_hist_2lep_MET->Fill(MET, m_weight);

  if (!metcut) return EL::StatusCode::SUCCESS;

  m_hist_2lep_PtZ->Fill(ZVec.Pt()/1000, m_weight);
  if (isMu) 
    m_hist_2lep_PtZmumu->Fill(ZVec.Pt()/1000, m_weight);
  else
    m_hist_2lep_PtZee->Fill(ZVec.Pt()/1000, m_weight);

  // primary vertex after cuts
  m_hist_2lep_ZPV->Fill(ZPV, m_weight);
  

  m_hist_2lep_NJets->Fill(nsigjets,m_weight);
  for (int i(0) ; i<nsigjets ; ++i) {
    m_hist_2lep_PtJets->Fill(sigJets.at(i)->pt()/1000.,m_weight);
    m_hist_2lep_EtaJets->Fill(sigJets.at(i)->eta(),m_weight);
  } 
  if (nsigjets>1) 
    m_hist_2lep_DPhiJets->Fill(acos(cos (sigJets.at(0)->phi() - sigJets.at(1)->phi())),m_weight);

  // signal jets
  for (int i(0) ; i<nsigjets ; ++i) {
    // fill btagger histograms
    bTaggerHists(1,sigJets.at(i));
    m_hist_2lep_PtJets->Fill(sigJets.at(i)->pt()/1000.,m_weight);
    m_hist_2lep_EtaJets->Fill(sigJets.at(i)->eta(),m_weight);
  } 


  return EL::StatusCode::SUCCESS;
}


EL::StatusCode AnalysisReader :: fill_VBFa()
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode AnalysisReader :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode AnalysisReader :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
  
  //  delete m_derivEffProvider;
  Info("finalize()", "Processed %i Events", m_eventCounter );

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode AnalysisReader :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  m_hist_mono_eff_MET->Divide(m_hist_mono_cutflow_mjWZ_MET,m_hist_mono_pre_MET,1,1,"B");

  double yields2_WZchan, yields2_Hchan;
  double yields2_All, yields2_pre, yields2_250, yields2_350,  yields2_fjet, yields2_el, yields2_mu, yields2_add, yields2_metOL, yields2_met, yields2_mjWZ, yields2_mjH, yields2_1subjet, yields2_2subjet, yields2_2btag;
  
  yields2_All =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("All"));
  yields2_All *= 5000;
  yields2_pre =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("pre-selection"));
  yields2_pre *= 5000;
  yields2_250 =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("met250"));
  yields2_250 *= 5000;
  yields2_350 =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("met350"));
  yields2_350 *= 5000;
  yields2_fjet =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("one fatjet"));
  yields2_fjet *= 5000;
  yields2_el =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("el veto"));
  yields2_el *= 5000;
  yields2_mu =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("mu veto"));
  yields2_mu *= 5000;
  yields2_add =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("addjetveto"));
  yields2_add *= 5000;
  yields2_metOL = m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("metjetOL"));
  yields2_metOL *= 5000;
  yields2_met =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("500met"));
  yields2_met *= 5000;
  yields2_mjWZ =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("mJ(WZ)"));
  yields2_mjWZ *= 5000;
  yields2_mjH =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("mJ(H)"));
  yields2_mjH *= 5000;
  yields2_1subjet =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("1subjet"));
  yields2_1subjet *= 5000;
  yields2_2subjet =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("2subjet"));
  yields2_2subjet *= 5000;
  yields2_2btag =  m_hist_mono_cutflow->GetBinContent(m_hist_mono_cutflow->GetXaxis()->FindBin("2btag"));
  yields2_2btag *= 5000;
  std::cout << "yields All/" <<  yields2_All << "/pre/" <<  yields2_pre << "/met250/" << yields2_250 << "/met350/" << yields2_350 << "/onefjet/" << yields2_fjet << "/el veto/" << yields2_el << "/mu veto/" << yields2_mu << "/addveto/" << yields2_add << "/metjetOL/" << yields2_metOL << "/500met/" << yields2_met << "/mJ(WZ)/" << yields2_mjWZ << "/mJ(H)/" << yields2_mjH << "/1subjet/" << yields2_1subjet << "/2subjet/" << yields2_2subjet << "/2btag/" << yields2_2btag << std::endl;
  yields2_WZchan = m_hist_mono_cutflow_mjWZ_MET->Integral();
  yields2_WZchan *= 5000;
  yields2_Hchan = m_hist_mono_cutflow_2btag_MET->Integral();
  yields2_Hchan *= 5000;

  std::cout << "yields(5fb^-1) in WZ channel :" <<  yields2_WZchan << std::endl;
  std::cout << "yields(5fb^-1) in H channel :" <<  yields2_Hchan << std::endl;

  int num_presel, num_MET250, num_MET350, num_fj, num_elveto, num_muveto, num_jetveto, num_metOL, num_MET500, num_mjWZ, num_mjH, num_1sj, num_2sj, num_2btag;
  num_presel = m_hist_mono_pre_MET->GetEntries();
  num_MET250 = m_hist_mono_cutflow_met250_MET->GetEntries();
  num_MET350 = m_hist_mono_cutflow_met350_MET->GetEntries();
  num_fj = m_hist_mono_cutflow_fjet_MET->GetEntries();
  num_elveto = m_hist_mono_cutflow_elveto_MET->GetEntries();
  num_muveto = m_hist_mono_cutflow_muveto_MET->GetEntries();
  num_jetveto = m_hist_mono_cutflow_jetveto_MET->GetEntries();
  num_metOL = m_hist_mono_cutflow_metjetOL_MET->GetEntries();
  num_MET500 = m_hist_mono_cutflow_met_MET->GetEntries();
  num_mjWZ = m_hist_mono_cutflow_mjWZ_MET->GetEntries();
  num_mjH = m_hist_mono_cutflow_mjH_MET->GetEntries();
  num_1sj = m_hist_mono_cutflow_1subjet_MET->GetEntries();
  num_2sj = m_hist_mono_cutflow_2subjet_MET->GetEntries();
  num_2btag = m_hist_mono_cutflow_2btag_MET->GetEntries();
  cout << "num_presel = " << num_presel << endl;
  cout << "num_MET250 = " << num_MET250 << endl;
  cout << "num_MET350 = " << num_MET350 << endl;
  cout << "num_fj = " << num_fj << endl;
  cout << "num_elveto = " << num_elveto << endl;
  cout << "num_muveto = " << num_muveto << endl;
  cout << "num_jetveto = " << num_jetveto << endl;
  cout << "num_metOL = " << num_metOL << endl;
  cout << "num_MET500 = " << num_MET500 << endl;
  cout << "num_mjWZ = " << num_mjWZ << endl;
  cout << "num_mjWZ = " << num_mjH << endl;
  cout << "num_1sj = " << num_1sj << endl;
  cout << "num_2sj = " << num_2sj << endl;
  cout << "num_2btag = " << num_2btag << endl;
 
  return EL::StatusCode::SUCCESS;
}
