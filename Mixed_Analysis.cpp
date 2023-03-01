#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

// include user defined histograms and auxiliary macros
#include "Histodef.cpp"
#include "Auxiliary.cpp"


// correctionlib
#include "correction.h"
using namespace std;
using correction::CorrectionSet;

#define MAX_ARRAY_SIZE 128
#define GEN_MAX_ARRAY_SIZE 1024

// function to calculate the weight for each event
// the weight is calculated as the product of luminosity and cross section of the process times the genWeight,
// LATER TO BE divided by the number of generated events OF ALL FILES OF THE DATASET(S)
double getWeight(double luminosity, double crossSection, Float_t genWeight, double SumWeights)
{
    return (luminosity * crossSection * genWeight); // / SumWeights;
}

void Mixed_Analysis(string inputFile, string ofile, double crossSection = -1, double IntLuminosity = 59.827879506, bool Signal = false)
{

    if (crossSection < 0. || IntLuminosity < 0.)
    {
        std::cout << "WARNING: crossection " << crossSection << " and Integrated luminosity " << IntLuminosity << endl;
    }

cout<<"Call completed!"<<endl;

    TFile *fin = TFile::Open(inputFile.c_str());
    TTree *trun = static_cast<TTree *>(fin->Get("Runs"));
    Long64_t genEventCount;
    Double_t genEventSumw;
    trun->SetBranchStatus("*", 0);
    trun->SetBranchStatus("genEventSumw", 1);
    trun->SetBranchStatus("genEventCount", 1);
    trun->SetBranchAddress("genEventSumw", &genEventSumw);
    trun->SetBranchAddress("genEventCount", &genEventCount);


    trun->GetEntry(0);

    TTree *tin = static_cast<TTree *>(fin->Get("Events"));

    // Set all branches to 0
    tin->SetBranchStatus("*", 0);
    // get the pt
    Float_t Tau_pt[MAX_ARRAY_SIZE], Electron_pt[MAX_ARRAY_SIZE], Jet_pt[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_pt", 1);
    tin->SetBranchAddress("Electron_pt", &Electron_pt);
    tin->SetBranchStatus("Tau_pt", 1);
    tin->SetBranchAddress("Tau_pt", &Tau_pt);
    tin->SetBranchStatus("Jet_pt", 1);
    tin->SetBranchAddress("Jet_pt", &Jet_pt);
    // get the number of Taus, electrons
    UInt_t nTau, nElectron;
    tin->SetBranchStatus("nElectron", 1);
    tin->SetBranchAddress("nElectron", &nElectron);
    tin->SetBranchStatus("nTau", 1);
    tin->SetBranchAddress("nTau", &nTau);
    // get the eta
    Float_t Tau_eta[MAX_ARRAY_SIZE], Electron_eta[MAX_ARRAY_SIZE], Jet_eta[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_eta", 1);
    tin->SetBranchAddress("Electron_eta", &Electron_eta);
    tin->SetBranchStatus("Tau_eta", 1);
    tin->SetBranchAddress("Tau_eta", &Tau_eta);
    tin->SetBranchStatus("Jet_eta", 1);
    tin->SetBranchAddress("Jet_eta", &Jet_eta);
    // get the phi
    Float_t Tau_phi[MAX_ARRAY_SIZE], Electron_phi[MAX_ARRAY_SIZE], Jet_phi[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_phi", 1);
    tin->SetBranchAddress("Electron_phi", &Electron_phi);
    tin->SetBranchStatus("Tau_phi", 1);
    tin->SetBranchAddress("Tau_phi", &Tau_phi);
    tin->SetBranchStatus("Jet_phi", 1);
    tin->SetBranchAddress("Jet_phi", &Jet_phi);
    // get the mass
    Float_t Tau_mass[MAX_ARRAY_SIZE], Electron_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];//, Tau_energy[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_mass", 1);
    tin->SetBranchAddress("Electron_mass", &Electron_mass);
    tin->SetBranchStatus("Tau_mass", 1);
    tin->SetBranchAddress("Tau_mass", &Tau_mass);
    /*tin->SetBranchStatus("Tau_energy", 1);
    tin->SetBranchAddress("Tau_energy", &Tau_energy);*/
    tin->SetBranchStatus("Jet_mass", 1);
    tin->SetBranchAddress("Jet_mass", &Jet_mass);

    // get gen quantities
    Int_t Tau_genPartIdx[MAX_ARRAY_SIZE], Electron_genPartIdx[MAX_ARRAY_SIZE];
    Int_t GenPart_pdgId[GEN_MAX_ARRAY_SIZE], GenPart_genPartIdxMother[GEN_MAX_ARRAY_SIZE], Jet_genJetIdx[MAX_ARRAY_SIZE],GenPart_statusFlags[GEN_MAX_ARRAY_SIZE];
    UChar_t Tau_genPartFlav[MAX_ARRAY_SIZE], Electron_genPartFlav[MAX_ARRAY_SIZE];
    Float_t Electron_ip3d[MAX_ARRAY_SIZE], Electron_sip3d[MAX_ARRAY_SIZE], Electron_dxy[MAX_ARRAY_SIZE], Electron_dz[MAX_ARRAY_SIZE],GenPart_pt[GEN_MAX_ARRAY_SIZE];
    UInt_t nGenPart;
    tin->SetBranchStatus("Electron_genPartIdx", 1);
    tin->SetBranchStatus("Electron_genPartFlav", 1);
    tin->SetBranchStatus("Tau_genPartIdx", 1);
    tin->SetBranchStatus("Tau_genPartFlav", 1);
    tin->SetBranchStatus("GenPart_pdgId", 1);
    tin->SetBranchStatus("GenPart_genPartIdxMother", 1);
    tin->SetBranchStatus("nGenPart", 1);
    tin->SetBranchStatus("Jet_genJetIdx",1);
    tin->SetBranchStatus("GenPart_pt", 1);
    tin->SetBranchStatus("GenPart_statusFlags", 1);
    tin->SetBranchAddress("nGenPart", &nGenPart);
    tin->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx);
    tin->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav);
    tin->SetBranchAddress("Tau_genPartIdx", &Tau_genPartIdx);
    tin->SetBranchAddress("Tau_genPartFlav", &Tau_genPartFlav);
    tin->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    tin->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    tin->SetBranchAddress("Jet_genJetIdx",&Jet_genJetIdx);
    tin->SetBranchAddress("GenPart_pt", &GenPart_pt);
    tin->SetBranchAddress("GenPart_statusFlags",&GenPart_statusFlags);

    tin->SetBranchStatus("Electron_ip3d", 1);
    tin->SetBranchStatus("Electron_sip3d", 1);
    tin->SetBranchStatus("Electron_dxy", 1);
    tin->SetBranchStatus("Electron_dz", 1);
    tin->SetBranchAddress("Electron_dz", &Electron_dz);
    tin->SetBranchAddress("Electron_ip3d", &Electron_ip3d);
    tin->SetBranchAddress("Electron_sip3d", &Electron_sip3d);
    tin->SetBranchAddress("Electron_dxy", &Electron_dxy);
    // collect the trigger information
    Bool_t HLT_Ele32_WPTight_Gsf;
    tin->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1);
    tin->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);

    // collect the triggger Ids
    Int_t Tau_charge[MAX_ARRAY_SIZE], Electron_charge[MAX_ARRAY_SIZE];
    Bool_t Electron_mvaFall17V2Iso_WP90[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Tau_charge", 1);
    tin->SetBranchStatus("Electron_charge", 1);
    tin->SetBranchStatus("Electron_mvaFall17V2Iso_WP90", 1);
    tin->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", &Electron_mvaFall17V2Iso_WP90);
    tin->SetBranchAddress("Tau_charge", &Tau_charge);
    tin->SetBranchAddress("Electron_charge", &Electron_charge);

    // Jet tagging and ID, FlavB is the recomended one, DeepB was used by Anup
    Float_t Jet_btagDeepFlavB[MAX_ARRAY_SIZE], Jet_btagDeepB[MAX_ARRAY_SIZE];
    UInt_t nJet;
    Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE],Jet_hadronFlavour[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Jet_btagDeepB", 1);
    tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
    tin->SetBranchStatus("nJet", 1);
    tin->SetBranchStatus("Jet_jetId", 1);
    tin->SetBranchStatus("Jet_puId", 1);
    tin->SetBranchStatus("Jet_hadronFlavour", 1);
    tin->SetBranchAddress("nJet", &nJet);
    tin->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);
    tin->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);   
    tin->SetBranchAddress("Jet_jetId", &Jet_jetId);
    tin->SetBranchAddress("Jet_puId", &Jet_puId);
    tin->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);

    UChar_t Tau_idDeepTau2017v2p1VSmu[MAX_ARRAY_SIZE],Tau_idDeepTau2017v2p1VSjet[MAX_ARRAY_SIZE], Tau_idDeepTau2017v2p1VSe[MAX_ARRAY_SIZE];
    Int_t Tau_decayMode[MAX_ARRAY_SIZE];

    tin->SetBranchStatus("Tau_idDeepTau2017v2p1VSmu", 1);
    tin->SetBranchStatus("Tau_idDeepTau2017v2p1VSjet", 1);
    tin->SetBranchStatus("Tau_idDeepTau2017v2p1VSe", 1);
    tin->SetBranchStatus("Tau_decayMode", 1);
    tin->SetBranchAddress("Tau_decayMode", &Tau_decayMode);
    tin->SetBranchAddress("Tau_idDeepTau2017v2p1VSmu", &Tau_idDeepTau2017v2p1VSmu);
    tin->SetBranchAddress("Tau_idDeepTau2017v2p1VSjet", &Tau_idDeepTau2017v2p1VSjet);   
    tin->SetBranchAddress("Tau_idDeepTau2017v2p1VSe", &Tau_idDeepTau2017v2p1VSe);
 
    //L1
    Float_t L1PreFiringWeight_Nom;
    tin->SetBranchStatus("L1PreFiringWeight_Nom", 1);
    tin->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom);

    // pu stuff
    Float_t N_pu_vertices;
    tin->SetBranchStatus("Pileup_nTrueInt", 1);
    tin->SetBranchAddress("Pileup_nTrueInt", &N_pu_vertices);

    // gen weight
    Float_t genWeight;
    tin->SetBranchStatus("genWeight", 1);
    tin->SetBranchAddress("genWeight", &genWeight);

    int non_matching_Tau = 0, non_matching_electron = 0;
    int n_dropped = 0;
    int trigger_dropped = 0;
    UInt_t nEv = tin->GetEntries();
    unsigned int n_events = nEv;
    TLorentzVector *Tau_p4 = new TLorentzVector();
    TLorentzVector *Electron_p4 = new TLorentzVector();
    TLorentzVector *MainBjet_p4 = new TLorentzVector();
    TLorentzVector *OppositeBjet_p4 = new TLorentzVector();
cout<<"Corr"<<endl;
    // open correctionfiles
    
    string Tau_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/tau.json.gz";
    string electron_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/electron.json.gz";
    string jets_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/jet_jmar.json";
    string b_tag_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/btagging.json.gz";
    string pileup_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/puWeights.json.gz";

    
    auto Tau_c_set = CorrectionSet::from_file(Tau_json);
    auto ele_c_set = CorrectionSet::from_file(electron_json);
    auto jet_c_set = CorrectionSet::from_file(jets_json);
    auto btag_c_set = CorrectionSet::from_file(b_tag_json);
    auto pu_c_set = CorrectionSet::from_file(pileup_json);

    auto Tau_Escale = Tau_c_set->at("tau_energy_scale");
    auto Tau_idvse = Tau_c_set->at("DeepTau2017v2p1VSe");
    auto Tau_idvsmu = Tau_c_set->at("DeepTau2017v2p1VSmu");
    auto Tau_idvsjet = Tau_c_set->at("DeepTau2017v2p1VSjet");
    //auto Tau_trig = Tau_c_set->at("");
    auto electron_id = ele_c_set->at("UL-Electron-ID-SF");
    auto jet_pu = jet_c_set->at("PUJetID_eff");
    auto b_tag = btag_c_set->at("deepJet_mujets");
    auto b_mistag= btag_c_set->at("deepJet_incl"); //only for light jets
    auto pu_correction = pu_c_set->at("Collisions18_UltraLegacy_goldenJSON");
    
    TFile *fecorr_trig = new TFile("/afs/cern.ch/user/g/gdamolin/public/Riccardo_egammaTriggerEfficiency_2018_20200422.root");
    TH2F * EleTrigHisto= static_cast<TH2F *>(fecorr_trig->Get("EGamma_SF2D"));

    TFile *fb_eff = new TFile("/afs/cern.ch/user/g/gdamolin/public/Beff_puLoose.root");
    TH2D * l_eff= static_cast<TH2D *>(fb_eff->Get("l_jets_tagged")); 
    TH2D * c_eff= static_cast<TH2D *>(fb_eff->Get("c_jets_tagged")); 
    TH2D * b_eff= static_cast<TH2D *>(fb_eff->Get("b_jets_tagged")); 
   
    // save the histograms in a new File
    // allow pt, inv mass, and eta to be stored in a Branch
    Float_t leading_lepton_pt, invMass, electron_eta, electron_pt, tau_eta, tau_pt;
    Float_t Tau_eta_from_W, Tau_pt_from_W, electron_eta_from_W, electron_pt_from_W;
    float Weight;
	cout<<"fout"<<endl;
    TFile *fout = new TFile(ofile.c_str(), "RECREATE");
    
    // create a new tree for the output
    TTree *tout = new TTree("tout", "tout");
    TTree *trun_out = new TTree("Run_out", "Run_out");
    // set the branches for the output tree
    tout->Branch("leading_lepton_pt", &leading_lepton_pt);
    tout->Branch("invMass", &invMass);
    tout->Branch("electron_eta", &electron_eta);
    tout->Branch("electron_pt", &electron_pt);
    tout->Branch("tau_eta", &tau_eta);
    tout->Branch("tau_pt", &tau_pt);
    tout->Branch("Weight", &Weight);

    int Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0, Nprongs=0;
    float dR_muE, dR_mujet, dR_ejet, dR_allJets, dR_lbJets, dR_mbJets, Apl_allJets, Apl_lbJets, Apl_mbJets, Phi_allJets, Phi_lbJets, Phi_mbJets, PTbjet,Acopl_etau;

 bool From2Taus=false, FromTau=false;

    tout->Branch("dR_mue", &dR_muE);
    tout->Branch("dR_mujet", &dR_mujet);
    tout->Branch("dR_ejet", &dR_ejet);
    tout->Branch("dR_allJets", &dR_allJets);
    tout->Branch("dR_lbJets", &dR_lbJets);
    tout->Branch("dR_mbJets", &dR_mbJets);
    tout->Branch("Apl_lbJets", &Apl_lbJets);
    tout->Branch("Apl_allJets", &Apl_allJets);
    tout->Branch("Apl_mbJets", &Apl_mbJets);
    tout->Branch("Phi_allJets", &Phi_allJets);
    tout->Branch("Phi_lbJets", &Phi_lbJets);
    tout->Branch("Phi_mbJets", &Phi_mbJets);
    tout->Branch("PTbjet", &PTbjet);
    tout->Branch("Nloose", &Nloose);
    tout->Branch("Nmedium", &Nmedium);
    tout->Branch("Ntight", &Ntight);
    tout->Branch("JetNotB", &JetsNotB);
    tout->Branch("Acopl_etau", &Acopl_etau);
    tout->Branch("Nprongs", &Nprongs);
    tout->Branch("FromTau", &FromTau);
    tout->Branch("From2Taus", &From2Taus);

    trun_out->Branch("genEventSumw", &genEventSumw);
    trun_out->Branch("IntLumi", &IntLuminosity);
    trun_out->Branch("xs", &crossSection);
    trun_out->Branch("nEv", &n_events);

    trun_out->Fill(); // we already called trun->GetEntry(0);

    /*size_t found = ofile.find_last_of("/");
    string oname=ofile.substr(found+1);
    string path=ofile.substr(0,found);
    string Tauname=path+"/Tau_"+oname;
    TFile *foutT = new TFile(Tauname.c_str(), "RECREATE");
    TTree *toutT = new TTree("toutT", "toutT");
    TTree *trun_outT = new TTree("Run_outT", "Run_outT");
    if (Signal) {
	toutT->Branch("leading_lepton_pt", &leading_lepton_pt);
    	toutT->Branch("invMass", &invMass);
    	toutT->Branch("electron_eta", &electron_eta);
   	toutT->Branch("electron_pt", &electron_pt);
   	toutT->Branch("tau_eta", &tau_eta);
    	toutT->Branch("tau_pt", &tau_pt);
    	toutT->Branch("Weight", &Weight);
	toutT->Branch("FromTau", &FromTau);
	toutT->Branch("From2Taus", &From2Taus);
	toutT->Branch("dR_mue", &dR_muE);
	toutT->Branch("dR_mujet", &dR_mujet);
	toutT->Branch("dR_ejet", &dR_ejet);
	    toutT->Branch("dR_allJets", &dR_allJets);
	    toutT->Branch("dR_lbJets", &dR_lbJets);
	    toutT->Branch("dR_mbJets", &dR_mbJets);
	    toutT->Branch("Apl_lbJets", &Apl_lbJets);
	    toutT->Branch("Apl_allJets", &Apl_allJets);
	    toutT->Branch("Apl_mbJets", &Apl_mbJets);
	    toutT->Branch("Phi_allJets", &Phi_allJets);
	    toutT->Branch("Phi_lbJets", &Phi_lbJets);
	    toutT->Branch("Phi_mbJets", &Phi_mbJets);
	    toutT->Branch("PTbjet", &PTbjet);
	    toutT->Branch("Nloose", &Nloose);
	    toutT->Branch("Nmedium", &Nmedium);
	    toutT->Branch("Ntight", &Ntight);
	    toutT->Branch("JetNotB", &JetsNotB);
	    toutT->Branch("Acopl_etau", &Acopl_etau);
	toutT->Branch("FromTau", &FromTau);
        toutT->Branch("From2Taus", &From2Taus);


	trun_outT->Branch("genEventSumw", &genEventSumw);
        trun_outT->Branch("IntLumi", &IntLuminosity);
        trun_outT->Branch("xs", &crossSection);
        trun_outT->Branch("nEv", &n_events);
	trun_outT->Fill();
	}*/
    fout->cd();
    
    //#pragma omp parallel for
    for (UInt_t i = 0; i <nEv; i++){
        tin->GetEntry(i); //cout<<"kill me"<<endl;
        if (i % 100000 == 0)
            cout << "Processing entry " << i << " of " << nEv << endl;
        // apply triggers
//cout<<"Skidding"<<endl;
        if (!( HLT_Ele32_WPTight_Gsf)){
            trigger_dropped++;
            continue;
        }

	bool OneProng=false, ThreeProng=false;
        Int_t Tau_idx = -1;
        for (UInt_t j = 0; j < nTau; j++){
		if (Tau_decayMode[j]>2 && Tau_decayMode[j]<10) continue;
             double ScaleE=Tau_Escale->evaluate({Tau_pt[j],abs(Tau_eta[j]),Tau_decayMode[j],Tau_genPartFlav[j],"DeepTau2017v2p1","nom"});
             Tau_p4->SetPtEtaPhiM(Tau_pt[j]*ScaleE, Tau_eta[j], Tau_phi[j], Tau_mass[j]*ScaleE);
            
            if ((Tau_p4->Pt()>22. && abs(Tau_eta[j])<2.3)&&(Tau_idDeepTau2017v2p1VSe[j]>=8&& Tau_idDeepTau2017v2p1VSmu[j]>=8 && Tau_idDeepTau2017v2p1VSjet[j]>=32)){ //Loose e- T mu T jet
		if (Tau_decayMode[j]<=2) {OneProng=true;}
		if (Tau_decayMode[j]>=10) {ThreeProng=true;}
		if (!(OneProng || ThreeProng)) {continue;}
                Tau_idx = j;
                break;
            }
        }
        if (Tau_idx==-1)  {
            n_dropped++;
            continue;
        }
        Weight = getWeight(IntLuminosity, crossSection, genWeight, genEventSumw);
	Weight *=  getTopPtWeight(GenPart_pdgId,GenPart_statusFlags,GenPart_pt,nGenPart);
	double Weight2=Weight;
	Weight*=L1PreFiringWeight_Nom;
        Weight *= pu_correction->evaluate({N_pu_vertices, "nominal"});
			
	Weight *=Tau_idvse->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Loose","nom"}); //Loose instead of VL
	int grr=Tau_genPartFlav[Tau_idx];
	cout<<"For eles, I got flavor "<< grr << " and I multiply by the weight "<< Tau_idvse->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Loose","nom"})<<endl;
	Weight *=Tau_idvsmu->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Tight","nom"});
	cout<<"For Mus, I got flavor "<< grr << " and I multiply by the weight "<< Tau_idvse->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Loose","nom"})<<endl;
	Weight *=Tau_idvsjet->evaluate({Tau_p4->Pt(), Tau_decayMode[Tau_idx],Tau_genPartFlav[Tau_idx],"Tight","nom","pt"});
  
        Int_t electron_idx = -1;
        for (UInt_t j = 0; j < nElectron; j++)
        {
            if ((Electron_pt[j] > 35 && abs(Electron_eta[j]) < 2.4 && Electron_mvaFall17V2Iso_WP90[j]  && abs(Electron_dxy[j])<0.2 && abs(Electron_dz[j])<0.5))
            {
		if((abs(Electron_eta[j])>1.44) && (abs(Electron_eta[j])<1.57)) {continue;}

                Electron_p4->SetPtEtaPhiM(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);

                if 	(Electron_p4->DeltaR(*Tau_p4) < 0.4)  {continue;}
                else	{electron_idx = j;   break;}
            }
        }//end ele for
        if (electron_idx==-1) {
            n_dropped++;
            continue;
        }
        /*if(Signal){
		bool secondistau=isFromTau(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Electron_genPartIdx[electron_idx]);
		if (secondistau){FromTau=true;}
		else {FromTau=false;}
		}*/

        Weight *= electron_id->evaluate({"2018", "sf", "wp90iso", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]}); 
        Weight *= electron_id->evaluate({"2018", "sf", "RecoAbove20", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]});
		
        if(HLT_Ele32_WPTight_Gsf) {
            //retrieve Histo
            int bin = EleTrigHisto->FindBin(Electron_eta[electron_idx],Electron_pt[electron_idx]);
            float temp= EleTrigHisto->GetBinContent(bin);
            Weight*=temp;
            }

        bool selection = ((Tau_idx > -1) && (electron_idx > -1));
        // check the seleected objects for opposite charge
        selection = selection && (Tau_charge[Tau_idx] * Electron_charge[electron_idx]) < 0;
        // the tight working point is 0.71, medium 0.2783, loose 0.0490
        Float_t jet_btag_deepFlav_wp = 0.2783;
        bool one_Bjet = false;
        int id_m_jet = -1;
	int njets=0;
        Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0;
	//vectors for applying b-tag corrections
	vector<int> njet_in_collection;
	vector<int> flavor;
	vector<bool> tagged;
	double t_weight=1.;
        for (size_t j = 0; j < nJet; j++)
        {
            if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6)){
	    TLorentzVector *Tjet_p4 = new TLorentzVector();
	    Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	    if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	    else {delete Tjet_p4;}
            //correction for pileupID
            int MC_pu = Jet_genJetIdx[j];
            float tempSF=1.,tempEff;
            //if is pileUpjet
            if (MC_pu<0 ) {
		tempSF=1.;
            	tempEff= 0;
            	}
            //if is truly a jet
            else { if (Jet_pt[j]<=50){
            	     tempSF= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"nom", "L"});
            	     tempEff= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"MCEff", "L"});
		    }
            	}
            bool passesPUID=(Jet_puId[j]>=4);
		
            if(!(Jet_pt[j]>50 || passesPUID ))	{t_weight*=(1-tempSF*tempEff)/(1-tempEff); }
            if((Jet_pt[j]>50 || passesPUID)) { 
             if(Jet_pt[j]<=50) t_weight*=tempSF; //else you are in pT>50 case: apply no sf
              //correction for b-tag
              njet_in_collection.push_back(j);
              flavor.push_back(abs(Jet_hadronFlavour[j]));
              tagged.push_back((Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp));
		njets++;
        
	      if (Jet_btagDeepFlavB[j] < 0.0490) JetsNotB++;
	      if (Jet_btagDeepFlavB[j] > 0.0490)
                   Nloose++;
              if (Jet_btagDeepFlavB[j] > 0.2783)
                   Nmedium++;
              if (Jet_btagDeepFlavB[j] > 0.71)
                   Ntight++;
            
              if (Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp){
                 if (!one_Bjet){
                      MainBjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
                      OppositeBjet_p4->SetPtEtaPhiM(Jet_pt[j], -1 * Jet_eta[j], InvertPhi(Jet_phi[j]), Jet_mass[j]);

                      if (MainBjet_p4->DeltaR(*Tau_p4) > 0.4 && MainBjet_p4->DeltaR(*Electron_p4) > 0.4){
                      	one_Bjet = true;
                      	id_m_jet = j;
			 }
                      }
               }
            }//end if(jetpt>50 !!puid==7)
          }//end kinematic if
        }//end for
          //corrections of jets already applied 

        Weight*=t_weight; 
            
	for(int jj=0;jj<flavor.size();jj++){
		int convflav=flavor[jj];
		if (flavor[jj]<4) convflav==0;
		if (!(convflav==0 || convflav==4 || convflav==5)) {cout<<"Something weird in the flavor of jet"<<endl;}
		if(tagged[jj]){
			if (convflav!=0) 
				Weight *= b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			else  Weight *= b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			continue;
			}

		//if not tagged
		if(!tagged[jj]) {
			double Eff=1.;
			double SF=1;
			if (convflav!=0) SF=b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			else SF=b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			
			//Get Eff
			if(convflav==0) {
				int bin =l_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=l_eff->GetBinContent(bin);
				}
			if(convflav==4) {
				int bin =c_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=c_eff->GetBinContent(bin);
				}
			if(convflav==5) {
				int bin =b_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=b_eff->GetBinContent(bin);
				}
			
			Weight*=(1-SF*Eff)/(1-Eff);
			}
		
		}
        //filling before jet selections
	
        selection = selection && (one_Bjet);
        if (!selection){
            n_dropped++;
            continue;
        }
	Acopl_etau=M_PI-(Electron_p4->DeltaPhi(*Tau_p4));
	if(OneProng){
		h_LooseJets_p1->Fill(Nloose, Weight);
		h_MediumJets_p1->Fill(Nmedium, Weight);
		h_TightJets_p1->Fill(Ntight, Weight);
		h_acopla_etau_p1->Fill(Acopl_etau,Weight);
		h_NJets_p1->Fill(njets,Weight);
		}
	if(ThreeProng){
		h_LooseJets_p3->Fill(Nloose, Weight);
		h_MediumJets_p3->Fill(Nmedium, Weight);
		h_TightJets_p3->Fill(Ntight, Weight);
		h_acopla_etau_p3->Fill(Acopl_etau,Weight);
		h_NJets_p3->Fill(njets,Weight);
		}

	if(OneProng) {Nprongs=1;}
	if(ThreeProng) {Nprongs=3;}
        PTbjet = MainBjet_p4->Pt();

        dR_mujet = Tau_p4->DeltaR(*MainBjet_p4);
        dR_ejet = Electron_p4->DeltaR(*MainBjet_p4);
        dR_muE = Tau_p4->DeltaR(*Electron_p4);

	// fill the tree
        tau_pt = Tau_p4->Pt();
        tau_eta = Tau_p4->Eta();
        electron_pt = Electron_pt[electron_idx];
        electron_eta = Electron_eta[electron_idx];
        // check whether Tau or electron is the leading one
        if (Tau_p4->Pt() > Electron_p4->Pt()) {leading_lepton_pt = Tau_p4->Pt();}
        else	{leading_lepton_pt = Electron_p4->Pt();}

        if(OneProng) {
		h_leading_lepton_pt_p1->Fill(leading_lepton_pt,Weight2);
		h_leading_lepton_pt_weighted_p1->Fill(leading_lepton_pt, Weight);
		h_Tau_pt_p1->Fill(tau_pt,Weight2);
		h_Tau_eta_p1->Fill(tau_eta,Weight2);
		h_Electron_pt_p1->Fill(electron_pt,Weight2);
		h_Electron_eta_p1->Fill(electron_eta,Weight2);
		h_Tau_pt_weighted_p1->Fill(tau_pt, Weight);
		h_Tau_eta_weighted_p1->Fill(tau_eta, Weight);
		h_Electron_pt_weighted_p1->Fill(electron_pt, Weight);
		h_Electron_eta_weighted_p1->Fill(electron_eta, Weight);
		b_pt_p1->Fill(MainBjet_p4->Pt(), Weight);
   		jethole_p1->Fill(MainBjet_p4->Eta(),MainBjet_p4->Phi(),Weight);
   		tauhole_p1->Fill(Tau_p4->Eta(),Tau_p4->Phi(),Weight);
   		ehole_p1->Fill(Electron_p4->Phi(),Electron_p4->Phi(),Weight);
		}

	 if(ThreeProng) {
		h_leading_lepton_pt_p3->Fill(leading_lepton_pt,Weight2);
		h_leading_lepton_pt_weighted_p3->Fill(leading_lepton_pt, Weight);
		h_Tau_pt_p3->Fill(tau_pt,Weight2);
		h_Tau_eta_p3->Fill(tau_eta,Weight2);
		h_Electron_pt_p3->Fill(electron_pt,Weight2);
		h_Electron_eta_p3->Fill(electron_eta,Weight2);
		h_Tau_pt_weighted_p3->Fill(tau_pt, Weight);
		h_Tau_eta_weighted_p3->Fill(tau_eta, Weight);
		h_Electron_pt_weighted_p3->Fill(electron_pt, Weight);
		h_Electron_eta_weighted_p3->Fill(electron_eta, Weight);
		b_pt_p3->Fill(MainBjet_p4->Pt(), Weight);
   		jethole_p3->Fill(MainBjet_p4->Eta(),MainBjet_p4->Phi(),Weight);
   		tauhole_p3->Fill(Tau_p4->Eta(),Tau_p4->Phi(),Weight);
   		ehole_p3->Fill(Electron_p4->Phi(),Electron_p4->Phi(),Weight);
		}
        // only for signal
        /*if (Signal)
        {
            // cross check which index the objects have that actually originate from the W
            size_t nTau_p4 = 0, nElectron_p4 = 0;
            for (UInt_t j = 0; j < nTau; j++)
            {
                // match the Tau to the PID of the W boson (PID=24)
                // printMCTree(nGenPart, GenPart_pdgId,GenPart_genPartIdxMother, Tau_genPartIdx[j]);
                if (isFromW(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Tau_genPartIdx[j]))
                {
		double ScaleE=1;
		if (Tau_decayMode[j]<3 && Tau_decayMode[j]>9 ) ScaleE=Tau_Escale->evaluate({Tau_pt[j],abs(Tau_eta[j]),Tau_decayMode[j],Tau_genPartFlav[j],"DeepTau2017v2p1","nom"});
                    Tau_pt_from_W = Tau_pt[j]*ScaleE;
                    Tau_eta_from_W = Tau_eta[j];
                    h_Tau_pt_from_W->Fill(Tau_pt_from_W);
                    h_Tau_eta_from_W->Fill(Tau_eta_from_W);
                    h_Tau_pt_weighted_from_W->Fill(Tau_pt_from_W, Weight);
                    h_Tau_eta_weighted_from_W->Fill(Tau_eta_from_W, Weight);
                    if (Tau_idx != j)
                        non_matching_Tau++;
                }
            }

            for (UInt_t j = 0; j < nElectron; j++)
            {
                if (isFromW(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Electron_genPartIdx[j]))
                {
                    electron_pt_from_W = Electron_pt[j];
                    electron_eta_from_W = Electron_eta[j];
                    h_Electron_pt_from_W->Fill(electron_pt_from_W);
                    h_Electron_eta_from_W->Fill(electron_eta_from_W);
                    h_Electron_pt_weighted_from_W->Fill(electron_pt_from_W, Weight);
                    h_Electron_eta_weighted_from_W->Fill(electron_eta_from_W, Weight);
                    if (electron_idx != j)
                        non_matching_electron++;
                }
            }
        }*/
        // END only for signal

        dR_allJets = 999, dR_lbJets = 999, dR_mbJets = 999;
        Apl_allJets = 1.1, Apl_lbJets = 1.1, Apl_mbJets = 1.1;
        for (size_t j = 0; j < nJet; j++)
        {  
          if (j == id_m_jet)
                continue;
	 TLorentzVector *Tjet_p4 = new TLorentzVector();
	 Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	 if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	 
         if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || (Jet_puId[j]>=4))){
            TLorentzVector *tempJet = Tjet_p4;
            double temp = OppositeBjet_p4->DeltaR(*tempJet);

            TVector3 A(tempJet->X(), tempJet->Y(), tempJet->Z());
            TVector3 B(MainBjet_p4->X(), MainBjet_p4->Y(), MainBjet_p4->Z());

            double tempApl = A.Dot(B) / (A.Mag() * B.Mag());

            if (temp < dR_allJets){ dR_allJets = temp;}
            if (tempApl < Apl_allJets){Apl_allJets = tempApl;}

            if (Jet_btagDeepFlavB[j] > 0.0490){
                if (temp < dR_lbJets) {dR_lbJets = temp;}
                if (tempApl < Apl_lbJets) {Apl_lbJets = tempApl;}
            }
            if (Jet_btagDeepFlavB[j] > 0.2783){
		if (temp < dR_mbJets) {dR_mbJets = temp; }
                if (tempApl < Apl_mbJets) {Apl_mbJets = tempApl;}
            }
            delete Tjet_p4;
         }//end if
        } //end for

        // dphi
        Phi_allJets = 999, Phi_lbJets = 999, Phi_mbJets = 999;
        for (size_t j = 0; j < nJet; j++){ 
          if (j == id_m_jet) continue;
	  TLorentzVector *Tjet_p4 = new TLorentzVector();
	  Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	  if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
          else{delete Tjet_p4;}
          if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || (Jet_puId[j]>=4))){
            double temp = Jet_phi[j] - OppositeBjet_p4->Phi();
            if (temp < -1 * M_PI) temp += 2 * M_PI;
            if (temp > M_PI) temp -= 2 * M_PI;
            if (temp < 0) temp *= (-1);

            if (temp < Phi_allJets)   {Phi_allJets = temp; }
            if ((Jet_btagDeepFlavB[j] > 0.0490) && (temp < Phi_lbJets))   {Phi_lbJets = temp; }
            if ((Jet_btagDeepFlavB[j] > 0.2783) && (temp < Phi_mbJets))   {Phi_mbJets = temp;}
          }//end if
        }//endfor
	if(OneProng){
		h_dR_allJets_p1->Fill(dR_allJets,Weight);
		h_dR_lbJets_p1->Fill(dR_lbJets,Weight);
		h_dR_mbJets_p1->Fill(dR_mbJets,Weight);
		h_Apl_allJets_p1->Fill(Apl_allJets,Weight);
		h_Apl_lbJets_p1->Fill(Apl_lbJets,Weight);
		h_Apl_mbJets_p1->Fill(Apl_mbJets,Weight);
		h_Phi_allJets_p1->Fill(Phi_allJets,Weight);
		h_Phi_lbJets_p1->Fill(Phi_lbJets,Weight);
		h_Phi_mbJets_p1->Fill(Phi_mbJets,Weight);
		}
	if(ThreeProng){
		h_dR_allJets_p3->Fill(dR_allJets,Weight);
		h_dR_lbJets_p3->Fill(dR_lbJets,Weight);
		h_dR_mbJets_p3->Fill(dR_mbJets,Weight);
		h_Apl_allJets_p3->Fill(Apl_allJets,Weight);
		h_Apl_lbJets_p3->Fill(Apl_lbJets,Weight);
		h_Apl_mbJets_p3->Fill(Apl_mbJets,Weight);
		h_Phi_allJets_p3->Fill(Phi_allJets,Weight);
		h_Phi_lbJets_p3->Fill(Phi_lbJets,Weight);
		h_Phi_mbJets_p3->Fill(Phi_mbJets,Weight);
		}

	
         h_e_3dsig->Fill(Electron_sip3d[electron_idx],Weight); 
         h_e_3d->Fill(Electron_ip3d[electron_idx],Weight);
         h_e_dxy->Fill(abs(Electron_dxy[electron_idx]),Weight);
	
        if (Tau_idx > -1 && electron_idx > -1){
            invMass = (*(Tau_p4) + *(Electron_p4)).M();
	    if(OneProng){
		    h_Tau_Electron_invariant_mass_p1->Fill(invMass,Weight2);
		    h_Tau_Electron_invariant_mass_weighted_p1->Fill(invMass,Weight);
		    }
	    if(ThreeProng){
		    h_Tau_Electron_invariant_mass_p3->Fill(invMass,Weight2);
		    h_Tau_Electron_invariant_mass_weighted_p3->Fill(invMass,Weight);
		    }
        }
     
        /*if(Signal && FromTau) {toutT->Fill();}
	else {tout->Fill();}*/
	tout->Fill();
    }

    delete fecorr_trig;

    std::cout << "non_matching_Tau = " << non_matching_Tau << endl;
    std::cout << "non_matching_electron = " << non_matching_electron << endl;

    std::cout << "NeV = " << nEv << endl;
    std::cout << "trigger dropped = " << trigger_dropped << endl;
    std::cout << "selections dropped = " << n_dropped << endl; //remember the cross trigger in Data

    std::cout << "Fraction of events discarded by trigger = " << (trigger_dropped * 1. / nEv) << endl;
    int Rem_trigger=nEv-trigger_dropped; //remember the cross trigger in Data
    std::cout << "Fraction of events removed by selections = " << (n_dropped * 1. / Rem_trigger) << endl;
    std::cout << "Final number of events "<< Rem_trigger - n_dropped<<endl;

    tout->Write();
    trun_out->Write();
   
    HistWrite();

    fout->Write();
    fout->Close();
   /* if (Signal) {
	cout<<"Saving Tau File!"<<endl;
	foutT->cd();
	toutT->Write();
	trun_outT->Write();
	foutT->Write();
	foutT->Close();
	}
   else cout<<"Not writing tree in Tau File"<<endl;*/
}

int main(int argc, char **argv)
{

    string inputFile = argv[1];
    string outputFile = argv[2];
    double crossSection = atof(argv[3]);
    double IntLuminosity = atof(argv[4]);
    string boolstr = argv[5];
    bool Signal = (boolstr == "true");

     HistIniz();

    Mixed_Analysis(inputFile, outputFile, crossSection, IntLuminosity, Signal);

    return 0;
}
