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
#include "Auxiliary.cc"

// correctionlib
#include "correction.h"
using namespace std;
using correction::CorrectionSet;

#define MAX_ARRAY_SIZE 128
#define GEN_MAX_ARRAY_SIZE 1024


double getWeight(double luminosity, double crossSection, Float_t genWeight, double SumWeights){
    return (luminosity * crossSection * genWeight);
}

void Mixed_Analysis(string inputFile, string ofile, double crossSection = -1, double IntLuminosity = 59.827879506, bool Data= false, bool systematics=false, string processname="A"){
    if (crossSection < 0. || IntLuminosity < 0.){
        std::cout << "WARNING: crossection " << crossSection << " and Integrated luminosity " << IntLuminosity << endl;
    }

cout<<"Call completed!"<<endl;
 if (Data) {systematics=false;}

    TFile *fin = TFile::Open(inputFile.c_str());
    TTree *trun;
    Long64_t genEventCount;
    Double_t genEventSumw;
    if(!Data){
	    trun = static_cast<TTree *>(fin->Get("Runs"));
	    trun->SetBranchStatus("*", 0);
	    trun->SetBranchStatus("genEventSumw", 1);
	    trun->SetBranchStatus("genEventCount", 1);
	    trun->SetBranchAddress("genEventSumw", &genEventSumw);
	    trun->SetBranchAddress("genEventCount", &genEventCount);

	    trun->GetEntry(0);
    }

    TTree *tin = static_cast<TTree *>(fin->Get("Events"));

    tin->SetBranchStatus("*", 0);
    Float_t Tau_pt[MAX_ARRAY_SIZE], Electron_pt[MAX_ARRAY_SIZE], Jet_pt[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_pt", 1);
    tin->SetBranchAddress("Electron_pt", &Electron_pt);
    tin->SetBranchStatus("Tau_pt", 1);
    tin->SetBranchAddress("Tau_pt", &Tau_pt);
    tin->SetBranchStatus("Jet_pt", 1);
    tin->SetBranchAddress("Jet_pt", &Jet_pt);
    UInt_t nTau, nElectron;
    tin->SetBranchStatus("nElectron", 1);
    tin->SetBranchAddress("nElectron", &nElectron);
    tin->SetBranchStatus("nTau", 1);
    tin->SetBranchAddress("nTau", &nTau);
    Float_t Tau_eta[MAX_ARRAY_SIZE], Electron_eta[MAX_ARRAY_SIZE], Jet_eta[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_eta", 1);
    tin->SetBranchAddress("Electron_eta", &Electron_eta);
    tin->SetBranchStatus("Tau_eta", 1);
    tin->SetBranchAddress("Tau_eta", &Tau_eta);
    tin->SetBranchStatus("Jet_eta", 1);
    tin->SetBranchAddress("Jet_eta", &Jet_eta);
    Float_t Tau_phi[MAX_ARRAY_SIZE], Electron_phi[MAX_ARRAY_SIZE], Jet_phi[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_phi", 1);
    tin->SetBranchAddress("Electron_phi", &Electron_phi);
    tin->SetBranchStatus("Tau_phi", 1);
    tin->SetBranchAddress("Tau_phi", &Tau_phi);
    tin->SetBranchStatus("Jet_phi", 1);
    tin->SetBranchAddress("Jet_phi", &Jet_phi);
    Float_t Tau_mass[MAX_ARRAY_SIZE], Electron_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_mass", 1);
    tin->SetBranchAddress("Electron_mass", &Electron_mass);
    tin->SetBranchStatus("Tau_mass", 1);
    tin->SetBranchAddress("Tau_mass", &Tau_mass);
    tin->SetBranchStatus("Jet_mass", 1);
    tin->SetBranchAddress("Jet_mass", &Jet_mass);

    // get gen quantities
    Int_t Tau_genPartIdx[MAX_ARRAY_SIZE], Electron_genPartIdx[MAX_ARRAY_SIZE],Jet_hadronFlavour[MAX_ARRAY_SIZE];
    Int_t GenPart_pdgId[GEN_MAX_ARRAY_SIZE], GenPart_genPartIdxMother[GEN_MAX_ARRAY_SIZE], GenPart_statusFlags[GEN_MAX_ARRAY_SIZE], Jet_genJetIdx[MAX_ARRAY_SIZE];
    UChar_t Tau_genPartFlav[MAX_ARRAY_SIZE], Electron_genPartFlav[MAX_ARRAY_SIZE];
    UInt_t nGenPart;
    Float_t GenPart_pt[GEN_MAX_ARRAY_SIZE];
    Float_t genWeight, N_pu_vertices, L1PreFiringWeight_Nom;
    if(!Data){
	    tin->SetBranchStatus("Electron_genPartIdx", 1);
	    tin->SetBranchStatus("Electron_genPartFlav", 1);
	    tin->SetBranchStatus("Tau_genPartIdx", 1);
	    tin->SetBranchStatus("Tau_genPartFlav", 1);
	    tin->SetBranchStatus("GenPart_pdgId", 1);
	    tin->SetBranchStatus("GenPart_genPartIdxMother", 1);
	    tin->SetBranchStatus("nGenPart", 1);
	    tin->SetBranchStatus("Jet_genJetIdx",1);
	    tin->SetBranchStatus("GenPart_pt",1);
	    tin->SetBranchStatus("GenPart_statusFlags",1);
	    tin->SetBranchStatus("Jet_hadronFlavour", 1);
	    tin->SetBranchAddress("nGenPart", &nGenPart);
	    tin->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx);
	    tin->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav);
	    tin->SetBranchAddress("Tau_genPartIdx", &Tau_genPartIdx);
	    tin->SetBranchAddress("Tau_genPartFlav", &Tau_genPartFlav);
	    tin->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
	    tin->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
	    tin->SetBranchAddress("Jet_genJetIdx",&Jet_genJetIdx);
	    tin->SetBranchAddress("GenPart_pt",&GenPart_pt);
	    tin->SetBranchAddress("GenPart_statusFlags",&GenPart_statusFlags);
	    tin->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
	    tin->SetBranchStatus("Pileup_nTrueInt", 1);
	    tin->SetBranchAddress("Pileup_nTrueInt", &N_pu_vertices);
	    tin->SetBranchStatus("genWeight", 1);
	    tin->SetBranchAddress("genWeight", &genWeight);
	    tin->SetBranchStatus("L1PreFiringWeight_Nom", 1);
	    tin->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom);
    }

    // collect the trigger information
    Bool_t HLT_Ele32_WPTight_Gsf;
    tin->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1);
    tin->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);

    // collect the triggger Ids
    Int_t Tau_charge[MAX_ARRAY_SIZE], Electron_charge[MAX_ARRAY_SIZE];
    Bool_t Electron_mvaFall17V2Iso_WP90[MAX_ARRAY_SIZE];
    Float_t  Electron_ip3d[MAX_ARRAY_SIZE], Electron_sip3d[MAX_ARRAY_SIZE], Electron_dxy[MAX_ARRAY_SIZE], Electron_dz[MAX_ARRAY_SIZE];

    tin->SetBranchStatus("Tau_charge", 1);
    tin->SetBranchStatus("Electron_charge", 1);
    tin->SetBranchStatus("Electron_mvaFall17V2Iso_WP90", 1);
    tin->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", &Electron_mvaFall17V2Iso_WP90);
    tin->SetBranchAddress("Tau_charge", &Tau_charge);
    tin->SetBranchAddress("Electron_charge", &Electron_charge);
    tin->SetBranchStatus("Electron_ip3d", 1);
    tin->SetBranchStatus("Electron_sip3d", 1);
    tin->SetBranchStatus("Electron_dxy", 1);
    tin->SetBranchStatus("Electron_dz", 1);
    tin->SetBranchAddress("Electron_dz", &Electron_dz);
    tin->SetBranchAddress("Electron_ip3d", &Electron_ip3d);
    tin->SetBranchAddress("Electron_sip3d", &Electron_sip3d);
    tin->SetBranchAddress("Electron_dxy", &Electron_dxy);

    // Jet tagging and ID, FlavB is the recomended one, DeepB was used by Anup
    Float_t Jet_btagDeepFlavB[MAX_ARRAY_SIZE], Jet_btagDeepB[MAX_ARRAY_SIZE];
    UInt_t nJet;
    Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Jet_btagDeepB", 1);
    tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
    tin->SetBranchStatus("nJet", 1);
    tin->SetBranchStatus("Jet_jetId", 1);
    tin->SetBranchStatus("Jet_puId", 1);
    tin->SetBranchAddress("nJet", &nJet);
    tin->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);
    tin->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);   
    tin->SetBranchAddress("Jet_jetId", &Jet_jetId);
    tin->SetBranchAddress("Jet_puId", &Jet_puId);

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
   

    int non_matching_Tau = 0, non_matching_electron = 0, evenottrigMatch=0,n_dropped = 0, trigger_dropped = 0;
    UInt_t nEv = tin->GetEntries();
    unsigned int n_events = nEv;
    TLorentzVector *Tau_p4 = new TLorentzVector();
    TLorentzVector *Electron_p4 = new TLorentzVector();
    TLorentzVector *MainBjet_p4 = new TLorentzVector();
    TLorentzVector *OppositeBjet_p4 = new TLorentzVector();

    // allow pt, inv mass, and eta to be stored in a Branch
    Float_t leading_lepton_pt, invMass, electron_eta, electron_pt, tau_eta, tau_pt;
    float Weight;

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
   
    TRandom3 * RndGen=new TRandom3();
     
    TFile *fout = new TFile(ofile.c_str(), "RECREATE");
    
    // create a new tree for the output
    TTree *tout = new TTree("tout", "tout");
    TTree *trun_out;
    
    // set the branches for the output tree
    tout->Branch("leading_lepton_pt", &leading_lepton_pt);
    tout->Branch("invMass", &invMass);
    tout->Branch("electron_eta", &electron_eta);
    tout->Branch("electron_pt", &electron_pt);
    tout->Branch("tau_eta", &tau_eta);
    tout->Branch("tau_pt", &tau_pt);

    int Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0;
    float PTbjet,Acopl_etau;

    tout->Branch("PTbjet", &PTbjet);
    tout->Branch("Nloose", &Nloose);
    tout->Branch("Nmedium", &Nmedium);
    tout->Branch("Ntight", &Ntight);
    tout->Branch("JetNotB", &JetsNotB);
    tout->Branch("Acopl_etau", &Acopl_etau);

    if(!Data){
	    trun_out = new TTree("Run_out", "Run_out");
	    trun_out->Branch("genEventSumw", &genEventSumw);
	    trun_out->Branch("IntLumi", &IntLuminosity);
	    trun_out->Branch("xs", &crossSection);
	    trun_out->Branch("nEv", &n_events);
	    trun_out->Fill();
    }
    
    vector<string> observables= {processname+"_Tau_pt",processname+"_Electron_pt",processname+"_B_pt",processname+"_Invariant_Mass",processname+"_Lepton_Acoplanarity",processname+"_Njets"}; 
    vector<string> systs= {"Tau_Id_mu","Tau_Id_Ele","Tau_Id_Jet","Ele_Reco","Ele_IdIso","Ele_trigger"};
    vector<string> shift={"","Up","Down"};
    
    
    TH1D* temp=NULL;
    vector<TH1D*>  Histos;
    if(systematics){ for(int i=0;i<observables.size()*(systs.size()*2+1);i++) Histos.push_back(temp); }
    else for(int i=0;i<observables.size();i++) Histos.push_back(temp);
    
    //get histos nominal values
    temp= new TH1D(observables[0].c_str(),observables[0].c_str(),44,0,220);
    Histos[0] = (TH1D*)temp->Clone();
    temp = new TH1D(observables[1].c_str(),observables[1].c_str(),40, 0, 200);
    Histos[1] = (TH1D*)temp->Clone();
    temp = new TH1D(observables[2].c_str(),observables[2].c_str(),40,25,425);
    Histos[2]= (TH1D*)temp->Clone(); 
    temp= new TH1D(observables[3].c_str(),observables[3].c_str(),40, 12, 412);
    Histos[3] = (TH1D*)temp->Clone();
    temp = new TH1D(observables[4].c_str(),observables[4].c_str(),40,0, 2*M_PI);
    Histos[4]= (TH1D*)temp->Clone();
    temp = new TH1D(observables[5].c_str(),observables[5].c_str(),12,0,12);
    Histos[5] = (TH1D*)temp->Clone();
    
    int auxindex=5;
    //here loop on systematics, clone the histo and do your things
    if(systematics){
    	for(int i=0; i<systs.size(); i++){
	    	 for(int j=0;j<observables.size();j++){
		    	 Histos[++auxindex]=cloneDims1d(Histos[j],(systs[i]+"Up").c_str());	
	    	}
		for(int j=0;j<observables.size();j++){
		    	 Histos[++auxindex]=cloneDims1d(Histos[j],(systs[i]+"Down").c_str());	
	    	}
    	} 
   }

   if(systematics && auxindex!=observables.size()*(systs.size()*2+1)-1) cout<<"Something may go terribly wrong here!"<< auxindex<<" vs "<<observables.size()*(systs.size()*2+1)-1 <<endl;

    #pragma omp parallel for
    for (UInt_t i = 0; i <nEv; i++){
        tin->GetEntry(i);
        if (i % 100000 == 0) std::cout << "Processing entry " << i << " of " << nEv << endl;
        
	// apply triggers
        if (!(HLT_Ele32_WPTight_Gsf)){trigger_dropped++; continue;}
        
        int Master_index=1; //the first entry is the nominal one!
	vector<double> VecWeights(systs.size()*2+1,1.);
	
	bool OneProng=false, ThreeProng=false; //Not used at the moment
        Int_t Tau_idx = -1;
        for (UInt_t j = 0; j < nTau; j++){
		if (Tau_decayMode[j]>2 && Tau_decayMode[j]<10) continue;
		if(!Data){
			     double ScaleE=Tau_Escale->evaluate({Tau_pt[j],abs(Tau_eta[j]),Tau_decayMode[j],Tau_genPartFlav[j],"DeepTau2017v2p1","nom"});
			     Tau_p4->SetPtEtaPhiM(Tau_pt[j]*ScaleE, Tau_eta[j], Tau_phi[j], Tau_mass[j]*ScaleE);
			    }
            if ((Tau_p4->Pt()>22. && abs(Tau_eta[j])<2.3)&&(Tau_idDeepTau2017v2p1VSe[j]>=8&& Tau_idDeepTau2017v2p1VSmu[j]>=8 && Tau_idDeepTau2017v2p1VSjet[j]>=32)){ //Loose e- T mu T jet
		if (Tau_decayMode[j]<=2) {OneProng=true;}
		if (Tau_decayMode[j]>=10) {ThreeProng=true;}
		if (!(OneProng || ThreeProng)) {continue;}
                Tau_idx = j;
                break;
            }
        }
        if (Tau_idx==-1)  { n_dropped++; continue;}	 
        Int_t electron_idx = -1;
        for (UInt_t j = 0; j < nElectron; j++){
            if ((Electron_pt[j] > 35 && abs(Electron_eta[j]) < 2.4 && Electron_mvaFall17V2Iso_WP90[j] && abs(Electron_dxy[j])<0.2 && abs(Electron_dz[j])<0.5)){
		if((abs(Electron_eta[j])>1.44) && (abs(Electron_eta[j])<1.57)) {continue;} //remove electrons in the acceptance break
                Electron_p4->SetPtEtaPhiM(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);
                if (Electron_p4->DeltaR(*Tau_p4) < 0.4) {continue;}
                else{electron_idx = j; break;}
            }
        }
        if (electron_idx==-1) {n_dropped++; continue;}

        bool selection = ((Tau_idx > -1) && (electron_idx > -1));
        selection = selection && ((Tau_charge[Tau_idx] * Electron_charge[electron_idx]) < 0);
        Float_t jet_btag_deepFlav_wp = 0.2783;
        bool one_Bjet = false;
        int id_m_jet = -1, njets=0;
        Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0;
	//vectors for applying b-tag corrections
	vector<int> njet_in_collection; vector<int> flavor; vector<bool> tagged;
	double t_weight=1.;
        for (size_t j = 0; j < nJet; j++) {
            if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6)){
	    TLorentzVector *Tjet_p4 = new TLorentzVector();
	    Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	    if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	    else {delete Tjet_p4;}
	
            //correction for pileupID
            float tempSF=1.,tempEff=0.;
            if(!Data){
		    int MC_pu = Jet_genJetIdx[j];
		    if (MC_pu<0 ) {tempSF=1.; tempEff= 0;}
		    //if is truly a jet
		    else { if (Jet_pt[j]<=50){
		    	     tempSF= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"nom", "L"});
		    	     tempEff= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"MCEff", "L"});
			    }
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
	      if (Jet_btagDeepFlavB[j] > 0.0490) Nloose++;
              if (Jet_btagDeepFlavB[j] > 0.2783) Nmedium++;
              if (Jet_btagDeepFlavB[j] > 0.71) Ntight++;
            
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

        selection = selection && (one_Bjet);
        if (!selection){ n_dropped++;  continue;}
         
        if(!Data){ 
		 
		Weight = getWeight(IntLuminosity, crossSection, genWeight, genEventSumw);
		Weight *=  getTopPtWeight(GenPart_pdgId,GenPart_statusFlags,GenPart_pt,nGenPart);
		Weight *= pu_correction->evaluate({N_pu_vertices, "nominal"});
		Weight*= L1PreFiringWeight_Nom;
		
		for(auto &k: VecWeights) k*=Weight;
		
				
		for(int j=0;j<VecWeights.size();j++){
			if(systematics){
				if(j==Master_index) {VecWeights[j]*= Tau_idvsmu->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Tight","up"}); continue;}
				if(j==Master_index+1) {VecWeights[j]*= Tau_idvsmu->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Tight","down"}); continue;}
				else VecWeights[j] *= Tau_idvsmu->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Tight","nom"});
			}
			else VecWeights[j] *= Tau_idvsmu->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Tight","nom"});
		} 
		Master_index+=2;
		for(int j=0;j<VecWeights.size();j++){
			if(systematics){
				if(j==Master_index) {VecWeights[j]*= Tau_idvse->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Loose","up"}); continue;}
				if(j==Master_index+1) {VecWeights[j]*= Tau_idvse->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Loose","down"}); continue;}
				else VecWeights[j]*= Tau_idvse->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Loose","nom"});
			}
			else VecWeights[j]*= Tau_idvse->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Loose","nom"});
		} 
		Master_index+=2;
		
		for(int j=0;j<VecWeights.size();j++){
			if(systematics){
				if(j==Master_index) {VecWeights[j]*= Tau_idvsjet->evaluate({Tau_p4->Pt(), Tau_decayMode[Tau_idx],Tau_genPartFlav[Tau_idx],"Tight","up","pt"}); continue;}
				if(j==Master_index+1) {VecWeights[j]*= Tau_idvsjet->evaluate({Tau_p4->Pt(), Tau_decayMode[Tau_idx],Tau_genPartFlav[Tau_idx],"Tight","down","pt"}); continue;}
				else VecWeights[j]*= Tau_idvsjet->evaluate({Tau_p4->Pt(), Tau_decayMode[Tau_idx],Tau_genPartFlav[Tau_idx],"Tight","nom","pt"});
			}
			else VecWeights[j]*=Tau_idvsjet->evaluate({Tau_p4->Pt(), Tau_decayMode[Tau_idx],Tau_genPartFlav[Tau_idx],"Tight","nom","pt"});
		} 
		Master_index+=2;
		
		for(int j=0;j<VecWeights.size();j++){
			if(systematics){
				if(j==Master_index) {VecWeights[j]*= electron_id->evaluate({"2018", "sfup", "RecoAbove20", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]}); continue;}
				if(j==Master_index+1) {VecWeights[j]*= electron_id->evaluate({"2018", "sfdown", "RecoAbove20", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]}); continue;}
				else VecWeights[j]*= electron_id->evaluate({"2018", "sf", "RecoAbove20", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]});
			}
			else VecWeights[j]*= electron_id->evaluate({"2018", "sf", "RecoAbove20", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]});
		} 
		Master_index+=2;
		for(int j=0;j<VecWeights.size();j++){
			if(systematics){
				if(j==Master_index) {VecWeights[j]*= electron_id->evaluate({"2018", "sfup", "wp90iso", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]}); continue;}
				if(j==Master_index+1) {VecWeights[j]*= electron_id->evaluate({"2018", "sfdown", "wp90iso", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]}); continue;}
				else VecWeights[j]*= electron_id->evaluate({"2018", "sf", "wp90iso", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]});
			}
			else VecWeights[j]*= electron_id->evaluate({"2018", "sf", "wp90iso", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]});
		} 
		Master_index+=2;
		
		if(HLT_Ele32_WPTight_Gsf) {
			 int bin = EleTrigHisto->FindBin(Electron_eta[electron_idx],Electron_pt[electron_idx]);
			 float temp= EleTrigHisto->GetBinContent(bin);
			 float downErr=EleTrigHisto->GetBinErrorLow(bin);
			 float upErr=EleTrigHisto->GetBinErrorUp(bin);
			 for(int j=0;j<VecWeights.size();j++){
				if(systematics){
					if(j==Master_index) {VecWeights[j]*=(temp+upErr); continue;}
					if(j==Master_index+1) {VecWeights[j]*= (temp-downErr); continue;}
					else VecWeights[j]*= temp;
				}
				else VecWeights[j]*= temp;
			} 
			Master_index+=2;
		}
		    
		
		for(auto &k: VecWeights) k*=t_weight; 
		double Weight2=1.;
		for(int jj=0;jj<flavor.size();jj++){
			int convflav=flavor[jj];
			if (flavor[jj]<4) convflav==0;
			if (!(convflav==0 || convflav==4 || convflav==5)) {cout<<"Something weird in the flavor of jet"<<endl;}
			if(tagged[jj]){
				if (convflav!=0) {Weight2 *= b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});}
				else  {Weight2 *= b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});}
				continue;}

			if(!tagged[jj]) {
				double Eff=1.;
				double SF=1.;
				if (convflav!=0) SF=b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
				else SF=b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
				
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
				Weight2*=(1-SF*Eff)/(1-Eff);
				}
			
			}
		for(auto &k: VecWeights) k*=Weight2;
	}
	
	tau_pt = Tau_p4->Pt();
        tau_eta = Tau_eta[Tau_idx];
        electron_pt = Electron_pt[electron_idx];
        electron_eta = Electron_eta[electron_idx];
        PTbjet = MainBjet_p4->Pt();
        
        if (Tau_idx > -1 && electron_idx > -1){
            invMass = (*(Tau_p4) + *(Electron_p4)).M();
        }
        Acopl_etau=M_PI-(Electron_p4->DeltaPhi(*Tau_p4));
	
	for(int k=0;k<VecWeights.size();k++){
		Histos[k*observables.size()] ->Fill(tau_pt,VecWeights[k]);
    		Histos[k*observables.size()+1] ->Fill(electron_pt,VecWeights[k]);
    		Histos[k*observables.size()+2] ->Fill(MainBjet_p4->Pt(),VecWeights[k]);
    		Histos[k*observables.size()+3] ->Fill(invMass,VecWeights[k]);
    		Histos[k*observables.size()+4] ->Fill(Acopl_etau,VecWeights[k]);
    		Histos[k*observables.size()+5] ->Fill(njets,VecWeights[k]);
		if(!systematics){break;}
		}

        tout->Fill();
    }

    delete fecorr_trig;


    std::cout << "NeV = " << nEv << endl;
    std::cout << "trigger dropped = " << trigger_dropped << endl;
    std::cout << "selections dropped = " << n_dropped << endl; //remember the cross trigger in Data

    std::cout << "Fraction of events discarded by trigger = " << (trigger_dropped * 1. / nEv) << endl;
    int Rem_trigger=nEv-trigger_dropped; //remember the cross trigger in Data
    std::cout << "Fraction of events removed by selections = " << (n_dropped * 1. / Rem_trigger) << endl;
    std::cout << "Final number of events "<< Rem_trigger - n_dropped<<endl;

    tout->Write();
    if(!Data) trun_out->Write();

    for(auto &k: Histos) HistWrite(k);

    fout->Close();
}

int main(int argc, char **argv){
    string inputFile = argv[1];
    string outputFile = argv[2];
    double crossSection = atof(argv[3]);
    double IntLuminosity = atof(argv[4]);
    string boolstr = argv[5];
    bool Data = (boolstr == "true" || boolstr == "True");
    string boolstr2 = argv[6];
    bool systematics = (boolstr2 == "true" || boolstr2 == "True");
    string processname = argv[7];

cout<<"Process "<<processname<<endl;

    Mixed_Analysis(inputFile, outputFile, crossSection, IntLuminosity, Data, systematics,processname);

    return 0;
}
