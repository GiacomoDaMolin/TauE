#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

// include user defined histograms and auxiliary macros
#include "Histodef.cpp"
#include "Auxiliary.cpp"
#include "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR.cc"

using namespace std;

#define MAX_ARRAY_SIZE 128

void DataAnalysis(string inputFile, string ofile, bool IsFirstDataSet)
{

    TFile *fin = TFile::Open(inputFile.c_str());
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
    Float_t Tau_mass[MAX_ARRAY_SIZE], Electron_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_mass", 1);
    tin->SetBranchAddress("Electron_mass", &Electron_mass);
    tin->SetBranchStatus("Tau_mass", 1);
    tin->SetBranchAddress("Tau_mass", &Tau_mass);
    tin->SetBranchStatus("Jet_mass", 1);
    tin->SetBranchAddress("Jet_mass", &Jet_mass);

    // collect the trigger information
    Bool_t  HLT_Ele32_WPTight_Gsf;
    tin->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1);
    tin->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);

    // collect the triggger Ids
    Int_t Tau_charge[MAX_ARRAY_SIZE], Electron_charge[MAX_ARRAY_SIZE];
    Bool_t Electron_mvaFall17V2Iso_WP90[MAX_ARRAY_SIZE];
    Float_t Electron_ip3d[MAX_ARRAY_SIZE], Electron_sip3d[MAX_ARRAY_SIZE], Electron_dxy[MAX_ARRAY_SIZE], Electron_dz[MAX_ARRAY_SIZE];
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

    // Jet tagging , FlavB is the recomennded one, DeepB was used by Anup
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


    int non_matching_Tau = 0, non_matching_electron = 0;
    int n_dropped = 0;
    int trigger_dropped = 0,crosstrigger=0;
    const auto nEv = tin->GetEntries();
    TLorentzVector *Tau_p4 = new TLorentzVector();
    TLorentzVector *Electron_p4 = new TLorentzVector();
    TLorentzVector *MainBjet_p4 = new TLorentzVector();
    TLorentzVector *OppositeBjet_p4 = new TLorentzVector();

       // allow pt, inv mass, and eta to be stored in a Branch
    Float_t leading_lepton_pt, invMass, electron_eta, electron_pt, tau_eta, tau_pt;
    Float_t Tau_eta_from_W, Tau_pt_from_W, electron_eta_from_W, electron_pt_from_W;
    TFile *fout =new TFile(ofile.c_str(),"RECREATE");
    // create a new tree for the output
    TTree *tout = new TTree("tout","tout");
    // set the branches for the output tree
    tout->Branch("leading_lepton_pt", &leading_lepton_pt);
    tout->Branch("invMass", &invMass);
    tout->Branch("electron_eta", &electron_eta);
    tout->Branch("electron_pt", &electron_pt);
    tout->Branch("tau_eta", &tau_eta);
    tout->Branch("tau_pt", &tau_pt);

    int Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0, Nprongs=0;
    float dR_muE, dR_mujet, dR_ejet, dR_allJets, dR_lbJets, dR_mbJets, Apl_allJets, Apl_lbJets, Apl_mbJets, Phi_allJets, Phi_lbJets, Phi_mbJets, PTbjet, Acopl_etau;

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

    #pragma omp parallel for
    for (UInt_t i = 0; i < nEv; i++)
    {
        tin->GetEntry(i);
        if (i % 100000 == 0)
            std::cout << "Processing entry " << i << " of " << nEv << std::endl;
        // apply triggers

        if (!(HLT_Ele32_WPTight_Gsf))
        { trigger_dropped++; continue; };

        // avoid cross triggers
        /*if (!IsFirstDataSet && HLT_Ele32_WPTight_Gsf && HLT_IsoMu24)
        {crosstrigger++;  continue;}*/

        // loop over the Taus and electrons and only keep the fist ones that pass the requirements
	bool OneProng=false, ThreeProng=false;
        Int_t Tau_idx = -1;
        for (UInt_t j = 0; j < nTau; j++){
            if ((Tau_pt[j]>22. && abs(Tau_eta[j])<2.3)&&(Tau_idDeepTau2017v2p1VSe[j]>=8 && Tau_idDeepTau2017v2p1VSmu[j]>=8 && Tau_idDeepTau2017v2p1VSjet[j]>=32)){ //Loose e- T mu T jet
		if (Tau_decayMode[j]<=2) {OneProng=true;}
		if (Tau_decayMode[j]>=10) {ThreeProng=true;}
		if (!(OneProng || ThreeProng)) {continue;}
                Tau_idx = j;
                Tau_p4->SetPtEtaPhiM(Tau_pt[j], Tau_eta[j], Tau_phi[j], Tau_mass[j]);
                break;
            }
        }
        if (Tau_idx==-1) {n_dropped++; continue; }

        Int_t electron_idx = -1;
        for (UInt_t j = 0; j < nElectron; j++){
            if ((Electron_pt[j]>35 && abs(Electron_eta[j])<2.4 && Electron_mvaFall17V2Iso_WP90[j] && abs(Electron_dxy[j])<0.2 && abs(Electron_dz[j])<0.5)){
		if((abs(Electron_eta[j])>1.44) && (abs(Electron_eta[j])<1.57)) {continue;}
                Electron_p4->SetPtEtaPhiM(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);
		if(Electron_p4->DeltaR(*Tau_p4)<0.4) {continue;}
                else {electron_idx = j; break;}
            }
        }
        if (electron_idx==-1) {
            n_dropped++;
            continue;
        }
        bool selection = ((Tau_idx > -1) && (electron_idx > -1));
        // check the seleected objects for opposite charge
        selection = selection && ((Tau_charge[Tau_idx] * Electron_charge[electron_idx]) < 0);
        
        // the tight working point is 0.71, medium 0.2783, loose 0.0490
        Float_t jet_btag_deepFlav_wp = 0.2783;
        bool one_Bjet = false;
        int id_m_jet=-1;
	int njets=0;
	Nloose=0, Nmedium=0, Ntight=0,JetsNotB=0;
        for (size_t j = 0; j < nJet; j++){
          if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || (Jet_puId[j]>=4)))
            {
	    TLorentzVector *Tjet_p4 = new TLorentzVector();
	    Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	    if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	    else {delete Tjet_p4;}
	    njets++;
	    if (Jet_btagDeepFlavB[j] < 0.0490) JetsNotB++;
	    if (Jet_btagDeepFlavB[j] > 0.0490)	Nloose++;
	    if (Jet_btagDeepFlavB[j] > 0.2783)	Nmedium++;
	    if (Jet_btagDeepFlavB[j] > 0.71)	Ntight++;
            if (Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp){
		if(!one_Bjet) {
			MainBjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
			OppositeBjet_p4->SetPtEtaPhiM(Jet_pt[j], -1*Jet_eta[j], InvertPhi(Jet_phi[j]), Jet_mass[j]);
				
			if ( MainBjet_p4->DeltaR(*Tau_p4)>0.4 && MainBjet_p4->DeltaR(*Electron_p4)>0.4){
                		one_Bjet = true; id_m_jet=j;
				}
                	}	
            } //end tag
	  }//end kinematic if
        }
        selection = selection && (one_Bjet);

	

        if (!selection){
            n_dropped++;
            continue;
        }
	Acopl_etau=M_PI-(Electron_p4->DeltaPhi(*Tau_p4));

	if(OneProng){
		h_LooseJets_p1->Fill(Nloose);
		h_MediumJets_p1->Fill(Nmedium);
		h_TightJets_p1->Fill(Ntight);
		h_acopla_etau_p1->Fill(Acopl_etau);
		h_NJets_p1->Fill(njets);
		}
        if(ThreeProng){
		h_LooseJets_p3->Fill(Nloose);
		h_MediumJets_p3->Fill(Nmedium);
		h_TightJets_p3->Fill(Ntight);
		h_acopla_etau_p3->Fill(Acopl_etau);
		h_NJets_p3->Fill(njets);
		}
	if(OneProng) {Nprongs=1;}
	if(ThreeProng) {Nprongs=3;}
        PTbjet=MainBjet_p4->Pt();

        dR_mujet=Tau_p4->DeltaR(*MainBjet_p4);
	dR_ejet=Electron_p4->DeltaR(*MainBjet_p4);
	dR_muE=Tau_p4->DeltaR(*Electron_p4);
	// fill the tree
        tau_pt = Tau_pt[Tau_idx];
        tau_eta = Tau_eta[Tau_idx];
        electron_pt = Electron_pt[electron_idx];
        electron_eta = Electron_eta[electron_idx];

        // check whether Tau or electron is the leading one
        if (Tau_p4->Pt() > Electron_p4->Pt())	{leading_lepton_pt = Tau_p4->Pt();}
	else					{leading_lepton_pt = Electron_p4->Pt();}

	if(OneProng){
		h_leading_lepton_pt_p1->Fill(leading_lepton_pt);
		h_leading_lepton_pt_weighted_p1->Fill(leading_lepton_pt);
		h_Tau_pt_p1->Fill(tau_pt);
		h_Tau_eta_p1->Fill(tau_eta);
		h_Electron_pt_p1->Fill(electron_pt);
		h_Electron_eta_p1->Fill(electron_eta);
		h_Tau_pt_weighted_p1->Fill(tau_pt);
		h_Tau_eta_weighted_p1->Fill(tau_eta);
		h_Electron_pt_weighted_p1->Fill(electron_pt);
		h_Electron_eta_weighted_p1->Fill(electron_eta);
		b_pt_p1->Fill(MainBjet_p4->Pt());
   		jethole_p1->Fill(MainBjet_p4->Eta(),MainBjet_p4->Phi());
   		tauhole_p1->Fill(Tau_p4->Eta(),Tau_p4->Phi());
   		ehole_p1->Fill(Electron_p4->Phi(),Electron_p4->Phi());
		}
	if(ThreeProng){
		h_leading_lepton_pt_p3->Fill(leading_lepton_pt);
		h_leading_lepton_pt_weighted_p3->Fill(leading_lepton_pt);
		h_Tau_pt_p3->Fill(tau_pt);
		h_Tau_eta_p3->Fill(tau_eta);
		h_Electron_pt_p3->Fill(electron_pt);
		h_Electron_eta_p3->Fill(electron_eta);
		h_Tau_pt_weighted_p3->Fill(tau_pt);
		h_Tau_eta_weighted_p3->Fill(tau_eta);
		h_Electron_pt_weighted_p3->Fill(electron_pt);
		h_Electron_eta_weighted_p3->Fill(electron_eta);
		b_pt_p3->Fill(MainBjet_p4->Pt());
   		jethole_p3->Fill(MainBjet_p4->Eta(),MainBjet_p4->Phi());
   		tauhole_p3->Fill(Tau_p4->Eta(),Tau_p4->Phi());
   		ehole_p3->Fill(Electron_p4->Phi(),Electron_p4->Phi());
		}

        
	dR_allJets=999, dR_lbJets=999, dR_mbJets=999;
	Apl_allJets=1.1,Apl_lbJets=1.1,Apl_mbJets=1.1;
	bool ok1=false,ok2=false ,ok3=false;

	for (size_t j = 0; j < nJet; j++){
		if (j==id_m_jet) continue;
		TLorentzVector *Tjet_p4 = new TLorentzVector();
	        Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	        if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	        
                if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || (Jet_puId[j]>=4))){
		 TLorentzVector *tempJet = Tjet_p4;
		 double temp=OppositeBjet_p4->DeltaR(*tempJet);

		 TVector3 A(tempJet->X(),tempJet->Y(),tempJet->Z());	
		 TVector3 B(MainBjet_p4->X(),MainBjet_p4->Y(),MainBjet_p4->Z());

		 double tempApl=A.Dot(B)/(A.Mag()*B.Mag());

		 if(temp<dR_allJets) {dR_allJets=temp; ok1=true;}
		 if(tempApl<Apl_allJets) {Apl_allJets=tempApl;}
		 if(Jet_btagDeepFlavB[j] > 0.0490){ ok2=true;
			if (temp<dR_lbJets){dR_lbJets=temp;}
			if (tempApl<Apl_lbJets) {Apl_lbJets=tempApl;}
			}
		 if(Jet_btagDeepFlavB[j] > 0.2783){ ok3=true;
			if (temp<dR_mbJets){dR_mbJets=temp;}
			if (tempApl<Apl_mbJets) {Apl_mbJets=tempApl;}
			}
  
			
		 delete Tjet_p4;	
		 }//end if
        	}//end for

//dphi
	Phi_allJets=999, Phi_lbJets=999, Phi_mbJets=999;
	for (size_t j = 0; j < nJet; j++){
		if(j==id_m_jet) continue;

		TLorentzVector *Tjet_p4 = new TLorentzVector();
	        Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	       	if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	        else {delete Tjet_p4;}
		if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || (Jet_puId[j]>=4))){
			double temp=Jet_phi[j]-OppositeBjet_p4->Phi();
			if (temp<-1*M_PI) temp+=2*M_PI;
			if (temp>M_PI) temp-=2*M_PI;
			if (temp<0) temp*=(-1);

			if(temp<Phi_allJets) {Phi_allJets=temp;}
			if((Jet_btagDeepFlavB[j] > 0.0490) && (temp<Phi_lbJets)) {Phi_lbJets=temp;}
			if((Jet_btagDeepFlavB[j] > 0.2783) && (temp<Phi_mbJets)) {Phi_mbJets=temp;}
		}//end if
        }//end for

	if(OneProng){
		h_dR_allJets_p1->Fill(dR_allJets);
		h_dR_lbJets_p1->Fill(dR_lbJets);
		h_dR_mbJets_p1->Fill(dR_mbJets);
		h_Apl_allJets_p1->Fill(Apl_allJets);
		h_Apl_lbJets_p1->Fill(Apl_lbJets);
		h_Apl_mbJets_p1->Fill(Apl_mbJets);
		h_Phi_allJets_p1->Fill(Phi_allJets);
		h_Phi_lbJets_p1->Fill(Phi_lbJets);
		h_Phi_mbJets_p1->Fill(Phi_mbJets);
		}
	if(ThreeProng){
		h_dR_allJets_p3->Fill(dR_allJets);
		h_dR_lbJets_p3->Fill(dR_lbJets);
		h_dR_mbJets_p3->Fill(dR_mbJets);
		h_Apl_allJets_p3->Fill(Apl_allJets);
		h_Apl_lbJets_p3->Fill(Apl_lbJets);
		h_Apl_mbJets_p3->Fill(Apl_mbJets);
		h_Phi_allJets_p3->Fill(Phi_allJets);
		h_Phi_lbJets_p3->Fill(Phi_lbJets);
		h_Phi_mbJets_p3->Fill(Phi_mbJets);
		}

         h_e_3dsig->Fill(Electron_sip3d[electron_idx]); 
         h_e_3d->Fill(Electron_ip3d[electron_idx]);
         h_e_dxy->Fill(abs(Electron_dxy[electron_idx]));

        if (Tau_idx > -1 && electron_idx > -1)
        {
            invMass = (*(Tau_p4) + *(Electron_p4)).M();
	    if(OneProng){
		    h_Tau_Electron_invariant_mass_p1->Fill(invMass);
		    h_Tau_Electron_invariant_mass_weighted_p1->Fill(invMass);
		    }
	    if(ThreeProng){
		    h_Tau_Electron_invariant_mass_p3->Fill(invMass);
		    h_Tau_Electron_invariant_mass_weighted_p3->Fill(invMass);
		    }
        }
	tout->Fill();
    }
    std::cout << "Total number of events: " << nEv << std::endl;
    int NumbEv=nEv;
    std::cout << "Removed because in another sample = " << crosstrigger << endl;
    std::cout << "trigger dropped = " << trigger_dropped << endl;
    std::cout << "selections dropped = " << n_dropped << endl; //remember the cross trigger in Data

    NumbEv-=crosstrigger;
    std::cout << "Fraction of events discarded by trigger = " << (trigger_dropped * 1. / NumbEv) << endl;
    int Rem_trigger=NumbEv-trigger_dropped; //remember the cross trigger in Data
    std::cout << "Fraction of events removed by selections = " << (n_dropped * 1. / Rem_trigger) << endl;
    std::cout << "Final number of events "<< Rem_trigger - n_dropped<<endl;

    // Write the histograms to the file
    HistWrite();

    fout->Write();
    fout->Close();
}

int main(int argc, char **argv)
{
    string inputFile = argv[1];
    string outputFile = argv[2];
    string boolstr=argv[3];
    bool IsFirstDataset= (boolstr=="true")||(boolstr=="True");
    if (IsFirstDataset) {std::cout<<"############## It is first dataset! ##################"<<std::endl;}
    if (!IsFirstDataset) {std::cout<<"@@@@@@@@@@@@@@ NOT first dataset! @@@@@@@@@@@@@@@@@@"<<std::endl;}
	
    HistIniz();

    DataAnalysis(inputFile, outputFile, IsFirstDataset);
}
