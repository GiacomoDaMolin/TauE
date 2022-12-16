#ifndef Histodef_cpp
#define Histodef_cpp
#include "TH1.h"

    // Define the histogram objects with TAU with 1 PRONG
    TH1D* h_Tau_pt_p1 = new TH1D("h_Tau_pt_p1","h_Tau_pt",50,0,200);
    TH1D* h_Tau_eta_p1 = new TH1D("h_Tau_eta_p1","h_Tau_eta",50,-5,5);
    TH1D* h_Electron_pt_p1 = new TH1D("h_Electron_pt_p1", "Electron_pt", 40, 0, 200);
    TH1D* h_Electron_eta_p1 = new TH1D("h_Electron_eta_p1", "Electron_eta", 50, -5, 5);
    // define the weighted histograms
    TH1D* h_Tau_pt_weighted_p1 = new TH1D("h_Tau_pt_weighted_p1","h_Tau_pt_weighted",50,0,200);
    TH1D* h_Tau_eta_weighted_p1 = new TH1D("h_Tau_eta_weighted_p1","h_Tau_eta_weighted",50,-5,5);
    TH1D* h_Electron_pt_weighted_p1 = new TH1D("h_Electron_pt_weighted_p1", "Electron_pt_weighted", 40, 0, 200);
    TH1D* h_Electron_eta_weighted_p1 = new TH1D("h_Electron_eta_weighted_p1", "Electron_eta_weighted", 50, -5, 5);

    TH1D* h_Tau_Electron_invariant_mass_p1 = new TH1D("Tau_Electron_invariant_mass_p1", "Tau_Electron_invariant_mass", 50, 12, 412);
    TH1D* h_Tau_Electron_invariant_mass_weighted_p1 = new TH1D("Tau_Electron_invariant_mass_weighted_p1", "Tau_Electron_invariant_mass_weighted", 50, 12, 412);

    TH1D* h_leading_lepton_pt_p1 = new TH1D("leading_lepton_pt_p1", "leading_lepton_pt", 45, 20, 200);
    TH1D* h_leading_lepton_pt_weighted_p1 = new TH1D("leading_lepton_pt_weighted_p1", "leading_lepton_pt_weighted", 45, 20, 200);
    // other histos (weighted only)
    TH1D *h_NJets_p1 = new TH1D("NJets_p1","NJets",15,0,15);
    TH1D *h_LooseJets_p1 = new TH1D("N_looseJets_p1","N_looseJets",12,0,12);
    TH1D *h_MediumJets_p1 = new TH1D("N_mediumJets_p1","N_mediumJets",8,0,8);
    TH1D *h_TightJets_p1 = new TH1D("N_tightJets_p1","N_tightJets",6,0,6);

	TH1D* h_dR_allJets_p1 = new TH1D("h_dR_allJets_p1","h_dR_allJets",50,0,4);
	TH1D* h_dR_lbJets_p1 = new TH1D("h_dR_lbJets_p1","h_dR_lbJets",50,0,4);
	TH1D* h_dR_mbJets_p1 = new TH1D("h_dR_mbJets_p1","h_dR_mbJets",50,0,4);

	TH1D* h_Apl_allJets_p1 = new TH1D("h_Apl_allJets_p1","h_Apl_allJets",20,-1,1);
	TH1D* h_Apl_lbJets_p1 = new TH1D("h_Apl_lbJets_p1","h_Apl_lbJets",20,-1,1);
	TH1D* h_Apl_mbJets_p1 = new TH1D("h_Apl_mbJets_p1","h_Apl_mbJets",20,-1,1);

	TH1D* h_Phi_allJets_p1 = new TH1D("h_Phi_allJets_p1","h_Phi_allJets",32,0,3.2);
	TH1D* h_Phi_lbJets_p1 = new TH1D("h_Phi_lbJets_p1","h_dR_lbJets",32,0,3.2);
	TH1D* h_Phi_mbJets_p1 = new TH1D("h_Phi_mbJets_p1","h_dR_mbJets",32,0,3.2);

	TH1D * h_acopla_etau_p1 =new TH1D("h_Acopla_etau_p1","h_Acopla_etau",30,0, 2*M_PI);
	
	// Define the histogram objects with TAU with 3 PRONGS
   TH1D* h_Tau_pt_p3 = new TH1D("h_Tau_pt_p3","h_Tau_pt",50,0,200);
    TH1D* h_Tau_eta_p3 = new TH1D("h_Tau_eta_p3","h_Tau_eta",50,-5,5);
    TH1D* h_Electron_pt_p3 = new TH1D("h_Electron_pt_p3", "Electron_pt", 40, 0, 200);
    TH1D* h_Electron_eta_p3 = new TH1D("h_Electron_eta_p3", "Electron_eta", 50, -5, 5);
    // define the weighted histograms
    TH1D* h_Tau_pt_weighted_p3 = new TH1D("h_Tau_pt_weighted_p3","h_Tau_pt_weighted",50,0,200);
    TH1D* h_Tau_eta_weighted_p3 = new TH1D("h_Tau_eta_weighted_p3","h_Tau_eta_weighted",50,-5,5);
    TH1D* h_Electron_pt_weighted_p3 = new TH1D("h_Electron_pt_weighted_p3", "Electron_pt_weighted", 40, 0, 200);
    TH1D* h_Electron_eta_weighted_p3 = new TH1D("h_Electron_eta_weighted_p3", "Electron_eta_weighted", 50, -5, 5);

    TH1D* h_Tau_Electron_invariant_mass_p3 = new TH1D("Tau_Electron_invariant_mass_p3", "Tau_Electron_invariant_mass", 50, 12, 412);
    TH1D* h_Tau_Electron_invariant_mass_weighted_p3 = new TH1D("Tau_Electron_invariant_mass_weighted_p3", "Tau_Electron_invariant_mass_weighted", 50, 12, 412);

    TH1D* h_leading_lepton_pt_p3 = new TH1D("leading_lepton_pt_p3", "leading_lepton_pt", 45, 20, 200);
    TH1D* h_leading_lepton_pt_weighted_p3 = new TH1D("leading_lepton_pt_weighted_p3", "leading_lepton_pt_weighted", 45, 20, 200);
    // other histos (weighted only)
    TH1D *h_NJets_p3 = new TH1D("NJets_p3","NJets",15,0,15);
    TH1D *h_LooseJets_p3 = new TH1D("N_looseJets_p3","N_looseJets",12,0,12);
    TH1D *h_MediumJets_p3 = new TH1D("N_mediumJets_p3","N_mediumJets",8,0,8);
    TH1D *h_TightJets_p3 = new TH1D("N_tightJets_p3","N_tightJets",6,0,6);

	TH1D* h_dR_allJets_p3 = new TH1D("h_dR_allJets_p3","h_dR_allJets",50,0,4);
	TH1D* h_dR_lbJets_p3 = new TH1D("h_dR_lbJets_p3","h_dR_lbJets",50,0,4);
	TH1D* h_dR_mbJets_p3 = new TH1D("h_dR_mbJets_p3","h_dR_mbJets",50,0,4);

	TH1D* h_Apl_allJets_p3 = new TH1D("h_Apl_allJets_p3","h_Apl_allJets",20,-1,1);
	TH1D* h_Apl_lbJets_p3 = new TH1D("h_Apl_lbJets_p3","h_Apl_lbJets",20,-1,1);
	TH1D* h_Apl_mbJets_p3 = new TH1D("h_Apl_mbJets_p3","h_Apl_mbJets",20,-1,1);

	TH1D* h_Phi_allJets_p3 = new TH1D("h_Phi_allJets_p3","h_Phi_allJets",32,0,3.2);
	TH1D* h_Phi_lbJets_p3 = new TH1D("h_Phi_lbJets_p3","h_dR_lbJets",32,0,3.2);
	TH1D* h_Phi_mbJets_p3 = new TH1D("h_Phi_mbJets_p3","h_dR_mbJets",32,0,3.2);

	TH1D * h_acopla_etau_p3 =new TH1D("h_Acopla_etau_p3","h_Acopla_etau",30,0, 2*M_PI);
	
	TH1D * h_Electron_pt_from_W= new TH1D("h_Electron_pt_from_W", "h_Electron_pt_from_W", 40, 0, 200);
        TH1D * h_Electron_eta_from_W= new TH1D("h_Electron_eta_from_W", "Electron_eta_from_W", 50, -5, 5);
        TH1D * h_Electron_pt_weighted_from_W= new TH1D("h_Electron_pt_from_W_weighted", "h_Electron_pt_from_W_weighted", 40, 0, 200);
        TH1D * h_Electron_eta_weighted_from_W= new TH1D("h_Electron_eta_from_W_weighted", "Electron_eta_from_W_weighted", 50, -5, 5);
        
	TH1D * h_Tau_pt_from_W= new TH1D("h_Tau_pt_p1","h_Tau_pt",50,0,200);
        TH1D * h_Tau_eta_from_W= new TH1D("h_Tau_eta_from_W", "Tau_eta_from_W", 50, -5, 5);
        TH1D * h_Tau_pt_weighted_from_W= new TH1D("h_Tau_pt_from_W_weighted", "h_Tau_pt_from_W_weighted", 50, 0, 200);
        TH1D * h_Tau_eta_weighted_from_W= new TH1D("h_Tau_eta_from_W_weighted", "Tau_eta_from_W_weighted", 50, -5, 5);
	

#endif // Histodef_cpp
