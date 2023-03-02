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
    Float_t Tau_pt[MAX_ARRAY_SIZE], Muon_pt[MAX_ARRAY_SIZE], Jet_pt[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_pt", 1);
    tin->SetBranchAddress("Muon_pt", &Muon_pt);
    tin->SetBranchStatus("Tau_pt", 1);
    tin->SetBranchAddress("Tau_pt", &Tau_pt);
    tin->SetBranchStatus("Jet_pt", 1);
    tin->SetBranchAddress("Jet_pt", &Jet_pt);
    // get the number of Taus, Muons
    UInt_t nTau, nMuon;
    tin->SetBranchStatus("nMuon", 1);
    tin->SetBranchAddress("nMuon", &nMuon);
    tin->SetBranchStatus("nTau", 1);
    tin->SetBranchAddress("nTau", &nTau);
    // get the eta
    Float_t Tau_eta[MAX_ARRAY_SIZE], Muon_eta[MAX_ARRAY_SIZE], Jet_eta[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_eta", 1);
    tin->SetBranchAddress("Muon_eta", &Muon_eta);
    tin->SetBranchStatus("Tau_eta", 1);
    tin->SetBranchAddress("Tau_eta", &Tau_eta);
    tin->SetBranchStatus("Jet_eta", 1);
    tin->SetBranchAddress("Jet_eta", &Jet_eta);
    // get the phi
    Float_t Tau_phi[MAX_ARRAY_SIZE], Muon_phi[MAX_ARRAY_SIZE], Jet_phi[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_phi", 1);
    tin->SetBranchAddress("Muon_phi", &Muon_phi);
    tin->SetBranchStatus("Tau_phi", 1);
    tin->SetBranchAddress("Tau_phi", &Tau_phi);
    tin->SetBranchStatus("Jet_phi", 1);
    tin->SetBranchAddress("Jet_phi", &Jet_phi);
    // get the mass
    Float_t Tau_mass[MAX_ARRAY_SIZE], Muon_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_mass", 1);
    tin->SetBranchAddress("Muon_mass", &Muon_mass);
    tin->SetBranchStatus("Tau_mass", 1);
    tin->SetBranchAddress("Tau_mass", &Tau_mass);
    tin->SetBranchStatus("Jet_mass", 1);
    tin->SetBranchAddress("Jet_mass", &Jet_mass);

    Bool_t HLT_IsoMu24;
    tin->SetBranchStatus("HLT_IsoMu24", 1);
    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);

    // collect the triggger Ids
    Int_t Tau_charge[MAX_ARRAY_SIZE], Muon_charge[MAX_ARRAY_SIZE];
    Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE],Muon_sip3d[MAX_ARRAY_SIZE];
    Bool_t Muon_tightId[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Tau_charge", 1);
    tin->SetBranchStatus("Muon_charge", 1);
    tin->SetBranchStatus("Muon_pfRelIso04_all", 1);
    tin->SetBranchStatus("Muon_tightId", 1);
    tin->SetBranchStatus("Muon_sip3d", 1);
    tin->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);
    tin->SetBranchAddress("Tau_charge", &Tau_charge);
    tin->SetBranchAddress("Muon_charge", &Muon_charge);
    tin->SetBranchAddress("Muon_tightId", &Muon_tightId);
    tin->SetBranchAddress("Muon_sip3d", &Muon_sip3d);

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


    int non_matching_Tau = 0, non_matching_Muon = 0;
    int n_dropped = 0;
    int trigger_dropped = 0,crosstrigger=0;
    const auto nEv = tin->GetEntries();
    TLorentzVector *Tau_p4 = new TLorentzVector();
    TLorentzVector *Muon_p4 = new TLorentzVector();

    RoccoR rc;
    rc.init("/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR2018UL.txt");


       // allow pt, inv mass, and eta to be stored in a Branch
    Float_t leading_lepton_pt, invMass, muon_eta, muon_pt, tau_eta, tau_pt;
    TFile *fout =new TFile(ofile.c_str(),"RECREATE");
    // create a new tree for the output
    TTree *tout = new TTree("tout","tout");
    // set the branches for the output tree
    tout->Branch("leading_lepton_pt", &leading_lepton_pt);
    tout->Branch("invMass", &invMass);
    tout->Branch("muon_eta", &muon_eta);
    tout->Branch("muon_pt", &muon_pt);
    tout->Branch("tau_eta", &tau_eta);
    tout->Branch("tau_pt", &tau_pt);

    int Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0, Nprongs=0;
    float Acopl_etau,dR_mue;
    
    tout->Branch("dR_mue", &dR_mue);
    tout->Branch("Nloose", &Nloose);
    tout->Branch("Acopl_etau", &Acopl_etau);
    tout->Branch("Nprongs", &Nprongs);

    #pragma omp parallel for
    for (UInt_t i = 0; i < nEv; i++)
    {
        tin->GetEntry(i);
        if (i % 100000 == 0)
            std::cout << "Processing entry " << i << " of " << nEv << std::endl;

        if (!(HLT_IsoMu24))        { trigger_dropped++; continue; }

        // loop over the Taus and Muons and only keep the fist ones that pass the requirements
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

        Int_t muon_idx = -1;
        for (UInt_t j = 0; j < nMuon; j++){
            if ((abs(Muon_eta[j])<2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15)){
		double scmDT=rc.kScaleDT(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j]);
       		Muon_pt[j]*= scmDT;
		if (Muon_pt[j]>27.){
		        muon_idx = j;
			Muon_p4->SetPtEtaPhiM(Muon_pt[muon_idx], Muon_eta[muon_idx], Muon_phi[muon_idx], Muon_mass[muon_idx]);
			break;
                }
            }
        }
        if (muon_idx==-1) {
            n_dropped++;
            continue;
        }
	if(Muon_p4->DeltaR(*Tau_p4)<0.4) {
		cout<<"DeltaR removed"<<endl; n_dropped++; continue;
		}
        bool selection = ((Tau_idx > -1) && (muon_idx > -1));
        // check the seleected objects for opposite charge
        selection = selection && (Tau_charge[Tau_idx] * Muon_charge[muon_idx]) < 0;
        
        // the tight working point is 0.71, medium 0.2783, loose 0.0490
        Float_t jet_btag_deepFlav_wp = 0.2783;
        bool one_Bjet = false;
        int id_m_jet=-1;
	int njets=0,nbjet=0;
	Nloose=0, Nmedium=0, Ntight=0,JetsNotB=0;
        for (size_t j = 0; j < nJet; j++){
          if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || (Jet_puId[j]>=4))){
	    TLorentzVector *Tjet_p4 = new TLorentzVector();
	    Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	    if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Muon_p4)<0.4)) {delete Tjet_p4; continue;}
	    else {delete Tjet_p4;}
	    njets++;
	    if (Jet_btagDeepFlavB[j] < 0.0490) JetsNotB++;
	    if (Jet_btagDeepFlavB[j] > 0.0490)	Nloose++;
	    if (Jet_btagDeepFlavB[j] > 0.2783)	Nmedium++;
	    if (Jet_btagDeepFlavB[j] > 0.71)	Ntight++;
            if (Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp){
                		one_Bjet = true; nbjet++;	
            } //end tag
	  }//end kinematic if
        }

        selection = selection && !(one_Bjet);

        if (!selection){
            n_dropped++;
            continue;
        }
	Acopl_etau=M_PI-(Muon_p4->DeltaPhi(*Tau_p4));

	if(OneProng){
		h_LooseJets_p1->Fill(Nloose);
		h_acopla_etau_p1->Fill(Acopl_etau);
		h_NJets_p1->Fill(njets);
		}
        if(ThreeProng){
		h_LooseJets_p3->Fill(Nloose);
		h_acopla_etau_p3->Fill(Acopl_etau);
		h_NJets_p3->Fill(njets);
		}
	if(OneProng) {Nprongs=1;}
	if(ThreeProng) {Nprongs=3;}

	dR_mue=Tau_p4->DeltaR(*Muon_p4);
	// fill the tree
        tau_pt = Tau_pt[Tau_idx];
        tau_eta = Tau_eta[Tau_idx];
        muon_pt = Muon_pt[muon_idx];
        muon_eta = Muon_eta[muon_idx];

        // check whether Tau or Muon is the leading one
        if (Tau_p4->Pt() > Muon_p4->Pt())	{leading_lepton_pt = Tau_p4->Pt();}
	else					{leading_lepton_pt = Muon_p4->Pt();}


	if(OneProng){
		h_leading_lepton_pt_p1->Fill(leading_lepton_pt);
		h_leading_lepton_pt_weighted_p1->Fill(leading_lepton_pt);
		h_Tau_pt_p1->Fill(tau_pt);
		h_Tau_eta_p1->Fill(tau_eta);
		h_Muon_pt_p1->Fill(muon_pt);
		h_Muon_eta_p1->Fill(muon_eta);
		h_Tau_pt_weighted_p1->Fill(tau_pt);
		h_Tau_eta_weighted_p1->Fill(tau_eta);
		h_Muon_pt_weighted_p1->Fill(muon_pt);
		h_Muon_eta_weighted_p1->Fill(muon_eta);
		h_NJets_p1->Fill(njets);
		h_Muon_sip3d_p1->Fill(Muon_sip3d[muon_idx]);
		tauhole_p1->Fill(Tau_p4->Eta(),Tau_p4->Phi());
		}
	if(ThreeProng){
		h_leading_lepton_pt_p3->Fill(leading_lepton_pt);
		h_leading_lepton_pt_weighted_p3->Fill(leading_lepton_pt);
		h_Tau_pt_p3->Fill(tau_pt);
		h_Tau_eta_p3->Fill(tau_eta);
		h_Muon_pt_p3->Fill(muon_pt);
		h_Muon_eta_p3->Fill(muon_eta);
		h_Tau_pt_weighted_p3->Fill(tau_pt);
		h_Tau_eta_weighted_p3->Fill(tau_eta);
		h_Muon_pt_weighted_p3->Fill(muon_pt);
		h_Muon_eta_weighted_p3->Fill(muon_eta);
		h_NJets_p3->Fill(njets);
		h_Muon_sip3d_p3->Fill(Muon_sip3d[muon_idx]);		
		tauhole_p3->Fill(Tau_p4->Eta(),Tau_p4->Phi());
		}

        if (Tau_idx > -1 && muon_idx > -1)
        {
            invMass = (*(Tau_p4) + *(Muon_p4)).M();
	    if(OneProng){
		    h_Tau_Muon_invariant_mass_p1->Fill(invMass);
		    h_Tau_Muon_invariant_mass_weighted_p1->Fill(invMass);
		    }
	    if(ThreeProng){
		    h_Tau_Muon_invariant_mass_p3->Fill(invMass);
		    h_Tau_Muon_invariant_mass_weighted_p3->Fill(invMass);
		    }
        }
	tout->Fill();
    }
    std::cout << "Total number of events: " << nEv << std::endl;
    int NumbEv=nEv;
    std::cout << "Removed because in another sample = " << crosstrigger << endl;
    std::cout << "trigger dropped = " << trigger_dropped << endl;
    std::cout << "selections dropped = " << n_dropped << endl; //remember the cross trigger in Data

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
