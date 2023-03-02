#ifndef Auxiliary_cpp
#define Auxiliary_cpp
#include <iostream>
#include"Histodef.cpp"
// inputs are: size of arrays GenID and GenPar
// the pointer at the start of the array of GenParticles_PDGID
// the pointer at the start of the array of GenParticles_Parent
// initialID: last index of the particle (ex: Muon_genPartIdx[j]=initialID)

// should be called from func as isFromW(MAX_ARRAY_SIZE,GenPart_pdgId,GenPart_genPartIdxMother,Muon_genPartIdx[j])
//arguments: (array of PdgID of Gen Particles (Int_t), array of GenParent of Gen Particles (Int_t, contains the index to the parent of the selected GenParticle
// InitialID  is the index of the starting muon in the GenParticles array)
bool isFromW(int size, Int_t *GenId, Int_t *GenParent, int initialID)
{
	if (initialID < 0)
	{
		return false;
	}
	// retrieve first PDG ID number
	int startPdg = GenId[initialID];
	int newID = initialID, newPdg = startPdg;
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg)
	{
		if (newID > size)
		{
			std::cout << "WARNING: index " << newID << " exceeding max size " << size << std::endl;
		}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		if (abs(newPdg) == 24)
			return true;
	}
	return false;
}

bool isFromTau(int size, Int_t *GenId, Int_t *GenParent, int initialID)
{
	if (initialID < 0)
	{
		return false;
	}
	// retrieve first PDG ID number
	int startPdg = GenId[initialID];
	int newID = initialID, newPdg = startPdg;
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg)
	{
		if (newID > size)
		{
			std::cout << "WARNING: index " << newID << " exceeding max size " << size << std::endl;
		}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		if (abs(newPdg) == 15)
			return true;
	}
	return false;
}

double InvertPhi(double phi){
 double invphi=phi+M_PI;
 if (invphi>M_PI){invphi=invphi-2*M_PI;}
 return invphi;
}

void HistIniz(){
// Define the histogram objects with TAU with 1 PRONG
     h_Tau_pt_p1->Sumw2();
     h_Tau_eta_p1->Sumw2();
     h_Muon_pt_p1->Sumw2();
     h_Muon_eta_p1->Sumw2();
     h_Tau_pt_weighted_p1->Sumw2();
     h_Tau_eta_weighted_p1->Sumw2();
     h_Muon_pt_weighted_p1->Sumw2();
     h_Muon_eta_weighted_p1->Sumw2();
     h_Tau_Muon_invariant_mass_p1->Sumw2();
     h_Tau_Muon_invariant_mass_weighted_p1->Sumw2();
     h_leading_lepton_pt_p1->Sumw2();
     h_leading_lepton_pt_weighted_p1->Sumw2();
h_NJets_p1->Sumw2();
h_LooseJets_p1->Sumw2();
h_MediumJets_p1->Sumw2();
h_TightJets_p1->Sumw2();
	 h_dR_allJets_p1->Sumw2();
	 h_dR_lbJets_p1->Sumw2();
	 h_dR_mbJets_p1->Sumw2();
	 h_Apl_allJets_p1->Sumw2();
	 h_Apl_lbJets_p1->Sumw2();
	 h_Apl_mbJets_p1->Sumw2();
	 h_Phi_allJets_p1->Sumw2();
	 h_Phi_lbJets_p1->Sumw2();
	 h_Phi_mbJets_p1->Sumw2();
	h_acopla_etau_p1->Sumw2();
	
	// Define the histogram objects with TAU with 3 PRONGS
    h_Tau_pt_p3->Sumw2();
     h_Tau_eta_p3->Sumw2();
     h_Muon_pt_p3->Sumw2();
     h_Muon_eta_p3->Sumw2();
     h_Tau_pt_weighted_p3->Sumw2();
     h_Tau_eta_weighted_p3->Sumw2();
     h_Muon_pt_weighted_p3->Sumw2();
     h_Muon_eta_weighted_p3->Sumw2();
     h_Tau_Muon_invariant_mass_p3->Sumw2();
     h_Tau_Muon_invariant_mass_weighted_p3->Sumw2();
     h_leading_lepton_pt_p3->Sumw2();
     h_leading_lepton_pt_weighted_p3->Sumw2();
h_NJets_p3->Sumw2();
h_LooseJets_p3->Sumw2();
h_MediumJets_p3->Sumw2();
h_TightJets_p3->Sumw2();
	 h_dR_allJets_p3->Sumw2();
	 h_dR_lbJets_p3->Sumw2();
	 h_dR_mbJets_p3->Sumw2();
	 h_Apl_allJets_p3->Sumw2();
	 h_Apl_lbJets_p3->Sumw2();
	 h_Apl_mbJets_p3->Sumw2();
	 h_Phi_allJets_p3->Sumw2();
	 h_Phi_lbJets_p3->Sumw2();
	 h_Phi_mbJets_p3->Sumw2();
h_acopla_etau_p3->Sumw2();

h_Muon_sip3d_p1->Sumw2();
h_Muon_sip3d_p3->Sumw2();

jethole_p1->Sumw2();
tauhole_p1->Sumw2();

jethole_p3->Sumw2();
tauhole_p3->Sumw2();
}

void HistWrite(){
h_Tau_pt_weighted_p1->SetBinContent(h_Tau_pt_weighted_p1->GetNbinsX(), h_Tau_pt_weighted_p1->GetBinContent(h_Tau_pt_weighted_p1->GetNbinsX()) + h_Tau_pt_weighted_p1->GetBinContent(h_Tau_pt_weighted_p1->GetNbinsX() + 1));
h_Tau_pt_weighted_p3->SetBinContent(h_Tau_pt_weighted_p3->GetNbinsX(), h_Tau_pt_weighted_p3->GetBinContent(h_Tau_pt_weighted_p3->GetNbinsX()) + h_Tau_pt_weighted_p3->GetBinContent(h_Tau_pt_weighted_p3->GetNbinsX() + 1));
h_Muon_pt_weighted_p1->SetBinContent(h_Muon_pt_weighted_p1->GetNbinsX(), h_Muon_pt_weighted_p1->GetBinContent(h_Muon_pt_weighted_p1->GetNbinsX()) + h_Muon_pt_weighted_p1->GetBinContent(h_Muon_pt_weighted_p1->GetNbinsX() + 1));
h_Muon_pt_weighted_p3->SetBinContent(h_Muon_pt_weighted_p3->GetNbinsX(), h_Muon_pt_weighted_p3->GetBinContent(h_Muon_pt_weighted_p3->GetNbinsX()) + h_Muon_pt_weighted_p3->GetBinContent(h_Muon_pt_weighted_p3->GetNbinsX() + 1));
h_Tau_Muon_invariant_mass_weighted_p1->SetBinContent(h_Tau_Muon_invariant_mass_weighted_p1->GetNbinsX(), h_Tau_Muon_invariant_mass_weighted_p1->GetBinContent(h_Tau_Muon_invariant_mass_weighted_p1->GetNbinsX()) + h_Tau_Muon_invariant_mass_weighted_p1->GetBinContent(h_Tau_Muon_invariant_mass_weighted_p1->GetNbinsX() + 1));
h_Tau_Muon_invariant_mass_weighted_p3->SetBinContent(h_Tau_Muon_invariant_mass_weighted_p3->GetNbinsX(), h_Tau_Muon_invariant_mass_weighted_p3->GetBinContent(h_Tau_Muon_invariant_mass_weighted_p3->GetNbinsX()) + h_Tau_Muon_invariant_mass_weighted_p3->GetBinContent(h_Tau_Muon_invariant_mass_weighted_p3->GetNbinsX() + 1));
// Define the histogram objects with TAU with 1 PRONG
     h_Tau_pt_p1->Write();
     h_Tau_eta_p1->Write();
     h_Muon_pt_p1->Write();
     h_Muon_eta_p1->Write();
     h_Tau_pt_weighted_p1->Write();
     h_Tau_eta_weighted_p1->Write();
     h_Muon_pt_weighted_p1->Write();
     h_Muon_eta_weighted_p1->Write();
     h_Tau_Muon_invariant_mass_p1->Write();
     h_Tau_Muon_invariant_mass_weighted_p1->Write();
     h_leading_lepton_pt_p1->Write();
     h_leading_lepton_pt_weighted_p1->Write();
h_NJets_p1->Write();
h_LooseJets_p1->Write();
h_MediumJets_p1->Write();
h_TightJets_p1->Write();
	 h_dR_allJets_p1->Write();
	 h_dR_lbJets_p1->Write();
	 h_dR_mbJets_p1->Write();
	 h_Apl_allJets_p1->Write();
	 h_Apl_lbJets_p1->Write();
	 h_Apl_mbJets_p1->Write();
	 h_Phi_allJets_p1->Write();
	 h_Phi_lbJets_p1->Write();
	 h_Phi_mbJets_p1->Write();
	h_acopla_etau_p1->Write();
	
	// Define the histogram objects with TAU with 3 PRONGS
    h_Tau_pt_p3->Write();
     h_Tau_eta_p3->Write();
     h_Muon_pt_p3->Write();
     h_Muon_eta_p3->Write();
     h_Tau_pt_weighted_p3->Write();
     h_Tau_eta_weighted_p3->Write();
     h_Muon_pt_weighted_p3->Write();
     h_Muon_eta_weighted_p3->Write();
     h_Tau_Muon_invariant_mass_p3->Write();
     h_Tau_Muon_invariant_mass_weighted_p3->Write();
     h_leading_lepton_pt_p3->Write();
     h_leading_lepton_pt_weighted_p3->Write();
h_NJets_p3->Write();
h_LooseJets_p3->Write();
h_MediumJets_p3->Write();
h_TightJets_p3->Write();
	 h_dR_allJets_p3->Write();
	 h_dR_lbJets_p3->Write();
	 h_dR_mbJets_p3->Write();
	 h_Apl_allJets_p3->Write();
	 h_Apl_lbJets_p3->Write();
	 h_Apl_mbJets_p3->Write();
	 h_Phi_allJets_p3->Write();
	 h_Phi_lbJets_p3->Write();
	 h_Phi_mbJets_p3->Write();
h_acopla_etau_p3->Write();

 h_Muon_sip3d_p1->Write();
 h_Muon_sip3d_p3->Write();

//jethole_p1->Write();
tauhole_p1->Write();

//jethole_p3->Write();
tauhole_p3->Write();

}
#endif // Auxiliary_cpp
