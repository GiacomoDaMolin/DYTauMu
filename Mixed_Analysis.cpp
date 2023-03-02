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
#include "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR.cc"


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

    if (crossSection < 0. || IntLuminosity < 0.){
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

    // get gen quantities
    Int_t Tau_genPartIdx[MAX_ARRAY_SIZE], Muon_genPartIdx[MAX_ARRAY_SIZE];
    Int_t GenPart_pdgId[GEN_MAX_ARRAY_SIZE], GenPart_genPartIdxMother[GEN_MAX_ARRAY_SIZE], Jet_genJetIdx[MAX_ARRAY_SIZE];
    Float_t GenPart_pt[GEN_MAX_ARRAY_SIZE];
    UChar_t Tau_genPartFlav[MAX_ARRAY_SIZE], Muon_genPartFlav[MAX_ARRAY_SIZE];
    UInt_t nGenPart;
    tin->SetBranchStatus("Muon_genPartIdx", 1);
    tin->SetBranchStatus("Muon_genPartFlav", 1);
    tin->SetBranchStatus("Tau_genPartIdx", 1);
    tin->SetBranchStatus("Tau_genPartFlav", 1);
    tin->SetBranchStatus("GenPart_pdgId", 1);
    tin->SetBranchStatus("GenPart_genPartIdxMother", 1);
    tin->SetBranchStatus("nGenPart", 1);
    tin->SetBranchStatus("Jet_genJetIdx",1);
    tin->SetBranchStatus("GenPart_pt",1);
    tin->SetBranchAddress("nGenPart", &nGenPart);
    tin->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx);
    tin->SetBranchAddress("Muon_genPartFlav", &Muon_genPartFlav);
    tin->SetBranchAddress("Tau_genPartIdx", &Tau_genPartIdx);
    tin->SetBranchAddress("Tau_genPartFlav", &Tau_genPartFlav);
    tin->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    tin->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    tin->SetBranchAddress("Jet_genJetIdx",&Jet_genJetIdx);
    tin->SetBranchAddress("GenPart_pt",&GenPart_pt);
    // collect the trigger information
     Bool_t HLT_IsoMu24;
    tin->SetBranchStatus("HLT_IsoMu24", 1);
    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);

    // collect the triggger Ids
    Int_t Tau_charge[MAX_ARRAY_SIZE];
    
    tin->SetBranchStatus("Tau_charge", 1);
    tin->SetBranchAddress("Tau_charge", &Tau_charge);

    // collect the triggger Ids
    Int_t Muon_charge[MAX_ARRAY_SIZE],Muon_nTrackerLayers[MAX_ARRAY_SIZE];
    Bool_t Muon_tightId[MAX_ARRAY_SIZE];
    Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE], Muon_sip3d[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_tightId", 1);
    tin->SetBranchStatus("Muon_charge", 1);
    tin->SetBranchStatus("Muon_pfRelIso04_all", 1);
    tin->SetBranchStatus("Muon_nTrackerLayers", 1);
    tin->SetBranchStatus("Muon_sip3d", 1);
    tin->SetBranchAddress("Muon_tightId", &Muon_tightId);
    tin->SetBranchAddress("Muon_charge", &Muon_charge);
    tin->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);
    tin->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers);
    tin->SetBranchAddress("Muon_sip3d", &Muon_sip3d);

    // Jet tagging and ID, FlavB is the recomended one, DeepB was used by Anup
    Float_t Jet_btagDeepFlavB[MAX_ARRAY_SIZE];
    UInt_t nJet;
    Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE],Jet_hadronFlavour[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
    tin->SetBranchStatus("nJet", 1);
    tin->SetBranchStatus("Jet_jetId", 1);
    tin->SetBranchStatus("Jet_puId", 1);
    tin->SetBranchStatus("Jet_hadronFlavour", 1);
    tin->SetBranchAddress("nJet", &nJet);
    tin->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);
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

    // pu stuff
    Float_t N_pu_vertices;
    tin->SetBranchStatus("Pileup_nTrueInt", 1);
    tin->SetBranchAddress("Pileup_nTrueInt", &N_pu_vertices);

    // gen weight
    Float_t genWeight;
    tin->SetBranchStatus("genWeight", 1);
    tin->SetBranchAddress("genWeight", &genWeight);

    //L1
    Float_t L1PreFiringWeight_Nom;
    tin->SetBranchStatus("L1PreFiringWeight_Nom", 1);
    tin->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom);

    int non_matching_Tau = 0, non_matching_Muon = 0;
    int n_dropped = 0;
    int trigger_dropped = 0;
    UInt_t nEv = tin->GetEntries();
    unsigned int n_events = nEv;
    TLorentzVector *Tau_p4 = new TLorentzVector();
    TLorentzVector *Muon_p4 = new TLorentzVector();
    // open correctionfiles
    
    string Tau_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/tau.json.gz";
    string muon_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/muon_Z.json.gz";
    string jets_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/jet_jmar.json";
    string b_tag_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/btagging.json.gz";
    string pileup_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/puWeights.json.gz";

    
    auto Tau_c_set = CorrectionSet::from_file(Tau_json);
    auto muon_c_set = CorrectionSet::from_file(muon_json);
    auto jet_c_set = CorrectionSet::from_file(jets_json);
    auto btag_c_set = CorrectionSet::from_file(b_tag_json);
    auto pu_c_set = CorrectionSet::from_file(pileup_json);

    auto Tau_Escale = Tau_c_set->at("tau_energy_scale");
    auto Tau_idvse = Tau_c_set->at("DeepTau2017v2p1VSe");
    auto Tau_idvsmu = Tau_c_set->at("DeepTau2017v2p1VSmu");
    auto Tau_idvsjet = Tau_c_set->at("DeepTau2017v2p1VSjet");
    //auto Tau_trig = Tau_c_set->at("");
    auto muon_trigger = muon_c_set->at("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight");
    auto muon_id = muon_c_set->at("NUM_TightID_DEN_genTracks");
    auto muon_iso = muon_c_set->at("NUM_TightRelIso_DEN_TightIDandIPCut");
    auto jet_pu = jet_c_set->at("PUJetID_eff");
    auto b_tag = btag_c_set->at("deepJet_mujets");
    auto b_mistag= btag_c_set->at("deepJet_incl"); //only for light jets
    auto pu_correction = pu_c_set->at("Collisions18_UltraLegacy_goldenJSON");

    TFile *fb_eff = new TFile("/afs/cern.ch/user/g/gdamolin/public/Beff_puLoose.root");
    TH2D * l_eff= static_cast<TH2D *>(fb_eff->Get("l_jets_tagged")); 
    TH2D * c_eff= static_cast<TH2D *>(fb_eff->Get("c_jets_tagged")); 
    TH2D * b_eff= static_cast<TH2D *>(fb_eff->Get("b_jets_tagged")); 
   
    // save the histograms in a new File
    // allow pt, inv mass, and eta to be stored in a Branch
    Float_t leading_lepton_pt, invMass, muon_eta, muon_pt, tau_eta, tau_pt;
    float Weight;
    TFile *fout = new TFile(ofile.c_str(), "RECREATE");
    
    // create a new tree for the output
    TTree *tout = new TTree("tout", "tout");
    TTree *trun_out = new TTree("Run_out", "Run_out");
    // set the branches for the output tree
    tout->Branch("leading_lepton_pt", &leading_lepton_pt);
    tout->Branch("invMass", &invMass);
    tout->Branch("muon_eta", &muon_eta);
    tout->Branch("muon_pt", &muon_pt);
    tout->Branch("tau_eta", &tau_eta);
    tout->Branch("tau_pt", &tau_pt);
    tout->Branch("Weight", &Weight);

    int Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0, Nprongs=0;
    float Acopl_etau, dR_mue;

 bool From2Taus=false, FromTau=false;

    tout->Branch("dR_mue", &dR_mue);
    tout->Branch("Acopl_etau", &Acopl_etau);
    tout->Branch("Nprongs", &Nprongs);
    tout->Branch("FromTau", &FromTau);
    tout->Branch("From2Taus", &From2Taus);

    trun_out->Branch("genEventSumw", &genEventSumw);
    trun_out->Branch("IntLumi", &IntLuminosity);
    trun_out->Branch("xs", &crossSection);
    trun_out->Branch("nEv", &n_events);

    trun_out->Fill(); // we already called trun->GetEntry(0);

    fout->cd();
    
    TRandom3 * RndGen=new TRandom3();
    RoccoR rc;
    rc.init("/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR2018UL.txt");
    #pragma omp parallel for
    for (UInt_t i = 0; i <nEv; i++){
        tin->GetEntry(i);
        if (i % 100000 == 0)
            cout << "Processing entry " << i << " of " << nEv << endl;

	if (!(HLT_IsoMu24))        { trigger_dropped++; continue; }

	bool OneProng=false, ThreeProng=false;
        Int_t Tau_idx = -1;
        for (UInt_t j = 0; j < nTau; j++){
		if (Tau_decayMode[j]>2 && Tau_decayMode[j]<10) continue;
             double ScaleE=Tau_Escale->evaluate({Tau_pt[j],abs(Tau_eta[j]),Tau_decayMode[j],Tau_genPartFlav[j],"DeepTau2017v2p1","nom"});
             Tau_p4->SetPtEtaPhiM(Tau_pt[j]*ScaleE, Tau_eta[j], Tau_phi[j], Tau_mass[j]*ScaleE);
            
            if ((Tau_p4->Pt()>22. && abs(Tau_eta[j])<2.3)&&(Tau_idDeepTau2017v2p1VSe[j]>=8 && Tau_idDeepTau2017v2p1VSmu[j]>=8 && Tau_idDeepTau2017v2p1VSjet[j]>=32)){ //Loose e- T mu T jet
		if (Tau_decayMode[j]<=2) {OneProng=true;}
		if (Tau_decayMode[j]>=10) {ThreeProng=true;}
		if (!(OneProng || ThreeProng)) {continue;}
                Tau_idx = j;
                break;
            }
        }
        if (Tau_idx==-1)  {n_dropped++; continue;}
        Weight = getWeight(IntLuminosity, crossSection, genWeight, genEventSumw);
	double Weight2=Weight;
        Weight *= pu_correction->evaluate({N_pu_vertices, "nominal"});
	Weight*= L1PreFiringWeight_Nom;
			
	Weight *=Tau_idvse->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Loose","nom"}); //Loose instead of VL
	Weight *=Tau_idvsmu->evaluate({abs(Tau_eta[Tau_idx]),Tau_genPartFlav[Tau_idx],"Tight","nom"});
	Weight *=Tau_idvsjet->evaluate({Tau_p4->Pt(),Tau_decayMode[Tau_idx],Tau_genPartFlav[Tau_idx],"Tight","nom","pt"});
  
        Int_t muon_idx = -1;
        for (UInt_t j = 0; j < nMuon; j++){ if ((abs(Muon_eta[j]) < 2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15)){
		  int NMCparticle=Muon_genPartIdx[j];
		  double scmMC;
		  if(NMCparticle>=0) {scmMC=rc.kSpreadMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],GenPart_pt[NMCparticle]);}
		  else {scmMC=rc.kSmearMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_nTrackerLayers[j],RndGen->Rndm());}

		  Muon_pt[j]*= scmMC;
		  if ( Muon_pt[j] > 27.){
		        muon_idx = j;
		        Muon_p4->SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
			break;
		    }
          }
        }
        if (muon_idx==-1)  {
            n_dropped++;
            continue;
        }
        /*if(Signal){
		bool secondistau=isFromTau(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Muon_genPartIdx[muon_idx]);
		if (secondistau){FromTau=true;}
		else {FromTau=false;}
		}*/
	if(HLT_IsoMu24) {Weight *= muon_trigger->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"});} 
        Weight *= muon_id->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"}); 
        Weight *= muon_iso->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"});
		

        bool selection = ((Tau_idx > -1) && (muon_idx > -1));

			//TODO implement jet veto
        selection = selection && (Tau_charge[Tau_idx] * Muon_charge[muon_idx]) < 0;
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
	    if((Tjet_p4->DeltaR(*Tau_p4)<0.4) || (Tjet_p4->DeltaR(*Muon_p4)<0.4)) {delete Tjet_p4; continue;}
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
              TLorentzVector *MainBjet_p4 = new TLorentzVector();
	      MainBjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	      if((MainBjet_p4->DeltaR(*Tau_p4)<0.4) || (MainBjet_p4->DeltaR(*Muon_p4)<0.4)) {delete MainBjet_p4; continue;}
	      else {delete MainBjet_p4; njets++;}
              njet_in_collection.push_back(j);
              flavor.push_back(abs(Jet_hadronFlavour[j]));
              tagged.push_back((Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp));
              if (Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp){
                 if (!one_Bjet){one_Bjet = true;}
               }
            }//end if(jetpt>50 !!puid==7)
          }//end kinematic if
        }//end for

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
	
        selection = selection && !(one_Bjet);
        if (!selection){
            n_dropped++;
            continue;
        }
	Acopl_etau=M_PI-(Muon_p4->DeltaPhi(*Tau_p4));
	if(OneProng){
		h_acopla_etau_p1->Fill(Acopl_etau,Weight);
		h_NJets_p1->Fill(njets,Weight);
		}
	if(ThreeProng){
		h_acopla_etau_p3->Fill(Acopl_etau,Weight);
		h_NJets_p3->Fill(njets,Weight);
		}

	if(OneProng) {Nprongs=1;}
	if(ThreeProng) {Nprongs=3;}

        dR_mue = Tau_p4->DeltaR(*Muon_p4);

	// fill the tree
        tau_pt = Tau_p4->Pt();
        tau_eta = Tau_p4->Eta();
        muon_pt = Muon_pt[muon_idx];
        muon_eta = Muon_eta[muon_idx];
        // check whether Tau or Muon is the leading one
        if (Tau_p4->Pt() > Muon_p4->Pt()) {leading_lepton_pt = Tau_p4->Pt();}
        else	{leading_lepton_pt = Muon_p4->Pt();}

        if(OneProng) {
		h_leading_lepton_pt_p1->Fill(leading_lepton_pt,Weight2);
		h_leading_lepton_pt_weighted_p1->Fill(leading_lepton_pt, Weight);
		h_Tau_pt_p1->Fill(tau_pt,Weight2);
		h_Tau_eta_p1->Fill(tau_eta,Weight2);
		h_Muon_pt_p1->Fill(muon_pt,Weight2);
		h_Muon_eta_p1->Fill(muon_eta,Weight2);
		h_Tau_pt_weighted_p1->Fill(tau_pt, Weight);
		h_Tau_eta_weighted_p1->Fill(tau_eta, Weight);
		h_Muon_pt_weighted_p1->Fill(muon_pt, Weight);
		h_Muon_eta_weighted_p1->Fill(muon_eta, Weight);
		h_NJets_p1->Fill(njets,Weight);
		h_Muon_sip3d_p1->Fill(Muon_sip3d[muon_idx], Weight);
		tauhole_p1->Fill(Tau_p4->Eta(),Tau_p4->Phi(),Weight);
		}

	 if(ThreeProng) {
		h_leading_lepton_pt_p3->Fill(leading_lepton_pt,Weight2);
		h_leading_lepton_pt_weighted_p3->Fill(leading_lepton_pt, Weight);
		h_Tau_pt_p3->Fill(tau_pt,Weight2);
		h_Tau_eta_p3->Fill(tau_eta,Weight2);
		h_Muon_pt_p3->Fill(muon_pt,Weight2);
		h_Muon_eta_p3->Fill(muon_eta,Weight2);
		h_Tau_pt_weighted_p3->Fill(tau_pt, Weight);
		h_Tau_eta_weighted_p3->Fill(tau_eta, Weight);
		h_Muon_pt_weighted_p3->Fill(muon_pt, Weight);
		h_Muon_eta_weighted_p3->Fill(muon_eta, Weight);
		h_NJets_p3->Fill(njets,Weight);
		h_Muon_sip3d_p3->Fill(Muon_sip3d[muon_idx], Weight);
		tauhole_p1->Fill(Tau_p4->Eta(),Tau_p4->Phi(),Weight);
		}

        if (Tau_idx > -1 && muon_idx > -1){
            invMass = (*(Tau_p4) + *(Muon_p4)).M();
	    if(OneProng){
		    h_Tau_Muon_invariant_mass_p1->Fill(invMass,Weight2);
		    h_Tau_Muon_invariant_mass_weighted_p1->Fill(invMass,Weight);
		    }
	    if(ThreeProng){
		    h_Tau_Muon_invariant_mass_p3->Fill(invMass,Weight2);
		    h_Tau_Muon_invariant_mass_weighted_p3->Fill(invMass,Weight);
		    }
        }
     
        tout->Fill();
    }


    std::cout << "non_matching_Tau = " << non_matching_Tau << endl;
    std::cout << "non_matching_Muon = " << non_matching_Muon << endl;

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
