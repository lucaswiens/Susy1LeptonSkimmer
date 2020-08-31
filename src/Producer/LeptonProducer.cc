#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/LeptonProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>

#include <cmath>
#include <random>

LeptonProducer::LeptonProducer(const int& era, const float& ptCut, const float& etaCut, const float& dxyCut, const float& dzCut, const float& sip3dCut, const float& isoCut, TTreeReader& reader):
	BaseProducer(&reader),
	era(era),
	ptCut(ptCut),
	etaCut(etaCut),
	dxyCut(dxyCut),
	dzCut(dzCut),
	sip3dCut(sip3dCut),
	isoCut(isoCut)
	{}

template <typename T>
void LeptonProducer::SortByIndex(T& var, std::vector<int> idx, unsigned int vectorSize) {
	T tmp(vectorSize);
	for (unsigned int i = 0; i < vectorSize; i++) {
		tmp.at(i) = var.at(idx[i]);
	}
	var = std::move(tmp);
}


void LeptonProducer::BeginJob(TTree* tree, bool &isData) {
	//Set data bool
	this->isData = isData;

	//Path to files containing Scale Factors TODO find correct files for 17,18
	TString muonIdSFFileLocation        = TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_ID.root");
	TString muonIsolationSFFileLocation = TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_Iso.root");
	TString muonTriggerSFFileLocation   = TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_Trg.root");
	TString electronGSFSFFileLocation   = TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/EGM2D_eleGSF.root");
	TString electronMVASFFileLocation   = TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/EGM2D_eleMVA90.root");

	muonIdSFFile        = TFile::Open(muonIdSFFileLocation, "READ");
	muonIsolationSFFile = TFile::Open(muonIsolationSFFileLocation, "READ");
	muonTriggerSFFile   = TFile::Open(muonTriggerSFFileLocation, "READ");

	electronGSFSFFile = TFile::Open(electronGSFSFFileLocation, "READ");
	electronMVASFFile = TFile::Open(electronMVASFFileLocation, "READ");

	muonIdSFHist        = static_cast<TH2F*>(muonIdSFFile->Get("MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio"));
	muonIsolationSFHist = static_cast<TH2F*>(muonIsolationSFFile->Get("LooseISO_MediumID_pt_eta/pt_abseta_ratio"));
	muonTriggerSFHist   = static_cast<TH2F*>(muonTriggerSFFile->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio"));

	// SF Histogram seems to be empty, calculate it yourself
	electronGSFSFHist = static_cast<TH2F*>(electronGSFSFFile->Get("EGamma_SF2D"));
	electronMVASFHist = static_cast<TH2F*>(electronMVASFFile->Get("EGamma_SF2D"));

	//Initiliaze TTreeReaderValues
	muonNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nMuon");
	muonPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_pt");
	muonEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_eta");
	muonPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_phi");
	muonMass = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_mass");
	muonCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_charge");
	muonDxy = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_dxy");
	muonDz = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_dz");
	muonPdgId = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_pdgId");
	muonSip3d = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_sip3d");
	muonMiniPFRelIsoAll = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_miniPFRelIso_all");

	muonIsPFCand = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_isPFcand");
	muonLooseId = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_looseId");
	muonMediumId = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_mediumId");
	muonTightId = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_tightId");//maybe not needed?

	electronNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nElectron");
	electronPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_pt");
	electronEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_eta");
	electronPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_phi");
	electronMass = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_mass");
	electronCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_charge");
	electronPdgId = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_pdgId");
	electronMiniPFRelIsoAll = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_miniPFRelIso_all");
	electronCutBased = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_cutBased");

	//Set TTreeReader for genpart and trigger obj from BaseProducer
	SetCollection(this->isData);

	//Set Branches of output tree
	tree->Branch("LeptonPt", &Pt);
	tree->Branch("LeptonEta", &Eta);
	tree->Branch("LeptonPhi", &Phi);
	tree->Branch("LeptonMass", &Mass);
	tree->Branch("LeptonMiniPFRelIsoAll", &MiniPFRelIsoAll);
	if (!isData) {
		tree->Branch("LeptonScaleFactor", &ScaleFactor);
	}

	tree->Branch("LeptonLooseId", &LooseId);
	tree->Branch("LeptonMediumId", &MediumId);
	tree->Branch("LeptonTightId", &TightId);

	tree->Branch("LeptonCharge", &Charge);
	tree->Branch("LeptonPdgId", &PdgId);
	tree->Branch("LeptonCutBased", &CutBased);
}

void LeptonProducer::Produce(CutFlow& cutflow, Susy1LeptonProduct *product) {
	//Initialize all variables as -999
	Pt.clear();
	Eta.clear();
	Phi.clear();
	Mass.clear();
	MiniPFRelIsoAll.clear();
	ScaleFactor.clear();

	LooseId.clear();
	MediumId.clear();
	TightId.clear();

	Charge.clear();
	PdgId.clear();

	CutBased.clear();

	nMuon = 0;
	nElectron = 0;
	nLepton = 0;

	nMuon = *muonNumber->Get();
	nElectron = *electronNumber->Get();
	nLepton = nMuon + nElectron;

	int muonCounter = 0, electronCounter = 0;
	if (nLepton > 0) {
		for (unsigned int i = 0; i < nMuon; i++) {
			const float& pt = muonPt->At(i);
			const float& eta = muonEta->At(i);
			const float& phi = muonPhi->At(i);
			const float& mass = muonMass->At(i);
			const float& charge = muonCharge->At(i);
			const float& dxy = muonDxy->At(i);
			const float& dz = muonDz->At(i);
			const float& sip3d = muonSip3d->At(i);
			const float& miniPFRelIsoAll = muonMiniPFRelIsoAll->At(i);

			const bool& isPFCand = muonIsPFCand->At(i);

			const int& pdgId = muonPdgId->At(i);

			if (pt > ptCut && abs(eta) < etaCut && dxy < dxyCut && dz < dzCut && sip3d < sip3dCut && miniPFRelIsoAll < isoCut && isPFCand) {
				if (!isData) {
					WeightCalculator* wc = new WeightCalculator;
					const float& idSF = wc->Get2DWeight(pt, eta, muonIdSFHist);
					const float& isolationSF = wc->Get2DWeight(pt, eta, muonIsolationSFHist);
					const float& triggerSF = wc->Get2DWeight(pt, eta, muonTriggerSFHist);
					ScaleFactor.push_back(idSF * isolationSF * triggerSF);
					delete wc;
				}

				Pt.push_back(pt); Eta.push_back(eta); Phi.push_back(phi); Mass.push_back(mass); Charge.push_back(charge); MiniPFRelIsoAll.push_back(miniPFRelIsoAll); PdgId.push_back(pdgId);
				muonCounter++;

				const bool& looseId = muonLooseId->At(i);
				const bool& mediumId = muonMediumId->At(i);
				const bool& tightId = muonTightId->At(i);
				if (tightId) { CutBased.push_back(4); TightId.push_back(true); MediumId.push_back(true); LooseId.push_back(true);}
				else if (mediumId) { CutBased.push_back(3); TightId.push_back(false); MediumId.push_back(true); LooseId.push_back(true);}
				else if (looseId) { CutBased.push_back(2); TightId.push_back(false); MediumId.push_back(false); LooseId.push_back(true);}
				else { CutBased.push_back(1); TightId.push_back(false); MediumId.push_back(false); LooseId.push_back(false);}
			}
		}
		for (unsigned int i = 0; i < nElectron; i++) {
			const float& pt = electronPt->At(i);
			const float& eta = electronEta->At(i);
			const float& phi = electronPhi->At(i);
			const float& mass = electronMass->At(i);
			const float& charge = muonCharge->At(i);
			const float& miniPFRelIsoAll = electronMiniPFRelIsoAll->At(i);

			const int & pdgId = electronPdgId->At(i);

			if (pt > ptCut && abs(eta) < etaCut && miniPFRelIsoAll < isoCut) {
				if (!isData) {
					WeightCalculator* wc = new WeightCalculator;
					const float& GSFSF = wc->Get2DWeight(pt, eta, electronGSFSFHist);
					const float& MVASF = wc->Get2DWeight(pt, eta, electronMVASFHist);
					ScaleFactor.push_back(GSFSF * MVASF);
					delete wc;
				}

				Pt.push_back(pt); Eta.push_back(eta); Phi.push_back(phi); Mass.push_back(mass); Charge.push_back(charge); MiniPFRelIsoAll.push_back(miniPFRelIsoAll); PdgId.push_back(pdgId);
				electronCounter++;

				int cutbased = electronCutBased->At(i);
				CutBased.push_back(cutbased);
				if (cutbased == 4) {TightId.push_back(true); MediumId.push_back(true); LooseId.push_back(true);}
				else if (cutbased == 3) {TightId.push_back(false); MediumId.push_back(true); LooseId.push_back(true);}
				else if (cutbased == 2) {TightId.push_back(false); MediumId.push_back(false); LooseId.push_back(true);}
				else {TightId.push_back(false); MediumId.push_back(false); LooseId.push_back(false);}

			}
		}

		//Overwrite number of leptons by excluding leptons that do not survive the cuts
		nMuon = muonCounter;
		nElectron = electronCounter;
		nLepton = Pt.size();

		//Sort Vectors according to leptonPt
		std::vector<int> idx(nLepton);
		std::iota(idx.begin(), idx.end(), 0);
		std::stable_sort(idx.begin(), idx.end(), [&](int i1, int i2) {return Pt[i1] > Pt[i2];});
		SortByIndex<std::vector<float>>(Pt, idx, nLepton);
		SortByIndex<std::vector<float>>(Eta, idx, nLepton);
		SortByIndex<std::vector<float>>(Phi, idx, nLepton);
		SortByIndex<std::vector<float>>(Mass, idx, nLepton);
		SortByIndex<std::vector<float>>(MiniPFRelIsoAll, idx, nLepton);

		SortByIndex<std::vector<bool>>(LooseId, idx, nLepton);
		SortByIndex<std::vector<bool>>(MediumId, idx, nLepton);
		SortByIndex<std::vector<bool>>(TightId, idx, nLepton);

		SortByIndex<std::vector<int>>(Charge, idx, nLepton);
		SortByIndex<std::vector<int>>(PdgId, idx, nLepton);

		SortByIndex<std::vector<unsigned int>>(CutBased, idx, nLepton);

		if (!isData) {
			SortByIndex<std::vector<float>>(ScaleFactor, idx, nLepton);
		}
	}

	if (Pt.size() != 0) {
		product->leptonPt = Pt.at(0);
		product->leptonPhi = Phi.at(0);
		product->leptonEta = Eta.at(0);
		product->leptonMass = Mass.at(0);
		product->leptonPdgId = PdgId.at(0);
		product->leptonCharge = Charge.at(0);

		std::string cutName("N_{#ell} = 1 (no ID req and Iso < 0.4)");
		cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
	} else {
		cutflow.passed = false;
	}
}

void LeptonProducer::EndJob(TFile* file) {
	muonIdSFFile->Close();
	muonIsolationSFFile->Close();
	muonTriggerSFFile->Close();

	electronGSFSFFile->Close();
	electronMVASFFile->Close();
}
