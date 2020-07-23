#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/LeptonProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

LeptonProducer::LeptonProducer(const int& era, const float& ptCut, const float& etaCut, const float& dxyCut, const float& dzCut, const float& sip3dCut, const float& isoCut, TTreeReader& reader):
	BaseProducer(&reader),
	era(era),
	ptCut(ptCut),
	etaCut(etaCut),
	dxyCut(etaCut),
	dzCut(etaCut),
	sip3dCut(etaCut),
	isoCut(etaCut)
	{}

void LeptonProducer::BeginJob(TTree* tree, bool &isData, const bool& isSyst){
	//Set data bool
	this->isData = isData;
	this->isSyst = isSyst;

	//Path to files containing Scale Factors TODO find correct files for 17,18
	muonIdScaleFactorFile = {
		{2016, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_ID.root")},
		{2017, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_ID.root")},
		{2018, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_ID.root")}
	};
	muonIsolationScaleFactorFile = {
		{2016, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_Iso.root")},
		{2017, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_Iso.root")},
		{2018, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_Iso.root")}
	};
	muonTriggerScaleFactorFile = {
		{2016, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_Trg.root")},
		{2017, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_Trg.root")},
		{2018, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/Mu_Trg.root")}
	};

	electronGSFScaleFactorFile = {
		{2016, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/EGM2D_eleGSF.root")},
		{2017, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/EGM2D_eleGSF.root")},
		{2018, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/EGM2D_eleGSF.root")}
	};
	electronMVAScaleFactorFile = {
		{2016, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/EGM2D_eleMVA90.root")},
		{2017, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/EGM2D_eleMVA90.root")},
		{2018, TString("$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/leptonSF/EGM2D_eleMVA90.root")}
	};

	TFile* muonIdSFFile        = TFile::Open(muonTriggerScaleFactorFile[era]);
	TFile* muonIsolationSFFile = TFile::Open(muonIsolationScaleFactorFile[era]);
	TFile* muonTriggerSFFile   = TFile::Open(muonTriggerScaleFactorFile[era]);

	TFile* electronGSFSFFile = TFile::Open(electronGSFScaleFactorFile[era]);
	TFile* electronMVASFFile   = TFile::Open(electronMVAScaleFactorFile[era]);

	// Possibly name difference for different years?
	muonIdSFHist        = (TH2F*) muonIdSFFile->Get("MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
	muonIsolationSFHist = (TH2F*) muonIsolationSFFile->Get("LooseISO_MediumID_pt_eta/pt_abseta_ratio");
	muonTriggerSFHist   = (TH2F*) muonTriggerSFFile->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio");

	// SF Histogram seems to be empty, calculate it yourself
	electronGSFSFHist = (TH2F*) electronGSFSFFile->Get("EGamma_SF2D");
	electronMVASFHist = (TH2F*) electronMVASFFile->Get("EGamma_SF2D");

	//Initiliaze TTreeReaderValues
	muonNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nMuon");
	muonPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_pt");
	muonEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_eta");
	muonPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_phi");
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
	electronPdgId = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_pdgId");
	electronMiniPFRelIsoAll = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_miniPFRelIso_all");
	electronCutBased = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_cutBased");



	//Set TTreeReader for genpart and trigger obj from BaseProducer
	SetCollection(this->isData);

	//Set Branches of output tree
	tree->Branch("nMuon", &nMuon);
	tree->Branch("nElectron", &nElectron);
	tree->Branch("nLepton", &nLepton);

	tree->Branch("LeptonPt", &Pt);
	tree->Branch("LeptonEta", &Eta);
	tree->Branch("LeptonPhi", &Phi);
	tree->Branch("LeptonMiniPFRelIsoAll", &MiniPFRelIsoAll);
	tree->Branch("LeptonSF", &ScaleFactor);
	//tree->Branch("Lepton", &);
}

void LeptonProducer::Produce(CutFlow& cutflow){
	//Initialize all variables as -999
	Pt = -999;
	Eta = -999;
	Phi = -999;
	MiniPFRelIsoAll = -999;
	Charge = -999;
	PdgId = -999;
	CutBased = -999;
	ScaleFactor = -999;

	nMuon = *muonNumber->Get();//muonPt->GetSize();
	nElectron = *electronNumber->Get(); //electronPt->GetSize();
	nLepton = nMuon + nElectron;

	// Create WeightCalculator class for easier acces to SF
	WeightCalculator* wc;

	if(nLepton == 1){
		if (nMuon == 1){
			const float& pt = muonPt->At(0);
			const float& eta = muonEta->At(0);
			const float& phi = muonPhi->At(0);
			const float& dxy = muonDxy->At(0);
			const float& dz = muonDz->At(0);
			const float& sip3d = muonSip3d->At(0);
			const float& miniPFRelIsoAll = muonMiniPFRelIsoAll->At(0);
			const bool& isPFCand = muonIsPFCand->At(0);

			if(pt > ptCut && abs(eta) < etaCut && dxy < dxyCut && dz < dzCut && sip3d < sip3dCut && miniPFRelIsoAll < isoCut && isPFCand){
				if(!isData){
					const float& idSF = wc->Get2DWeight(pt, eta, muonIdSFHist);
					const float& isolationSF = wc->Get2DWeight(pt, eta, muonIsolationSFHist);
					const float& triggerSF = wc->Get2DWeight(pt, eta, muonTriggerSFHist);
					ScaleFactor = idSF * isolationSF * triggerSF;
				}
				Pt = pt; Eta = eta; Phi = phi; MiniPFRelIsoAll = miniPFRelIsoAll;
				const bool& looseId = muonLooseId->At(0);
				const bool& mediumId = muonMediumId->At(0);
				const bool& tightId = muonTightId->At(0);
				if(tightId){ CutBased = 4;}
				else if (mediumId){ CutBased = 3;}
				else if (looseId){ CutBased = 2;}
				else { CutBased = 1;}
			}
		} else if (nElectron == 1){
			const float& pt = electronPt->At(0);
			const float& eta = electronEta->At(0);
			const float& phi = electronPhi->At(0);
			const float& miniPFRelIsoAll = electronMiniPFRelIsoAll->At(0);
			if(pt > ptCut && abs(eta) < etaCut && miniPFRelIsoAll < isoCut){
				if(!isData){
					const float& GSFSF = wc->Get2DWeight(pt, eta, electronGSFSFHist);
					const float& MVASF = wc->Get2DWeight(pt, eta, electronMVASFHist);
					ScaleFactor = GSFSF * MVASF;
				}
				Pt = pt; Eta = eta; Phi = phi; MiniPFRelIsoAll = miniPFRelIsoAll;
				CutBased = electronCutBased->At(0);
			}
		}
	}

	if (Pt != 999){
		std::string cutName("N_{#ell} = 1 (no ID req and Iso < 0.4)");
		cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
	} else {
		cutflow.passed = false;
	}
}

void LeptonProducer::EndJob(TFile* file){
}
