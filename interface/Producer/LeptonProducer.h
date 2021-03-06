#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

class LeptonProducer : public BaseProducer {
	private:
		//Check if it is data or MC
		bool isData, doSystematics;
		//bool isSignal // or maybe make a fastsim flag

		//Scale Factor Files
		TFile *muonIdSFFile;
		TFile *muonIsolationSFFile;
		TFile *muonTriggerSFFile;

		TFile *electronGSFSFFile;
		TFile *electronMVASFFile;
		//Histograms containing scalefactors
		TH2F *muonTriggerSFHist;
		TH2F *muonIsolationSFHist;
		TH2F *muonIdSFHist;

		TH2F *electronMVASFHist;
		TH2F *electronGSFSFHist;

		//Cut Variables
		int era;
		float ptCut, etaCut, dxyCut, dzCut, sip3dCut, isoCut;

		//Vector for the output variables
		std::vector<float> Pt, Eta, Phi, Mass, MiniPFRelIsoAll, dileptonMass,
			ScaleFactorId, ScaleFactorTrigger, ScaleFactorIsolation, ScaleFactorGSF, ScaleFactorMVA,
			ScaleFactorIdUp, ScaleFactorTriggerUp, ScaleFactorIsolationUp, ScaleFactorGSFUp, ScaleFactorMVAUp,
			ScaleFactorIdDown, ScaleFactorIsolationDown, ScaleFactorTriggerDown, ScaleFactorGSFDown, ScaleFactorMVADown;
		std::vector<bool> LooseId, MediumId, TightId;
		std::vector<int> Charge, PdgId;
		std::vector<unsigned int> CutBased;
		unsigned int nMuon, nElectron, nLepton;//, CutBased;

		//TTreeReader Values for NANO AOD analysis
		std::unique_ptr<TTreeReaderValue<unsigned int>> muonNumber, electronNumber;
		std::unique_ptr<TTreeReaderArray<float>> muonPt, muonEta, muonPhi, muonMass, muonIso, muonDxy, muonDz, muonSip3d, muonMiniPFRelIsoAll;
		std::unique_ptr<TTreeReaderArray<int>> muonPdgId, muonCharge;
		std::unique_ptr<TTreeReaderArray<bool>> muonLooseId, muonMediumId, muonTightId, muonIsPFCand;
		std::unique_ptr<TTreeReaderArray<float>> electronPt, electronEta, electronPhi, electronMass, electronMiniPFRelIsoAll;
		std::unique_ptr<TTreeReaderArray<int>> electronPdgId, electronCutBased, electronCharge;

		template <typename T>
		void SortByIndex(T &var, std::vector<int> idx, unsigned int vectorSize);
	public:
		LeptonProducer(const int &era, const float &ptCut, const float &etaCut, const float &dxyCut, const float &dzCut, const float &sip3dCut, const float &isoCut, TTreeReader &reader);

		void BeginJob(TTree *tree, bool &isData, bool &doSystematics);
		void Produce(CutFlow &cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile *file);
};
