#ifndef LEPTONPRODUCER_H
#define LEPTONPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

class LeptonProducer: public BaseProducer {
	private:
		//Check if it is data or MC
		bool isData;
		//bool isSignal // or maybe make a fastsim flag

		//Lepton Scale Factor files
		std::map<int, TString> muonTriggerScaleFactorFile;
		std::map<int, TString> muonIsolationScaleFactorFile;
		std::map<int, TString> muonIdScaleFactorFile;
		std::map<int, TString> electronMVAScaleFactorFile;
		std::map<int, TString> electronGSFScaleFactorFile;

		TH2F* muonTriggerSFHist;
		TH2F* muonIsolationSFHist;
		TH2F* muonIdSFHist;

		TH2F* electronMVASFHist;
		TH2F* electronGSFSFHist;

		//Cut Variables
		int era;
		float ptCut, etaCut, dxyCut, dzCut, sip3dCut, isoCut;

		//Vector for the output variables
		float Pt, Eta, Phi, MiniPFRelIsoAll, ScaleFactor;
		bool LooseId, MediumId, TightId, IsPFCand;
		unsigned int nMuon, nElectron, nLepton, PdgId, CutBased;
		int Charge;

		//TTreeReader Values for NANO AOD analysis
		std::unique_ptr<TTreeReaderValue<unsigned int>> muonNumber, electronNumber;
		std::unique_ptr<TTreeReaderArray<float>> muonPt, muonEta, muonPhi, muonIso, muonDxy, muonDz, muonSip3d, muonMiniPFRelIsoAll;
		std::unique_ptr<TTreeReaderArray<int>> muonPdgId, muonCharge;
		std::unique_ptr<TTreeReaderArray<bool>> muonLooseId, muonMediumId, muonTightId, muonIsPFCand;
		std::unique_ptr<TTreeReaderArray<float>> electronPt, electronEta, electronPhi, electronMiniPFRelIsoAll;
		std::unique_ptr<TTreeReaderArray<int>> electronPdgId, electronCutBased, electronCharge;

	public:
		LeptonProducer(const int& era, const float& ptCut, const float& etaCut, const float& dxyCut, const float& dzCut, const float& sip3dCut, const float& isoCut, TTreeReader& reader);

		void BeginJob(TTree* tree, bool& isData, const bool& isSyst=false);
		void Produce(CutFlow& cutflow);
		void EndJob(TFile* file);
};

#endif
