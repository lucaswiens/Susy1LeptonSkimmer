#ifndef SUSY1LEPTONPRODUCT_H
#define SUSY1LEPTONPRODUCT_H

#include <TTree.h>
#include <TFile.h>

class Susy1LeptonProduct {
	private:
		int era;
		bool preVFP, isData, isUp;
		std::string eraSelector, sampleName;
		char runPeriod;
		double xSection, luminosity;
		//std::shared_ptr<TTree> metaData;
	public:
		Susy1LeptonProduct(const int &era, const bool &isData, const std::string &sampleName, const char &runPeriod, const double &xSection, TFile &outputFile);
		void RegisterOutput(std::vector<std::shared_ptr<TTree>> outputTrees);
		int GetEra() {return era;}
		std::string GetEraSelector() { return eraSelector;}
		bool GetIsPreVFP() {return preVFP;}
		bool GetIsData() {return isData;}
		char GetRunPeriod() {return runPeriod;}

		// Max value for static arrays, should use assert to enforce nObject < nMax
		static const std::size_t nMax = 20;

		// Lepton Information
		int nLepton, nGoodLepton, nVetoLepton;

		// Muons Inforamtion
		int nMuon, nGoodMuon, nVetoMuon, nAntiSelectedMuon;
		std::array<double, nMax> muonPt, muonEta, muonPhi, muonMass, muonMiniIso, muonDxy, muonDz, muonSip3d,
			muonLooseIsoSf, muonLooseIsoSfUp, muonLooseIsoSfDown,
			muonTightIsoSf, muonTightIsoSfUp, muonTightIsoSfDown,
			muonLooseSf, muonLooseSfUp, muonLooseSfDown,
			muonMediumSf, muonMediumSfUp, muonMediumSfDown,
			muonTightSf, muonTightSfUp, muonTightSfDown,
			muonTriggerSf, muonTriggerSfUp, muonTriggerSfDown;
		std::array<int, nMax> muonPdgId, muonCharge, muonCutBasedId, muonGenMatchedIndex;
		std::array<bool, nMax> muonTightId, muonMediumId, muonLooseId,
			muonMvaId, muonIsPfCand,
			muonIsGood, muonIsVeto, muonIsAntiSelected;
		std::vector<double> muonPtVector;

		int nElectron, nGoodElectron, nVetoElectron, nAntiSelectedElectron;
		std::array<double, nMax> electronPt, electronEta, electronPhi, electronMass,
			electronDxy, electronDz,
			electronECorr,
			electronMiniIso, electronIso03, electronIso04, electronRelJetIso,
			electronEnergyScaleUp, electronEnergyScaleDown,
			electronEnergySigmaUp, electronEnergySigmaDown,
			electronRecoSf, electronRecoSfUp, electronRecoSfDown,
			electronVetoSf, electronVetoSfUp, electronVetoSfDown,
			electronLooseSf, electronLooseSfUp, electronLooseSfDown,
			electronMediumSf, electronMediumSfUp, electronMediumSfDown,
			electronTightSf, electronTightSfUp, electronTightSfDown,
			electronMediumMvaSf, electronMediumMvaSfUp, electronMediumMvaSfDown,
			electronTightMvaSf, electronTightMvaSfUp, electronTightMvaSfDown;
		std::array<int, nMax> electronCharge, electronCutBasedId, electronNLostHits;
		std::array<bool, nMax> electronLooseMvaId, electronMediumMvaId, electronTightMvaId,
			electronTightId, electronMediumId, electronLooseId, electronVetoId,
			electronIsGood, electronIsVeto, electronIsAntiSelected,
			electronConvVeto;

		int nJet,
			nDeepCsvBTag, nDeepJetBTag, jetId;
		double rho, metPt, metPhi, wBosonMinMass, wBosonMinMassPt, wBosonBestMass, wBosonBestMassPt, topBestMass, topBestMassPt;
		std::array<double, nMax> jetPt, jetEta, jetPhi, jetMass,
			jetArea, jetRawFactor,
			jetDeepCsv, jetDeepJet;
		std::array<bool, nMax> jetDeepCsvTightId, jetDeepCsvMediumId, jetDeepCsvLooseId,
			jetDeepJetTightId, jetDeepJetMediumId, jetDeepJetLooseId;

		int nFatJet;
		std::array<int, nMax> fatJetId;
		std::array<double, nMax> fatJetMass, fatJetPt, fatJetEta, fatJetPhi,
			fatJetArea, fatJetRawFactor,
			fatJetDeepTagMDTvsQCD, fatJetDeepTagMDWvsQCD,
			fatJetDeepTagTvsQCD, fatJetDeepTagWvsQCD;

};

#endif

