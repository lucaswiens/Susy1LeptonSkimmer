#ifndef SUSY1LEPTONPRODUCT_H
#define SUSY1LEPTONPRODUCT_H

#include <TTree.h>
#include <TFile.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

class Susy1LeptonProduct {
	private:
		int era;
		bool preVFP, isData, isFastSim;
		std::string eraSelector, sampleName;
		char runPeriod;
		float xSection, luminosity;

		TTree metaData;
	public:
		Susy1LeptonProduct(const int &era, const bool &isData, const bool &isFastSim, const std::string &sampleName, const char &runPeriod, const float &xSection, const pt::ptree &configTree, TFile &outputFile);
		void RegisterTrigger(const std::vector<std::string> &triggerNames, const std::vector<std::string> &metFilterNames, const std::vector<std::shared_ptr<TTree>> &outputTrees);
		void RegisterMetFilter(const std::vector<std::string> &filterNames, const std::vector<std::shared_ptr<TTree>> &outputTrees);
		void RegisterOutput(std::vector<std::shared_ptr<TTree>> outputTrees, const pt::ptree &configTree);
		void WriteMetaData(TFile &outputFile);
		int GetEra() {return era;}
		std::string GetEraSelector() { return eraSelector;}
		bool GetIsPreVFP() {return preVFP;}
		bool GetIsData() {return isData;}
		bool GetIsFastSim() {return isFastSim;}
		char GetRunPeriod() {return runPeriod;}

		// Max value for static arrays, should use assert to enforce nObject < nMax, so you know when to increase nMax
		static const int nMax = 50;

		// Lepton Information
		int nLepton, nGoodLepton, nVetoLepton;

		// Muons Inforamtion
		int nMuon, nGoodMuon, nVetoMuon, nAntiSelectedMuon;
		std::array<float, nMax> muonPt, muonEta, muonPhi, muonMass, muonMiniIso, muonDxy, muonDz, muonSip3d,
			muonLooseIsoSf, muonLooseIsoSfUp, muonLooseIsoSfDown,
			muonMediumIsoSf, muonMediumIsoSfUp, muonMediumIsoSfDown,
			muonLooseSf, muonLooseSfUp, muonLooseSfDown,
			muonMediumSf, muonMediumSfUp, muonMediumSfDown,
			muonTightSf, muonTightSfUp, muonTightSfDown,
			muonTriggerSf, muonTriggerSfUp, muonTriggerSfDown,
			muonLooseFastSf, muonLooseFastSfUp, muonLooseFastSfDown,
			muonMediumFastSf, muonMediumFastSfUp, muonMediumFastSfDown;
		std::array<int, nMax> muonPdgId, muonCharge, muonCutBasedId, muonGenMatchedIndex;
		std::array<bool, nMax> muonTightId, muonMediumId, muonLooseId,
			muonMvaId, muonIsPfCand,
			muonIsGood, muonIsVeto, muonIsAntiSelected;

		// Electron Inforamtion
		int nElectron, nGoodElectron, nVetoElectron, nAntiSelectedElectron;
		std::array<float, nMax> electronPt, electronEta, electronPhi, electronMass,
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
			electronTightMvaSf, electronTightMvaSfUp, electronTightMvaSfDown,
			electronVetoFastSf, electronVetoFastSfUp, electronVetoFastSfDown,
			electronTightFastSf, electronTightFastSfUp, electronTightFastSfDown,
			electronVetoMvaFastSf, electronVetoMvaFastSfUp, electronVetoMvaFastSfDown,
			electronTightMvaFastSf,electronTightMvaFastSfUp,electronTightMvaFastSfDown;
		std::array<int, nMax> electronCharge, electronCutBasedId, electronNLostHits;
		std::array<short, nMax> electronLooseMvaId, electronMediumMvaId, electronTightMvaId,
			electronTightId, electronMediumId, electronLooseId, electronVetoId,
			electronIsGood, electronIsVeto, electronIsAntiSelected,
			electronConvVeto;

		// Jet Inforamtion
		int nJet,
			nDeepJetLooseBTag, nDeepJetMediumBTag, nDeepJetTightBTag,
			jetId;
		float rho, metPt, metPtJerUp, metPtJerDown, metPhi;
		std::vector<float> metPtJecUp, metPtJecDown;
		std::array<int, nMax> jetPartFlav,
			jetDeepJetId;
		std::array<float, nMax> jetPt, jetPtJerUp, jetPtJerDown,
			jetMass, jetMassJerUp, jetMassJerDown,
			jetEta, jetPhi,
			jetArea,
			jetDeepJetLooseSf, jetDeepJetMediumSf, jetDeepJetTightSf,
			jetDeepJetLooseLightSf, jetDeepJetMediumLightSf, jetDeepJetTightLightSf,
			jetDeepJet;
		std::vector<std::array<float, nMax>> jetPtJecUp, jetPtJecDown,
			jetMassJecUp, jetMassJecDown;
		std::array<bool, nMax> jetDeepJetLooseId, jetDeepJetMediumId, jetDeepJetTightId;
		std::vector<std::array<float, nMax>>
			jetDeepJetLooseSfUp, jetDeepJetLooseSfDown,
			jetDeepJetMediumSfUp, jetDeepJetMediumSfDown,
			jetDeepJetTightSfUp, jetDeepJetTightSfDown,
			jetDeepJetLooseLightSfUp, jetDeepJetLooseLightSfDown,
			jetDeepJetMediumLightSfUp, jetDeepJetMediumLightSfDown,
			jetDeepJetTightLightSfUp, jetDeepJetTightLightSfDown;


		// FatJet Inforamtion
		int nFatJet;
		std::array<int, nMax> fatJetId;
		std::array<float, nMax> fatJetPt, fatJetPtJerUp, fatJetPtJerDown,
			fatJetMass, fatJetMassJerUp, fatJetMassJerDown,
			fatJetEta, fatJetPhi,
			fatJetArea,
			fatJetDeepTagMDTvsQCD, fatJetDeepTagMDWvsQCD,
			fatJetDeepTagTvsQCD, fatJetDeepTagWvsQCD;
		std::vector<std::array<float, nMax>> fatJetPtJecUp, fatJetPtJecDown, fatJetMassJecUp, fatJetMassJecDown;


		// High Level Inforamtion
		float HT, LT, LP, deltaPhi, absoluteDeltaPhi, wBosonMt;
		bool isSignalRegion;

		// IsolatedTrack Information
		int nIsoTrack;
		bool isoTrackVeto;
		std::array<float, nMax> isoTrackMt2, isoTrackPt;
		std::array<bool, nMax> isoTrackIsHadronicDecay;
		std::array<int, nMax> isoTrackPdgId;

		// PileUp Weights
		int nPdfWeight, nScaleWeight;
		float nTrueInt,
			preFire, preFireUp, preFireDown,
			pileUpWeight, pileUpWeightUp, pileUpWeightDown;
		std::array<float, 103> pdfWeight; // size is fixed
		std::array<float, 9> scaleWeight; // size is fixed

		// Trigger Informaion
		//std::vector<bool> triggerValues, metTriggerValues, metFilterValues;
		std::vector<short> triggerValues, metTriggerValues, metFilterValues;
		bool hltEleOr, hltMuOr, hltMetOr;

		// Gen Level Information
		int nGenPart, nGenLepton, nGenTau, nGenLeptonFromTau, nGenMatchedW, nGenNeutrino,
			leptonDecayChannelFlag;
		float genWeight;
		bool isDiLeptonEvent, isHadTauEvent, leptonsInAcceptance;
		std::array<float, nMax> genDeltaPhiLepWSum, genDeltaPhiLepWDirect, genWSumMass, genWDirectMass, genMTLepNu, genNeutrinoPt;
		std::array<int, nMax> grandMotherPdgId, genTauGrandMotherPdgId, genTauMotherPdgId, genLepGrandMotherPdgId, genLepMotherPdgId;
};

#endif;

