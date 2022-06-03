#ifndef SUSY1LEPTONPRODUCT_H
#define SUSY1LEPTONPRODUCT_H

#include <TTree.h>
#include <TFile.h>

class Susy1LeptonProduct {
	private:
		int era;
		bool preVFP, isData;
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

		// Max value for static arrays, should use assert to enforce nObject < nMax
		static const std::size_t nMax = 20;

		// Lepton Information
		int nLepton;

		// Muons Inforamtion
		int nMuon, nGoodMuon, nVetoMuon, nAntiSelectedMuon;
		std::array<double, nMax> muonPt, muonRoccorPt, muonRoccorPtUp, muonRoccorPtDown, muonEta, muonPhi, muonMass, muonMiniIso, muonDxy, muonDz, muonSip3d;
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
			electronEnergySigmaUp, electronEnergySigmaDown;
		std::array<int, nMax> electronCharge, electronCutBasedId, electronNLostHits;
		std::array<bool, nMax> electronLooseMvaId, electronMediumMvaId, electronTightMvaId,
			electronTightId, electronMediumId, electronLooseId, electronVetoId,
			electronLoose, electronMedium, electronTight,
			electronIsGood, electronIsVeto, electronIsAntiSelected,
			electronConvVeto;

		// FIXME Old Info Delete Later
		float leptonPt = -999;
		float leptonEta = -999;
		float leptonPhi = -999;
		float leptonMass = -999;
		int leptonPdgId = -999;
		int leptonCharge = -999;

		int nJet = 0;
		int nLooseCSVBTagJet = 0;
		int nMediumCSVBTagJet = 0;
		int nTightCSVBTagJet = 0;
		int nLooseDFBTagJet = 0;
		int nMediumDFBTagJet = 0;
		int nTightDFBTagJet = 0;
		std::vector<float> jetPt;
		std::vector<float> jetPhi;
		std::vector<float> jetEta;
		std::vector<float> jetMass;
		std::vector<bool> jetMediumCSVBTag;
		std::vector<bool> jetMediumDFBTag;
		//Required for Signal Region Flag

		float metPt = -999;
		float metPhi = -999;

		int nGenPart = -999;
		int nGenWBoson = -999;
		std::vector<float> genWBosonPt;
		std::vector<float> genWBosonEta;
		std::vector<float> genWBosonPhi;
		std::vector<float> genWBosonMass;

		int nGenTop = -999;
		std::vector<float> genTopPt;
		std::vector<float> genTopEta;
		std::vector<float> genTopPhi;
		std::vector<float> genTopMass;
};

#endif

