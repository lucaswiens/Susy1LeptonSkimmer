#ifndef DATAREADER_H
#define DATAREADER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Utility.h>

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TRandom.h>

class DataReader {
	private:
		std::shared_ptr<TFile> inputFile;
		std::shared_ptr<TTree> inputTree;

		// Event Number
		int entry;

		// MetaData Leafs
		TLeaf *processNameLeaf;
		// Muon Leafs
		TLeaf *nMuonLeaf,
			*muonPtLeaf, *muonEtaLeaf, *muonPhiLeaf, *muonMassLeaf, *muonIsoLeaf, *muonDxyLeaf, *muonDzLeaf, *muonSip3dLeaf, *muonMiniIsoLeaf,
			*muonPdgIdLeaf, *muonChargeLeaf, *muonNTrackerLayersLeaf,
			*muonLooseIdLeaf, *muonMediumIdLeaf, *muonTightIdLeaf, *muonMvaIdLeaf, *muonIsPfCandLeaf;

		// Electron Leafs
		TLeaf *nElectronLeaf,
			*electronPtLeaf, *electronEtaLeaf, *electronPhiLeaf, *electronMassLeaf,
			*electronDxyLeaf, *electronDzLeaf,
			*electronChargeLeaf,
			*electronECorrLeaf,
			*electronMiniIsoLeaf, *electronIso03Leaf, *electronIso04Leaf, *electronRelJetIsoLeaf,
			*electronCutBasedIdLeaf, *electronLooseMvaIdLeaf, *electronMediumMvaIdLeaf, *electronTightMvaIdLeaf,
			*electronConvVetoLeaf,
			*electronNLostHitsLeaf,
			*electronEnergyScaleUpLeaf, *electronEnergyScaleDownLeaf,
			*electronEnergySigmaUpLeaf, *electronEnergySigmaDownLeaf;

		// Jet Leafs
		TLeaf *nJetLeaf, *nFatJetLeaf,
			//*jetFlavourLeaf,
			*jetIdLeaf,
			*jetMassLeaf, *jetPtLeaf, *jetEtaLeaf, *jetPhiLeaf,
			*jetAreaLeaf, *jetRawFactorLeaf,
			*jetDeepCsvLeaf, *jetDeepJetLeaf,
			*rhoLeaf,
			*metPhiLeaf, *metPtLeaf;
			//*fatJet_deepTagMD_H4qvsQCDLeaf, *fatJet_deepTagMD_HbbvsQCDLeaf, *fatJet_deepTagMD_TvsQCDLeaf, *fatJet_deepTagMD_WvsQCDLeaf,
			//*fatJet_deepTagMD_ZHbbvsQCDLeaf, *fatJet_deepTagMD_ZHccvsQCDLeaf, *fatJet_deepTagMD_ZbbvsQCDLeaf, *fatJet_deepTagMD_ZvsQCDLeaf,
			//*fatJet_deepTagMD_bbvsLightLeaf, *fatJet_deepTagMD_ccvsLightLeaf,
			//*fatJet_deepTag_HLeaf, *fatJet_deepTag_QCDLeaf, *fatJet_deepTag_QCDothersLeaf,
			//*fatJet_deepTag_TvsQCDLeaf, *fatJet_deepTag_WvsQCDLeaf, *fatJet_deepTag_ZvsQCDLeaf;

		// GenJet for smearing
		TLeaf *nGenJetLeaf, *nGenFatJetLeaf,
			*genJetPtLeaf, *genJetEtaLeaf, *genJetPhiLeaf,
			*genFatJetPtLeaf, *genFatJetEtaLeaf, *genFatJetPhiLeaf;



		// Generator Particle Leafs
		TLeaf *nGenPartLeaf,
			*genPDGLeaf, *genMotherIndexLeaf,
			*genPtLeaf, *genEtaLeaf, *genPhiLeaf, *genMassLeaf;

	public:
		DataReader(const std::string &fileName, const std::string &treeName);

		// Event entry information
		int GetEntries(){return inputTree->GetEntries();}
		void SetEntry(const int &entry){this->entry = entry;}

		// Muon
		void ReadMuonEntry();
		void GetMuonValues(const int &index);
		int nMuon;
		double muonPt, muonEta, muonPhi, muonMass, muonIso, muonDxy, muonDz, muonSip3d, muonMiniIso, muonRandomNumber;
		int muonPdgId, muonCharge, muonMvaId, muonNTrackerLayers;
		bool muonLooseId, muonMediumId, muonTightId, muonIsPfCand;
		std::map<char, bool> muonIdMap;

		// Electron
		void ReadElectronEntry();
		void GetElectronValues(const int &index);
		int nElectron;
		double electronPt, electronEta, electronPhi, electronMass,
			electronDxy, electronDz,
			electronECorr,
			electronMiniIso, electronIso03, electronIso04, electronRelJetIso,
			electronEnergyScaleUp, electronEnergyScaleDown,
			electronEnergySigmaUp, electronEnergySigmaDown;
		int electronCharge, electronCutBasedId, electronConvVeto, electronNLostHits;
		bool electronLooseMvaId, electronMediumMvaId, electronTightMvaId;
		std::map<char, bool> electronIdMap;

		// Jet
		void ReadJetEntry();
		void GetJetValues(const int &index);
		int nJet, nFatJet,
			jetId;
			//jetFlavour;
		double jetMass, jetPt, jetEta, jetPhi,
			jetArea, jetRawFactor,
			jetDeepCsv, jetDeepJet,
			rho,
			metPhi, metPt;

		// Gen Jet for Smearing
		void ReadGenJetEntry();
		void GetGenJetValues(const int &index);
		int nGenJet;
		double genJetPt, genJetEta, genJetPhi;
		void ReadGenFatJetEntry();
		void GetGenFatJetValues(const int &index);
		int nGenFatJet;
		double genFatJetPt, genFatJetEta, genFatJetPhi;


		// Generator Particle
		void ReadGenEntry();
		void GetGenValues(const int &index);
		int nGenPart;
		int genPDG, genMotherIndex;
		double genPt, genEta, genPhi, genMass;

		int GetGenMatchedIndex(const double &recoPt, const double &recoPhi, const double &recoEta, const int& recoPDG, const double &deltaRCut, const double &deltaPtCut);
		int LastGenCopy(const int& index);
		std::vector<int> alreadyMatchedIndex;

};

#endif

