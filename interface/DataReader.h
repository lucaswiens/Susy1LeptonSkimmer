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
		std::size_t entry;

		// MetaData Leafs
		TLeaf *processNameLeaf;
		// Muons Leafs
		TLeaf *nMuonLeaf,
			*muonPtLeaf, *muonEtaLeaf, *muonPhiLeaf, *muonMassLeaf, *muonIsoLeaf, *muonDxyLeaf, *muonDzLeaf, *muonSip3dLeaf, *muonMiniIsoLeaf,
			*muonPdgIdLeaf, *muonChargeLeaf, *muonNTrackerLayersLeaf,
			*muonLooseIdLeaf, *muonMediumIdLeaf, *muonTightIdLeaf, *muonMvaIdLeaf, *muonIsPfCandLeaf;

		// Electrons Leafs
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

		// Generator Particle Leafs
		TLeaf *nGenPartLeaf,
			*genPDGLeaf, *genMotherIndexLeaf,
			*genPtLeaf, *genEtaLeaf, *genPhiLeaf, *genMassLeaf;
	public:
		DataReader(const std::string &fileName, const std::string &treeName);

		// Event entry information
		std::size_t GetEntries(){return inputTree->GetEntries();}
		void SetEntry(const std::size_t &entry){this->entry = entry;}

		// Muon
		void ReadMuonEntry();
		void GetMuonValues(const std::size_t &index);
		int nMuon;
		double muonPt, muonEta, muonPhi, muonMass, muonIso, muonDxy, muonDz, muonSip3d, muonMiniIso, muonRandomNumber;
		int muonPdgId, muonCharge, muonMvaId, muonNTrackerLayers;
		bool muonLooseId, muonMediumId, muonTightId, muonIsPfCand;

		// Electron
		void ReadElectronEntry();
		void GetElectronValues(const std::size_t &index);
		int nElectron;
		double electronPt, electronEta, electronPhi, electronMass,
			electronDxy, electronDz,
			electronECorr,
			electronMiniIso, electronIso03, electronIso04, electronRelJetIso,
			electronEnergyScaleUp, electronEnergyScaleDown,
			electronEnergySigmaUp, electronEnergySigmaDown;
		int electronCharge, electronCutBasedId, electronConvVeto, electronNLostHits;
		bool electronLooseMvaId, electronMediumMvaId, electronTightMvaId;

		// Generator Particle
		void ReadGenPartEntry();
		void GetGenPartValues(const std::size_t &index);
		int nGenPart;
		int genPDG, genMotherIndex;
		double genPt, genEta, genPhi, genMass;

		int GetGenMatchedIndex(const double &recoPt, const double &recoPhi, const double &recoEta, const int& recoPDG, const double &deltaRCut, const double &deltaPtCut);
		int LastGenCopy(const int& index);
		std::vector<int> alreadyMatchedIndex;

};

#endif

