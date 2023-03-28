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

		// Event Number and meta information
		int entry;
		bool isData, isFastSim;

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
			*electronPdgIdLeaf, *electronChargeLeaf,
			*electronECorrLeaf,
			*electronMiniIsoLeaf, *electronIso03Leaf, *electronIso04Leaf, *electronRelJetIsoLeaf,
			*electronCutBasedIdLeaf, *electronLooseMvaIdLeaf, *electronMediumMvaIdLeaf, *electronTightMvaIdLeaf,
			*electronConvVetoLeaf,
			*electronNLostHitsLeaf,
			*electronEnergyScaleUpLeaf, *electronEnergyScaleDownLeaf,
			*electronEnergySigmaUpLeaf, *electronEnergySigmaDownLeaf;

		// Jet Leafs
		TLeaf *nJetLeaf,
			*jetPartFlavLeaf,
			*jetIdLeaf,
			*jetMassLeaf, *jetPtLeaf, *jetEtaLeaf, *jetPhiLeaf,
			*jetAreaLeaf, *jetRawFactorLeaf,
			*jetDeepJetLeaf,
			*rhoLeaf,
			*metPhiLeaf, *metPtLeaf, *rawMetPhiLeaf, *rawMetPtLeaf, *caloMetPtLeaf;

		// Fat Jet Leafs
		TLeaf *nFatJetLeaf,
			*fatJetIdLeaf,
			*fatJetMassLeaf, *fatJetPtLeaf, *fatJetEtaLeaf, *fatJetPhiLeaf,
			*fatJetAreaLeaf, *fatJetRawFactorLeaf,
			*fatJetDeepTagMDTvsQCDLeaf, *fatJetDeepTagMDWvsQCDLeaf,
			*fatJetDeepTagTvsQCDLeaf, *fatJetDeepTagWvsQCDLeaf;

		// IsoTrack Leafs
		TLeaf *nIsoTrackLeaf,
			*isoTrackPdgIdLeaf, *isoTrackChargeLeaf,
			*isoTrackPtLeaf, *isoTrackEtaLeaf, *isoTrackPhiLeaf;

		// GenJet for smearing
		TLeaf *nGenJetLeaf, *nGenFatJetLeaf,
			*genJetPtLeaf, *genJetEtaLeaf, *genJetPhiLeaf,
			*genFatJetPtLeaf, *genFatJetEtaLeaf, *genFatJetPhiLeaf;

		// PileUp Leafs
		TLeaf *nTrueIntLeaf,
			*nPdfWeightLeaf, *pdfWeightLeaf,
			*nScaleWeightLeaf, *scaleWeightLeaf,
			*preFireLeaf, *preFireUpLeaf, *preFireDownLeaf;

		// Trigger and METFilter Leafs
		std::vector<TLeaf*> triggerLeafs, metTriggerLeafs, metFilterLeafs;

		// Generator Particle Leafs
		TLeaf *nGenPartLeaf,
			*genPdgIdLeaf, *genMotherIndexLeaf,
			*genPtLeaf, *genEtaLeaf, *genPhiLeaf, *genMassLeaf,
			*genMetPtLeaf, *genMetPhiLeaf,
			*genWeightLeaf,
			*genStatusLeaf, *genStatusFlagsLeaf;

		// FastSim Leafs
		std::vector<TLeaf*> genModel;

	public:
		DataReader(const std::string &fileName, const std::string &treeName, const bool &isData, const bool &isFastSim);

		// Event entry information
		int GetEntries(){return inputTree->GetEntries();}
		void SetEntry(const int &entry){this->entry = entry;}

		// Muon
		void ReadMuonEntry();
		void GetMuonValues(const int &index);
		int nMuon;
		float muonPt, muonEta, muonPhi, muonMass, muonIso, muonDxy, muonDz, muonSip3d, muonMiniIso, muonRandomNumber;
		int muonPdgId, muonCharge, muonMvaId, muonNTrackerLayers;
		bool muonLooseId, muonMediumId, muonTightId, muonIsPfCand;
		std::map<char, bool> muonIdMap;

		// Electron
		void ReadElectronEntry();
		void GetElectronValues(const int &index);
		int nElectron;
		float electronPt, electronEta, electronPhi, electronMass,
			electronDxy, electronDz,
			electronECorr,
			electronMiniIso, electronIso03, electronIso04, electronRelJetIso,
			electronEnergyScaleUp, electronEnergyScaleDown,
			electronEnergySigmaUp, electronEnergySigmaDown;
		int electronPdgId, electronCharge, electronCutBasedId, electronConvVeto, electronNLostHits;
		bool electronLooseMvaId, electronMediumMvaId, electronTightMvaId;
		std::map<char, bool> electronIdMap;

		// Jet
		void ReadJetEntry();
		void GetJetValues(const int &index);
		int nJet,
			jetId,
			jetPartFlav;
		float jetMass, jetPt, jetEta, jetPhi,
			jetArea, jetRawFactor,
			jetDeepJet,
			rho,
			metPhi, metPt, rawMetPhi, rawMetPt, caloMetPt;

		void ReadFatJetEntry();
		void GetFatJetValues(const int &index);
		int nFatJet,
			fatJetId;
		float fatJetMass, fatJetPt, fatJetEta, fatJetPhi,
			fatJetArea, fatJetRawFactor,
			fatJetDeepTagMDTvsQCD, fatJetDeepTagMDWvsQCD,
			fatJetDeepTagTvsQCD, fatJetDeepTagWvsQCD;

		// Isolated Tracks
		void ReadIsoTrackEntry();
		void GetIsoTrackValues(const int &index);
		int nIsoTrack,
			isoTrackCharge,
			isoTrackPdgId;
		float isoTrackPt, isoTrackEta, isoTrackPhi;

		// Gen Jet for Smearing
		void ReadGenJetEntry();
		void GetGenJetValues(const int &index);
		int nGenJet;
		float genJetPt, genJetEta, genJetPhi;
		void ReadGenFatJetEntry();
		void GetGenFatJetValues(const int &index);
		int nGenFatJet;
		float genFatJetPt, genFatJetEta, genFatJetPhi;

		// PileUp
		void ReadPileUpEntry();
		void GetPileUpValues();
		int nPdfWeight, nScaleWeight;
		float nTrueInt,
			preFire, preFireUp, preFireDown;
		std::array<float, 103> pdfWeight; // size is fixed
		std::array<float, 9> scaleWeight; // size is fixed

		// Trigger and MET Filter
		std::vector<short> triggerValues, metTriggerValues, metFilterValues;
		std::vector<std::string> triggerNames, metTriggerNames, metFilterNames;
		void RegisterTrigger();
		void ReadTrigger();
		void GetTrigger();
		void RegisterMetFilter(const std::vector<std::string> &metFilterNames);
		void ReadMetFilter();
		void GetMetFilter();

		// Generator Particle
		void ReadGenEntry();
		void GetGenValues(const int &index);
		int nGenPart,
			genPdgId, genMotherIndex,
			genStatus, genStatusFlags;
		float genPt, genEta, genPhi, genMass,
			genMetPt, genMetPhi,
			genWeight;

		int GetGenMatchedIndex(const float &recoPt, const float &recoPhi, const float &recoEta, const int &recoPDG, const float &deltaRCut, const float &deltaPtCut);
		int LastGenCopy(const int &index);
		std::vector<int> alreadyMatchedIndex;

		std::pair<int, int> GetGenModel();
};

#endif

