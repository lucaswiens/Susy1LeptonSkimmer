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
		TLeaf *nJetLeaf,
			//*jetFlavourLeaf,
			*jetIdLeaf,
			*jetMassLeaf, *jetPtLeaf, *jetEtaLeaf, *jetPhiLeaf,
			*jetAreaLeaf, *jetRawFactorLeaf,
			*jetDeepCsvLeaf, *jetDeepJetLeaf,
			*rhoLeaf,
			*metPhiLeaf, *metPtLeaf;
		// Fat Jet Leafs
		TLeaf *nFatJetLeaf,
			*fatJetIdLeaf,
			*fatJetMassLeaf, *fatJetPtLeaf, *fatJetEtaLeaf, *fatJetPhiLeaf,
			*fatJetAreaLeaf, *fatJetRawFactorLeaf,
			*fatJetDeepTagMDTvsQCDLeaf, *fatJetDeepTagMDWvsQCDLeaf,
			*fatJetDeepTagTvsQCDLeaf, *fatJetDeepTagWvsQCDLeaf;

		// IsoTrack Leafs
		TLeaf *nIsoTrackLeaf,
			*isoTrackPdgIdLeaf,
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

		// Generator Particle Leafs
		TLeaf *nGenPartLeaf,
			*genPdgIdLeaf, *genMotherIndexLeaf,
			*genPtLeaf, *genEtaLeaf, *genPhiLeaf, *genMassLeaf,
			*genStatusLeaf, *genStatusFlagsLeaf;

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
		int nJet,
			jetId;
			//jetFlavour;
		double jetMass, jetPt, jetEta, jetPhi,
			jetArea, jetRawFactor,
			jetDeepCsv, jetDeepJet,
			rho,
			metPhi, metPt;

		void ReadFatJetEntry();
		void GetFatJetValues(const int &index);
		int nFatJet,
			fatJetId;
		double fatJetMass, fatJetPt, fatJetEta, fatJetPhi,
			fatJetArea, fatJetRawFactor,
			fatJetDeepTagMDTvsQCD, fatJetDeepTagMDWvsQCD,
			fatJetDeepTagTvsQCD, fatJetDeepTagWvsQCD;

		// Isolated Tracks
		void ReadIsoTrackEntry();
		void GetIsoTrackValues(const int &index);
		int nIsoTrack,
			isoTrackPdgId;
		double isoTrackPt, isoTrackEta, isoTrackPhi;

		// Gen Jet for Smearing
		void ReadGenJetEntry();
		void GetGenJetValues(const int &index);
		int nGenJet;
		double genJetPt, genJetEta, genJetPhi;
		void ReadGenFatJetEntry();
		void GetGenFatJetValues(const int &index);
		int nGenFatJet;
		double genFatJetPt, genFatJetEta, genFatJetPhi;

		// PileUp
		void ReadPileUpEntry();
		void GetPileUpValues();
		int nPdfWeight, nScaleWeight;
		double nTrueInt,
			preFire, preFireUp, preFireDown;
		std::array<double, 103> pdfWeight; // size is fixed
		std::array<double, 9> scaleWeight; // size is fixed

		// Generator Particle
		void ReadGenEntry();
		void GetGenValues(const int &index);
		int nGenPart,
			genPdgId, genMotherIndex,
			genStatus, genStatusFlags;
		double genPt, genEta, genPhi, genMass;

		int GetGenMatchedIndex(const double &recoPt, const double &recoPhi, const double &recoEta, const int& recoPDG, const double &deltaRCut, const double &deltaPtCut);
		int LastGenCopy(const int& index);
		std::vector<int> alreadyMatchedIndex;

};

#endif

