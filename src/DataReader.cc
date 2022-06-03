#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/DataReader.h>
#include<iostream>

DataReader::DataReader(const std::string &fileName, const std::string &treeName) {
	inputFile = std::shared_ptr<TFile>(TFile::Open(fileName.c_str(), "READ"));
	inputTree.reset(static_cast<TTree*>(inputFile->Get(treeName.c_str())));

	// Muons
	nMuonLeaf               = inputTree->GetLeaf("nMuon");
	muonPtLeaf              = inputTree->GetLeaf("Muon_pt");
	muonEtaLeaf             = inputTree->GetLeaf("Muon_eta");
	muonPhiLeaf             = inputTree->GetLeaf("Muon_phi");
	muonMassLeaf            = inputTree->GetLeaf("Muon_mass");
	muonDxyLeaf             = inputTree->GetLeaf("Muon_dxy");
	muonDzLeaf              = inputTree->GetLeaf("Muon_dz");
	muonSip3dLeaf           = inputTree->GetLeaf("Muon_sip3d");
	muonMiniIsoLeaf = inputTree->GetLeaf("Muon_miniPFRelIso_all");
	muonPdgIdLeaf           = inputTree->GetLeaf("Muon_pdgId");
	muonChargeLeaf          = inputTree->GetLeaf("Muon_charge");
	muonLooseIdLeaf         = inputTree->GetLeaf("Muon_looseId");
	muonMediumIdLeaf        = inputTree->GetLeaf("Muon_mediumId");
	muonTightIdLeaf         = inputTree->GetLeaf("Muon_tightId");
	muonMvaIdLeaf           = inputTree->GetLeaf("Muon_mvaId");
	muonIsPfCandLeaf        = inputTree->GetLeaf("Muon_isPFcand");
	muonNTrackerLayersLeaf  = inputTree->GetLeaf("Muon_nTrackerLayers");

	// Electrons
	nElectronLeaf               = inputTree->GetLeaf("nElectron");
	electronPtLeaf              = inputTree->GetLeaf("Electron_pt");
	electronEtaLeaf             = inputTree->GetLeaf("Electron_eta");
	electronPhiLeaf             = inputTree->GetLeaf("Electron_phi");
	electronMassLeaf            = inputTree->GetLeaf("Electron_mass");
	electronDxyLeaf             = inputTree->GetLeaf("Electron_dxy");
	electronDzLeaf              = inputTree->GetLeaf("Electron_dz");
	electronChargeLeaf          = inputTree->GetLeaf("Electron_charge");
	electronECorrLeaf           = inputTree->GetLeaf("Electron_eCorr");
	electronMiniIsoLeaf         = inputTree->GetLeaf("Electron_miniPFRelIso_all");
	//electronIso03Leaf           = inputTree->GetLeaf("Electron_pfRelIso03_all");
	//electronIso04Leaf           = inputTree->GetLeaf("Electron_pfRelIso04_all");
	electronRelJetIsoLeaf       = inputTree->GetLeaf("Electron_jetRelIso");
	electronCutBasedIdLeaf      = inputTree->GetLeaf("Electron_cutBased");
	electronLooseMvaIdLeaf      = inputTree->GetLeaf("Electron_mvaFall17V2Iso_WPL");
	electronMediumMvaIdLeaf     = inputTree->GetLeaf("Electron_mvaFall17V2Iso_WP90"); // 90% signal efficiency; https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Training_Details_and_Working_Poi
	electronTightMvaIdLeaf      = inputTree->GetLeaf("Electron_mvaFall17V2Iso_WP80"); // 80% signal efficiency
	electronConvVetoLeaf        = inputTree->GetLeaf("Electron_convVeto");
	electronNLostHitsLeaf       = inputTree->GetLeaf("Electron_lostHits");
	electronEnergyScaleUpLeaf   = inputTree->GetLeaf("Electron_dEscaleUp");
	electronEnergyScaleDownLeaf = inputTree->GetLeaf("Electron_dEscaleDown");
	electronEnergySigmaUpLeaf   = inputTree->GetLeaf("Electron_dEsigmaUp");
	electronEnergySigmaDownLeaf = inputTree->GetLeaf("Electron_dEsigmaDown");

	//Gen part related
	nGenPartLeaf       = inputTree->GetLeaf("nGenPart");
	genPDGLeaf         = inputTree->GetLeaf("GenPart_pdgId");
	genMotherIndexLeaf = inputTree->GetLeaf("GenPart_genPartIdxMother");
	genPtLeaf          = inputTree->GetLeaf("GenPart_pt");
	genPhiLeaf         = inputTree->GetLeaf("GenPart_phi");
	genEtaLeaf         = inputTree->GetLeaf("GenPart_eta");
	genMassLeaf        = inputTree->GetLeaf("GenPart_mass");
}

void DataReader::ReadMuonEntry() {
	if(nMuonLeaf->GetBranch()->GetReadEntry() == entry) { return;}

	muonRandomNumber = gRandom->Rndm();

	nMuonLeaf->GetBranch()->GetEntry(entry);
	muonPtLeaf->GetBranch()->GetEntry(entry);
	muonEtaLeaf->GetBranch()->GetEntry(entry);
	muonPhiLeaf->GetBranch()->GetEntry(entry);
	muonMassLeaf->GetBranch()->GetEntry(entry);
	muonDxyLeaf->GetBranch()->GetEntry(entry);
	muonDzLeaf->GetBranch()->GetEntry(entry);
	muonSip3dLeaf->GetBranch()->GetEntry(entry);
	muonMiniIsoLeaf->GetBranch()->GetEntry(entry);
	muonPdgIdLeaf->GetBranch()->GetEntry(entry);
	muonChargeLeaf->GetBranch()->GetEntry(entry);
	muonLooseIdLeaf->GetBranch()->GetEntry(entry);
	muonMediumIdLeaf->GetBranch()->GetEntry(entry);
	muonTightIdLeaf->GetBranch()->GetEntry(entry);
	muonMvaIdLeaf->GetBranch()->GetEntry(entry);
	muonIsPfCandLeaf->GetBranch()->GetEntry(entry);

	nMuon = nMuonLeaf->GetValue();
}

void DataReader::GetMuonValues(const std::size_t &index) {
	muonPt             = muonPtLeaf->GetValue(index);
	muonEta            = muonEtaLeaf->GetValue(index);
	muonPhi            = muonPhiLeaf->GetValue(index);
	muonMass           = muonMassLeaf->GetValue(index);
	muonDxy            = muonDxyLeaf->GetValue(index);
	muonDz             = muonDzLeaf->GetValue(index);
	muonSip3d          = muonSip3dLeaf->GetValue(index);
	muonMiniIso        = muonMiniIsoLeaf->GetValue(index);
	muonPdgId          = muonPdgIdLeaf->GetValue(index);
	muonCharge         = muonChargeLeaf->GetValue(index);
	muonLooseId        = muonLooseIdLeaf->GetValue(index);
	muonMediumId       = muonMediumIdLeaf->GetValue(index);
	muonTightId        = muonTightIdLeaf->GetValue(index);
	muonMvaId          = muonMvaIdLeaf->GetValue(index);
	muonIsPfCand       = muonIsPfCandLeaf->GetValue(index);
	muonNTrackerLayers = muonNTrackerLayersLeaf->GetValue(index);
	//std::cout << "T = " << muonTightId << "; M = " << muonMediumId << "; L" << muonLooseId << std::endl;
}


void DataReader::ReadElectronEntry() {
	if(nElectronLeaf->GetBranch()->GetReadEntry() == entry) return;
	nElectronLeaf->GetBranch()->GetEntry(entry);
	electronPtLeaf->GetBranch()->GetEntry(entry);
	electronEtaLeaf->GetBranch()->GetEntry(entry);
	electronPhiLeaf->GetBranch()->GetEntry(entry);
	electronMassLeaf->GetBranch()->GetEntry(entry);
	electronDxyLeaf->GetBranch()->GetEntry(entry);
	electronDzLeaf->GetBranch()->GetEntry(entry);
	electronChargeLeaf->GetBranch()->GetEntry(entry);
	electronECorrLeaf->GetBranch()->GetEntry(entry);
	electronMiniIsoLeaf->GetBranch()->GetEntry(entry);
	//electronIso03Leaf->GetBranch()->GetEntry(entry);
	//electronIso04Leaf->GetBranch()->GetEntry(entry);
	electronRelJetIsoLeaf->GetBranch()->GetEntry(entry);
	electronCutBasedIdLeaf->GetBranch()->GetEntry(entry);
	electronLooseMvaIdLeaf->GetBranch()->GetEntry(entry);
	electronMediumMvaIdLeaf->GetBranch()->GetEntry(entry);
	electronTightMvaIdLeaf->GetBranch()->GetEntry(entry);
	electronConvVetoLeaf->GetBranch()->GetEntry(entry);
	electronNLostHitsLeaf->GetBranch()->GetEntry(entry);
	electronEnergyScaleUpLeaf->GetBranch()->GetEntry(entry);
	electronEnergyScaleDownLeaf->GetBranch()->GetEntry(entry);
	electronEnergySigmaUpLeaf->GetBranch()->GetEntry(entry);
	electronEnergySigmaDownLeaf->GetBranch()->GetEntry(entry);

	nElectron = nElectronLeaf->GetValue();
}

void DataReader::GetElectronValues(const std::size_t &index) {
	nElectron               = nElectronLeaf->GetValue(index);
	electronPt              = electronPtLeaf->GetValue(index);
	electronEta             = electronEtaLeaf->GetValue(index);
	electronPhi             = electronPhiLeaf->GetValue(index);
	electronMass            = electronMassLeaf->GetValue(index);
	electronDxy             = electronDxyLeaf->GetValue(index);
	electronDz              = electronDzLeaf->GetValue(index);
	electronCharge          = electronChargeLeaf->GetValue(index);
	electronECorr           = electronECorrLeaf->GetValue(index);
	electronMiniIso         = electronMiniIsoLeaf->GetValue(index);
	//electronIso03           = electronIso03Leaf->GetValue(index);
	//electronIso04           = electronIso04Leaf->GetValue(index);
	electronRelJetIso       = electronRelJetIsoLeaf->GetValue(index);
	electronCutBasedId      = electronCutBasedIdLeaf->GetValue(index);
	electronLooseMvaId      = electronLooseMvaIdLeaf->GetValue(index);
	electronMediumMvaId     = electronMediumMvaIdLeaf->GetValue(index);
	electronTightMvaId      = electronTightMvaIdLeaf->GetValue(index);
	electronConvVeto        = electronConvVetoLeaf->GetValue(index);
	electronNLostHits       = electronNLostHitsLeaf->GetValue(index);
	electronEnergyScaleUp   = electronEnergyScaleUpLeaf->GetValue(index);
	electronEnergyScaleDown = electronEnergyScaleDownLeaf->GetValue(index);
	electronEnergySigmaUp   = electronEnergySigmaUpLeaf->GetValue(index);
	electronEnergySigmaDown = electronEnergySigmaDownLeaf->GetValue(index);
}

void DataReader::ReadGenPartEntry() {
	if(nGenPartLeaf->GetBranch()->GetReadEntry() == entry) { return;}
	nGenPartLeaf->GetBranch()->GetEntry(entry);
	genPDGLeaf->GetBranch()->GetEntry(entry);
	genMotherIndexLeaf->GetBranch()->GetEntry(entry);
	genPtLeaf->GetBranch()->GetEntry(entry);
	genPhiLeaf->GetBranch()->GetEntry(entry);
	genEtaLeaf->GetBranch()->GetEntry(entry);
	genMassLeaf->GetBranch()->GetEntry(entry);

	nGenPart = nGenPartLeaf->GetValue();
}

void DataReader::GetGenPartValues(const std::size_t &index) {
	genPDG         = genPDGLeaf->GetValue(index);
	genMotherIndex = genMotherIndexLeaf->GetValue(index);
	genPt          = genPtLeaf->GetValue(index);
	genPhi         = genPhiLeaf->GetValue(index);
	genEta         = genEtaLeaf->GetValue(index);
	genMass        = genMassLeaf->GetValue(index);
}

int DataReader::LastGenCopy(const int& index){
	GetGenPartValues(index);

	int partIndex = index, motherIndex = genMotherIndex;
	int partPDG = genPDG;

	while(true){
		GetGenPartValues(motherIndex);

		if(partPDG == genPDG){
			partIndex = motherIndex;
			motherIndex = genMotherIndex;
		}

		else break;
	}

	return partIndex;
}

int DataReader::GetGenMatchedIndex(const double &recoPt, const double &recoPhi, const double &recoEta, const int& recoPDG, const double &deltaRCut, const double &deltaPtCut){
	int genIndex = -999;
	double deltaR,
		deltaPt,
		deltaRMin = std::numeric_limits<double>::max(),
		deltaPtMin = std::numeric_limits<double>::max();

	ReadGenPartEntry();
	for(int iGen = 0; iGen < nGenPart; iGen++){
		GetGenPartValues(iGen);

		deltaR = Utility::DeltaR(recoEta, recoPhi, genEta, genPhi);
		deltaPt = std::abs(recoPt - genPt) / recoPt;

		if(deltaR > deltaRCut || deltaPt > deltaPtCut) continue;

		if(deltaR < deltaRMin && deltaPt < deltaPtMin && recoPDG == std::abs(genPDG)){
			int index = LastGenCopy(iGen);
			if(std::find(alreadyMatchedIndex.begin(), alreadyMatchedIndex.end(), index) != alreadyMatchedIndex.end()) continue;

			genIndex = index;
			deltaRMin = deltaR;
			deltaPtMin = deltaPt;
		}
	}

	return genIndex;
}
