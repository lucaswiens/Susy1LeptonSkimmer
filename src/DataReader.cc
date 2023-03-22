#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/DataReader.h>
#include<iostream>

DataReader::DataReader(const std::string &fileName, const std::string &treeName, const bool &isData, const bool &isFastSim) : isData(isData), isFastSim(isFastSim) {
	std::cout << std::endl << "Starting to open:" << std::endl << fileName << std::endl;
	inputFile = std::shared_ptr<TFile>(TFile::Open(fileName.c_str(), "READ"));
	std::cout << "File Opened!" << std::endl;
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
	muonMiniIsoLeaf         = inputTree->GetLeaf("Muon_miniPFRelIso_all");
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
	electronPdgIdLeaf           = inputTree->GetLeaf("Electron_pdgId");
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

	if (!isFastSim) {
		electronEnergyScaleUpLeaf   = inputTree->GetLeaf("Electron_dEscaleUp");
		electronEnergyScaleDownLeaf = inputTree->GetLeaf("Electron_dEscaleDown");
		electronEnergySigmaUpLeaf   = inputTree->GetLeaf("Electron_dEsigmaUp");
		electronEnergySigmaDownLeaf = inputTree->GetLeaf("Electron_dEsigmaDown");
	}

	// Jets
	nJetLeaf        = inputTree->GetLeaf("nJet");
	jetMassLeaf     = inputTree->GetLeaf("Jet_mass");
	jetPtLeaf       = inputTree->GetLeaf("Jet_pt");
	jetEtaLeaf      = inputTree->GetLeaf("Jet_eta");
	jetPhiLeaf      = inputTree->GetLeaf("Jet_phi");
	jetAreaLeaf     = inputTree->GetLeaf("Jet_area");
	jetDeepJetLeaf  = inputTree->GetLeaf("Jet_btagDeepFlavB");
	if (!isData) {
		jetPartFlavLeaf = inputTree->GetLeaf("Jet_partonFlavour");
	}
	jetRawFactorLeaf   = inputTree->GetLeaf("Jet_rawFactor");
	jetIdLeaf       = inputTree->GetLeaf("Jet_jetId");
	//jetPUIDLeaf     = inputTree->GetLeaf("Jet_puId");
	rhoLeaf         = inputTree->GetLeaf("fixedGridRhoFastjetAll");

	// Fat Jet
	nFatJetLeaf               = inputTree->GetLeaf("nFatJet");
	fatJetMassLeaf            = inputTree->GetLeaf("FatJet_mass");
	fatJetPtLeaf              = inputTree->GetLeaf("FatJet_pt");
	fatJetEtaLeaf             = inputTree->GetLeaf("FatJet_eta");
	fatJetPhiLeaf             = inputTree->GetLeaf("FatJet_phi");
	fatJetAreaLeaf            = inputTree->GetLeaf("FatJet_area");
	fatJetRawFactorLeaf       = inputTree->GetLeaf("FatJet_rawFactor");
	fatJetIdLeaf              = inputTree->GetLeaf("FatJet_jetId");
	fatJetDeepTagMDTvsQCDLeaf = inputTree->GetLeaf("FatJet_deepTagMD_TvsQCD");
	fatJetDeepTagMDWvsQCDLeaf = inputTree->GetLeaf("FatJet_deepTagMD_WvsQCD");
	fatJetDeepTagTvsQCDLeaf   = inputTree->GetLeaf("FatJet_deepTag_TvsQCD");
	fatJetDeepTagWvsQCDLeaf   = inputTree->GetLeaf("FatJet_deepTag_WvsQCD");

	// Isolated Tracks
	nIsoTrackLeaf     = inputTree->GetLeaf("nIsoTrack");
	isoTrackPdgIdLeaf = inputTree->GetLeaf("IsoTrack_pdgId");
	isoTrackPtLeaf    = inputTree->GetLeaf("IsoTrack_pt");
	isoTrackEtaLeaf   = inputTree->GetLeaf("IsoTrack_eta");
	isoTrackPhiLeaf   = inputTree->GetLeaf("IsoTrack_phi");

	// Gen Jet for Smearing
	nGenJetLeaf = inputTree->GetLeaf("nGenJet");
	genJetPtLeaf = inputTree->GetLeaf("GenJet_pt");
	genJetEtaLeaf = inputTree->GetLeaf("GenJet_eta");
	genJetPhiLeaf = inputTree->GetLeaf("GenJet_phi");

	nGenFatJetLeaf = inputTree->GetLeaf("nGenJetAK8");
	genFatJetPtLeaf = inputTree->GetLeaf("GenJetAK8_pt");
	genFatJetEtaLeaf = inputTree->GetLeaf("GenJetAK8_eta");
	genFatJetPhiLeaf = inputTree->GetLeaf("GenJetAK8_phi");

	// MET
	metPtLeaf  = inputTree->GetLeaf("MET_pt");
	metPhiLeaf = inputTree->GetLeaf("MET_phi");
	rawMetPtLeaf  = inputTree->GetLeaf("RawMET_pt");
	rawMetPhiLeaf = inputTree->GetLeaf("RawMET_phi");
	caloMetPtLeaf = inputTree->GetLeaf("CaloMET_pt");

	// PileUp
	nTrueIntLeaf     = inputTree->GetLeaf("Pileup_nTrueInt");
	nPdfWeightLeaf   = inputTree->GetLeaf("nLHEPdfWeight");
	pdfWeightLeaf    = inputTree->GetLeaf("LHEPdfWeight");
	nScaleWeightLeaf = inputTree->GetLeaf("nLHEScaleWeight");
	scaleWeightLeaf  = inputTree->GetLeaf("LHEScaleWeight");
	if (!isFastSim) {
		preFireLeaf      = inputTree->GetLeaf("L1PreFiringWeight_Nom");
		preFireUpLeaf    = inputTree->GetLeaf("L1PreFiringWeight_Up");
		preFireDownLeaf  = inputTree->GetLeaf("L1PreFiringWeight_Dn");
	}

	//Gen part related
	nGenPartLeaf       = inputTree->GetLeaf("nGenPart");
	genPdgIdLeaf       = inputTree->GetLeaf("GenPart_pdgId");
	genMotherIndexLeaf = inputTree->GetLeaf("GenPart_genPartIdxMother");
	genStatusLeaf      = inputTree->GetLeaf("GenPart_status");
	genStatusFlagsLeaf = inputTree->GetLeaf("GenPart_statusFlags");
	genPtLeaf          = inputTree->GetLeaf("GenPart_pt");
	genPhiLeaf         = inputTree->GetLeaf("GenPart_phi");
	genEtaLeaf         = inputTree->GetLeaf("GenPart_eta");
	genMassLeaf        = inputTree->GetLeaf("GenPart_mass");
	genWeightLeaf      = inputTree->GetLeaf("genWeight");
	genMetPtLeaf      = inputTree->GetLeaf("GenMET_pt");
	genMetPhiLeaf      = inputTree->GetLeaf("GenMET_phi");

	if (isFastSim) {
		for (TBranch *branch : TRangeDynCast<TBranch>(inputTree->GetListOfBranches())) {
			const std::string &branchName = branch->GetName();
			if (branchName.find("GenModel") != std::string::npos) {
				genModel.push_back(inputTree->GetLeaf(branch->GetName()));
			}
		}
	}
}

void DataReader::ReadMuonEntry() {
	if (nMuonLeaf->GetBranch()->GetReadEntry() == entry) { return;}

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
	muonNTrackerLayersLeaf->GetBranch()->GetEntry(entry);

	nMuon = nMuonLeaf->GetValue();
}

void DataReader::GetMuonValues(const int &index) {
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

	muonIdMap = {
		{'T', muonTightId},
		{'M', muonMediumId},
		{'L', muonLooseId}
	};
}


void DataReader::ReadElectronEntry() {
	if (nElectronLeaf->GetBranch()->GetReadEntry() == entry) return;
	nElectronLeaf->GetBranch()->GetEntry(entry);
	electronPtLeaf->GetBranch()->GetEntry(entry);
	electronEtaLeaf->GetBranch()->GetEntry(entry);
	electronPhiLeaf->GetBranch()->GetEntry(entry);
	electronMassLeaf->GetBranch()->GetEntry(entry);
	electronDxyLeaf->GetBranch()->GetEntry(entry);
	electronDzLeaf->GetBranch()->GetEntry(entry);
	electronChargeLeaf->GetBranch()->GetEntry(entry);
	electronPdgIdLeaf->GetBranch()->GetEntry(entry);
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
	if (!isFastSim) {
		electronEnergyScaleUpLeaf->GetBranch()->GetEntry(entry);
		electronEnergyScaleDownLeaf->GetBranch()->GetEntry(entry);
		electronEnergySigmaUpLeaf->GetBranch()->GetEntry(entry);
		electronEnergySigmaDownLeaf->GetBranch()->GetEntry(entry);
	}

	nElectron = nElectronLeaf->GetValue();
}

void DataReader::GetElectronValues(const int &index) {
	nElectron               = nElectronLeaf->GetValue(index);
	electronPt              = electronPtLeaf->GetValue(index);
	electronEta             = electronEtaLeaf->GetValue(index);
	electronPhi             = electronPhiLeaf->GetValue(index);
	electronMass            = electronMassLeaf->GetValue(index);
	electronDxy             = electronDxyLeaf->GetValue(index);
	electronDz              = electronDzLeaf->GetValue(index);
	electronCharge          = electronChargeLeaf->GetValue(index);
	electronPdgId          = electronPdgIdLeaf->GetValue(index);
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
	if (!isFastSim) {
		electronEnergyScaleUp   = electronEnergyScaleUpLeaf->GetValue(index);
		electronEnergyScaleDown = electronEnergyScaleDownLeaf->GetValue(index);
		electronEnergySigmaUp   = electronEnergySigmaUpLeaf->GetValue(index);
		electronEnergySigmaDown = electronEnergySigmaDownLeaf->GetValue(index);
	}

	electronIdMap = {
		{'T', electronCutBasedId >= 4},
		{'M', electronCutBasedId >= 3},
		{'L', electronCutBasedId >= 2},
		{'V', electronCutBasedId >= 1}
	};
}

void DataReader::ReadJetEntry() {
	if (nJetLeaf->GetBranch()->GetReadEntry() == entry) { return;}
	nJetLeaf->GetBranch()->GetEntry(entry);
	nJet = nJetLeaf->GetValue();

	// Rho
	rhoLeaf->GetBranch()->GetEntry(entry);
	rho = rhoLeaf->GetValue();
	// MET
	metPtLeaf->GetBranch()->GetEntry(entry);
	metPt  = metPtLeaf->GetValue();

	metPhiLeaf->GetBranch()->GetEntry(entry);
	metPhi = metPhiLeaf->GetValue();

	rawMetPtLeaf->GetBranch()->GetEntry(entry);
	rawMetPt  = rawMetPtLeaf->GetValue();

	rawMetPhiLeaf->GetBranch()->GetEntry(entry);
	rawMetPhi = rawMetPhiLeaf->GetValue();

	caloMetPtLeaf->GetBranch()->GetEntry(entry);
	caloMetPt = caloMetPtLeaf->GetValue();

	jetMassLeaf->GetBranch()->GetEntry(entry);
	jetPtLeaf->GetBranch()->GetEntry(entry);
	jetEtaLeaf->GetBranch()->GetEntry(entry);
	jetPhiLeaf->GetBranch()->GetEntry(entry);
	jetAreaLeaf->GetBranch()->GetEntry(entry);
	jetDeepJetLeaf->GetBranch()->GetEntry(entry);
	if (!isData) {
		jetPartFlavLeaf->GetBranch()->GetEntry(entry);
	}
	jetRawFactorLeaf->GetBranch()->GetEntry(entry);
	jetIdLeaf->GetBranch()->GetEntry(entry);
	//jetPUIDLeaf->GetBranch()->GetEntry(entry);
}

void DataReader::GetJetValues(const int &index) {
	jetMass = jetMassLeaf->GetValue(index);
	jetPt = jetPtLeaf->GetValue(index);
	jetEta = jetEtaLeaf->GetValue(index);
	jetPhi = jetPhiLeaf->GetValue(index);
	jetArea = jetAreaLeaf->GetValue(index);
	jetDeepJet = jetDeepJetLeaf->GetValue(index);
	if (!isData) {
		jetPartFlav = jetPartFlavLeaf->GetValue(index);
	}
	jetRawFactor = jetRawFactorLeaf->GetValue(index);
	jetId = jetIdLeaf->GetValue(index);
	//jetPUID = jetPUIDLeaf->GetValue(index);
}

void DataReader::ReadFatJetEntry() {
	if (nFatJetLeaf->GetBranch()->GetReadEntry() == entry) { return;}
	nFatJetLeaf->GetBranch()->GetEntry(entry);
	nFatJet = nFatJetLeaf->GetValue();

	fatJetMassLeaf->GetBranch()->GetEntry(entry);
	fatJetPtLeaf->GetBranch()->GetEntry(entry);
	fatJetEtaLeaf->GetBranch()->GetEntry(entry);
	fatJetPhiLeaf->GetBranch()->GetEntry(entry);
	fatJetAreaLeaf->GetBranch()->GetEntry(entry);
	fatJetRawFactorLeaf->GetBranch()->GetEntry(entry);
	fatJetIdLeaf->GetBranch()->GetEntry(entry);
	fatJetDeepTagMDTvsQCDLeaf->GetBranch()->GetEntry(entry);
	fatJetDeepTagMDWvsQCDLeaf->GetBranch()->GetEntry(entry);
	fatJetDeepTagTvsQCDLeaf->GetBranch()->GetEntry(entry);
	fatJetDeepTagWvsQCDLeaf->GetBranch()->GetEntry(entry);
}

void DataReader::GetFatJetValues(const int &index) {
	fatJetMass = fatJetMassLeaf->GetValue(index);
	fatJetPt = fatJetPtLeaf->GetValue(index);
	fatJetEta = fatJetEtaLeaf->GetValue(index);
	fatJetPhi = fatJetPhiLeaf->GetValue(index);
	fatJetArea = fatJetAreaLeaf->GetValue(index);
	fatJetRawFactor = fatJetRawFactorLeaf->GetValue(index);
	fatJetId = fatJetIdLeaf->GetValue(index);
	fatJetDeepTagMDTvsQCD = fatJetDeepTagMDTvsQCDLeaf->GetValue(index);
	fatJetDeepTagMDWvsQCD = fatJetDeepTagMDWvsQCDLeaf->GetValue(index);
	fatJetDeepTagTvsQCD   = fatJetDeepTagTvsQCDLeaf->GetValue(index);
	fatJetDeepTagWvsQCD   = fatJetDeepTagWvsQCDLeaf->GetValue(index);
}

void DataReader::ReadIsoTrackEntry() {
	if (nIsoTrackLeaf->GetBranch()->GetReadEntry() == entry) { return;}
	nIsoTrackLeaf->GetBranch()->GetEntry(entry);
	nIsoTrack = nIsoTrackLeaf->GetValue();

	isoTrackPdgIdLeaf->GetBranch()->GetEntry(entry);
	isoTrackPtLeaf->GetBranch()->GetEntry(entry);
	isoTrackEtaLeaf->GetBranch()->GetEntry(entry);
	isoTrackPhiLeaf->GetBranch()->GetEntry(entry);
}

void DataReader::GetIsoTrackValues(const int &index) {
	isoTrackPdgId = isoTrackPdgIdLeaf->GetValue(index);
	isoTrackPt = isoTrackPtLeaf->GetValue(index);
	isoTrackEta = isoTrackEtaLeaf->GetValue(index);
	isoTrackPhi = isoTrackPhiLeaf->GetValue(index);
}


void DataReader::ReadPileUpEntry() {
	if(nPdfWeightLeaf->GetBranch()->GetReadEntry() == entry) { return;}
	nTrueIntLeaf->GetBranch()->GetEntry(entry);
	nTrueInt = nTrueIntLeaf->GetValue();

	nPdfWeightLeaf->GetBranch()->GetEntry(entry);
	nPdfWeight = nPdfWeightLeaf->GetValue();

	nScaleWeightLeaf->GetBranch()->GetEntry(entry);
	nScaleWeight = nScaleWeightLeaf->GetValue();

	pdfWeightLeaf->GetBranch()->GetEntry(entry);
	scaleWeightLeaf->GetBranch()->GetEntry(entry);
	if (!isFastSim) {
		preFireLeaf->GetBranch()->GetEntry(entry);
		preFireUpLeaf->GetBranch()->GetEntry(entry);
		preFireDownLeaf->GetBranch()->GetEntry(entry);
	}
}

void DataReader::GetPileUpValues() {
	for (int iPdf = 0; iPdf < nPdfWeight; iPdf++) {
		pdfWeight[iPdf] = pdfWeightLeaf->GetValue(iPdf);
	}

	for (int iScale = 0; iScale < nScaleWeight; iScale++) {
		scaleWeight [iScale] = scaleWeightLeaf->GetValue(iScale);
	}

	if (!isFastSim) {
		preFire     = preFireLeaf->GetValue();
		preFireUp   = preFireUpLeaf->GetValue();
		preFireDown = preFireDownLeaf->GetValue();
	}
}

void DataReader::RegisterTrigger() {
	std::vector<int> triggerRemovalIndices;
	for (int iTrigger = 0; iTrigger < triggerNames.size(); iTrigger++) {
		const std::string &name = triggerNames.at(iTrigger);
		TLeaf *trigger = inputTree->GetLeaf(name.c_str());
		if (trigger == nullptr) {
			std::cout << "Not Found! Skipping Trigger: " << name << std::endl;
			triggerRemovalIndices.push_back(iTrigger);
			continue;
		}
		std::cout << "Register Trigger: " << name << std::endl;
		triggerLeafs.push_back(trigger);
		triggerValues.push_back(true);
	}

	for (const std::string &name : metTriggerNames) {
		TLeaf *trigger = inputTree->GetLeaf(name.c_str());
		if (trigger == nullptr) {
			std::cout << "Not Found! Skipping Trigger: " << name << std::endl;
			continue;
		}
		std::cout << "Register Trigger: " << name << std::endl;
		metTriggerLeafs.push_back(trigger);
		metTriggerValues.push_back(true);
	}

	// Remove trigger names if they are not available to get correct evaluation of the HLT_EleOr and HLT_MuonOr
	for (int iRemove = triggerRemovalIndices.size() - 1; iRemove >= 0; iRemove--) {
		const int &removeIndex = triggerRemovalIndices.at(iRemove);
		std::cout << "Removing " << triggerNames.at(removeIndex) << " at index = " << removeIndex << std::endl;
		triggerNames.erase(triggerNames.begin() + removeIndex);
	}
}

void DataReader::ReadTrigger() {
	for(TLeaf *triggerLeaf : triggerLeafs) {
		triggerLeaf->GetBranch()->GetEntry(entry);
	}
	for(TLeaf *triggerLeaf : metTriggerLeafs) {
		triggerLeaf->GetBranch()->GetEntry(entry);
	}
}

void DataReader::GetTrigger() {
	for(int i = 0; i < triggerLeafs.size(); ++i) {
		triggerValues[i] = triggerLeafs[i]->GetValue();
	}
	for(int i = 0; i < metTriggerLeafs.size(); ++i) {
		metTriggerValues[i] = metTriggerLeafs[i]->GetValue();
	}
}

void DataReader::RegisterMetFilter(const std::vector<std::string>& metFilterNames) {
	for(const std::string& filterName : metFilterNames) {
		if (filterName == "Flag_eeBadScFilter" && isData) { continue;}
		TLeaf *filter = inputTree->GetLeaf(filterName.c_str());
		if(filter == nullptr) {
			std::cout << "MET filter not found: '" + filterName + "', continue without it!" << std::endl;
			continue;
		}
		metFilterLeafs.push_back(filter);
		metFilterValues.push_back(true);
	}
}

void DataReader::ReadMetFilter() {
	for(TLeaf *filterLeaf : metFilterLeafs) {
		filterLeaf->GetBranch()->GetEntry(entry);
	}
}

void DataReader::GetMetFilter() {
	for(int i = 0; i < metFilterLeafs.size(); ++i) {
		metFilterValues[i] = metFilterLeafs[i]->GetValue();
	}
}

void DataReader::ReadGenJetEntry() {
	if(nGenJetLeaf->GetBranch()->GetReadEntry() == entry) return;
	nGenJetLeaf->GetBranch()->GetEntry(entry);
	nGenJet = nGenJetLeaf->GetValue();

	genJetPtLeaf->GetBranch()->GetEntry(entry);
	genJetEtaLeaf->GetBranch()->GetEntry(entry);
	genJetPhiLeaf->GetBranch()->GetEntry(entry);

	genMetPtLeaf->GetBranch()->GetEntry(entry);
	genMetPhiLeaf->GetBranch()->GetEntry(entry);
	genMetPt = genMetPtLeaf->GetValue();
	genMetPhi = genMetPhiLeaf->GetValue();
}

void DataReader::GetGenJetValues(const int &index) {
	genJetPt = genJetPtLeaf->GetValue(index);
	genJetEta = genJetEtaLeaf->GetValue(index);
	genJetPhi = genJetPhiLeaf->GetValue(index);
}


void DataReader::ReadGenFatJetEntry() {
	if(nGenFatJetLeaf->GetBranch()->GetReadEntry() == entry) return;
	nGenFatJetLeaf->GetBranch()->GetEntry(entry);
	nGenFatJet = nGenFatJetLeaf->GetValue();

	genFatJetPtLeaf->GetBranch()->GetEntry(entry);
	genFatJetEtaLeaf->GetBranch()->GetEntry(entry);
	genFatJetPhiLeaf->GetBranch()->GetEntry(entry);
}

void DataReader::GetGenFatJetValues(const int &index) {
	genFatJetPt = genFatJetPtLeaf->GetValue(index);
	genFatJetEta = genFatJetEtaLeaf->GetValue(index);
	genFatJetPhi = genFatJetPhiLeaf->GetValue(index);
}

void DataReader::ReadGenEntry() {
	if (nGenPartLeaf->GetBranch()->GetReadEntry() == entry) { return;}
	nGenPartLeaf->GetBranch()->GetEntry(entry);
	genPdgIdLeaf->GetBranch()->GetEntry(entry);
	genMotherIndexLeaf->GetBranch()->GetEntry(entry);
	genStatusLeaf->GetBranch()->GetEntry(entry);
	genStatusFlagsLeaf->GetBranch()->GetEntry(entry);
	genPtLeaf->GetBranch()->GetEntry(entry);
	genPhiLeaf->GetBranch()->GetEntry(entry);
	genEtaLeaf->GetBranch()->GetEntry(entry);
	genMassLeaf->GetBranch()->GetEntry(entry);
	genWeightLeaf->GetBranch()->GetEntry(entry);

	nGenPart = nGenPartLeaf->GetValue();
	genWeight = genWeightLeaf->GetValue();
}

void DataReader::GetGenValues(const int &index) {
	genPdgId       = genPdgIdLeaf->GetValue(index);
	genMotherIndex = genMotherIndexLeaf->GetValue(index);
	genStatus      = genStatusLeaf->GetValue(index);
	genStatusFlags = genStatusFlagsLeaf->GetValue(index);
	genPt          = genPtLeaf->GetValue(index);
	genPhi         = genPhiLeaf->GetValue(index);
	genEta         = genEtaLeaf->GetValue(index);
	genMass        = genMassLeaf->GetValue(index);
}

std::pair<int, int> DataReader::GetGenModel() {
	//TODO
}

int DataReader::LastGenCopy(const int& index) {
	GetGenValues(index);

	int partIndex = index, motherIndex = genMotherIndex;
	int partPDG = genPdgId;

	while(true) {
		GetGenValues(motherIndex);

		if (partPDG == genPdgId) {
			partIndex = motherIndex;
			motherIndex = genMotherIndex;
		} else {
			break;
		}
	}

	return partIndex;
}

int DataReader::GetGenMatchedIndex(const float &recoPt, const float &recoPhi, const float &recoEta, const int& recoPDG, const float &deltaRCut, const float &deltaPtCut) {
	int genIndex = -999;
	float deltaR,
		deltaPt,
		deltaRMin = std::numeric_limits<float>::max(),
		deltaPtMin = std::numeric_limits<float>::max();

	ReadGenEntry();
	for(int iGen = 0; iGen < nGenPart; iGen++) {
		GetGenValues(iGen);

		deltaR = Utility::DeltaR(recoEta, recoPhi, genEta, genPhi);
		deltaPt = std::abs(recoPt - genPt) / recoPt;

		if (deltaR > deltaRCut || deltaPt > deltaPtCut) continue;

		if (deltaR < deltaRMin && deltaPt < deltaPtMin && recoPDG == std::abs(genPdgId)) {
			int index = LastGenCopy(iGen);
			if (std::find(alreadyMatchedIndex.begin(), alreadyMatchedIndex.end(), index) != alreadyMatchedIndex.end()) continue;

			genIndex = index;
			deltaRMin = deltaR;
			deltaPtMin = deltaPt;
		}
	}

	return genIndex;
}
