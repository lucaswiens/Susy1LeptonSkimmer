#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>
#include<iostream>

Susy1LeptonProduct::Susy1LeptonProduct(const int &era, const bool &isData, const bool &isFastSim, const std::string &sampleName, const char &runPeriod, const float &xSection, const pt::ptree &configTree, TFile &outputFile) :
	era(era),
	isData(isData),
	isFastSim(isFastSim),
	sampleName(sampleName),
	runPeriod(runPeriod),
	xSection(xSection) {
		std::string outputFileName = outputFile.GetName();
		//std::string &primaryDataset;
		if (era == 2016) {
			if (outputFileName.find("UL16NanoAODAPVv") != std::string::npos) {
				// See Eras: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDatasetsUL2016
				this->preVFP = true; // pre-VFP aka HIPM aka APV because let's have 3 names for the exact same thing
				eraSelector = std::to_string(era) + "PreVFP";
			} else {
				this->preVFP = false; // aka postVFP, no-HIPM, no-APV which is just default track construction
				eraSelector = std::to_string(era) + "PostVFP";
			}
		} else {
			eraSelector = std::to_string(era);
			this->preVFP = false; // aka postVFP, no-HIPM, no-APV which is just default track construction
		}

		// if clause has to be shorter than length of str but > 0
		if (sampleName.find("SingleEle") < 20) {
				this->primaryDataset = "isSingleElectron";
				//datasetDecider = "PD_SingleEle";
				//std::cout << primaryDataset << sampleName.find("SingleEle") << std::endl;
			}
		else if (sampleName.find("SingleMu") < 20) {
				this->primaryDataset = "isSingleMuon";
				//datasetDecider = "PD_SingleMu";
				//std::cout << primaryDataset << sampleName.find("SingleMu") << std::endl;
			}
		else if (sampleName.find("MET_") < 20) {
				this->primaryDataset = "isMet";
				//std::cout << primaryDataset << std::endl;
			}


		for (const std::pair<std::string, boost::property_tree::ptree> jecSyst : configTree.get_child("Producer.Jet.JECSystematic")) {
			jetPtJecUp.push_back(std::array<float, nMax>());
			jetPtJecDown.push_back(std::array<float, nMax>());
			jetMassJecUp.push_back(std::array<float, nMax>());
			jetMassJecDown.push_back(std::array<float, nMax>());
		}

		for (const std::pair<std::string, boost::property_tree::ptree> jecSyst : configTree.get_child("Producer.Jet.JECSystematic")) {
			fatJetPtJecUp.push_back(std::array<float, nMax>());
			fatJetPtJecDown.push_back(std::array<float, nMax>());
			fatJetMassJecUp.push_back(std::array<float, nMax>());
			fatJetMassJecDown.push_back(std::array<float, nMax>());
		}

		for (const std::pair<std::string, boost::property_tree::ptree> j : configTree.get_child("Producer.Jet.BTagSystematic")) {
			jetDeepJetLooseSfUp.push_back(std::array<float, nMax>());
			jetDeepJetLooseSfDown.push_back(std::array<float, nMax>());
			jetDeepJetMediumSfUp.push_back(std::array<float, nMax>());
			jetDeepJetMediumSfDown.push_back(std::array<float, nMax>());
			jetDeepJetTightSfUp.push_back(std::array<float, nMax>());
			jetDeepJetTightSfDown.push_back(std::array<float, nMax>());
		}
		for (const std::pair<std::string, boost::property_tree::ptree> j : configTree.get_child("Producer.Jet.LightTagSystematic")) {
			jetDeepJetLooseLightSfUp.push_back(std::array<float, nMax>());
			jetDeepJetLooseLightSfDown.push_back(std::array<float, nMax>());
			jetDeepJetMediumLightSfUp.push_back(std::array<float, nMax>());
			jetDeepJetMediumLightSfDown.push_back(std::array<float, nMax>());
			jetDeepJetTightLightSfUp.push_back(std::array<float, nMax>());
			jetDeepJetTightLightSfDown.push_back(std::array<float, nMax>());
		}
}

void Susy1LeptonProduct::RegisterTrigger (const std::vector<std::string> &triggerNames,const std::vector<std::string> &metTriggerNames, const std::vector<std::shared_ptr<TTree>> &outputTrees) {
	/*##################################################################################################
	#   The std::vector<bool> is special and cannot be used for getting the address of its elements    #
	#   But you can use std::vector<short> and store it as a boolian in the output tree                #
	#   https://en.cppreference.com/w/cpp/container/vector_bool                                        #
	##################################################################################################*/
	triggerValues = std::vector<short>(triggerNames.size(), true);
	metTriggerValues = std::vector<short>(metTriggerNames.size(), true);

	for (const std::shared_ptr<TTree>& tree: outputTrees) {
		for (int iTrigger = 0; iTrigger < triggerNames.size(); iTrigger++) {
			tree->Branch(triggerNames[iTrigger].c_str (), &triggerValues[iTrigger], (triggerNames[iTrigger] + "/O").c_str ());
		}
		for (int iTrigger = 0; iTrigger < metTriggerNames.size(); iTrigger++) {
			tree->Branch(metTriggerNames[iTrigger].c_str (), &metTriggerValues[iTrigger], (metTriggerNames[iTrigger] + "/O").c_str ());
		}
	}

	for (const std::shared_ptr<TTree> &tree : outputTrees) {
		tree->Branch("HLT_EleOr", &hltEleOr);
		tree->Branch("HLT_MuonOr", &hltMuOr);
		tree->Branch("HLT_MetOr", &hltMetOr);
	}
}

void Susy1LeptonProduct::RegisterMetFilter (const std::vector<std::string> &metFilterNames, const std::vector<std::shared_ptr<TTree>> &outputTrees) {
	/*##################################################################################################
	#   The std::vector<bool> is special and cannot be used for getting the address of its elements    #
	#   But you can use std::vector<short> and store it as a boolian in the output tree                #
	#   https://en.cppreference.com/w/cpp/container/vector_bool                                        #
	##################################################################################################*/
	metFilterValues = std::vector<short>(metFilterNames.size(), true);

	//for (const std::shared_ptr<TTree>& tree: outputTrees) {
	//	for (int iFilter = 0; iFilter < metFilterNames.size(); iFilter++) {
	//		tree->Branch(metFilterNames[iFilter].c_str (), &metFilterValues[iFilter], (metFilterNames[iFilter] + "/O").c_str ());
	//	}
	//}
}

void Susy1LeptonProduct::RegisterOutput(std::vector<std::shared_ptr<TTree>> outputTrees, const pt::ptree &configTree) {
	/*#########################################################################################
	#   Set Branches of output tree                                                           #
	#   See the following link for type definitions                                           #
	#   https://root.cern.ch/doc/master/classTBranch.html#ac0412c423e6c8388b42247e0410cf822   #
	#########################################################################################*/
	for (const std::shared_ptr<TTree> &tree : outputTrees) {
		tree->Branch("nMuon", &nMuon);
		tree->Branch("nGoodMuon", &nGoodMuon);
		tree->Branch("nVetoMuon", &nVetoMuon);
		tree->Branch("nAntiSelectedMuon", &nAntiSelectedMuon);

		tree->Branch("MuonPt", muonPt.data(), "MuonPt[nMuon]/F");
		tree->Branch("MuonEta", muonEta.data(), "MuonEta[nMuon]/F");
		tree->Branch("MuonPhi", muonPhi.data(), "MuonPhi[nMuon]/F");
		tree->Branch("MuonMass", muonMass.data(), "MuonMass[nMuon]/F");
		tree->Branch("MuonMiniIso", muonMiniIso.data(), "MuonMiniIso[nMuon]/F");

		tree->Branch("MuonLooseId", muonLooseId.data(), "MuonLooseId[nMuon]/O");
		tree->Branch("MuonMediumId", muonMediumId.data(), "MuonMediumId[nMuon]/O");
		tree->Branch("MuonTightId", muonTightId.data(), "MuonTightId[nMuon]/O");
		tree->Branch("MuonIsGood", muonIsGood.data(), "MuonIsGood[nMuon]/O");
		tree->Branch("MuonIsVeto", muonIsVeto.data(), "MuonIsVeto[nMuon]/O");
		tree->Branch("MuonIsAntiSelected", muonIsAntiSelected.data(), "MuonIsAntiSelected[nMuon]/O");

		if (!isData) {
			tree->Branch("MuonLooseIsoSf", muonLooseIsoSf.data(), "MuonLooseIsoSf[nMuon]/F");
			tree->Branch("MuonLooseIsoSfUp", muonLooseIsoSfUp.data(), "MuonLooseIsoSfUp[nMuon]/F");
			tree->Branch("MuonLooseIsoSfDown", muonLooseIsoSfDown.data(), "MuonLooseIsoSfDown[nMuon]/F");
			tree->Branch("MuonMediumIsoSf", muonMediumIsoSf.data(), "MuonMediumIsoSf[nMuon]/F");
			tree->Branch("MuonMediumIsoSfUp", muonMediumIsoSfUp.data(), "MuonMediumIsoSfUp[nMuon]/F");
			tree->Branch("MuonMediumIsoSfDown", muonMediumIsoSfDown.data(), "MuonMediumIsoSfDown[nMuon]/F");
			tree->Branch("MuonLooseSf", muonLooseSf.data(), "MuonLooseSf[nMuon]/F");
			tree->Branch("MuonLooseSfUp", muonLooseSfUp.data(), "MuonLooseSfUp[nMuon]/F");
			tree->Branch("MuonLooseSfDown", muonLooseSfDown.data(), "MuonLooseSfDown[nMuon]/F");
			tree->Branch("MuonMediumSf", muonMediumSf.data(), "MuonMediumSf[nMuon]/F");
			tree->Branch("MuonMediumSfUp", muonMediumSfUp.data(), "MuonMediumSfUp[nMuon]/F");
			tree->Branch("MuonMediumSfDown", muonMediumSfDown.data(), "MuonMediumSfDown[nMuon]/F");
			tree->Branch("MuonTightSf", muonTightSf.data(), "MuonTightSf[nMuon]/F");
			tree->Branch("MuonTightSfUp", muonTightSfUp.data(), "MuonTightSfUp[nMuon]/F");
			tree->Branch("MuonTightSfDown", muonTightSfDown.data(), "MuonTightSfDown[nMuon]/F");
			tree->Branch("MuonTriggerSf", muonTriggerSf.data(), "MuonTriggerSf[nMuon]/F");
			tree->Branch("MuonTriggerSfUp", muonTriggerSfUp.data(), "MuonTriggerSfUp[nMuon]/F");
			tree->Branch("MuonTriggerSfDown", muonTriggerSfDown.data(), "MuonTriggerSfDown[nMuon]/F");

			if (isFastSim) {
				tree->Branch("MuonLooseFastSf", muonLooseFastSf.data(), "MuonLooseFastSf[nMuon]/F");
				tree->Branch("MuonLooseFastSfUp", muonLooseFastSfUp.data(), "MuonLooseFastSfUp[nMuon]/F");
				tree->Branch("MuonLooseFastSfDown", muonLooseFastSfDown.data(), "MuonLooseFastSfDown[nMuon]/F");

				tree->Branch("MuonMediumFastSf", muonLooseFastSf.data(), "MuonLooseFastSf[nMuon]/F");
				tree->Branch("MuonMediumFastSfUp", muonLooseFastSfUp.data(), "MuonLooseFastSfUp[nMuon]/F");
				tree->Branch("MuonMediumFastSfDown", muonLooseFastSfDown.data(), "MuonLooseFastSfDown[nMuon]/F");

				tree->Branch("ElectronVetoFastSf", electronVetoFastSf.data(), "ElectronVetoFastSf[nElectron]/F");
				tree->Branch("ElectronVetoFastSfUp", electronVetoFastSfUp.data(), "ElectronVetoFastSfUp[nElectron]/F");
				tree->Branch("ElectronVetoFastSfDown", electronVetoFastSfDown.data(), "ElectronVetoFastSfDown[nElectron]/F");

				tree->Branch("ElectronTightFastSf", electronTightFastSf.data(), "ElectronTightFastSf[nElectron]/F");
				tree->Branch("ElectronTightFastSfUp", electronTightFastSfUp.data(), "ElectronTightFastSfUp[nElectron]/F");
				tree->Branch("ElectronTightFastSfDown", electronTightFastSfDown.data(), "ElectronTightFastSfDown[nElectron]/F");

				tree->Branch("ElectronVetoMVAFastSf", electronVetoMvaFastSf.data(), "ElectronVetoMVAFastSf[nElectron]/F");
				tree->Branch("ElectronVetoMVAFastSfUp", electronVetoMvaFastSfUp.data(), "ElectronVetoMVAFastSfUp[nElectron]/F");
				tree->Branch("ElectronVetoMVAFastSfDown", electronVetoMvaFastSfDown.data(), "ElectronVetoMVAFastSfDown[nElectron]/F");

				tree->Branch("ElectronTightMVAFastSf", electronTightMvaFastSf.data(), "ElectronTightMVAFastSf[nElectron]/F");
				tree->Branch("ElectronTightMVAFastSfUp", electronTightMvaFastSfUp.data(), "ElectronTightMVAFastSfUp[nElectron]/F");
				tree->Branch("ElectronTightMVAFastSfDown", electronTightMvaFastSfDown.data(), "ElectronTightMVAFastSfDown[nElectron]/F");
			}
		}

		tree->Branch("MuonCharge", muonCharge.data(), "MuonCharge[nMuon]/I");
		tree->Branch("MuonPdgId", muonPdgId.data(), "MuonPdgId[nMuon]/I");
		tree->Branch("MuonGenMatchedIndex", muonMass.data(), "MuonGenMatchedIndex[nMuon]/I");
		tree->Branch("MuonCutBasedId", muonCutBasedId.data(), "MuonCutBasedId[nMuon]/I");

		tree->Branch("nElectron", &nElectron);
		tree->Branch("nGoodElectron", &nGoodElectron);
		tree->Branch("nVetoElectron", &nVetoElectron);
		tree->Branch("nAntiSelectedElectron", &nAntiSelectedElectron);

		tree->Branch("ElectronPt", electronPt.data(), "ElectronPt[nElectron]/F");
		tree->Branch("ElectronEta", electronEta.data(), "ElectronEta[nElectron]/F");
		tree->Branch("ElectronPhi", electronPhi.data(), "ElectronPhi[nElectron]/F");
		tree->Branch("ElectronMass", electronMass.data(), "ElectronMass,[nElectron]/F");
		tree->Branch("ElectronDxy", electronDxy.data(), "ElectronDxy[nElectron]/F");
		tree->Branch("ElectronDz", electronDz.data(), "ElectronDz[nElectron]/F");
		tree->Branch("ElectronECorr", electronECorr.data(), "ElectronECorr,[nElectron]/F");
		tree->Branch("ElectronMiniIso", electronMiniIso.data(), "ElectronMiniIso[nElectron]/F");
		//tree->Branch("ElectronIso03", electronIso03.data(), "ElectronIso03[nMuon]/F");
		//tree->Branch("ElectronIso04", electronIso04.data(), "ElectronIso04[nMuon]/F");
		tree->Branch("ElectronRelJetIso", electronRelJetIso.data(), "ElectronRelJetIso,[nElectron]/F");
		if (isFastSim) {
			tree->Branch("ElectronEnergyScaleUp", electronEnergyScaleUp.data(), "ElectronEnergyScaleUp[nElectron]/F");
			tree->Branch("ElectronEnergyScaleDown", electronEnergyScaleDown.data(), "ElectronEnergyScaleDown,[nElectron]/F");
			tree->Branch("ElectronEnergySigmaUp", electronEnergySigmaUp.data(), "ElectronEnergySigmaUp[nElectron]/F");
			tree->Branch("ElectronEnergySigmaDown", electronEnergySigmaDown.data(), "ElectronEnergySigmaDown[nElectron]/F");
		}

		tree->Branch("ElectronLooseMvaId", electronLooseMvaId.data(), "ElectronLooseMvaId[nElectron]/O");
		tree->Branch("ElectronMediumMvaId", electronMediumMvaId.data(), "ElectronMediumMvaId[nElectron]/O");
		tree->Branch("ElectronTightMvaId", electronTightMvaId.data(), "ElectronTightMvaId[nElectron]/O");

		tree->Branch("ElectronIsGood", muonIsGood.data(), "ElectronIsGood[nElectron]/O");
		tree->Branch("ElectronIsVeto", muonIsVeto.data(), "ElectronIsVeto[nElectron]/O");

		tree->Branch("ElectronCharge", electronCharge.data(), "ElectronCharge[nElectron]/I");
		tree->Branch("ElectronPdgId", electronPdgId.data(), "ElectronPdgId[nElectron]/I");
		tree->Branch("ElectronCutBasedId", electronCutBasedId.data(), "ElectronCutBasedId[nElectron]/I");
		//tree->Branch("ElectronNLostHits", electronNLostHits.data(), "ElectronNLostHits[nMuon]/I");

		tree->Branch("ElectronLooseId", electronLooseId.data(), "ElectronLooseId[nElectron]/O");
		tree->Branch("ElectronMediumId", electronMediumId.data(), "ElectronMediumId[nElectron]/O");
		tree->Branch("ElectronTightId", electronTightId.data(), "ElectronTightId[nElectron]/O");
		tree->Branch("ElectronVetoId", electronVetoId.data(), "ElectronVetoId[nElectron]/O");
		//tree->Branch("ElectronConvVeto", electronConvVeto.data(), "ElectronConvVeto[nMuon]/O");

		if (!isData) {
			tree->Branch("ElectronRecoSf", electronRecoSf.data(), "ElectronRecoSf[nElectron]/F");
			tree->Branch("ElectronLooseSf", electronLooseSf.data(), "ElectronLooseSf[nElectron]/F");
			tree->Branch("ElectronMediumSf", electronMediumSf.data(), "ElectronMediumSf[nElectron]/F");
			tree->Branch("ElectronTightSf", electronTightSf.data(), "ElectronTightSf[nElectron]/F");
			tree->Branch("ElectronMediumMvaSf", electronMediumMvaSf.data(), "ElectronMediumMvaSf[nElectron]/F");
			tree->Branch("ElectronTightMvaSf", electronTightMvaSf.data(), "ElectronTightMvaSf[nElectron]/F");
			tree->Branch("ElectronRecoSfUp", electronRecoSfUp.data(), "ElectronRecoSfUp[nElectron]/F");
			tree->Branch("ElectronRecoSfDown", electronRecoSfDown.data(), "ElectronRecoSfDown[nElectron]/F");
			tree->Branch("ElectronLooseSfUp", electronLooseSfUp.data(), "ElectronLooseSfUp[nElectron]/F");
			tree->Branch("ElectronLooseSfDown", electronLooseSfDown.data(), "ElectronLooseSfDown[nElectron]/F");
			tree->Branch("ElectronMediumSfUp", electronMediumSfUp.data(), "ElectronMediumSfUp[nElectron]/F");
			tree->Branch("ElectronMediumSfDown", electronMediumSfDown.data(), "ElectronMediumSfDown[nElectron]/F");
			tree->Branch("ElectronTightSfUp", electronTightSfUp.data(), "ElectronTightSfUp[nElectron]/F");
			tree->Branch("ElectronTightSfDown", electronTightSfDown.data(), "ElectronTightSfDown[nElectron]/F");
			tree->Branch("ElectronMediumMvaSfUp", electronMediumMvaSfUp.data(), "ElectronMediumMvaSfUp[nElectron]/F");
			tree->Branch("ElectronMediumMvaSfDown", electronMediumMvaSfDown.data(), "ElectronMediumMvaSfDown[nElectron]/F");
			tree->Branch("ElectronTightMvaSfUp", electronTightMvaSfUp.data(), "ElectronTightMvaSfUp[nElectron]/F");
			tree->Branch("ElectronTightMvaSfDown", electronTightMvaSfDown.data(), "ElectronTightMvaSfDown[nElectron]/F");
		}

		tree->Branch("MetPt", &metPt);
		tree->Branch("CaloMET_pt", &CaloMET_pt);
		tree->Branch("MetPhi", &metPhi);
		tree->Branch("nJet", &nJet);
		tree->Branch("JetPt", jetPt.data(), "JetPt[nJet]/F");
		tree->Branch("JetPt_JERUp", jetPtJerUp.data(), "JetPt_JERUp[nJet]/F");
		tree->Branch("JetPt_JERDown", jetPtJerDown.data(), "JetPt_JERDown[nJet]/F");
		tree->Branch("JetEta", jetEta.data(), "JetEta[nJet]/F");
		tree->Branch("JetPhi", jetPhi.data(), "JetPhi[nJet]/F");
		tree->Branch("JetMass", jetMass.data(), "JetMass,[nJet]/F");
		tree->Branch("JetMass_JERUp", jetMassJerUp.data(), "JetMass_JERUp[nJet]/F");
		tree->Branch("JetMass_JERDown", jetMassJerDown.data(), "JetMass_JERDown[nJet]/F");
		//tree->Branch("JetDeepJet", jetDeepJet.data(), "JetDeepJet,[nJet]/F");
		tree->Branch("JetDeepJetLooseId", jetDeepJetLooseId.data(), "JetDeepJetLooseId;[nJet]/O");
		tree->Branch("JetDeepJetMediumId", jetDeepJetMediumId.data(), "JetDeepJetMediumId[nJet]/O");
		tree->Branch("JetDeepJetTightId", jetDeepJetTightId.data(), "JetDeepJetTightId[nJet]/O");
		tree->Branch("nDeepJetLooseBTag", &nDeepJetLooseBTag);
		tree->Branch("nDeepJetMediumBTag", &nDeepJetMediumBTag);
		tree->Branch("nDeepJetTightBTag", &nDeepJetTightBTag);
		tree->Branch("JetId", &jetId);

		if (!isData) {
			tree->Branch("JetDeepJetLooseSf", jetDeepJetLooseSf.data(), "JetDeepJetLooseSf[nJet]/F");
			tree->Branch("JetDeepJetMediumSf", jetDeepJetMediumSf.data(), "JetDeepJetMediumSf[nJet]/F");
			tree->Branch("JetDeepJetTightSf", jetDeepJetTightSf.data(), "JetDeepJetTightSf[nJet]/F");

			if (isFastSim) {
				tree->Branch("JetDeepJetLooseFastSf", jetDeepJetLooseFastSf.data(), "JetDeepJetLooseFastSf[nJet]/F");
				tree->Branch("JetDeepJetMediumFastSf", jetDeepJetMediumFastSf.data(), "JetDeepJetMediumFastSf[nJet]/F");
				tree->Branch("JetDeepJetTightFastSf", jetDeepJetTightFastSf.data(), "JetDeepJetTightFastSf[nJet]/F");
			}
		}

		tree->Branch("nFatJet", &nFatJet);
		tree->Branch("FatJetPt", fatJetPt.data(), "FatJetPt[nFatJet]/F");

		tree->Branch("FatJetEta", fatJetEta.data(), "FatJetEta[nFatJet]/F");
		tree->Branch("FatJetPhi", fatJetPhi.data(), "FatJetPhi[nFatJet]/F");
		tree->Branch("FatJetMass", fatJetMass.data(), "FatJetMass[nFatJet]/F");
		tree->Branch("FatJetArea", fatJetArea.data(), "FatJetArea[nFatJet]/F");
		//tree->Branch("FatJetRawFactor", fatJetRawFactor.data(), "FatJetRawFactor[nFatJet]/F");
		tree->Branch("FatJetId", fatJetId.data(), "FatJetId[nFatJet]/F");
		tree->Branch("FatJetDeepTagMDTvsQCD", fatJetDeepTagMDTvsQCD.data(), "FatJetDeepTagMDTvsQCD[nFatJet]/F");
		tree->Branch("FatJetDeepTagMDWvsQCD", fatJetDeepTagMDWvsQCD.data(), "FatJetDeepTagMDWvsQCD[nFatJet]/F");
		tree->Branch("FatJetDeepTagTvsQCD", fatJetDeepTagTvsQCD.data(), "FatJetDeepTagTvsQCD[nFatJet]/F");
		tree->Branch("FatJetDeepTagWvsQCD", fatJetDeepTagWvsQCD.data(), "FatJetDeepTagWvsQCD[nFatJet]/F");

		tree->Branch("HT", &HT);
		tree->Branch("LT", &LT);
		tree->Branch("LP", &LP);
		tree->Branch("DeltaPhi", &deltaPhi);
		tree->Branch("WBosonMt", &wBosonMt);

		tree->Branch("nIsoTrack", &nIsoTrack);
		tree->Branch("IsoTrackVeto", &isoTrackVeto);
		tree->Branch("IsoTrackPdgId", isoTrackPdgId.data(), "IsoTrackPdgId[nIsoTrack]/I");
		tree->Branch("IsoTrackPt", isoTrackPt.data(), "IsoTrackPt[nIsoTrack]/F");
		tree->Branch("IsoTrackMt2", isoTrackMt2.data(), "IsoTrackMt2[nIsoTrack]/F");
		tree->Branch("IsoTrackIsHadronicDecay", isoTrackIsHadronicDecay.data(), "IsoTrackIsHadronicDecay[nIsoTrack]/O");

		if (!isData) {
			tree->Branch("nTrueInt", &nTrueInt);
			tree->Branch("PileUpWeight", &pileUpWeight);
			tree->Branch("PileUpWeightUp", &pileUpWeightUp);
			tree->Branch("PileUpWeightDown", &pileUpWeightDown);

			if (!isFastSim) {
				tree->Branch("PreFireWeight", &preFire);
				tree->Branch("PreFireWeightUp", &preFireUp);
				tree->Branch("PreFireWeightDown", &preFireDown);
			}
			tree->Branch("nPdfWeight", &nPdfWeight);
			tree->Branch("LHEPdfWeight", pdfWeight.data(), "pdfWeight[nPdfWeight]/F");
			tree->Branch("nScaleWeight", &nScaleWeight);
			tree->Branch("LHEScaleWeight", scaleWeight.data(), "scaleWeight[nScaleWeight]/F");
			tree->Branch("GenWeight", &genWeight);

			tree->Branch("nGenLepton", &nGenLepton);
			tree->Branch("GrandMotherPdgId", grandMotherPdgId.data(), "GrandMotherPdgId[nGenLepton]/I");
			tree->Branch("GenLepGrandMotherPdgId", genLepGrandMotherPdgId.data(), "GenLepGrandMotherPdgId[nGenLepton]/I");
			tree->Branch("GenLepMotherPdgId", genLepMotherPdgId.data(), "GenLepMotherPdgId[nGenLepton]/I");

			tree->Branch("nGenTau", &nGenTau);
			tree->Branch("GenTauGrandMotherPdgId", genTauGrandMotherPdgId.data(), "GenTauGrandMotherPdgId[nGenTau]/I");
			tree->Branch("GenTauMotherPdgId", genTauMotherPdgId.data(), "GenTauMotherPdgId[nGenTau]/I");

			tree->Branch("nGenLeptonFromTau", &nGenLeptonFromTau);
			tree->Branch("LeptonDecayChannelFlag", &leptonDecayChannelFlag);
			tree->Branch("IsDiLeptonEvent", &isDiLeptonEvent);
			tree->Branch("IsHadTauEvent", &isHadTauEvent);
			tree->Branch("LeptonsInAcceptance", &leptonsInAcceptance);

			tree->Branch("nGenMatchedW", &nGenMatchedW);
			tree->Branch("GenDeltaPhiLepWSum", genDeltaPhiLepWSum.data(), "GenDeltaPhiLepWSum[nGenMatchedW]/F");
			tree->Branch("GenDeltaPhiLepWDirect", genDeltaPhiLepWDirect.data(), "GenDeltaPhiLepWDirect[nGenMatchedW]/F");
			tree->Branch("GenWSumMass", genWSumMass.data(), "GenWSumMass[nGenMatchedW]/F");
			tree->Branch("GenWDirectMass", genWDirectMass.data(), "GenWDirectMass[nGenMatchedW]/F");
			tree->Branch("GenMTLepNu", genMTLepNu.data(), "GenMTLepNu[nGenMatchedW]/F");

			tree->Branch("nGenNeutrino", &nGenNeutrino);
			tree->Branch("GenNeutrinoPt", genNeutrinoPt.data(), "GenNeutrinoPt[nGenNeutrino]/F");

			int iSyst = 0;
			for (const std::pair<std::string, boost::property_tree::ptree> jecSyst : configTree.get_child("Producer.Jet.JECSystematic")) {
				std::string syst = jecSyst.second.get_value<std::string>();
				tree->Branch(("JetPt_" + syst + "Up").c_str (), jetPtJecUp[iSyst].data(), ("JetPt_" + syst + "Up[nJet]/F").c_str ());
				tree->Branch(("JetPt_" + syst + "Down").c_str (), jetPtJecDown[iSyst].data(), ("JetPt_" + syst + "Down[nJet]/F").c_str ());
				iSyst++;
			}

			iSyst = 0;
			for (const std::pair<std::string, boost::property_tree::ptree> jecSyst : configTree.get_child("Producer.Jet.JECSystematic")) {
				std::string syst = jecSyst.second.get_value<std::string>();
				tree->Branch(("FatJetPt_" + syst + "Up").c_str (), fatJetPtJecUp[iSyst].data(), ("FatJetPt_" + syst + "Up[nFatJet]/F").c_str ());
				tree->Branch(("FatJetPt_" + syst + "Down").c_str (), fatJetPtJecDown[iSyst].data(), ("FatJetPt_" + syst + "Down[nFatJet]/F").c_str ());
				iSyst++;
			}

			for (const std::pair<std::string, boost::property_tree::ptree> j : configTree.get_child("Producer.Jet.BTagSystematic")) {
				std::string bSyst = j.second.get_value<std::string>();
				tree->Branch(("JetDeepJetLooseSf_" + bSyst + "Up").c_str (), jetDeepJetLooseSfUp[iSyst].data(), ("JetDeepJetLooseSf_" + bSyst + "Up[nJet]/F").c_str ());
				tree->Branch(("JetDeepJetLooseSf_" + bSyst + "Down").c_str (), jetDeepJetLooseSfDown[iSyst].data(), ("JetDeepJetLooseSf_" + bSyst + "Down[nJet]/F").c_str ());

				tree->Branch(("JetDeepJetMediumSf_" + bSyst + "Up").c_str (), jetDeepJetMediumSfUp[iSyst].data(), ("JetDeepJetMediumSf_" + bSyst + "Up[nJet]/F").c_str ());
				tree->Branch(("JetDeepJetMediumSf_" + bSyst + "Down").c_str (), jetDeepJetMediumSfDown[iSyst].data(), ("JetDeepJetMediumSf_" + bSyst + "Down[nJet]/F").c_str ());

				tree->Branch(("JetDeepJetTightSf_" + bSyst + "Up").c_str (), jetDeepJetTightSfUp[iSyst].data(), ("JetDeepJetTightSf_" + bSyst + "Up[nJet]/F").c_str ());
				tree->Branch(("JetDeepJetTightSf_" + bSyst + "Down").c_str (), jetDeepJetTightSfDown[iSyst].data(), ("JetDeepJetTightSf_" + bSyst + "Down[nJet]/F").c_str ());
				iSyst++;
			}

			iSyst = 0;
			for (const std::pair<std::string, boost::property_tree::ptree> j : configTree.get_child("Producer.Jet.LightTagSystematic")) {
				std::string bSyst = j.second.get_value<std::string>();
				tree->Branch(("JetDeepJetLooseLightSf_" + bSyst + "Up").c_str (), jetDeepJetLooseLightSfUp[iSyst].data(), ("JetDeepJetLooseLightSf_" + bSyst + "Up[nJet]/F").c_str ());
				tree->Branch(("JetDeepJetLooseLightSf_" + bSyst + "Down").c_str (), jetDeepJetLooseLightSfDown[iSyst].data(), ("JetDeepJetLooseLightSf_" + bSyst + "Down[nJet]/F").c_str ());

				tree->Branch(("JetDeepJetMediumLightSf_" + bSyst + "Up").c_str (), jetDeepJetMediumLightSfUp[iSyst].data(), ("JetDeepJetMediumLightSf_" + bSyst + "Up[nJet]/F").c_str ());
				tree->Branch(("JetDeepJetMediumLightSf_" + bSyst + "Down").c_str (), jetDeepJetMediumLightSfDown[iSyst].data(), ("JetDeepJetMediumLightSf_" + bSyst + "Down[nJet]/F").c_str ());

				tree->Branch(("JetDeepJetTightLightSf_" + bSyst + "Up").c_str (), jetDeepJetTightLightSfUp[iSyst].data(), ("JetDeepJetTightLightSf_" + bSyst + "Up[nJet]/F").c_str ());
				tree->Branch(("JetDeepJetTightLightSf_" + bSyst + "Down").c_str (), jetDeepJetTightLightSfDown[iSyst].data(), ("JetDeepJetTightLightSf_" + bSyst + "Down[nJet]/F").c_str ());
				iSyst++;
			}

			if (isFastSim) {
				tree->Branch("JetDeepJetLooseFastSfUp", jetDeepJetLooseFastSfUp.data(), "JetDeepJetLooseFastSfUp[nJet]/F");
				tree->Branch("JetDeepJetMediumFastSfUp", jetDeepJetMediumFastSfUp.data(), "JetDeepJetMediumFastSfUp[nJet]/F");
				tree->Branch("JetDeepJetTightFastSfUp", jetDeepJetTightFastSfUp.data(), "JetDeepJetTightFastSfUp[nJet]/F");

				tree->Branch("JetDeepJetLooseFastSfDown", jetDeepJetLooseFastSfDown.data(), "JetDeepJetLooseFastSfDown[nJet]/F");
				tree->Branch("JetDeepJetMediumFastSfDown", jetDeepJetMediumFastSfDown.data(), "JetDeepJetMediumFastSfDown[nJet]/F");
				tree->Branch("JetDeepJetTightFastSfDown", jetDeepJetTightFastSfDown.data(), "JetDeepJetTightFastSfDown[nJet]/F");
			}

		}

		if (isFastSim) {
			tree->Branch("susyXSectionNLO", &susyXSectionNLO);
			tree->Branch("susyXSectionNLLO", &susyXSectionNLLO);
			tree->Branch("susyXSectionNLOUp", &susyXSectionNLOUp);
			tree->Branch("susyXSectionNLLOUp", &susyXSectionNLLOUp);
			tree->Branch("susyXSectionNLODown", &susyXSectionNLODown);
			tree->Branch("susyXSectionNLLODown", &susyXSectionNLLODown);
			tree->Branch("mStop", &stopMass);
			tree->Branch("mGluino", &gluinoMass);
			tree->Branch("mNeutralino", &neutralinoMass);
			tree->Branch("mChargino", &charginoMass);
		}
	}
}

void Susy1LeptonProduct::WriteMetaData(TFile &outputFile) {
	//TTree metaData("MetaData", "MetaData");
	metaData.SetName("MetaData");
	metaData.SetTitle("MetaData");
	metaData.Branch("Era", &era);
	metaData.Branch("PreVFP", &preVFP);
	metaData.Branch("IsData", &isData);
	metaData.Branch("IsFastSim", &isFastSim);
	metaData.Branch("SampleName", &sampleName);
	metaData.Branch("primaryDataset", &primaryDataset);
	if (isData) {
		metaData.Branch("runPeriod", &runPeriod);
	} else {
		metaData.Branch("xSection", &xSection);
		if (era == 2016) {
			if (preVFP) {
				luminosity = 19.5; //fb
			} else {
				luminosity = 16.5; //fb
			}
		} else if (era == 2017) {
			luminosity = 41.48; //fb
		} else if (era == 2018) {
			luminosity = 59.83; //fb
		}
		metaData.Branch("Luminosity", &luminosity);
	}
	//std::cout << "Filling MetaData primaryDataset with" <<  primaryDataset <<  std::endl; // datasetDecider <<
	metaData.SetDirectory(&outputFile);
	metaData.Fill();
	metaData.Write(0, TObject::kOverwrite);
}
