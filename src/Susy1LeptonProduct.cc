#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

Susy1LeptonProduct::Susy1LeptonProduct(const int &era, const bool &isData, const std::string &sampleName, const char &runPeriod, const double &xSection, TFile &outputFile) :
	era(era),
	isData(isData),
	sampleName(sampleName),
	runPeriod(runPeriod),
	xSection(xSection)
	{
		std::string outputFileName = outputFile.GetName();
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

		//metaData = std::make_shared<TTree>("MetaData", "MetaData");
		TTree metaData("MetaData", "MetaData");
		metaData.Branch("era", &this->era);
		metaData.Branch("preVFP", &this->preVFP);
		metaData.Branch("isData", &this->isData);
		metaData.Branch("sampleName", &this->sampleName);
		if (isData) {
			metaData.Branch("runPeriod", &this->runPeriod);
		} else {
			metaData.Branch("xSection", &this->xSection);
			if (era == 2016) {
				if (this->preVFP) {
					luminosity = 19.5; //fb
				} else {
					luminosity = 16.5; //fb
				}
			} else if (era == 2017) {
				luminosity = 41.48; //fb
			} else if (era == 2018) {
				luminosity = 59.83; //fb
			}
			metaData.Branch("luminosity", &this->luminosity);
		}

		metaData.SetDirectory(&outputFile);
		metaData.Fill();
		metaData.Write(0, TObject::kOverwrite);
}


void Susy1LeptonProduct::RegisterOutput(std::vector<std::shared_ptr<TTree>> outputTrees) {
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

		tree->Branch("MuonPt", muonPt.data(), "MuonPt[nMuon]/D");
		tree->Branch("MuonPtVector", &muonPtVector);
		tree->Branch("MuonEta", muonEta.data(), "MuonEta[nMuon]/D");
		tree->Branch("MuonPhi", muonPhi.data(), "MuonPhi[nMuon]/D");
		tree->Branch("MuonMass", muonMass.data(), "MuonMass[nMuon]/D");
		tree->Branch("MuonMiniIso", muonMiniIso.data(), "MuonMiniIso[nMuon]/D");

		tree->Branch("MuonLooseId", muonLooseId.data(), "MuonLooseId[nMuon]/O");
		tree->Branch("MuonMediumId", muonMediumId.data(), "MuonMediumId[nMuon]/O");
		tree->Branch("MuonTightId", muonTightId.data(), "MuonTightId[nMuon]/O");
		tree->Branch("MuonIsGood", muonIsGood.data(), "MuonIsGood[nMuon]/O");
		tree->Branch("MuonIsVeto", muonIsVeto.data(), "MuonIsVeto[nMuon]/O");
		tree->Branch("MuonIsAntiSelected", muonIsAntiSelected.data(), "MuonIsAntiSelected[nMuon]/O");

		tree->Branch("MuonCharge", muonCharge.data(), "MuonCharge[nMuon]/I");
		tree->Branch("MuonPdgId", muonPdgId.data(), "MuonPdgId[nMuon]/I");
		tree->Branch("MuonGenMatchedIndex", muonMass.data(), "MuonGenMatchedIndex[nMuon]/I");
		tree->Branch("MuonCutBasedId", muonCutBasedId.data(), "MuonCutBasedId[nMuon]/I");

		if (isData) {
			tree->Branch("MuonLooseIsoSf", muonLooseIsoSf.data(), "MuonLooseIsoSf[nMuon]/D");
			tree->Branch("MuonLooseIsoSfUp", muonLooseIsoSfUp.data(), "MuonLooseIsoSfUp[nMuon]/D");
			tree->Branch("MuonLooseIsoSfDown", muonLooseIsoSfDown.data(), "MuonLooseIsoSfDown[nMuon]/D");
			tree->Branch("MuonTightIsoSf", muonTightIsoSf.data(), "MuonTightIsoSf[nMuon]/D");
			tree->Branch("MuonTightIsoSfUp", muonTightIsoSfUp.data(), "MuonTightIsoSfUp[nMuon]/D");
			tree->Branch("MuonTightIsoSfDown", muonTightIsoSfDown.data(), "MuonTightIsoSfDown[nMuon]/D");
			tree->Branch("MuonLooseSf", muonLooseSf.data(), "MuonLooseSf[nMuon]/D");
			tree->Branch("MuonLooseSfUp", muonLooseSfUp.data(), "MuonLooseSfUp[nMuon]/D");
			tree->Branch("MuonLooseSfDown", muonLooseSfDown.data(), "MuonLooseSfDown[nMuon]/D");
			tree->Branch("MuonMediumSf", muonMediumSf.data(), "MuonMediumSf[nMuon]/D");
			tree->Branch("MuonMediumSfUp", muonMediumSfUp.data(), "MuonMediumSfUp[nMuon]/D");
			tree->Branch("MuonMediumSfDown", muonMediumSfDown.data(), "MuonMediumSfDown[nMuon]/D");
			tree->Branch("MuonTightSf", muonTightSf.data(), "MuonTightSf[nMuon]/D");
			tree->Branch("MuonTightSfUp", muonTightSfUp.data(), "MuonTightSfUp[nMuon]/D");
			tree->Branch("MuonTightSfDown", muonTightSfDown.data(), "MuonTightSfDown[nMuon]/D");
			tree->Branch("MuonTriggerSf", muonTriggerSf.data(), "MuonTriggerSf[nMuon]/D");
			tree->Branch("MuonTriggerSfUp", muonTriggerSfUp.data(), "MuonTriggerSfUp[nMuon]/D");
			tree->Branch("MuonTriggerSfDown", muonTriggerSfDown.data(), "MuonTriggerSfDown[nMuon]/D");
		}

		tree->Branch("nElectron", &nElectron);
		tree->Branch("nGoodElectron", &nGoodElectron);
		tree->Branch("nVetoElectron", &nVetoElectron);
		tree->Branch("nAntiSelectedElectron", &nAntiSelectedElectron);

		tree->Branch("ElectronPt", electronPt.data(), "ElectronPt[nMuon]/D");
		tree->Branch("ElectronEta", electronEta.data(), "ElectronEta[nMuon]/D");
		tree->Branch("ElectronPhi", electronPhi.data(), "ElectronPhi[nMuon]/D");
		tree->Branch("ElectronMass,", electronMass.data(), "ElectronMass,[nMuon]/D");
		tree->Branch("ElectronDxy", electronDxy.data(), "ElectronDxy[nMuon]/D");
		tree->Branch("ElectronDz", electronDz.data(), "ElectronDz[nMuon]/D");
		tree->Branch("ElectronECorr,", electronECorr.data(), "ElectronECorr,[nMuon]/D");
		tree->Branch("ElectronMiniIso", electronMiniIso.data(), "ElectronMiniIso[nMuon]/D");
		//tree->Branch("ElectronIso03", electronIso03.data(), "ElectronIso03[nMuon]/D");
		//tree->Branch("ElectronIso04", electronIso04.data(), "ElectronIso04[nMuon]/D");
		tree->Branch("ElectronRelJetIso,", electronRelJetIso.data(), "ElectronRelJetIso,[nMuon]/D");
		tree->Branch("ElectronEnergyScaleUp", electronEnergyScaleUp.data(), "ElectronEnergyScaleUp[nMuon]/D");
		tree->Branch("ElectronEnergyScaleDown,", electronEnergyScaleDown.data(), "ElectronEnergyScaleDown,[nMuon]/D");
		tree->Branch("ElectronEnergySigmaUp", electronEnergySigmaUp.data(), "ElectronEnergySigmaUp[nMuon]/D");
		tree->Branch("ElectronEnergySigmaDown", electronEnergySigmaDown.data(), "ElectronEnergySigmaDown[nMuon]/D");

		tree->Branch("ElectronLooseMvaId", electronLooseMvaId.data(), "ElectronLooseMvaId[nMuon]/O");
		tree->Branch("ElectronMediumMvaId", electronMediumMvaId.data(), "ElectronMediumMvaId[nMuon]/O");
		tree->Branch("ElectronTightMvaId", electronTightMvaId.data(), "ElectronTightMvaId[nMuon]/O");

		tree->Branch("ElectronCharge", electronCharge.data(), "ElectronCharge[nMuon]/I");
		tree->Branch("ElectronCutBasedId", electronCutBasedId.data(), "ElectronCutBasedId[nMuon]/I");
		//tree->Branch("ElectronNLostHits", electronNLostHits.data(), "ElectronNLostHits[nMuon]/I");

		tree->Branch("ElectronLooseId", electronLooseId.data(), "ElectronLooseId[nMuon]/O");
		tree->Branch("ElectronMediumId", electronMediumId.data(), "ElectronMediumId[nMuon]/O");
		tree->Branch("ElectronTightId", electronTightId.data(), "ElectronTightId[nMuon]/O");
		tree->Branch("ElectronVetoId", electronVetoId.data(), "ElectronVetoId[nMuon]/O");
		//tree->Branch("ElectronConvVeto", electronConvVeto.data(), "ElectronConvVeto[nMuon]/O");

		if (isData) {
			tree->Branch("ElectronRecoSf", electronRecoSf.data(), "ElectronRecoSf[nMuon]/D");
			tree->Branch("ElectronLooseSf", electronLooseSf.data(), "ElectronLooseSf[nMuon]/D");
			tree->Branch("ElectronMediumSf", electronMediumSf.data(), "ElectronMediumSf[nMuon]/D");
			tree->Branch("ElectronTightSf", electronTightSf.data(), "ElectronTightSf[nMuon]/D");
			tree->Branch("ElectronMediumMvaSf", electronMediumMvaSf.data(), "ElectronMediumMvaSf[nMuon]/D");
			tree->Branch("ElectronTightMvaSf", electronTightMvaSf.data(), "ElectronTightMvaSf[nMuon]/D");
			tree->Branch("ElectronRecoSfUp", electronRecoSfUp.data(), "ElectronRecoSfUp[nMuon]/D");
			tree->Branch("ElectronRecoSfDown", electronRecoSfDown.data(), "ElectronRecoSfDown[nMuon]/D");
			tree->Branch("ElectronLooseSfUp", electronLooseSfUp.data(), "ElectronLooseSfUp[nMuon]/D");
			tree->Branch("ElectronLooseSfDown", electronLooseSfDown.data(), "ElectronLooseSfDown[nMuon]/D");
			tree->Branch("ElectronMediumSfUp", electronMediumSfUp.data(), "ElectronMediumSfUp[nMuon]/D");
			tree->Branch("ElectronMediumSfDown", electronMediumSfDown.data(), "ElectronMediumSfDown[nMuon]/D");
			tree->Branch("ElectronTightSfUp", electronTightSfUp.data(), "ElectronTightSfUp[nMuon]/D");
			tree->Branch("ElectronTightSfDown", electronTightSfDown.data(), "ElectronTightSfDown[nMuon]/D");
			tree->Branch("ElectronMediumMvaSfUp", electronMediumMvaSfUp.data(), "ElectronMediumMvaSfUp[nMuon]/D");
			tree->Branch("ElectronMediumMvaSfDown", electronMediumMvaSfDown.data(), "ElectronMediumMvaSfDown[nMuon]/D");
			tree->Branch("ElectronTightMvaSfUp", electronTightMvaSfUp.data(), "ElectronTightMvaSfUp[nMuon]/D");
			tree->Branch("ElectronTightMvaSfDown", electronTightMvaSfDown.data(), "ElectronTightMvaSfDown[nMuon]/D");
		}

		tree->Branch("nJet", &nJet);
		tree->Branch("nDeepCsvLooseBTag", &nDeepCsvLooseBTag);
		tree->Branch("nDeepCsvMediumBTag", &nDeepCsvMediumBTag);
		tree->Branch("nDeepCsvTightBTag", &nDeepCsvTightBTag);
		tree->Branch("nDeepJetLooseBTag", &nDeepJetLooseBTag);
		tree->Branch("nDeepJetMediumBTag", &nDeepJetMediumBTag);
		tree->Branch("nDeepJetTightBTag", &nDeepJetTightBTag);

		tree->Branch("JetId", &jetId);
		tree->Branch("MetPt", &metPt);
		tree->Branch("MetPhi", &metPhi);
		tree->Branch("JetPt", jetPt.data(), "JetPt[nJet]/D");
		tree->Branch("JetEta", jetEta.data(), "JetEta[nJet]/D");
		tree->Branch("JetPhi", jetPhi.data(), "JetPhi[nJet]/D");
		tree->Branch("JetMass", jetMass.data(), "JetMass,[nJet]/D");
		tree->Branch("JetArea", jetArea.data(), "JetArea[nJet]/D");
		tree->Branch("JetRawFactor", jetRawFactor.data(), "JetRawFactor,[nJet]/D");
		tree->Branch("JetDeepCSV", jetDeepCsv.data(), "JetDeepCSV[nJet]/D");
		tree->Branch("JetDeepJet", jetDeepJet.data(), "JetDeepJet,[nJet]/D");
		tree->Branch("JetDeepCSVLooseId", jetDeepCsvLooseId.data(), "JetDeepCSVLooseId,[nJet]/O");
		tree->Branch("JetDeepCSVMediumId", jetDeepCsvMediumId.data(), "JetDeepCSVMediumId[nJet]/O");
		tree->Branch("JetDeepCSVTightId", jetDeepCsvTightId.data(), "JetDeepCSVTightId[nJet]/O");
		tree->Branch("JetDeepJetLooseId", jetDeepJetLooseId.data(), "JetDeepJetLooseId;[nJet]/O");
		tree->Branch("JetDeepJetMediumId", jetDeepJetMediumId.data(), "JetDeepJetMediumId[nJet]/O");
		tree->Branch("JetDeepJetTightId", jetDeepJetTightId.data(), "JetDeepJetTightId[nJet]/O");
		tree->Branch("wBosonMinMass", &wBosonMinMass);
		tree->Branch("wBosonMinMassPt", &wBosonMinMassPt);
		tree->Branch("wBosonBestMass", &wBosonBestMass);
		tree->Branch("wBosonBestMassPt", &wBosonBestMassPt);
		tree->Branch("topBestMass", &topBestMass);
		tree->Branch("topBestMassPt", &topBestMassPt);

		tree->Branch("nFatJet", &nFatJet);
		tree->Branch("FatJetMass", fatJetMass.data(), "FatJetMass[nFatJet]/D");
		tree->Branch("FatJetPt", fatJetPt.data(), "FatJetPt[nFatJet]/D");
		tree->Branch("FatJetEta", fatJetEta.data(), "FatJetEta[nFatJet]/D");
		tree->Branch("FatJetPhi", fatJetPhi.data(), "FatJetPhi[nFatJet]/D");
		tree->Branch("FatJetArea", fatJetArea.data(), "FatJetArea[nFatJet]/D");
		tree->Branch("FatJetRawFactor", fatJetRawFactor.data(), "FatJetRawFactor[nFatJet]/D");
		tree->Branch("FatJetId", fatJetId.data(), "FatJetId[nFatJet]/D");
		tree->Branch("FatJetDeepTagMDTvsQCD", fatJetDeepTagMDTvsQCD.data(), "FatJetDeepTagMDTvsQCD[nFatJet]/D");
		tree->Branch("FatJetDeepTagMDWvsQCD", fatJetDeepTagMDWvsQCD.data(), "FatJetDeepTagMDWvsQCD[nFatJet]/D");
		tree->Branch("FatJetDeepTagTvsQCD", fatJetDeepTagTvsQCD.data(), "FatJetDeepTagTvsQCD[nFatJet]/D");
		tree->Branch("FatJetDeepTagWvsQCD", fatJetDeepTagWvsQCD.data(), "FatJetDeepTagWvsQCD[nFatJet]/D");

		tree->Branch("HT", &HT);
		tree->Branch("LT", &LT);
		tree->Branch("LP", &LP);
		tree->Branch("DeltaPhi", &deltaPhi);
		tree->Branch("AbsoluteDeltaPhi", &absoluteDeltaPhi);
		tree->Branch("WBosonMt", &wBosonMt);

		tree->Branch("IsSignalRegion", &isSignalRegion);
		tree->Branch("nIsoTrack", &nIsoTrack);
		tree->Branch("IsoTrackVeto", &isoTrackVeto);
		tree->Branch("IsoTrackPdgId", isoTrackPdgId.data(), "IsoTrackPdgId[nIsoTrack]/I");
		tree->Branch("IsoTrackPt", isoTrackPt.data(), "IsoTrackPt[nIsoTrack]/D");
		tree->Branch("IsoTrackMt2", isoTrackMt2.data(), "IsoTrackMt2[nIsoTrack]/D");
		tree->Branch("IsoTrackIsHadronicDecay", isoTrackIsHadronicDecay.data(), "IsoTrackIsHadronicDecay[nIsoTrack]/O");
	}
}
