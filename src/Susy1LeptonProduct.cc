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
				eraSelector = std::to_string(era * 10);
			} else {
				this->preVFP = false; // aka postVFP, no-HIPM, no-APV which is just default track construction
				eraSelector = std::to_string(era);
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
	// Set Branches of output tree
	for (const std::shared_ptr<TTree> &tree : outputTrees) {
		tree->Branch("nMuon", &nMuon);
		tree->Branch("nGoodMuon", &nGoodMuon);
		tree->Branch("nVetoMuon", &nVetoMuon);
		tree->Branch("nAntiSelectedMuon", &nAntiSelectedMuon);

		tree->Branch("MuonPt", muonPt.data(), "MuonPt[nMuon]/D");
		tree->Branch("MuonRoccorPt", muonRoccorPt.data(), "MuonRoccorPt[nMuon]/D");
		tree->Branch("MuonRoccorPtUp", muonRoccorPtUp.data(), "MuonRoccorPtUp[nMuon]/D");
		tree->Branch("MuonRoccorPtDown", muonRoccorPtDown.data(), "MuonRoccorPtDown[nMuon]/D");
		tree->Branch("MuonPtVector", &muonPtVector);
		tree->Branch("MuonEta", muonEta.data(), "MuonEta[nMuon]/D");
		tree->Branch("MuonPhi", muonPhi.data(), "MuonPhi[nMuon]/D");
		tree->Branch("MuonMass", muonMass.data(), "MuonMass[nMuon]/D");
		tree->Branch("MuonMiniIso", muonMiniIso.data(), "MuonMiniIso[nMuon]/D");

		tree->Branch("MuonLooseId", muonLooseId.data(), "MuonLooseId[nMuon]/B");
		tree->Branch("MuonMediumId", muonMediumId.data(), "MuonMediumId[nMuon]/B");
		tree->Branch("MuonTightId", muonTightId.data(), "MuonTightId[nMuon]/B");
		tree->Branch("MuonIsGood", muonIsGood.data(), "MuonIsGood[nMuon]/B");
		tree->Branch("MuonIsVeto", muonIsVeto.data(), "MuonIsVeto[nMuon]/B");
		tree->Branch("MuonIsAntiSelected", muonIsAntiSelected.data(), "MuonIsAntiSelected[nMuon]/B");

		tree->Branch("MuonCharge", muonCharge.data(), "MuonCharge[nMuon]/I");
		tree->Branch("MuonPdgId", muonPdgId.data(), "MuonPdgId[nMuon]/I");
		tree->Branch("MuonGenMatchedIndex", muonMass.data(), "MuonGenMatchedIndex[nMuon]/I");
		tree->Branch("MuonCutBasedId", muonCutBasedId.data(), "MuonCutBasedId[nMuon]/I");

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

		tree->Branch("ElectronLooseMvaId", electronLooseMvaId.data(), "ElectronLooseMvaId[nMuon]/B");
		tree->Branch("ElectronMediumMvaId", electronMediumMvaId.data(), "ElectronMediumMvaId[nMuon]/B");
		tree->Branch("ElectronTightMvaId", electronTightMvaId.data(), "ElectronTightMvaId[nMuon]/B");

		tree->Branch("ElectronCharge", electronCharge.data(), "ElectronCharge[nMuon]/I");
		tree->Branch("ElectronCutBasedId", electronCutBasedId.data(), "ElectronCutBasedId[nMuon]/I");
		//tree->Branch("ElectronNLostHits", electronNLostHits.data(), "ElectronNLostHits[nMuon]/I");

		tree->Branch("ElectronTightId", electronTightId.data(), "ElectronTightId[nMuon]/B");
		tree->Branch("ElectronMediumId", electronMediumId.data(), "ElectronMediumId[nMuon]/B");
		tree->Branch("ElectronLooseId", electronLooseId.data(), "ElectronLooseId[nMuon]/B");
		tree->Branch("ElectronVetoId", electronVetoId.data(), "ElectronVetoId[nMuon]/B");
		//tree->Branch("ElectronConvVeto", electronConvVeto.data(), "ElectronConvVeto[nMuon]/B");
	}
}
