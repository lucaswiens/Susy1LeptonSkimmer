
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/ScaleFactorProducer.h>

ScaleFactorProducer::ScaleFactorProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {
	TH1::AddDirectory(false);
	std::string cmsswBase = std::getenv("CMSSW_BASE");

	std::cout << "electronSf: " << cmsswBase + "/src" + scaleFactorTree.get<std::string>("Electron." + eraSelector + ".JSON") << std::endl;
	std::cout << "muonSf:     " << cmsswBase + "/src" + scaleFactorTree.get<std::string>("Muon." + eraSelector + "ScaleFactor.JSON") << std::endl;
	std::cout << "bTagSf:     " << cmsswBase + "/src" + scaleFactorTree.get<std::string>("Jet.BTag." + eraSelector) << std::endl;

	electronSf = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Electron." + eraSelector + ".JSON"));
	muonSf = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Muon." + eraSelector + "ScaleFactor.JSON"));

	electronEraAlias = scaleFactorTree.get<std::string>("Electron." + eraSelector + ".eraAlias");
	muonEraAlias = scaleFactorTree.get<std::string>("Muon." + eraSelector + "ScaleFactor.eraAlias");
	muonTriggName = scaleFactorTree.get<std::string>("Muon." + eraSelector + "ScaleFactor.triggerName");

}

void ScaleFactorProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	float electronPt, muonEta, muonPt;

	//Loop over selected electrons
	for(int iElectron = 0; iElectron < product.nElectron; iElectron++){
		//electronPt = product.electronPt[iElectron] >= 30 ? product.electronPt[iElectron] : 30;
		electronPt = product.electronPt[iElectron];

		std::string recoPtThreshold = product.electronPt[iElectron] >= 20 ? "RecoAbove20" : "RecoBelow20";
		product.electronRecoSf[iElectron]      = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", recoPtThreshold, product.electronEta[iElectron], electronPt});
		product.electronVetoSf[iElectron]     = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "Veto", product.electronEta[iElectron], electronPt});
		product.electronMediumSf[iElectron]    = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "Medium", product.electronEta[iElectron], electronPt});
		product.electronTightSf[iElectron]     = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "Tight", product.electronEta[iElectron], electronPt});
		product.electronMediumMvaSf[iElectron] = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "wp90iso", product.electronEta[iElectron], electronPt});
		product.electronTightMvaSf[iElectron]  = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "wp80iso", product.electronEta[iElectron], electronPt});

		product.electronRecoSfUp[iElectron]   = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", recoPtThreshold, product.electronEta[iElectron], electronPt});
		product.electronVetoSfUp[iElectron]     = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "Veto", product.electronEta[iElectron], electronPt});
		product.electronMediumSfUp[iElectron]    = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "Medium", product.electronEta[iElectron], electronPt});
		product.electronTightSfUp[iElectron]     = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "Tight", product.electronEta[iElectron], electronPt});
		product.electronMediumMvaSfUp[iElectron] = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "wp90iso", product.electronEta[iElectron], electronPt});
		product.electronTightMvaSfUp[iElectron]  = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "wp80iso", product.electronEta[iElectron], electronPt});

		product.electronRecoSfDown[iElectron]      = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfdown", recoPtThreshold, product.electronEta[iElectron], electronPt});
		product.electronVetoSfDown[iElectron]     = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfdown", "Veto", product.electronEta[iElectron], electronPt});
		product.electronMediumSfDown[iElectron]    = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfdown", "Medium", product.electronEta[iElectron], electronPt});
		product.electronTightSfDown[iElectron]     = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfdown", "Tight", product.electronEta[iElectron], electronPt});
		product.electronMediumMvaSfDown[iElectron] = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfdown", "wp90iso", product.electronEta[iElectron], electronPt});
		product.electronTightMvaSfDown[iElectron]  = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "wp80iso", product.electronEta[iElectron], electronPt});
	}

	//Loop over selected muons
	for(int iMuon = 0; iMuon < product.nMuon; iMuon++){
		muonPt  = product.muonPt[iMuon] >= 30 ? product.muonPt[iMuon] : 30;
		muonEta = std::abs(product.muonEta[iMuon]);

		//Scale factors
		product.muonLooseIsoSf[iMuon] = muonSf->at("NUM_LooseRelIso_DEN_LooseID")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonTightIsoSf[iMuon] = muonSf->at("NUM_TightRelIso_DEN_TightIDandIPCut")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonLooseSf[iMuon]    = muonSf->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonMediumSf[iMuon]   = muonSf->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonTightSf[iMuon]    = muonSf->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonTriggerSf[iMuon]  = muonSf->at(muonTriggName)->evaluate({muonEraAlias, muonEta, muonPt, "sf"});

		product.muonLooseIsoSfUp[iMuon] = muonSf->at("NUM_LooseRelIso_DEN_LooseID")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonTightIsoSfUp[iMuon] = muonSf->at("NUM_TightRelIso_DEN_TightIDandIPCut")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonLooseSfUp[iMuon]    = muonSf->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonMediumSfUp[iMuon]   = muonSf->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonTightSfUp[iMuon]    = muonSf->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonTriggerSfUp[iMuon]  = muonSf->at(muonTriggName)->evaluate({muonEraAlias, muonEta, muonPt, "systup"});

		product.muonLooseIsoSfDown[iMuon] = muonSf->at("NUM_LooseRelIso_DEN_LooseID")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonTightIsoSfDown[iMuon] = muonSf->at("NUM_TightRelIso_DEN_TightIDandIPCut")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonLooseSfDown[iMuon]    = muonSf->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonMediumSfDown[iMuon]   = muonSf->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonTightSfDown[iMuon]    = muonSf->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonTriggerSfDown[iMuon]  = muonSf->at(muonTriggName)->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
	}
}

void ScaleFactorProducer::EndJob(TFile &) {}
