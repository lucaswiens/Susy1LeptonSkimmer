
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/ScaleFactorProducer.h>

ScaleFactorProducer::ScaleFactorProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {
	//TH1::AddDirectory(false);
	std::string cmsswBase = std::getenv("CMSSW_BASE");


	electronSf = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Electron." + eraSelector + ".JSON"));
	muonSf = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Muon." + eraSelector + ".ScaleFactor.JSON"));
	bTagSf = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.BTag." + eraSelector));

	electronEraAlias = scaleFactorTree.get<std::string>("Electron." + eraSelector + ".eraAlias");
	muonEraAlias = scaleFactorTree.get<std::string>("Muon." + eraSelector + ".ScaleFactor.eraAlias");
	muonTriggName = scaleFactorTree.get<std::string>("Muon." + eraSelector + ".ScaleFactor.triggerName");

	bTagSyst = Utility::GetVector<std::string>(configTree, "Producer.Jet.BTagSystematics");
}

void ScaleFactorProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	double electronPt, muonEta, muonPt;

	//Loop over selected electrons
	for (int iElectron = 0; iElectron < product.nElectron; iElectron++) {
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
	for (int iMuon = 0; iMuon < product.nMuon; iMuon++) {
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

	//Btag efficiency and Sf
	for (int iJet = 0; iJet < product.nJet; iJet++) {
		int flav = (std::abs(product.jetPartFlav[iJet]) == 4 || std::abs(product.jetPartFlav[iJet]) == 5) ? std::abs(product.jetPartFlav[iJet]) : 0;
		std::string postFix = flav != 0 ? "_comb" : "_incl";

		product.jetDeepCSVLooseSf[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({"central", "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		product.jetDeepCSVMediumSf[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({"central", "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		product.jetDeepCSVTightSf[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({"central", "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});

		product.jetDeepJetLooseSf[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({"central", "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		product.jetDeepJetMediumSf[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({"central", "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		product.jetDeepJetTightSf[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({"central", "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});

		for (int bSyst = 0; bSyst < bTagSyst.size(); ++bSyst) {
			std::string shiftUp = flav != 0 ? "up_" + bTagSyst.at(bSyst) : "central",
				shiftDown = flav != 0 ? "down_" + bTagSyst.at(bSyst) : "central";

			product.jetDeepCSVLooseSfUp.at(bSyst)[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({shiftUp, "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepCSVLooseSfDown.at(bSyst)[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({shiftDown, "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepCSVMediumSfUp.at(bSyst)[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({shiftUp, "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepCSVMediumSfDown.at(bSyst)[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({shiftDown, "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepCSVTightSfUp.at(bSyst)[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({shiftUp, "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepCSVTightSfDown.at(bSyst)[iJet] = bTagSf->at("deepCSV" + postFix)->evaluate({shiftDown, "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});

			product.jetDeepJetLooseSfUp.at(bSyst)[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetLooseSfDown.at(bSyst)[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetMediumSfUp.at(bSyst)[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetMediumSfDown.at(bSyst)[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetTightSfUp.at(bSyst)[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetTightSfDown.at(bSyst)[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		}
	}
}

void ScaleFactorProducer::EndJob(TFile &) {}
