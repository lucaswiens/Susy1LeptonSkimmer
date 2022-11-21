#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/ScaleFactorProducer.h>

ScaleFactorProducer::ScaleFactorProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector, TFile &outputFile) {
	Name = "ScaleFactorProducer";
	//TH1::AddDirectory(false);
	std::string cmsswBase = std::getenv("CMSSW_BASE");

	//Read information needed
	//run = configTree.get<std::string>("run");
	//era = configTree.get<std::string>("era");
	//isData = run != "MC";

	electronSf = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Electron." + eraSelector + ".JSON"));
	muonSf = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Muon." + eraSelector + ".ScaleFactor.JSON"));
	bTagSf = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.BTag." + eraSelector));

	electronEraAlias = scaleFactorTree.get<std::string>("Electron." + eraSelector + ".EraAlias");
	muonEraAlias = scaleFactorTree.get<std::string>("Muon." + eraSelector + ".ScaleFactor.EraAlias");
	muonTriggName = scaleFactorTree.get<std::string>("Muon." + eraSelector + ".ScaleFactor.TriggerName");

	bTagSyst = Utility::GetVector<std::string>(configTree, "Producer.Jet.BTagSystematic"); bTagSystLight = Utility::GetVector<std::string>(configTree, "Producer.Jet.LightTagSystematic");
	////Set histograms
	float ptCut = configTree.get<float>("Producer.Jet.Pt");
	float etaCut = configTree.get<float>("Producer.Jet.Eta");

	std::vector<float> etaBins = {-etaCut, -1.4, 1.4, etaCut};
	std::vector<float> ptBins;

	ptBins = {ptCut, 50, 70, 90, 200};

	bTotal = std::make_shared<TH2D>("nTrueB", "TotalB", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
	cTotal = std::make_shared<TH2D>("nTrueC", "TotalC", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
	lightTotal = std::make_shared<TH2D>("nTrueLight", "TotalLight", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

	bTagEffBLooseDeepJet = std::make_shared<TH2D>("nLooseBbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
	bTagEffBMediumDeepJet = std::make_shared<TH2D>("nMediumBbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
	bTagEffBTightDeepJet = std::make_shared<TH2D>("nTightBbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

	bTagEffCLooseDeepJet = std::make_shared<TH2D>("nLooseCbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
	bTagEffCMediumDeepJet = std::make_shared<TH2D>("nMediumCbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
	bTagEffCTightDeepJet = std::make_shared<TH2D>("nTightCbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

	bTagEffLightLooseDeepJet = std::make_shared<TH2D>("nLooseLightbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
	bTagEffLightMediumDeepJet = std::make_shared<TH2D>("nMediumLightbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
	bTagEffLightTightDeepJet = std::make_shared<TH2D>("nTightLightbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
}

void ScaleFactorProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	float electronPt, muonEta, muonPt;

	//Loop over selected electrons
	for (int iElectron = 0; iElectron < product.nElectron; iElectron++) {
		//electronPt = product.electronPt[iElectron] >= 30 ? product.electronPt[iElectron] : 30;
		electronPt = product.electronPt[iElectron];

		std::string recoPtThreshold            = product.electronPt[iElectron] >= 20 ? "RecoAbove20" : "RecoBelow20";
		product.electronRecoSf[iElectron]      = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", recoPtThreshold, product.electronEta[iElectron], electronPt});
		product.electronVetoSf[iElectron]      = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "Veto", product.electronEta[iElectron], electronPt});
		product.electronMediumSf[iElectron]    = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "Medium", product.electronEta[iElectron], electronPt});
		product.electronTightSf[iElectron]     = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "Tight", product.electronEta[iElectron], electronPt});
		product.electronMediumMvaSf[iElectron] = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "wp90iso", product.electronEta[iElectron], electronPt});
		product.electronTightMvaSf[iElectron]  = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sf", "wp80iso", product.electronEta[iElectron], electronPt});

		product.electronRecoSfUp[iElectron]      = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", recoPtThreshold, product.electronEta[iElectron], electronPt});
		product.electronVetoSfUp[iElectron]      = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "Veto", product.electronEta[iElectron], electronPt});
		product.electronMediumSfUp[iElectron]    = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "Medium", product.electronEta[iElectron], electronPt});
		product.electronTightSfUp[iElectron]     = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "Tight", product.electronEta[iElectron], electronPt});
		product.electronMediumMvaSfUp[iElectron] = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "wp90iso", product.electronEta[iElectron], electronPt});
		product.electronTightMvaSfUp[iElectron]  = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfup", "wp80iso", product.electronEta[iElectron], electronPt});

		product.electronRecoSfDown[iElectron]      = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfdown", recoPtThreshold, product.electronEta[iElectron], electronPt});
		product.electronVetoSfDown[iElectron]      = electronSf->at("UL-Electron-ID-SF")->evaluate({electronEraAlias, "sfdown", "Veto", product.electronEta[iElectron], electronPt});
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
		product.muonMediumIsoSf[iMuon] = muonSf->at("NUM_LooseRelIso_DEN_MediumID")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonLooseSf[iMuon]    = muonSf->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonMediumSf[iMuon]   = muonSf->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonTightSf[iMuon]    = muonSf->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "sf"});
		product.muonTriggerSf[iMuon]  = muonSf->at(muonTriggName)->evaluate({muonEraAlias, muonEta, muonPt, "sf"});

		product.muonLooseIsoSfUp[iMuon] = muonSf->at("NUM_LooseRelIso_DEN_LooseID")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonMediumIsoSfUp[iMuon] = muonSf->at("NUM_LooseRelIso_DEN_MediumID")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonLooseSfUp[iMuon]    = muonSf->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonMediumSfUp[iMuon]   = muonSf->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonTightSfUp[iMuon]    = muonSf->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systup"});
		product.muonTriggerSfUp[iMuon]  = muonSf->at(muonTriggName)->evaluate({muonEraAlias, muonEta, muonPt, "systup"});

		product.muonLooseIsoSfDown[iMuon] = muonSf->at("NUM_LooseRelIso_DEN_LooseID")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonMediumIsoSfDown[iMuon] = muonSf->at("NUM_LooseRelIso_DEN_MediumID")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonLooseSfDown[iMuon]    = muonSf->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonMediumSfDown[iMuon]   = muonSf->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonTightSfDown[iMuon]    = muonSf->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
		product.muonTriggerSfDown[iMuon]  = muonSf->at(muonTriggName)->evaluate({muonEraAlias, muonEta, muonPt, "systdown"});
	}

	/*######################################################################################################
	#   Btag efficiency and Sf                                                                             #
	#   https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1a_Event_reweighting_using_scale        #
	#   Due to large file sizes, only store the systematic variations of the recommended DeepJet bTagger   #
	######################################################################################################*/
	for (int iJet = 0; iJet < product.nJet; iJet++) {
		// btag Sf
		int flav = (std::abs(product.jetPartFlav[iJet]) == 4 || std::abs(product.jetPartFlav[iJet]) == 5) ? std::abs(product.jetPartFlav[iJet]) : 0;
		std::string postFix = flav != 0 ? "_comb" : "_incl";

		product.jetDeepJetLooseSf[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({"central", "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		product.jetDeepJetMediumSf[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({"central", "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		product.jetDeepJetTightSf[iJet] = bTagSf->at("deepJet" + postFix)->evaluate({"central", "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});

		for (int bSyst = 0; bSyst < bTagSyst.size(); ++bSyst) {
			std::string shiftUp = flav != 0 ? "up_" + bTagSyst[bSyst] : "central",
				shiftDown = flav != 0 ? "down_" + bTagSyst[bSyst] : "central";
			product.jetDeepJetLooseSfUp[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetLooseSfDown[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetMediumSfUp[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetMediumSfDown[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetTightSfUp[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetTightSfDown[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		}

		for (int bSyst = 0; bSyst < bTagSystLight.size(); ++bSyst) {
			std::string shiftUp = flav == 0 ? "up_" + bTagSyst[bSyst] : "central",
				shiftDown = flav == 0 ? "down_" + bTagSyst[bSyst] : "central";
			product.jetDeepJetLooseLightSfUp[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetLooseLightSfDown[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "L", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetMediumLightSfUp[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetMediumLightSfDown[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "M", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetTightLightSfUp[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftUp, "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
			product.jetDeepJetTightLightSfDown[bSyst][iJet] = bTagSf->at("deepJet" + postFix)->evaluate({shiftDown, "T", flav, std::abs(product.jetEta[iJet]), product.jetPt[iJet]});
		}

		if (std::abs(product.jetPartFlav[iJet]) == 5) {
			bTotal->Fill(product.jetPt[iJet], product.jetEta[iJet]);

			if (product.jetDeepJetId[iJet] >= 3) {
				bTagEffBLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffBMediumDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffBTightDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			} else if (product.jetDeepJetId[iJet] >= 2) {
				bTagEffBLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffBMediumDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			} else if (product.jetDeepJetId[iJet] >= 1) {
				bTagEffBLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			}
		} else if (std::abs(product.jetPartFlav[iJet]) == 4) {
			cTotal->Fill(product.jetPt[iJet], product.jetEta[iJet]);

			if (product.jetDeepJetId[iJet] >= 3) {
				bTagEffCLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffCMediumDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffCTightDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			} else if (product.jetDeepJetId[iJet] >= 2) {
				bTagEffCLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffCMediumDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			} else if (product.jetDeepJetId[iJet] >= 1) {
				bTagEffCLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			}
		} else {
			lightTotal->Fill(product.jetPt[iJet], product.jetEta[iJet]);

			if (product.jetDeepJetId[iJet] >= 3) {
				bTagEffLightLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffLightMediumDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffLightTightDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			} else if (product.jetDeepJetId[iJet] >= 2) {
				bTagEffLightLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
				bTagEffLightMediumDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			} else if (product.jetDeepJetId[iJet] >= 1) {
				bTagEffLightLooseDeepJet->Fill(product.jetPt[iJet], product.jetEta[iJet]);
			}
		}
	}
}

void ScaleFactorProducer::EndJob(TFile &file) {
	bTotal->Write();
	cTotal->Write();
	lightTotal->Write();

	bTagEffBLooseDeepJet->Write();
	bTagEffBMediumDeepJet->Write();
	bTagEffBTightDeepJet->Write();

	bTagEffCLooseDeepJet->Write();
	bTagEffCMediumDeepJet->Write();
	bTagEffCTightDeepJet->Write();

	bTagEffLightLooseDeepJet->Write();
	bTagEffLightMediumDeepJet->Write();
	bTagEffLightTightDeepJet->Write();
}
