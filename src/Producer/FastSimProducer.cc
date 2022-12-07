#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/FastSimProducer.h>

FastSimProducer::FastSimProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector, TFile &outputFile) {
	Name = "FastSimProducer";
	std::string cmsswBase = std::getenv("CMSSW_BASE");

	stopPdgId       = configTree.get<int>("Producer.FastSim.StopPdgId");
	gluinoPdgId     = configTree.get<int>("Producer.FastSim.GluinoPdgId");
	neutralinoPdgId = configTree.get<int>("Producer.FastSim.NeutralinoPdgId");
	charginoPdgId   = configTree.get<int>("Producer.FastSim.CharginoPdgId");

	// Open json file with cross section
	pt::read_json(std::string(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/fastsim/" + configTree.get<std::string>("Producer.FastSim.CrossSection")).c_str(), xSectionTree);

	// Open root files stored in root files
	electronSfFile = (TFile*) TFile::Open((cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/fastsim/" + eraSelector + "/leptonScaleFactor/" + scaleFactorTree.get<std::string>("Electron." + eraSelector + ".FastSim.Path")).c_str());
	muonSfFile     = (TFile*) TFile::Open((cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/fastsim/" + eraSelector + "/leptonScaleFactor/" + scaleFactorTree.get<std::string>("Muon." + eraSelector + ".ScaleFactor.FastSim.Path")).c_str());

	/*************************************************************************************************
	*   Electron Fastsim Scale Factors                                                               *
	*   2016 https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Data_leading_order_FullSim_MC_co   *
	*   2017 https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Electrons_FullSim_FastSim_AN1      *
	*   2018 https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_compari   *
	*************************************************************************************************/
	electronVetoSf     = std::make_shared<TH2F>(*(TH2F*)electronSfFile->Get("CutBasedVetoNoIso94XV2_sf"));
	electronTightSf    = std::make_shared<TH2F>(*(TH2F*)electronSfFile->Get("CutBasedTightNoIso94XV2_sf"));
	electronVetoMvaSf = std::make_shared<TH2F>(*(TH2F*)electronSfFile->Get("MVAVLooseTightIP2DMini4_sf")); // not sure if this is correct
	electronTightMvaSf = std::make_shared<TH2F>(*(TH2F*)electronSfFile->Get("ConvIHit0_sf"));

	/*************************************************************************************************
	*   Muon Scale Factors                                                                           *
	*   2016 https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_com_AN1   *
	*   2017 https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Muons_FullSim_FastSim_AN1          *
	*   2018 https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Muons_FullSim_FastSim              *
	*************************************************************************************************/
	muonLooseSf      = std::make_shared<TH2F>(*(TH2F*)muonSfFile->Get("miniIso04_LooseId_sf"));
	muonMediumSf      = std::make_shared<TH2F>(*(TH2F*)muonSfFile->Get("miniIso02_MediumId_sf"));

	/***********************************************************************************************************
	*   bTag scale factor csv files                                                                            *
	*   https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy#FastSim_correction_factors   *
	*   https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X#FastSim_correction_factors          *
	*   https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X#FastSim_correction_factors         *
	/***********************************************************************************************************
	***********************************************************************************************************/
	btagCsvReader = std::make_unique<BTagCSVReader>(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/fastsim/" + eraSelector + "/bTagScaleFactor/" + scaleFactorTree.get<std::string>("Jet.BTag.FastSim.BTag." + eraSelector), true);
}

void FastSimProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	// Lepton Scale Factors between FullSim and FastSim
	for (int iMuon = 0; iMuon < product.nMuon; iMuon++) {
		const float &muonPt  = product.muonPt[iMuon];
		const float &muonEta = std::abs(product.muonEta[iMuon]);
		product.muonLooseFastSf[iMuon]  = Utility::Get2DWeight(muonEta, muonPt, static_cast<TH2F*>(muonLooseSf.get())); // function needs raw pointer, shared_ptr()->get() returns it
		product.muonMediumFastSf[iMuon] = Utility::Get2DWeight(muonEta, muonPt, static_cast<TH2F*>(muonMediumSf.get()));

		const float &looseFastSfErr  = Utility::Get2DWeightErr(muonEta, muonPt, static_cast<TH2F*>(muonLooseSf.get()));
		const float &mediumFastSfErr = Utility::Get2DWeightErr(muonEta, muonPt, static_cast<TH2F*>(muonMediumSf.get()));

		product.muonLooseFastSfUp[iMuon]  = product.muonLooseFastSf[iMuon]  + looseFastSfErr;
		product.muonMediumFastSfUp[iMuon] = product.muonMediumFastSf[iMuon] + mediumFastSfErr;

		product.muonLooseFastSfDown[iMuon]  = product.muonLooseFastSf[iMuon]  - looseFastSfErr;
		product.muonMediumFastSfDown[iMuon] = product.muonMediumFastSf[iMuon] - mediumFastSfErr;
	}

	for (int iElectron = 0; iElectron < product.nElectron; iElectron++) {
		const float &electronPt  = product.electronPt[iElectron];
		const float &electronEta = std::abs(product.electronEta[iElectron]);
		product.electronVetoFastSf[iElectron]     = Utility::Get2DWeight(electronEta, electronPt, static_cast<TH2F*>(electronVetoSf.get()));
		product.electronTightFastSf[iElectron]    = Utility::Get2DWeight(electronEta, electronPt, static_cast<TH2F*>(electronTightSf.get()));
		product.electronVetoMvaFastSf[iElectron]  = Utility::Get2DWeight(electronEta, electronPt, static_cast<TH2F*>(electronVetoMvaSf.get()));
		product.electronTightMvaFastSf[iElectron] = Utility::Get2DWeight(electronEta, electronPt, static_cast<TH2F*>(electronTightMvaSf.get()));

		const float &vetoFastSfErr     = Utility::Get2DWeightErr(electronEta, electronPt, static_cast<TH2F*>(electronVetoSf.get()));
		const float &tightFastSfErr    = Utility::Get2DWeightErr(electronEta, electronPt, static_cast<TH2F*>(electronTightSf.get()));
		const float &vetoMvaFastSfErr  = Utility::Get2DWeightErr(electronEta, electronPt, static_cast<TH2F*>(electronVetoMvaSf.get()));
		const float &tightMvaFastSfErr = Utility::Get2DWeightErr(electronEta, electronPt, static_cast<TH2F*>(electronTightMvaSf.get()));

		product.electronVetoFastSfUp[iElectron]    = product.electronVetoFastSf[iElectron]     + vetoFastSfErr;
		product.electronTightFastSfUp[iElectron]   = product.electronTightFastSf[iElectron]    + tightFastSfErr;
		product.electronVetoMvaFastSfUp[iElectron] = product.electronVetoMvaFastSf[iElectron]  + vetoMvaFastSfErr;
		product.electronTightMvaFastSfUp[iElectron]= product.electronTightMvaFastSf[iElectron] + tightMvaFastSfErr;

		product.electronVetoFastSfDown[iElectron]    = product.electronVetoFastSf[iElectron]     - vetoFastSfErr;
		product.electronTightFastSfDown[iElectron]   = product.electronTightFastSf[iElectron]    - tightFastSfErr;
		product.electronVetoMvaFastSfDown[iElectron] = product.electronVetoMvaFastSf[iElectron]  - vetoMvaFastSfErr;
		product.electronTightMvaFastSfDown[iElectron]= product.electronTightMvaFastSf[iElectron] - tightMvaFastSfErr;
	}

	for (int iJet = 0; iJet < product.nJet; iJet++) {
		product.jetDeepJetLooseFastSf[iJet]      = btagCsvReader->Get(product.jetPt.at(iJet), 'L');
		product.jetDeepJetMediumFastSf[iJet]     = btagCsvReader->Get(product.jetPt.at(iJet), 'M');
		product.jetDeepJetTightFastSf[iJet]      = btagCsvReader->Get(product.jetPt.at(iJet), 'T');

		product.jetDeepJetLooseFastSfUp[iJet]    = btagCsvReader->GetUp(product.jetPt.at(iJet), 'L');
		product.jetDeepJetMediumFastSfUp[iJet]   = btagCsvReader->GetUp(product.jetPt.at(iJet), 'M');
		product.jetDeepJetTightFastSfUp[iJet]    = btagCsvReader->GetUp(product.jetPt.at(iJet), 'T');

		product.jetDeepJetLooseFastSfDown[iJet]  = btagCsvReader->GetDown(product.jetPt.at(iJet), 'L');
		product.jetDeepJetMediumFastSfDown[iJet] = btagCsvReader->GetDown(product.jetPt.at(iJet), 'M');
		product.jetDeepJetTightFastSfDown[iJet]  = btagCsvReader->GetDown(product.jetPt.at(iJet), 'T');
	}

	dataReader.ReadGenEntry();
	product.stopMass = -999;
	product.gluinoMass = -999;
	product.neutralinoMass = -999;
	product.charginoMass = -999;

	product.susyXSectionNLO    = -999;
	product.susyXSectionNLLO   = -999;
	product.susyXSectionNLOUp  = -999;
	product.susyXSectionNLLOUp = -999;
	product.susyXSectionNLODown  = -999;
	product.susyXSectionNLLODown = -999;

	for (int iGen = 0; iGen < dataReader.nGenPart; iGen++) {
		dataReader.GetGenValues(iGen);
		if (std::abs(dataReader.genPdgId) == stopPdgId) {
			product.stopMass = dataReader.genMass;
		} else if (std::abs(dataReader.genPdgId) == gluinoPdgId) {
			product.gluinoMass = dataReader.genMass;
		} else if (std::abs(dataReader.genPdgId) == neutralinoPdgId) {
			product.neutralinoMass = dataReader.genMass;
		} else if (std::abs(dataReader.genPdgId) == charginoPdgId) {
			product.charginoMass = dataReader.genMass;
		}

		// Check if the mass exists in the cross section json file... Sometimes the masses are set to strange unexpected values, one might needs to correct this later or discard those events
		boost::optional< pt::ptree& > child = xSectionTree.get_child_optional("NLO+NLL Cross section [pb]." + std::to_string((int)product.gluinoMass));
		if (!child) { continue;}

		product.susyXSectionNLO                  = xSectionTree.get<float>("NLO+NLL Cross section [pb]." + std::to_string((int)product.gluinoMass));
		const float &susyXSectionNLOUncertainty  = xSectionTree.get<float>("Uncertainty (NLO + NLL) [%]." + std::to_string((int)product.gluinoMass));
		product.susyXSectionNLLO                 = xSectionTree.get<float>("NNLO+NNLL Cross section [pb]." + std::to_string((int)product.gluinoMass));
		const float &susyXSectionNLLOUncertainty = xSectionTree.get<float>("Uncertainty (NNLOapprox + NNLL) [%]." + std::to_string((int)product.gluinoMass));

		product.susyXSectionNLOUp  = product.susyXSectionNLO  + susyXSectionNLOUncertainty;
		product.susyXSectionNLLOUp = product.susyXSectionNLLO + susyXSectionNLLOUncertainty;

		product.susyXSectionNLODown  = product.susyXSectionNLO  - susyXSectionNLOUncertainty;
		product.susyXSectionNLLODown = product.susyXSectionNLLO - susyXSectionNLLOUncertainty;
	}
}

void FastSimProducer::EndJob(TFile &file) {
	electronSfFile->Close();
	muonSfFile->Close();
}
