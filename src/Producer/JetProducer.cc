#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/JetProducer.h>

JetProducer::JetProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product) {
	Name = "JetProducer";
	std::string cmsswBase = std::getenv("CMSSW_BASE");

	/*################################################################################################
	#   https://github.com/cms-nanoAOD/correctionlib                                                 #
	#   https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master                         #
	#   Instructions:                                                                                #
	#   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite   #
	################################################################################################*/
	era = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".Era");
	if (product.GetIsData()) { // Run Period and JEC Version for Data
		dataType = "Run" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".RunPeriod." + product.GetRunPeriod()) + "_" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".DATA");
	} else { // JEC Version for MC
		dataType  = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".MC");
	}
	jerVersion       = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".JER");
	ak4Algorithm     = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK4.Algorithm");
	ak8Algorithm     = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK8.Algorithm");
	ak4CorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK4.JSON"));
	ak8CorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK8.JSON"));
	jmeCorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".JME.JSON"));

	/**************************************************************
	*   Fastsim JEC corrections                                   *
	*   This is not available in the nice and easy json format.   *
	*   So just use the old JEC classes                           *
	**************************************************************/
	if (product.GetIsFastSim()) {
		std::vector<JetCorrectorParameters> ak4Vec, ak8Vec;
		for (std::string fileName: Utility::GetVector<std::string>(scaleFactorTree, "Jet.JERC." + product.GetEraSelector() + ".FastSim")) {
			ak4Vec.push_back(JetCorrectorParameters(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/JERC/" + product.GetEraSelector() + "/" + fileName + ak4Algorithm + ".txt"));
			ak8Vec.push_back(JetCorrectorParameters(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/JERC/" + product.GetEraSelector() + "/" + fileName + ak8Algorithm + ".txt"));
		}
		fastJetCorrectorAk4 = std::make_shared<FactorizedJetCorrector>(ak4Vec);
		fastJetCorrectorAk8 = std::make_shared<FactorizedJetCorrector>(ak8Vec);
	}


	// Set object to get JEC uncertainty
	if (product.GetIsData()) {
		jecSystematics = std::vector<std::string>{};
	} else if (product.GetIsFastSim()) {
		jecSystematics = {"Total"}; // hardcoded for now because only one kind of uncertainty exists for this.
		//!product.GetIsData() ? Utility::GetVector<std::string>(configTree, "Producer.Jet.JECSystematic") : std::vector<std::string>{};
	} else {
		jecSystematics = Utility::GetVector<std::string>(configTree, "Producer.Jet.JECSystematic");
	}

	const std::string &jecAk4SystematicPath = cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/JERC/" + product.GetEraSelector() + "/" + scaleFactorTree.get<std::string>("Jet.JECSystematic." + (product.GetIsFastSim() ? "FastSim." : std::string()) + product.GetEraSelector()) + ak4Algorithm + ".txt";
	const std::string &jecAk8SystematicPath = cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/JERC/" + product.GetEraSelector() + "/" + scaleFactorTree.get<std::string>("Jet.JECSystematic." + (product.GetIsFastSim() ? "FastSim." : std::string()) + product.GetEraSelector()) + ak8Algorithm + ".txt";

	//fastAk4CorrectionUncertainty, fastAk8CorrectionUncertainty;
	if (product.GetIsFastSim()) {
		ak4CorrectionUncertainty.insert({"Total", std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(jecAk4SystematicPath, ""))});
		ak8CorrectionUncertainty.insert({"Total", std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(jecAk8SystematicPath, ""))});
	} else {
		for (std::string systematic : jecSystematics) {
			ak4CorrectionUncertainty.insert({systematic, std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(jecAk4SystematicPath, systematic))});
			ak8CorrectionUncertainty.insert({systematic, std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(jecAk8SystematicPath, systematic))});
		}
	}

	// Smearing for fat jets not available in jsonpog -> use old school smearing TODO update when jsonpog makes them available
	if (product.GetEra() == 2017) { // TODO this is a bad hardcoded fix but good enough for now
		resolutionAk8 = JME::JetResolution(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/JERC/" + product.GetEraSelector() + "/" + era + "_JRV3_MC_PtResolution_" + ak8Algorithm + ".txt");
		resolutionSfAk8 = JME::JetResolutionScaleFactor(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/JERC/" + product.GetEraSelector() + "/" + era + "_JRV3_MC_SF_" + ak8Algorithm + ".txt");
	} else {
		std::cout << era + "_" + jerVersion + "_PtResolution_" + ak8Algorithm << std::endl;
		resolutionAk8 = JME::JetResolution(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/JERC/" + product.GetEraSelector() + "/" + era + "_" + jerVersion + "_PtResolution_" + ak8Algorithm + ".txt");
		resolutionSfAk8 = JME::JetResolutionScaleFactor(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/JERC/" + product.GetEraSelector() + "/" + era + "_" + jerVersion + "_SF_" + ak8Algorithm + ".txt");
	}

	/*################################################################################
	#   Definition of Working Points come from                                       #
	#   https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP    #
	#   https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP   #
	#   https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL17      #
	#   https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP   #
	################################################################################*/
	deepCsvBTagMap = {
		{'L', configTree.get<float>("Producer.Jet.DeepCSV." + product.GetEraSelector() + ".Loose")},
		{'M', configTree.get<float>("Producer.Jet.DeepCSV." + product.GetEraSelector() + ".Medium")},
		{'T', configTree.get<float>("Producer.Jet.DeepCSV." + product.GetEraSelector() + ".Tight")},
	};

	deepJetBTagMap = {
		{'L', configTree.get<float>("Producer.Jet.DeepJet." + product.GetEraSelector() + ".Loose")},
		{'M', configTree.get<float>("Producer.Jet.DeepJet." + product.GetEraSelector() + ".Medium")},
		{'T', configTree.get<float>("Producer.Jet.DeepJet." + product.GetEraSelector() + ".Tight")},
	};

	deepAk8TopTagMap = {
		{'L', configTree.get<float>("Producer.Jet.DeepAK8.Top." + product.GetEraSelector() + ".Loose")},
		{'M', configTree.get<float>("Producer.Jet.DeepAK8.Top." + product.GetEraSelector() + ".Medium")},
		{'T', configTree.get<float>("Producer.Jet.DeepAK8.Top." + product.GetEraSelector() + ".Tight")},
		{'t', configTree.get<float>("Producer.Jet.DeepAK8.Top." + product.GetEraSelector() + ".VeryTight")},
	};

	deepAk8TopMDTagMap = {
		{'L', configTree.get<float>("Producer.Jet.DeepAK8.TopMassDecorrelated." + product.GetEraSelector() + ".Loose")},
		{'M', configTree.get<float>("Producer.Jet.DeepAK8.TopMassDecorrelated." + product.GetEraSelector() + ".Medium")},
		{'T', configTree.get<float>("Producer.Jet.DeepAK8.TopMassDecorrelated." + product.GetEraSelector() + ".Tight")},
		{'t', configTree.get<float>("Producer.Jet.DeepAK8.TopMassDecorrelated." + product.GetEraSelector() + ".VeryTight")},
	};

	deepAk8WTagMap = {
		{'l', configTree.get<float>("Producer.Jet.DeepAK8.W." + product.GetEraSelector() + ".VeryLoose")},
		{'L', configTree.get<float>("Producer.Jet.DeepAK8.W." + product.GetEraSelector() + ".Loose")},
		{'M', configTree.get<float>("Producer.Jet.DeepAK8.W." + product.GetEraSelector() + ".Medium")},
		{'T', configTree.get<float>("Producer.Jet.DeepAK8.W." + product.GetEraSelector() + ".Tight")},
	};

	deepAk8WMDTagMap = {
		{'l', configTree.get<float>("Producer.Jet.DeepAK8.WMassDecorrelated." + product.GetEraSelector() + ".VeryLoose")},
		{'L', configTree.get<float>("Producer.Jet.DeepAK8.WMassDecorrelated." + product.GetEraSelector() + ".Loose")},
		{'M', configTree.get<float>("Producer.Jet.DeepAK8.WMassDecorrelated." + product.GetEraSelector() + ".Medium")},
		{'T', configTree.get<float>("Producer.Jet.DeepAK8.WMassDecorrelated." + product.GetEraSelector() + ".Tight")},
	};


	jetPtCut  = configTree.get<float>("Producer.Jet.Pt");
	jetEtaCut = configTree.get<float>("Producer.Jet.Eta");
	std::cout << std::endl <<
		"The following cuts are applied to Jets:"   << std::endl <<
		"|Eta| < " << jetEtaCut << std::endl <<
		"|Pt|  > " << jetPtCut << std::endl;
}

void JetProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	product.jetPt.fill(0); // This is better for sorting
	product.fatJetPt.fill(0); // This is better for sorting
	dataReader.ReadJetEntry();
	assert(dataReader.nJet <= product.nMax);
	int jetCounter = 0;

	float metPx = dataReader.metPt * std::cos(dataReader.metPhi),
		metPy = dataReader.metPt * std::sin(dataReader.metPhi);

	for (int iJet = 0; iJet < dataReader.nJet; iJet++) {
		dataReader.GetJetValues(iJet);
		if (!product.GetIsData()) { dataReader.ReadGenEntry();}
		const float &jetPtRaw = dataReader.jetPt * (1 - dataReader.jetRawFactor);
		const float &jetMassRaw = dataReader.jetMass * (1 - dataReader.jetRawFactor);

		// Jet Energy Correction (JEC)
		const float &correctionFactor = CorrectEnergy(dataReader, jetPtRaw, true, product.GetIsFastSim());
		// Jet Eergry Resolution (JER)
		std::map<char, float> smearFactor = {{'N', 1.0}, {'U', 1.0}, {'D', 1.0}};
		if (!product.GetIsData()) {
			smearFactor = JetProducer::SmearEnergy(dataReader, correctionFactor * jetPtRaw, true);
		}

		// Get Raw JetPt
		const float &jetPtCorrected        = jetPtRaw * correctionFactor * smearFactor.at('N');
		const float &jetPtCorrectedJerUp   = jetPtRaw * correctionFactor * smearFactor.at('U');
		const float &jetPtCorrectedJerDown = jetPtRaw * correctionFactor * smearFactor.at('D');
		const float &jetMassCorrected      = jetMassRaw * correctionFactor * smearFactor.at('N');

		// If any variation of jetPt passes the cut, keep it
		bool passesJetPtCut = jetPtCorrected >= jetPtCut || jetPtCorrectedJerUp >= jetPtCut || jetPtCorrectedJerDown >= jetPtCut;

		std::vector<float> jetJecUp, jetJecDown;
		if (!product.GetIsData()) {
			for (int iJec = 0; iJec < jecSystematics.size(); iJec++) {
				std::string systematic = jecSystematics.at(iJec);
				ak4CorrectionUncertainty.at(systematic)->setJetPt(correctionFactor * jetPtRaw);
				ak4CorrectionUncertainty.at(systematic)->setJetEta(dataReader.jetEta);
				ak4CorrectionUncertainty.at(systematic)->setJetEta(dataReader.jetPhi);
				const float &correctionFactorUp = correctionFactor * (1 + ak4CorrectionUncertainty.at(systematic)->getUncertainty(true));
				const float &smearFactorJecUp = product.GetIsData() ? 1.0 : JetProducer::SmearEnergy(dataReader, jetPtRaw * correctionFactorUp, true).at('N');
				jetJecUp.push_back(correctionFactorUp * smearFactorJecUp);

				ak4CorrectionUncertainty.at(systematic)->setJetPt(correctionFactor * jetPtRaw);
				ak4CorrectionUncertainty.at(systematic)->setJetEta(dataReader.jetEta);
				ak4CorrectionUncertainty.at(systematic)->setJetEta(dataReader.jetPhi);
				const float &correctionFactorDown = correctionFactor * (1 - ak4CorrectionUncertainty.at(systematic)->getUncertainty(false));
				const float &smearFactorJecDown = product.GetIsData() ? 1.0 : JetProducer::SmearEnergy(dataReader, jetPtRaw * correctionFactorDown, true).at('N');
				jetJecDown.push_back(correctionFactorDown * smearFactorJecDown);

				// if any variation of jetPt (JEC or JER) passes the cut, keep the Jet to avoid redundant skimming for each variation
				if (jetPtRaw * jetJecUp.at(iJec) >= jetPtCut || jetPtRaw * jetJecDown.at(iJec) >= jetPtCut) {
					passesJetPtCut = true;
				}
			}
		}

		if (!passesJetPtCut || std::abs(dataReader.jetEta) > jetEtaCut) { continue;}

		metPx += dataReader.jetPt * std::cos(dataReader.jetPhi) - jetPtCorrected * std::cos(dataReader.jetPhi);
		metPy += dataReader.jetPt * std::sin(dataReader.jetPhi) - jetPtCorrected * std::sin(dataReader.jetPhi);

		product.jetPt[jetCounter]   = jetPtCorrected;
		product.jetEta[jetCounter]  = dataReader.jetEta;
		product.jetPhi[jetCounter]  = dataReader.jetPhi;
		product.jetMass[jetCounter] = dataReader.jetMass * correctionFactor * smearFactor.at('N');

		if (!product.GetIsData()) {
			product.jetPtJerUp[jetCounter]     = jetPtCorrectedJerUp;
			product.jetPtJerDown[jetCounter]   = jetPtCorrectedJerDown;
			product.jetMassJerUp[jetCounter]   = jetMassRaw * correctionFactor * smearFactor.at('U');
			product.jetMassJerDown[jetCounter] = jetMassRaw * correctionFactor * smearFactor.at('D');
			for (int iJec = 0; iJec < jecSystematics.size(); iJec++) {
				product.jetPtJecUp.at(iJec)[jetCounter] = jetPtRaw * jetJecUp.at(iJec);
				product.jetPtJecDown.at(iJec)[jetCounter] = jetPtRaw * jetJecDown.at(iJec);
				product.jetMassJecUp.at(iJec)[jetCounter] = jetMassRaw * jetJecUp.at(iJec);
				product.jetMassJecDown.at(iJec)[jetCounter] = jetMassRaw * jetJecDown.at(iJec);
			}

			product.jetPartFlav[jetCounter] = dataReader.jetPartFlav;
		}

		product.jetDeepJetLooseId[jetCounter]  = dataReader.jetDeepJet > deepJetBTagMap.at('L');
		product.jetDeepJetMediumId[jetCounter] = dataReader.jetDeepJet > deepJetBTagMap.at('M');
		product.jetDeepJetTightId[jetCounter]  = dataReader.jetDeepJet > deepJetBTagMap.at('T');

		product.jetDeepJetId[jetCounter] = product.jetDeepJetLooseId[jetCounter] + product.jetDeepJetMediumId[jetCounter] + product.jetDeepJetTightId[jetCounter];

		jetCounter++;
	}

	// Sort all Jet related variables according to nominal JetPt
	std::vector<int> indices(jetCounter);
	std::iota(indices.begin(), indices.end(), 0);
	std::stable_sort(indices.begin(), indices.end(), [&](int i1, int i2) {return product.jetPt[i1] > product.jetPt[i2];});
	for (std::array<float, product.nMax> *jetVariable : {&product.jetPt, &product.jetEta, &product.jetPhi, &product.jetMass}) {
		Utility::SortByIndex(*jetVariable, indices, jetCounter);
	}
	for (std::array<bool, product.nMax> *jetVariable : {&product.jetDeepJetLooseId, &product.jetDeepJetMediumId, &product.jetDeepJetTightId}) {
		Utility::SortByIndex(*jetVariable, indices, jetCounter);
	}

	/*############################################################
	#   Jet Cleaning                                             #
	#   Remove jets from collection that match to good Leptons   #
	############################################################*/
	float deltaRMin = std::numeric_limits<float>::max();
	std::vector<int> muonIndices, electronIndices, jetRemovalIndices;
	int nearestMuonIndex = -999, nearestElectronIndex = -999;
	for (int iMuon = 0; iMuon < product.nMuon; iMuon++) {
		if (!product.muonIsGood[iMuon]) { continue;}
		int nearestJetIndex = -999;
		for (int iJet = 0; iJet < jetCounter; iJet++) {
			if (std::find(jetRemovalIndices.begin(), jetRemovalIndices.end(),iJet)!=jetRemovalIndices.end()) { continue;} // don't match if a jet is already marked for removal
			float deltaR = Utility::DeltaR(product.jetEta[iJet], product.jetPhi[iJet], product.muonEta[iMuon], product.muonPhi[iMuon]);
			if (deltaR < deltaRMin && deltaR < 0.4) {
				deltaRMin = deltaR;
				nearestMuonIndex = iMuon;
				nearestJetIndex = iJet;
			}
		}

		if (nearestJetIndex > 0) {
			muonIndices.push_back(iMuon);
			jetRemovalIndices.push_back(nearestJetIndex);
		}
	}

	for (int iElectron = 0; iElectron < product.nElectron; iElectron++) {
		if (!product.electronIsGood[iElectron]) { continue;}
		int nearestJetIndex = -999;
		for (int iJet = 0; iJet < jetCounter; iJet++) {
			if (std::find(jetRemovalIndices.begin(), jetRemovalIndices.end(),iJet)!=jetRemovalIndices.end()) { continue;}
			float deltaR = Utility::DeltaR(product.jetEta[iJet], product.jetPhi[iJet], product.electronEta[iElectron], product.electronPhi[iElectron]);
			if (deltaR < deltaRMin && deltaR < 0.4) {
				deltaRMin = deltaR;
				nearestElectronIndex = iElectron;
				nearestJetIndex = iJet;
			}
		}

		if (nearestJetIndex > 0) {
			jetRemovalIndices.push_back(nearestJetIndex);
			electronIndices.push_back(iElectron);
		}
	}

	int removeCounter = 0; // Remove the element corresponding to the index
	for (int iRemove : jetRemovalIndices) {
		for (std::array<float, product.nMax> *jetVariable : {&product.jetPt, &product.jetEta, &product.jetPhi, &product.jetMass}) {
			Utility::RemoveByIndex(jetVariable, iRemove - removeCounter, jetCounter);
		}
		for (std::array<bool, product.nMax> *jetVariable : {&product.jetDeepJetLooseId, &product.jetDeepJetMediumId, &product.jetDeepJetTightId}) {
			Utility::RemoveByIndex(jetVariable, iRemove - removeCounter, jetCounter);
		}
		removeCounter++; // Keep Tack of already deleted jets
	}
	jetCounter -= removeCounter;

	product.nJet = jetCounter;

	product.metPt  = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
	product.metPhi = std::atan2(metPy, metPx);
	product.caloMetPt = dataReader.caloMetPt;

	int fatJetCounter = 0;
	dataReader.ReadFatJetEntry();
	assert(dataReader.nFatJet <= product.nMax);
	for (int iFatJet = 0; iFatJet < dataReader.nFatJet; iFatJet++) {
		dataReader.GetFatJetValues(iFatJet);
		const float &fatJetPtRaw = dataReader.fatJetPt * (1 - dataReader.fatJetRawFactor);
		const float &fatJetMassRaw = dataReader.fatJetMass * (1 - dataReader.fatJetRawFactor);

		float correctionFactor = CorrectEnergy(dataReader, fatJetPtRaw, false, product.GetIsFastSim());
		std::map<char, float> smearFactor = {{'N', 1.0}, {'U', 1.0}, {'D', 1.0}};
		if (!product.GetIsData()) {
			smearFactor = JetProducer::SmearFatEnergy(dataReader, correctionFactor * fatJetPtRaw);
		}

		std::vector<float> fatJetJecUp, fatJetJecDown;
		for (int iJec = 0; iJec < jecSystematics.size(); iJec++) { // Store the correction factors
			std::string systematic = jecSystematics.at(iJec);
			ak8CorrectionUncertainty.at(systematic)->setJetPt(correctionFactor * fatJetPtRaw);
			ak8CorrectionUncertainty.at(systematic)->setJetEta(dataReader.fatJetEta);
			ak8CorrectionUncertainty.at(systematic)->setJetEta(dataReader.fatJetPhi);
			const float &correctionFactorUp = correctionFactor * (1 + ak8CorrectionUncertainty.at(systematic)->getUncertainty(true));
			const float &smearFactorJecUp = product.GetIsData() ? 1.0 : JetProducer::SmearFatEnergy(dataReader, fatJetPtRaw * correctionFactorUp).at('N');
			fatJetJecUp.push_back(correctionFactorUp * smearFactorJecUp);

			ak8CorrectionUncertainty.at(systematic)->setJetPt(correctionFactor * fatJetPtRaw);
			ak8CorrectionUncertainty.at(systematic)->setJetEta(dataReader.fatJetEta);
			ak8CorrectionUncertainty.at(systematic)->setJetEta(dataReader.fatJetPhi);
			const float &correctionFactorDown = correctionFactor * (1 - ak8CorrectionUncertainty.at(systematic)->getUncertainty(false));
			const float &smearFactorJecDown = product.GetIsData() ? 1.0 : JetProducer::SmearFatEnergy(dataReader, fatJetPtRaw * correctionFactorDown).at('N');
			fatJetJecDown.push_back(correctionFactorDown * smearFactorJecDown);
		}

		product.fatJetPt[fatJetCounter]              = fatJetPtRaw * correctionFactor * smearFactor.at('N');
		product.fatJetEta[fatJetCounter]             = dataReader.fatJetEta;
		product.fatJetPhi[fatJetCounter]             = dataReader.fatJetPhi;
		product.fatJetMass[fatJetCounter]            = dataReader.fatJetMass * correctionFactor * smearFactor.at('N');
		product.fatJetArea[fatJetCounter]            = dataReader.fatJetArea;
		product.fatJetId[fatJetCounter]              = dataReader.fatJetId;

		product.fatJetDeepTagMDTvsQCD[fatJetCounter] = dataReader.fatJetDeepTagMDTvsQCD;
		product.fatJetDeepTagMDWvsQCD[fatJetCounter] = dataReader.fatJetDeepTagMDWvsQCD;
		product.fatJetDeepTagTvsQCD[fatJetCounter]   = dataReader.fatJetDeepTagTvsQCD;
		product.fatJetDeepTagWvsQCD[fatJetCounter]   = dataReader.fatJetDeepTagWvsQCD;

		product.fatJetDeepAk8TopLooseId[fatJetCounter]       = dataReader.fatJetDeepTagTvsQCD > deepAk8TopTagMap.at('L');
		product.fatJetDeepAk8TopMediumId[fatJetCounter]      = dataReader.fatJetDeepTagTvsQCD > deepAk8TopTagMap.at('M');
		product.fatJetDeepAk8TopTightId[fatJetCounter]       = dataReader.fatJetDeepTagTvsQCD > deepAk8TopTagMap.at('T');
		product.fatJetDeepAk8TopVeryTightId[fatJetCounter]   = dataReader.fatJetDeepTagTvsQCD > deepAk8TopTagMap.at('t');

		product.fatJetDeepAk8TopMDLooseId[fatJetCounter]     = dataReader.fatJetDeepTagTvsQCD > deepAk8TopMDTagMap.at('L');
		product.fatJetDeepAk8TopMDMediumId[fatJetCounter]    = dataReader.fatJetDeepTagTvsQCD > deepAk8TopMDTagMap.at('M');
		product.fatJetDeepAk8TopMDTightId[fatJetCounter]     = dataReader.fatJetDeepTagTvsQCD > deepAk8TopMDTagMap.at('T');
		product.fatJetDeepAk8TopMDVeryTightId[fatJetCounter] = dataReader.fatJetDeepTagTvsQCD > deepAk8TopMDTagMap.at('t');

		product.fatJetDeepAk8WVeryLooseId[fatJetCounter]     = dataReader.fatJetDeepTagTvsQCD > deepAk8WTagMap.at('l');
		product.fatJetDeepAk8WLooseId[fatJetCounter]         = dataReader.fatJetDeepTagTvsQCD > deepAk8WTagMap.at('L');
		product.fatJetDeepAk8WMediumId[fatJetCounter]        = dataReader.fatJetDeepTagTvsQCD > deepAk8WTagMap.at('M');
		product.fatJetDeepAk8WTightId[fatJetCounter]         = dataReader.fatJetDeepTagTvsQCD > deepAk8WTagMap.at('T');

		product.fatJetDeepAk8WMDVeryLooseId[fatJetCounter]   = dataReader.fatJetDeepTagTvsQCD > deepAk8WMDTagMap.at('l');
		product.fatJetDeepAk8WMDLooseId[fatJetCounter]       = dataReader.fatJetDeepTagTvsQCD > deepAk8WMDTagMap.at('L');
		product.fatJetDeepAk8WMDMediumId[fatJetCounter]      = dataReader.fatJetDeepTagTvsQCD > deepAk8WMDTagMap.at('M');
		product.fatJetDeepAk8WMDTightId[fatJetCounter]       = dataReader.fatJetDeepTagTvsQCD > deepAk8WMDTagMap.at('T');

		for (int iJec = 0; iJec < jecSystematics.size(); iJec++) {
			product.fatJetPtJecUp.at(iJec)[fatJetCounter] = fatJetPtRaw * fatJetJecUp.at(iJec);
			product.fatJetPtJecDown.at(iJec)[fatJetCounter] = fatJetPtRaw * fatJetJecDown.at(iJec);
			product.fatJetMassJecUp.at(iJec)[fatJetCounter] = fatJetMassRaw * fatJetJecUp.at(iJec);
			product.fatJetMassJecDown.at(iJec)[fatJetCounter] = fatJetMassRaw * fatJetJecDown.at(iJec);
		}

		fatJetCounter++;
	}
	product.nFatJet = fatJetCounter;
}

float JetProducer::CorrectEnergy(DataReader &dataReader, const float &jetPtRaw, const bool &isAk4, const bool &isFastSim) {
	const float &jetEta = isAk4 ? dataReader.jetEta : dataReader.fatJetEta;
	const float &jetArea = isAk4 ? dataReader.jetArea : dataReader.fatJetArea;
	const std::string &jetAlgorithm = isAk4 ? ak4Algorithm : ak8Algorithm;
	if (isFastSim) {
		std::shared_ptr<FactorizedJetCorrector>& corr = isAk4 ? fastJetCorrectorAk4 : fastJetCorrectorAk8;
		corr->setJetPt(jetPtRaw);
		corr->setJetEta(jetEta);
		corr->setJetPhi(isAk4 ? dataReader.jetPhi : dataReader.fatJetPhi);
		corr->setRho(dataReader.rho);
		corr->setJetA(jetArea);
		return corr->getCorrection();
	} else {
		const std::unique_ptr<correction::CorrectionSet> &jetCorrectionSet = isAk4 ? ak4CorrectionSet : ak8CorrectionSet;
		return jetCorrectionSet->at(era + "_" + dataType + "_L1FastJet_" + jetAlgorithm)->evaluate({dataReader.jetArea, dataReader.jetEta, jetPtRaw, dataReader.rho}) *
			jetCorrectionSet->at(era + "_" + dataType + "_L2Relative_"   + jetAlgorithm)->evaluate({dataReader.jetEta, jetPtRaw}) *
			jetCorrectionSet->at(era + "_" + dataType + "_L3Absolute_"   + jetAlgorithm)->evaluate({dataReader.jetEta, jetPtRaw}) *
			jetCorrectionSet->at(era + "_" + dataType + "_L2L3Residual_" + jetAlgorithm)->evaluate({dataReader.jetEta, jetPtRaw});
	}
}

/*#######################################################################################################################
#   https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures                                      #
#   https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263   #
#   !!!Warning!!! - json format not available for AK8, use old school smearing with SmearFatEnergy                      *
#######################################################################################################################*/
std::map<char, float> JetProducer::SmearEnergy(DataReader &dataReader, const float &jetPtCorrected, const bool &isAk4) {
	const float &coneSize = isAk4 ? 0.2 : 0.4;
	const float &jetEta = isAk4 ? dataReader.jetEta : dataReader.fatJetEta;
	const float &jetPhi = isAk4 ? dataReader.jetPhi : dataReader.fatJetPhi;
	const std::unique_ptr<correction::CorrectionSet> &jetCorrectionSet = isAk4 ? ak4CorrectionSet : ak8CorrectionSet;

	float resolution       = jetCorrectionSet->at(era + "_" + jerVersion + "_PtResolution_" + (isAk4 ? ak4Algorithm : ak8Algorithm))->evaluate({jetEta, jetPtCorrected, dataReader.rho});
	float resolutionSf     = jetCorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + (isAk4 ? ak4Algorithm : ak8Algorithm))->evaluate({jetEta, "nom"});
	float resolutionSfUp   = jetCorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + (isAk4 ? ak4Algorithm : ak8Algorithm))->evaluate({jetEta, "up"});
	float resolutionSfDown = jetCorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + (isAk4 ? ak4Algorithm : ak8Algorithm))->evaluate({jetEta, "down"});

	/*#############################################################################################
	#   Do not use dataReader.GetMatchedIndex here                                                #
	#   That would not just match to GenJets but GenParts, thus be less efficient                 #
	#   Also the deltaPt criterion is different as it does not take the resolution into account   #
	#############################################################################################*/
	bool isMatched = false;
	float deltaR, deltaPt, matchedGenJetPt,
		deltaRMin = std::numeric_limits<float>::max(),
		deltaPtMin = std::numeric_limits<float>::max();

	if (isAk4) {
		dataReader.ReadGenJetEntry();
	} else {
		dataReader.ReadGenFatJetEntry();
	}

	int nGenJet = isAk4 ? dataReader.nGenJet : dataReader.nGenFatJet;
	for (int iGen = 0; iGen < nGenJet; iGen++) {
		if (isAk4) {
			dataReader.GetGenJetValues(iGen);
		} else {
			dataReader.GetGenFatJetValues(iGen);
		}

		float genJetPt = isAk4 ? dataReader.genJetPt : dataReader.genFatJetPt;
		float genJetEta = isAk4 ? dataReader.genJetEta : dataReader.genFatJetEta;
		float genJetPhi = isAk4 ? dataReader.genJetPhi : dataReader.genFatJetPhi;
		deltaR  = Utility::DeltaR(jetEta, jetPhi, genJetEta, genJetPhi);
		deltaPt = std::abs(genJetPt - jetPtCorrected) / jetPtCorrected;

		if (deltaR > deltaRMin) { continue;}

		if (deltaR < coneSize && deltaPt < 3 * resolution * jetPtCorrected) {
			matchedGenJetPt = genJetPt;
			deltaRMin = deltaR;
			isMatched = true;
		}
	}

	float smearFactor = 1.0, smearFactorUp = 1.0, smearFactorDown = 1.0;
	if (isMatched) {
		smearFactor     = 1. + (resolutionSf - 1)     * (jetPtCorrected - matchedGenJetPt) / jetPtCorrected;
		smearFactorUp   = 1. + (resolutionSfUp - 1)   * (jetPtCorrected - matchedGenJetPt) / jetPtCorrected;
		smearFactorDown = 1. + (resolutionSfDown - 1) * (jetPtCorrected - matchedGenJetPt) / jetPtCorrected;
	} else {
		if (resolutionSf > 1.) {
			std::default_random_engine generator; std::normal_distribution<> gaus(0, resolution * std::sqrt(resolutionSf * resolutionSf - 1));
			smearFactor = 1. + gaus(generator);
		}
		if (resolutionSfUp > 1.) {
			std::default_random_engine generator; std::normal_distribution<> gaus(0, resolution * std::sqrt(resolutionSfUp * resolutionSfUp - 1));
			smearFactorUp = 1. + gaus(generator);
		}
		if (resolutionSfDown > 1.) {
			std::default_random_engine generator; std::normal_distribution<> gaus(0, resolution * std::sqrt(resolutionSfDown * resolutionSfDown - 1));
			smearFactorDown = 1. + gaus(generator);
		}
	}

	/*#####################################################
	#   Negative or too small smearFactor.                #
	#   We would change direction of the jet              #
	#   and this is not what we want.                     #
	#   Recompute the smearing factor in order to have:   #
	#   jet.energy() == MIN_JET_ENERGY                    #
	#####################################################*/
	if (jetPtCorrected * smearFactor < 1e-2)     { smearFactor     = 1e-2 / jetPtCorrected;}
	if (jetPtCorrected * smearFactorUp < 1e-2)   { smearFactorUp   = 1e-2 / jetPtCorrected;}
	if (jetPtCorrected * smearFactorDown < 1e-2) { smearFactorDown = 1e-2 / jetPtCorrected;}

	/******************
	*   N : Nominal   *
	*   U : Up        *
	*   D : Down      *
	******************/
	//std::cout << "{{'N', " << smearFactor << "}, {'U', " << smearFactorUp << "}, {'D', " << smearFactorDown << "}}" << std::endl;
	return {{'N', smearFactor}, {'U', smearFactorUp}, {'D', smearFactorDown}};
}

// The energy smearing is not available for fatjets in the current release, so smear it using the old procedure
std::map<char, float> JetProducer::SmearFatEnergy(DataReader &dataReader, const float &jetPtCorrected) {
	float coneSize = 0.4;
	float jetEta = dataReader.fatJetEta;
	float jetPhi = dataReader.fatJetPhi;

	jetParameter.setJetPt(jetPtCorrected).setJetEta(jetEta).setRho(dataReader.rho);

	float resolution       = resolutionAk8.getResolution(jetParameter);
	float resolutionSf     = resolutionSfAk8.getScaleFactor(jetParameter);
	float resolutionSfUp   = resolutionSfAk8.getScaleFactor(jetParameter, Variation::UP);
	float resolutionSfDown = resolutionSfAk8.getScaleFactor(jetParameter, Variation::DOWN);

	/*#############################################################################################
	#   Do not use dataReader.GetMatchedIndex here                                                #
	#   That would not just match to GenJets but GenParts, thus be less efficient                 #
	#   Also the deltaPt criterion is different as it does not take the resolution into account   #
	#############################################################################################*/
	bool isMatched = false;
	float deltaR, deltaPt, matchedGenJetPt,
		deltaRMin = std::numeric_limits<float>::max(),
		deltaPtMin = std::numeric_limits<float>::max();

	dataReader.ReadGenFatJetEntry();

	int nGenJet = dataReader.nGenFatJet;
	for(int iGen = 0; iGen < nGenJet; iGen++) {
		dataReader.GetGenFatJetValues(iGen);

		float genJetPt = dataReader.genFatJetPt;
		float genJetEta = dataReader.genFatJetEta;
		float genJetPhi = dataReader.genFatJetPhi;
		deltaR  = Utility::DeltaR(jetEta, jetPhi, genJetEta, genJetPhi);
		deltaPt = std::abs(genJetPt - jetPtCorrected) / jetPtCorrected;

		if (deltaR > deltaRMin) { continue;}

		if (deltaR < coneSize && deltaPt < 3 * resolution * jetPtCorrected) {
			matchedGenJetPt = genJetPt;
			deltaRMin = deltaR;
			isMatched = true;
		}
	}

	float smearFactor = 1.0, smearFactorUp = 1.0, smearFactorDown = 1.0;
	if (isMatched) {
		smearFactor     = 1. + (resolutionSf - 1)     * (jetPtCorrected - matchedGenJetPt) / jetPtCorrected;
		smearFactorUp   = 1. + (resolutionSfUp - 1)   * (jetPtCorrected - matchedGenJetPt) / jetPtCorrected;
		smearFactorDown = 1. + (resolutionSfDown - 1) * (jetPtCorrected - matchedGenJetPt) / jetPtCorrected;
	} else {
		if (resolutionSf > 1.) {
			std::default_random_engine generator; std::normal_distribution<> gaus(0, resolution * std::sqrt(resolutionSf * resolutionSf - 1));
			smearFactor = 1. + gaus(generator);
		}
		if (resolutionSfUp > 1.) {
			std::default_random_engine generator; std::normal_distribution<> gaus(0, resolution * std::sqrt(resolutionSfUp * resolutionSfUp - 1));
			smearFactorUp = 1. + gaus(generator);
		}
		if (resolutionSfDown > 1.) {
			std::default_random_engine generator; std::normal_distribution<> gaus(0, resolution * std::sqrt(resolutionSfDown * resolutionSfDown - 1));
			smearFactorDown = 1. + gaus(generator);
		}
	}

	/*#####################################################
	#   Negative or too small smearFactor.                #
	#   We would change direction of the jet              #
	#   and this is not what we want.                     #
	#   Recompute the smearing factor in order to have:   #
	#   jet.energy() == MIN_JET_ENERGY                    #
	#####################################################*/
	if (jetPtCorrected * smearFactor < 1e-2)     { smearFactor     = 1e-2 / jetPtCorrected;}
	if (jetPtCorrected * smearFactorUp < 1e-2)   { smearFactorUp   = 1e-2 / jetPtCorrected;}
	if (jetPtCorrected * smearFactorDown < 1e-2) { smearFactorDown = 1e-2 / jetPtCorrected;}

	/******************
	*   N : Nominal   *
	*   U : Up        *
	*   D : Down      *
	******************/
	//std::cout << "{{'N', " << smearFactor << "}, {'U', " << smearFactorUp << "}, {'D', " << smearFactorDown << "}}" << std::endl;
	return {{'N', smearFactor}, {'U', smearFactorUp}, {'D', smearFactorDown}};

}

void JetProducer::EndJob(TFile &file) {}
