#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/JetProducer.h>

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>

#include <cmath>
#include <random>

JetProducer::JetProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product) {
	std::string cmsswBase = std::getenv("CMSSW_BASE");

	/*################################################################################################
	#   https://github.com/cms-nanoAOD/correctionlib                                                 #
	#   https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master                         #
	#   Instructions:                                                                                #
	#   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite   #
	################################################################################################*/
	era              = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".Era");
	dataType         = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + (product.GetIsData()? ".DATA" : ".MC"));
	runPeriod        = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".RunPeriod." + "M");
	jerVersion       = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".JER");
	ak4Algorithm     = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK4.Algorithm");
	ak8Algorithm     = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK8.Algorithm");
	ak4CorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK4.JSON"));
	ak8CorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK8.JSON"));
	jmeCorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".JME.JSON"));

	// Set object to get JEC uncertainty
	bool isJECSysteamtic = true; //FIXME maybe make this part of the product?
	if (isJECSystematic) {
		jetCorrectionUncertainty = std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JECUNC." + product.GetEraSelector()), "Total"));
	}


	/*################################################################################
	#   Definition of Working Points come from                                       #
	#   https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP    #
	#   https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP   #
	#   https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL17      #
	#   https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP   #
	################################################################################*/
	deepCsvBTagMap = {
		{'L', configTree.get<double>("Producer.Jet.DeepCSV." + product.GetEraSelector() + ".Loose")},
		{'M', configTree.get<double>("Producer.Jet.DeepCSV." + product.GetEraSelector() + ".Medium")},
		{'T', configTree.get<double>("Producer.Jet.DeepCSV." + product.GetEraSelector() + ".Tight")},
	};

	deepJetBTagMap = {
		{'L', configTree.get<double>("Producer.Jet.DeepJet." + product.GetEraSelector() + ".Loose")},
		{'M', configTree.get<double>("Producer.Jet.DeepJet." + product.GetEraSelector() + ".Medium")},
		{'T', configTree.get<double>("Producer.Jet.DeepJet." + product.GetEraSelector() + ".Tight")},
	};

	jetPtCut  = configTree.get<double>("Producer.Jet.Pt");
	jetEtaCut = configTree.get<double>("Producer.Jet.Eta");
	std::cout << std::endl <<
		"The following cuts are applied to Jets:"   << std::endl <<
		"|Eta| < " << jetEtaCut << std::endl <<
		"|Pt|  > " << jetPtCut << std::endl;
}

/*#######################################################################################################################
#   https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures                                      #
#   https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263   #
#######################################################################################################################*/
double JetProducer::SmearEnergy(DataReader &dataReader, const double &jetPtCorrected, const bool &isAk4) {
	double coneSize = isAk4 ? 0.2 : 0.4;
	double jetEta = isAk4 ? dataReader.jetEta : dataReader.fatJetEta;
	double jetPhi = isAk4 ? dataReader.jetPhi : dataReader.fatJetPhi;
	double resolution = ak4CorrectionSet->at(era + "_" + jerVersion + "_PtResolution_" + ak4Algorithm)->evaluate({jetEta, jetPtCorrected, dataReader.rho});
	double resolutionSF;
	if (isJERSystematic) {
		resolutionSF = isUp ? ak4CorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + ak4Algorithm)->evaluate({jetEta, "up"}) : ak4CorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + ak4Algorithm)->evaluate({jetEta, "down"});
	} else {
		resolutionSF = ak4CorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + ak4Algorithm)->evaluate({jetEta, "nom"});
	}

	/*#############################################################################################
	#   Do not use dataReader.GetMatchedIndex here                                                #
	#   That would not just match to GenJets but GenParts, thus be less efficient                 #
	#   Also the deltaPt criterion is different as it does not take the resolution into account   #
	#############################################################################################*/
	bool isMatched = false;
	double deltaR, deltaPt, matchedGenJetPt,
		deltaRMin = std::numeric_limits<double>::max(),
		deltaPtMin = std::numeric_limits<double>::max();

	if (isAk4) {
		dataReader.ReadGenJetEntry();
	} else {
		dataReader.ReadGenFatJetEntry();
	}

	int nGenJet = isAk4 ? dataReader.nGenJet : dataReader.nGenFatJet;
	for(int iGen = 0; iGen < nGenJet; iGen++) {
		if (isAk4) {
			dataReader.GetGenJetValues(iGen);
		} else {
			dataReader.GetGenFatJetValues(iGen);
		}

		double genJetPt = isAk4 ? dataReader.genJetPt : dataReader.genFatJetPt;
		double genJetEta = isAk4 ? dataReader.genJetEta : dataReader.genFatJetEta;
		double genJetPhi = isAk4 ? dataReader.genJetPhi : dataReader.genFatJetPhi;
		deltaR  = Utility::DeltaR(jetEta, jetPhi, genJetEta, genJetPhi);
		deltaPt = std::abs(genJetPt - jetPtCorrected) / jetPtCorrected;

		if (deltaR > deltaRMin) { continue;}

		if (deltaR < coneSize && deltaPt < 3 * resolution * jetPtCorrected) {
			matchedGenJetPt = genJetPt;
			deltaRMin = deltaR;
			isMatched = true;
		}
	}

	double smearFactor = 1.0;
	if (isMatched) {
		smearFactor = 1. + (resolutionSF - 1) * (jetPtCorrected - matchedGenJetPt) / jetPtCorrected;
	} else if (resolutionSF > 1.) {
		std::default_random_engine generator; std::normal_distribution<> gaus(0, resolution * std::sqrt(resolutionSF * resolutionSF - 1));
		smearFactor = 1. + gaus(generator);
	}

	/*#####################################################
	#   Negative or too small smearFactor.                #
	#   We would change direction of the jet              #
	#   and this is not what we want.                     #
	#   Recompute the smearing factor in order to have:   #
	#   jet.energy() == MIN_JET_ENERGY                    #
	#####################################################*/
	if (jetPtCorrected * smearFactor < 1e-2) { smearFactor = 1e-2 / jetPtCorrected;}

	return smearFactor;
}

void JetProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	product.jetPt.fill(0); // This might not be needed but it makes testing a lot easier
	dataReader.ReadJetEntry();
	assert(dataReader.nJet < product.nMax);
	int jetCounter = 0;

	double metPx = dataReader.metPt * std::cos(dataReader.metPhi),
		metPy = dataReader.metPt * std::sin(dataReader.metPhi);

	for (int iJet = 0; iJet < dataReader.nJet; iJet++) {
		dataReader.GetJetValues(iJet);
		if (!product.GetIsData()) { dataReader.ReadGenEntry();}

		// FIXME Missing runPeriod for data
		double correctionFactor = ak4CorrectionSet->at(era + "_" + dataType + "_L1FastJet_" + ak4Algorithm)->evaluate({dataReader.jetArea, dataReader.jetEta, dataReader.jetPt, dataReader.rho}) *
			ak4CorrectionSet->at(era + "_" + dataType + "_L2Relative_"   + ak4Algorithm)->evaluate({dataReader.jetEta, dataReader.jetPt}) *
			ak4CorrectionSet->at(era + "_" + dataType + "_L3Absolute_"   + ak4Algorithm)->evaluate({dataReader.jetEta, dataReader.jetPt}) *
			ak4CorrectionSet->at(era + "_" + dataType + "_L2L3Residual_" + ak4Algorithm)->evaluate({dataReader.jetEta, dataReader.jetPt});

		//if (isJECSystematic) {
		bool isUp = true; //FIXME
		if (false) {
			jetCorrectionUncertainty->setJetPt(correctionFactor * dataReader.jetPt);
			jetCorrectionUncertainty->setJetEta(dataReader.jetEta);
			const double uncertainty = jetCorrectionUncertainty->getUncertainty(isUp);
			correctionFactor *= isUp ? (1 + uncertainty) : (1 - uncertainty);
		}

		// Smearing TODO(maybe?) only implemented for AK4
		const double &smearFactor = product.GetIsData() ? 1.0 : JetProducer::SmearEnergy(dataReader, correctionFactor * dataReader.jetPt, true);

		const double &jetPtCorrected = dataReader.jetPt * correctionFactor * smearFactor;
		if (jetPtCorrected > jetPtCut && std::abs(dataReader.jetEta) > jetEtaCut) { continue;}

		metPx += dataReader.jetPt * std::cos(dataReader.jetPhi) - jetPtCorrected * std::cos(dataReader.jetPhi);
		metPy += dataReader.jetPt * std::sin(dataReader.jetPhi) - jetPtCorrected * std::sin(dataReader.jetPhi);

		product.jetPt[jetCounter]   = jetPtCorrected;
		product.jetEta[jetCounter]  = dataReader.jetEta;
		product.jetPhi[jetCounter]  = dataReader.jetPhi;
		product.jetMass[jetCounter] = dataReader.jetMass * correctionFactor * smearFactor;

		product.jetDeepCsvLooseId[jetCounter]  = dataReader.jetDeepCsv > deepCsvBTagMap.at('L');
		product.jetDeepCsvMediumId[jetCounter] = dataReader.jetDeepCsv > deepCsvBTagMap.at('M');
		product.jetDeepCsvTightId[jetCounter]  = dataReader.jetDeepCsv > deepCsvBTagMap.at('T');
		product.jetDeepJetLooseId[jetCounter]  = dataReader.jetDeepJet > deepJetBTagMap.at('L');
		product.jetDeepJetMediumId[jetCounter] = dataReader.jetDeepJet > deepJetBTagMap.at('M');
		product.jetDeepJetTightId[jetCounter]  = dataReader.jetDeepJet > deepJetBTagMap.at('T');

		jetCounter++;
	}

	std::vector<int> indices(jetCounter);
	std::iota(indices.begin(), indices.end(), 0);
	std::stable_sort(indices.begin(), indices.end(), [&](int i1, int i2) {return product.jetPt[i1] > product.jetPt[i2];});
	for (std::array<double, product.nMax> *jetVariable : {&product.jetPt, &product.jetEta, &product.jetPhi, &product.jetMass}) {
		SortByIndex(*jetVariable, indices, jetCounter);
	}
	for (std::array<bool, product.nMax> *jetVariable : {&product.jetDeepCsvLooseId, &product.jetDeepCsvMediumId, &product.jetDeepCsvTightId, &product.jetDeepJetLooseId, &product.jetDeepJetMediumId, &product.jetDeepJetTightId}) {
		SortByIndex(*jetVariable, indices, jetCounter);
	}

	/*############################################################
	#   Jet Cleaning                                             #
	#   Remove jets from collection that match to good Leptons   #
	############################################################*/
	double deltaRMin = std::numeric_limits<double>::max();
	std::vector<int> muonIndices, electronIndices, jetRemovalIndices;;
	int nearestMuonIndex = -999, nearestElectronIndex = -999;
	for (int iMuon = 0; iMuon < product.nMuon; iMuon++) {
		if (!product.muonIsGood[iMuon]) { continue;}
		if (std::find(muonIndices.begin(), muonIndices.end(), iMuon)!=muonIndices.end()) { continue;}
		int nearestJetIndex = -999;
		for (int iJet = 0; iJet < jetCounter; iJet++) {
			if (std::find(jetRemovalIndices.begin(), jetRemovalIndices.end(),iJet)!=jetRemovalIndices.end()) { continue;}
			double deltaR = Utility::DeltaR(product.jetEta[iJet], product.jetPhi[iJet], product.muonEta[iMuon], product.muonPhi[iMuon]);
			if (deltaR < deltaRMin && deltaR < 0.4) { // FIXME change 0.4 do configurable param
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
		if (std::find(electronIndices.begin(), electronIndices.end(),iElectron)!=electronIndices.end()) { continue;}
		int nearestJetIndex = -999;
		for (int iJet = 0; iJet < jetCounter; iJet++) {
			if (std::find(jetRemovalIndices.begin(), jetRemovalIndices.end(),iJet)!=jetRemovalIndices.end()) { continue;}
			double deltaR = Utility::DeltaR(product.jetEta[iJet], product.jetPhi[iJet], product.electronEta[iElectron], product.electronPhi[iElectron]);
			if (deltaR < deltaRMin && deltaR < 0.4) { // FIXME change 0.4 do configurable param
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

	int removeCounter = 0; // Since the index positions shift after removing one entry, one has to adjust for this
	for (int iRemove : jetRemovalIndices) {
		for (std::array<double, product.nMax> *jetVariable : {&product.jetPt, &product.jetEta, &product.jetPhi, &product.jetMass}) {
			RemoveByIndex(jetVariable, iRemove - removeCounter, jetCounter);
		}
		for (std::array<bool, product.nMax> *jetVariable : {&product.jetDeepCsvLooseId, &product.jetDeepCsvMediumId, &product.jetDeepCsvTightId, &product.jetDeepJetLooseId, &product.jetDeepJetMediumId, &product.jetDeepJetTightId}) {
			RemoveByIndex(jetVariable, iRemove - removeCounter, jetCounter);
		}
		removeCounter++;
	}
	jetCounter -= removeCounter;

	product.nJet = jetCounter;

	product.metPt  = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
	product.metPhi = std::atan2(metPy, metPx);

	// This is not strictly needed but was done in previous analysis FIXME Maybe delete this? Or Save index information
	product.wBosonMinMass    =  999,
	product.wBosonMinMassPt  = -999,
	product.wBosonBestMass   = -999,
	product.wBosonBestMassPt = -999,
	product.topBestMass      = -999,
	product.topBestMassPt    = -999;
	for (int iJet1 = 0; iJet1 < product.nJet; iJet1++) {
		ROOT::Math::PtEtaPhiMVector jet1P4 = ROOT::Math::PtEtaPhiMVector(product.jetPt.at(iJet1), product.jetEta.at(iJet1), product.jetPhi.at(iJet1), product.jetMass.at(iJet1));
		if (!product.jetDeepCsvMediumId.at(iJet1)) {
			for (int iJet2 = iJet1 + 1; iJet2 < product.nJet; iJet2++) {
				ROOT::Math::PtEtaPhiMVector jet2P4 = ROOT::Math::PtEtaPhiMVector(product.jetPt.at(iJet2), product.jetEta.at(iJet2), product.jetPhi.at(iJet2), product.jetMass.at(iJet2));
				ROOT::Math::PtEtaPhiMVector diJetP4 = jet1P4 + jet2P4;
				double diJetMass = diJetP4.M();
				if (diJetMass > 30 && diJetMass < product.wBosonMinMass) {
					product.wBosonMinMass = diJetMass;
					product.wBosonMinMassPt = diJetP4.Pt();
				}

				if (std::abs(diJetMass - pdgWBosonMass) < std::abs(product.wBosonBestMass - pdgWBosonMass)) {
					product.wBosonBestMass = diJetMass;
					product.wBosonBestMassPt = diJetP4.Pt();
				}

				for (int iBJet = 0; iBJet < product.nJet; iBJet++) {
					if (product.jetDeepCsvMediumId.at(iBJet) &&
							(Utility::DeltaR(product.jetEta.at(iBJet), product.jetPhi.at(iBJet), jet1P4.Eta(), jet1P4.Phi()) < 0.1 ||
							 Utility::DeltaR(product.jetEta.at(iBJet), product.jetPhi.at(iBJet), jet2P4.Eta(), jet2P4.Phi()) < 0.1)
					) {
						ROOT::Math::PtEtaPhiMVector bJetP4 = ROOT::Math::PtEtaPhiMVector(product.jetPt.at(iJet2), product.jetEta.at(iJet2), product.jetPhi.at(iJet2), product.jetMass.at(iJet2));
						ROOT::Math::PtEtaPhiMVector topP4 = diJetP4 + bJetP4;
						double topMass = topP4.M();
						if (std::abs(topMass - pdgTopMass) < std::abs(product.topBestMass - pdgTopMass)) {
							product.topBestMass = topMass;
							product.topBestMassPt = topP4.Pt();
						}
					}
				}
			}
		}
	}

	int fatJetCounter = 0;
	dataReader.ReadFatJetEntry();
	assert(dataReader.nFatJet < product.nMax);
	for (int iFatJet = 0; iFatJet < dataReader.nFatJet; iFatJet++) {
		dataReader.GetFatJetValues(iFatJet);

		double correctionFactor = ak8CorrectionSet->at(era + "_" + dataType + "_L1FastJet_" + ak8Algorithm)->evaluate({dataReader.fatJetArea, dataReader.fatJetEta, dataReader.fatJetPt, dataReader.rho}) *
			ak8CorrectionSet->at(era + "_" + dataType + "_L2Relative_"   + ak8Algorithm)->evaluate({dataReader.fatJetEta, dataReader.fatJetPt}) *
			ak8CorrectionSet->at(era + "_" + dataType + "_L3Absolute_"   + ak8Algorithm)->evaluate({dataReader.fatJetEta, dataReader.fatJetPt}) *
			ak8CorrectionSet->at(era + "_" + dataType + "_L2L3Residual_" + ak8Algorithm)->evaluate({dataReader.fatJetEta, dataReader.fatJetPt});

		const double &smearFactor = product.GetIsData() ? 1.0 : JetProducer::SmearEnergy(dataReader, correctionFactor * dataReader.fatJetPt, false);

		product.fatJetMass[fatJetCounter]            = dataReader.fatJetMass * correctionFactor * smearFactor;
		product.fatJetPt[fatJetCounter]              = dataReader.fatJetPt * correctionFactor * smearFactor;
		product.fatJetEta[fatJetCounter]             = dataReader.fatJetEta;
		product.fatJetPhi[fatJetCounter]             = dataReader.fatJetPhi;
		product.fatJetArea[fatJetCounter]            = dataReader.fatJetArea;
		product.fatJetRawFactor[fatJetCounter]       = dataReader.fatJetRawFactor;
		product.fatJetId[fatJetCounter]              = dataReader.fatJetId;
		product.fatJetDeepTagMDTvsQCD[fatJetCounter] = dataReader.fatJetDeepTagMDTvsQCD;
		product.fatJetDeepTagMDWvsQCD[fatJetCounter] = dataReader.fatJetDeepTagMDWvsQCD;
		product.fatJetDeepTagTvsQCD[fatJetCounter]   = dataReader.fatJetDeepTagTvsQCD;
		product.fatJetDeepTagWvsQCD[fatJetCounter]   = dataReader.fatJetDeepTagWvsQCD;
		fatJetCounter++;
	}
	product.nFatJet = fatJetCounter;
}

template <typename T>
void JetProducer::RemoveByIndex(T &array, int removeIndex, int arraySize) {
	std::copy(array->begin() + removeIndex + 1, // copy everything starting here
		array->begin() + arraySize,         // and ending here, not including it,
		array->begin() + removeIndex        // to this destination
	);
}

template <typename T>
void JetProducer::SortByIndex(T &array, std::vector<int> indices, int size) {
	T tmp = array;
	for (std::size_t i = 0; i < size; i++) {
		array[i] = tmp[indices[i]];
	}
}

void JetProducer::EndJob(TFile &file) {}
