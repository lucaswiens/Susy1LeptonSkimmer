#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/JetProducer.h>

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>

#include <cmath>
#include <random>

JetProducer::JetProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product) {
	std::string cmsswBase = std::getenv("CMSSW_BASE");

	era          = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".Era");
	dataType     = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + (product.GetIsData()? ".DATA" : ".MC"));
	runPeriod    = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".RunPeriod." + "M");
	jerVersion   = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".JER");
	ak4Algorithm = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK4.Algorithm");
	ak8Algorithm = scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK8.Algorithm");
	ak4CorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK4.JSON"));
	ak8CorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".AK8.JSON"));
	jmeCorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JERC." + product.GetEraSelector() + ".JME.JSON"));

	// Set object to get JEC uncertainty
	// TODO JSON format not avialable yet
	bool isJECSysteamtic = true; //FIXME maybe make this part of the product?
	if (isJECSystematic) {
		jetCorrectionUncertainty = std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Jet.JECUNC." + product.GetEraSelector()), "Total"));
	}

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

	jetPtCut    = configTree.get<double>("Producer.Jet.Pt");
	jetEtaCut   = configTree.get<double>("Producer.Jet.Eta");
	std::cout << std::endl <<
		"The following cuts are applied to Jets:"   << std::endl <<
		"|Eta| < " << jetEtaCut << std::endl <<
		"|Pt|  > " << jetPtCut << std::endl;
}

// For Testing TODO Delete this
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
/*
double JetProducer::CorrectEnergy(const double &pt, const double &eta, const double &rho, const double &area) {
	jetCorrector->setJetPt(pt);
	jetCorrector->setJetEta(eta);
	jetCorrector->setRho(rho);
	jetCorrector->setJetA(area);
	return jetCorrector->getCorrection();
}
*/


// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
double JetProducer::SmearEnergy(DataReader &dataReader, const double &jetPtCorrected, const double &coneSize) {
	double resolution = ak4CorrectionSet->at(era + "_" + jerVersion + "_PtResolution_" + ak4Algorithm)->evaluate({dataReader.jetEta, jetPtCorrected, dataReader.rho});
	double resolutionSF;
	if (isJERSystematic) {
		resolutionSF = isUp ? ak4CorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + ak4Algorithm)->evaluate({dataReader.jetEta, "up"}) : ak4CorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + ak4Algorithm)->evaluate({dataReader.jetEta, "down"});
	} else {
		resolutionSF = ak4CorrectionSet->at(era + "_" + jerVersion + "_ScaleFactor_" + ak4Algorithm)->evaluate({dataReader.jetEta, "nom"});
	}

	/*#############################################################################################
	#   Do not use dataReader.GetMatchedIndex here                                                #
	#   That would not just match to GenJets but GenParts, thus be less efficient                 #
	#   Also the deltaPt criterion is different as it does not take the resolution into account   #
	#############################################################################################*/
	bool isMatched = false;
	double deltaR, deltaPt, genJetPt,
		deltaRMin = std::numeric_limits<double>::max(),
		deltaPtMin = std::numeric_limits<double>::max();
	dataReader.ReadGenJetEntry();
	for(int iGen = 0; iGen < dataReader.nGenJet; iGen++) {
		dataReader.GetGenJetValues(iGen);

		deltaR  = Utility::DeltaR(dataReader.jetEta, dataReader.jetPhi, dataReader.genJetEta, dataReader.genJetPhi);
		deltaPt = std::abs(dataReader.genJetPt - jetPtCorrected) / jetPtCorrected;

		if (deltaR > deltaRMin) { continue;}

		if (deltaR < coneSize / 2 && deltaPt < 3 * resolution * jetPtCorrected) {
			genJetPt = dataReader.genJetPt;
			deltaRMin = deltaR;
			isMatched = true;
		}
	}

	double smearFactor = 1.0;
	if (isMatched) {
		smearFactor = 1. + (resolutionSF - 1) * (jetPtCorrected - dataReader.genJetPt) / jetPtCorrected;
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
	int jetCounter = 0, fatJetCounter =0,
		looseCSVBTagCounter = 0, mediumCSVBTagCounter = 0, tightCSVBTagCounter = 0,
		looseDeepJetBTagCounter = 0, mediumDeepJetBTagCounter = 0, tightDeepJetBTagCounter = 0;

	double metPx = dataReader.metPt * std::cos(dataReader.metPhi),
		metPy = dataReader.metPt * std::sin(dataReader.metPhi);

	for (int iJet = 0; iJet < dataReader.nJet; iJet++) {
		dataReader.GetJetValues(iJet);
		if (!product.GetIsData()) { dataReader.ReadGenEntry();}

		// Missing runPeriod for data
		double correctionFactor = ak4CorrectionSet->at(era + "_" + dataType + "_L1FastJet_"    + ak4Algorithm)->evaluate({dataReader.jetArea, dataReader.jetEta, dataReader.jetPt, dataReader.rho}) *
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
		const double &smearFactor = product.GetIsData() ? 1.0 : JetProducer::SmearEnergy(dataReader, correctionFactor * dataReader.jetPt, 0.4); // FIXE maybe do something more elegant for coneSize argument

		const double &jetPtCorrected = dataReader.jetPt * correctionFactor * smearFactor;
		if (jetPtCorrected > jetPtCut && abs(dataReader.jetEta) > jetEtaCut) { continue;}

		metPx += jetPtCorrected * std::cos(dataReader.jetPhi) - jetPtCorrected * std::cos(dataReader.jetPhi);
		metPy += jetPtCorrected * std::sin(dataReader.jetPhi) - jetPtCorrected * std::sin(dataReader.jetPhi);

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

	//std::cout << "\nPt(Before Sorting):\n"; for (int iJet = 0; iJet < jetCounter; iJet++) { std::cout << iJet << "(" << product.jetPt[iJet] << ", " << product.jetEta[iJet] << "); "; } std::cout << std::endl;

	std::vector<int> indices(jetCounter);
	std::iota(indices.begin(), indices.end(), 0);
	std::stable_sort(indices.begin(), indices.end(), [&](int i1, int i2) {return product.jetPt[i1] > product.jetPt[i2];});
	for (std::array<double, 20> *jetVariable : {&product.jetPt, &product.jetEta, &product.jetPhi, &product.jetMass}) {
		SortByIndex(*jetVariable, indices, jetCounter);
	}
	for (std::array<bool, 20> *jetVariable : {&product.jetDeepCsvLooseId, &product.jetDeepCsvMediumId, &product.jetDeepCsvTightId, &product.jetDeepJetLooseId, &product.jetDeepJetMediumId, &product.jetDeepJetTightId}) {
		SortByIndex(*jetVariable, indices, jetCounter);
	}
	//std::cout << "Pt(After Sorting):\n"; for (int iJet = 0; iJet < jetCounter; iJet++) { std::cout << iJet << "(" << product.jetPt[iJet] << ", " << product.jetEta[iJet] << "); "; } std::cout << std::endl;

	// Jet Cleaning
	double deltaRMin = std::numeric_limits<double>::max();
	std::vector<int> muonIndices, electronIndices, jetRemovalIndices;;
	int nearestMuonIndex = -999, nearestElectronIndex = -999;
	for (int iMuon = 0; iMuon < product.nMuon; iMuon++) {
		if (!product.muonIsGood[iMuon]) { continue;}
		if (std::find(muonIndices.begin(), muonIndices.end(),iMuon)!=muonIndices.end()) { continue;}
		int nearestJetIndex = -999;
		for (int iJet = 0; iJet < jetCounter; iJet++) {
			if (std::find(jetRemovalIndices.begin(), jetRemovalIndices.end(),iJet)!=jetRemovalIndices.end()) { continue;}
			double deltaR = Utility::DeltaR(product.jetEta[iJet], product.jetPhi[iJet], product.muonEta[iMuon], product.muonPhi[iMuon]);
			if (deltaR < deltaRMin && deltaR < 0.4) { // FIXME change 0.4 do configureable param
				deltaRMin = deltaR;
				nearestMuonIndex = iMuon;
				nearestJetIndex = iJet;
			}
		}

		if (nearestJetIndex > 0) {
			muonIndices.push_back(iMuon);
			std::cout << "Found Jet to be removed: " << nearestJetIndex << std::endl;
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
			if (deltaR < deltaRMin && deltaR < 0.4) { // FIXME change 0.4 do configureable param
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

	int removeCounter = 0;
	for (int iRemove : jetRemovalIndices) {
		for (std::array<double, 20> *jetVariable : {&product.jetPt, &product.jetEta, &product.jetPhi, &product.jetMass}) {
			RemoveByIndex(jetVariable, iRemove - removeCounter, jetCounter);
		}
		for (std::array<bool, 20> *jetVariable : {&product.jetDeepCsvLooseId, &product.jetDeepCsvMediumId, &product.jetDeepCsvTightId, &product.jetDeepJetLooseId, &product.jetDeepJetMediumId, &product.jetDeepJetTightId}) {
			RemoveByIndex(jetVariable, iRemove - removeCounter, jetCounter);
		}
		removeCounter++;
	}
	jetCounter -= removeCounter;

	//std::cout << "Pt(After Cleaning):\n"; for (int iJet = 0; iJet < jetCounter; iJet++) { std::cout << iJet << "(" << product.jetPt[iJet] << ", " << product.jetEta[iJet] << "); "; } std::cout << std::endl;

	product.nJet = jetCounter;

	for (int iFatJet = 0; iFatJet < dataReader.nFatJet; iFatJet++) {
		//FatJetDeepTagMD_H4qvsQCD.push_back(fatJetDeepTagMDH4qvsQCD->At(iFatJet));
		//FatJetDeepTagMD_HbbvsQCD.push_back(fatJetDeepTagMDHbbvsQCD->At(iFatJet));
		//FatJetDeepTagMD_TvsQCD.push_back(fatJetDeepTagMDTvsQCD->At(iFatJet));
		//FatJetDeepTagMD_WvsQCD.push_back(fatJetDeepTagMDWvsQCD->At(iFatJet));
		//FatJetDeepTagMD_ZHbbvsQCD.push_back(fatJetDeepTagMDZHbbvsQCD->At(iFatJet));
		//FatJetDeepTagMD_ZHccvsQCD.push_back(fatJetDeepTagMDZHccvsQCD->At(iFatJet));
		//FatJetDeepTagMD_ZbbvsQCD.push_back(fatJetDeepTagMDZbbvsQCD->At(iFatJet));
		//FatJetDeepTagMD_ZvsQCD.push_back(fatJetDeepTagMDZvsQCD->At(iFatJet));
		//FatJetDeepTagMD_bbvsLiFatJetght.push_back(fatJetDeepTagMDBbvsLiFatJetght->At(iFatJet));
		//FatJetDeepTagMD_ccvsLiFatJetght.push_back(fatJetDeepTagMDCcvsLiFatJetght->At(i));
		//FatJetDeepTag_H.push_back(fatJetDeepTagH->At(i));
		//FatJetDeepTag_QCD.push_back(fatJetDeepTagQCD->At(iFatJet));
		//FatJetDeepTag_QCDothers.push_back(fatJetDeepTagQCDothers->At(iFatJet));
		//FatJetDeepTag_TvsQCD.push_back(fatJetDeepTagTvsQCD->At(iFatJet));
		//FatJetDeepTag_WvsQCD.push_back(fatJetDeepTagWvsQCD->At(iFatJet));
		//FatJetDeepTag_ZvsQCD.push_back(fatJetDeepTagZvsQCD->At(iFatJet));
	}

	product.metPt  = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
	product.metPhi = std::atan2(metPy, metPx);

	product.wBosonMinMass    = 999,
	product.wBosonMinMassPt  = -999,
	product.wBosonBestMass   = -999,
	product.wBosonBestMassPt = -999,
	product.topBestMass      = -999,
	product.topBestMassPt    = -999;
	for (int iJet1 = 0; iJet1 < product.nJet; iJet1++) {
		ROOT::Math::PtEtaPhiMVector jet1P4 = ROOT::Math::PtEtaPhiMVector(product.jetPt.at(iJet1), product.jetEta.at(iJet1), product.jetPhi.at(iJet1), product.jetMass.at(iJet1));
		if (!product.jetMediumCSVBTag.at(iJet1)) {
			for (int iJet2 = iJet1 + 1; iJet2 < nJet; iJet2++) {
				ROOT::Math::PtEtaPhiMVector jet2P4 = ROOT::Math::PtEtaPhiMVector(product.jetPt.at(iJet2), product.jetEta.at(iJet2), product.jetPhi.at(iJet2), product.jetMass.at(iJet2));
				ROOT::Math::PtEtaPhiMVector diJetP4 = jet1P4 + jet2P4;
				double diJetMass = diJetP4.M();
				if (diJetMass > 30 && diJetMass < wBosonMinMass) { wBosonMinMass = diJetMass; wBosonMinMassPt = diJetP4.Pt();}
				if (abs(diJetMass - pdgWBosonMass) < abs(wBosonBestMass - pdgWBosonMass)) {
					wBosonBestMass = diJetMass; wBosonBestMassPt = diJetP4.Pt();
				}
				for (int b = 0; b < nJet; b++) {
					if (JetMediumCSVBTag.at(b) && (Utility::DeltaR(JetEta.at(b), JetPhi.at(b), jet1P4.Eta(), jet1P4.Phi()) < 0.1 || Utility::DeltaR(JetEta.at(b), JetPhi.at(b), jet2P4.Eta(), jet2P4.Phi()) < 0.1)) {
						ROOT::Math::PtEtaPhiMVector bJetP4 = ROOT::Math::PtEtaPhiMVector(JetPt.at(j), JetEta.at(j), JetPhi.at(j), JetMass.at(j));
						ROOT::Math::PtEtaPhiMVector topP4 = diJetP4 + bJetP4;
						double topMass = topP4.M();
						if (abs(topMass - pdgTopMass) < abs(topBestMass - pdgTopMass)) { topBestMass = topMass; topBestMassPt = topP4.Pt();}
					}
				}
			}
		}
	}
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
