#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/JetProducer.h>

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>

#include <cmath>
#include <random>

JetProducer::JetProducer(const int &era, const float &ptCut, const float &etaCut, const float &deltaRCut, const bool &preVFP, const char &runPeriod):
	preVFP(preVFP),
	runPeriod(runPeriod),
	era(era),
	ptCut(ptCut),
	etaCut(etaCut),
	deltaRCut(deltaRCut)
	{
		//if (era == 2016 && preVFP) {
		//	eraSelector = std::to_string(era * 10);
		//} else {
		//	eraSelector = std::to_string(era);
		//}
	}

template <typename T>
void JetProducer::SortByIndex(T &var, std::vector<int> idx, unsigned int vectorSize) {
	T tmp(vectorSize);
	for (unsigned int i = 0; i < vectorSize; i++) {
		tmp.at(i) = var.at(idx[i]);
	}
	var = std::move(tmp);
}

void JetProducer::SetCorrector(const char &runPeriod) {
	std::vector<JetCorrectorParameters> corrVec;

	for (std::string fileName : isData? jecData[eraSelector] : jecMC[eraSelector]) {

		if (fileName.find("@") != std::string::npos) {
			fileName.replace(fileName.find("@"), 1, runEras[eraSelector][runPeriod]);
		}
		corrVec.push_back(JetCorrectorParameters(fileName));
	}

	jetCorrector = new FactorizedJetCorrector(corrVec);
}

// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
float JetProducer::CorrectEnergy(const float &pt, const float &eta, const float &rho, const float &area) {
	jetCorrector->setJetPt(pt);
	jetCorrector->setJetEta(eta);
	jetCorrector->setRho(rho);
	jetCorrector->setJetA(area);
	return jetCorrector->getCorrection();
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
float JetProducer::SmearEnergy(const float &pt, const float &eta, const float &phi, const float &rho, const float &coneSize) {
	jetParameter.setJetPt(pt).setJetEta(eta).setRho(rho);

	float reso = resolution.getResolution(jetParameter);
	float resoSF;
	if (isJERSystematic) {
		resoSF = isUp ? resolution_sf.getScaleFactor(jetParameter, Variation::UP) : resolution_sf.getScaleFactor(jetParameter, Variation::DOWN);
	} else {
		resoSF = resolution_sf.getScaleFactor(jetParameter);
	}

	float smearFactor  = 1.;
	float dR;
	float genPt, genPhi, genEta, genMass;
	//unsigned int size;
	bool isMatched = false;

	unsigned int size = 999;//genJetPt->GetSize();

	//Loop over all gen jets and find match
	for (unsigned int i = 0; i < size; i++) {
		genPt = 999;//genJetPt->At(i);
		genPhi = -999;//genJetPhi->At(i);
		genEta = -999;//genJetEta->At(i);
		genMass = -999; //genJetMass->At(i);

		dR = std::sqrt(std::pow(phi - genPhi, 2) + std::pow(eta - genEta, 2));

		//Check if jet and gen jet are matched
		if (dR < coneSize/2. && abs(pt - genPt) < 3. * reso * pt) {
			genJet = ROOT::Math::PtEtaPhiMVector(genPt, genEta, genPhi, genMass);
			isMatched = true;
			break;
		}
	}

	//If you found gen matched
	if (isMatched) {
		smearFactor = 1. + (resoSF-1) * (pt - genPt) / pt;
	} else if (resoSF > 1.) {
		std::default_random_engine generator;
		std::normal_distribution<> gaus(0, reso * std::sqrt(resoSF * resoSF - 1));
		smearFactor = 1. + gaus(generator);
	}

	//Check if direction of jet not changed
	if (pt * smearFactor < 1e-2) { smearFactor = 1e-2 / pt;}

	return smearFactor;
}

//void JetProducer::BeginJob(std::shared_ptr<TTree> tree, bool &isData, bool &doSystematics) {
//	//Set data bool
//	this->isData = isData;
//	this->doSystematics = doSystematics;
//	isUp = ((std::string)tree->GetName()).find("Up") != std::string::npos ? true : false;
//	isJERSystematic = ((std::string)tree->GetName()).find("JER") != std::string::npos ? true : false;
//	isJECSystematic = ((std::string)tree->GetName()).find("JEC") != std::string::npos ? true : false;
//
//	// Path to Correction files
//	std::string cmsswBase = std::string(std::getenv("CMSSW_BASE"));
//	std::string dataFilePath = cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/";
//	// TODO Fastsim uses its own files.. Still needs to be added
//	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
//	// https://github.com/cms-jet/JECDatabase
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	jecMC = {
//		{2016, {dataFilePath + "JEC/Summer19UL16/Summer19UL16_V7_MC/Summer19UL16_V7_MC_L1FastJet_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16/Summer19UL16_V7_MC/Summer19UL16_V7_MC_L2Relative_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16/Summer19UL16_V7_MC/Summer19UL16_V7_MC_L3Absolute_AK4PFchs.txt"}
//		},
//		{20160, {dataFilePath + "JEC/Summer19UL16APV/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16APV/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC_L2Relative_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16APV/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC_L3Absolute_AK4PFchs.txt"}
//		},
//		{2017, {dataFilePath + "JEC/Summer19UL17/Summer19UL17_V5_MC/Summer19UL17_V5_MC_L1FastJet_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL17/Summer19UL17_V5_MC/Summer19UL17_V5_MC_L2Relative_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL17/Summer19UL17_V5_MC/Summer19UL17_V5_MC_L3Absolute_AK4PFchs.txt"}
//		},
//		{2018, {dataFilePath + "JEC/Summer19UL18/Summer19UL18_V5_MC/Summer19UL18_V5_MC_L1FastJet_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL18/Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2Relative_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL18/Summer19UL18_V5_MC/Summer19UL18_V5_MC_L3Absolute_AK4PFchs.txt"}
//		}
//	};
//
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	jecData = {
//		{2016, {dataFilePath + "JEC/Summer19UL16/Summer19UL16_V7_DATA/Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16/Summer19UL16_V7_DATA/Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16/Summer19UL16_V7_DATA/Summer19UL16_RunFGH_V7_DATA_L3Absolute_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16/Summer19UL16_V7_DATA/Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs.txt"}
//		},
//		{20160, {dataFilePath + "JEC/Summer19UL16APV/Summer19UL16APV_V7_DATA/Summer19UL16APV_Run@_V7_DATA_L1FastJet_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16APV/Summer19UL16APV_V7_DATA/Summer19UL16APV_Run@_V7_DATA_L2Relative_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16APV/Summer19UL16APV_V7_DATA/Summer19UL16APV_Run@_V7_DATA_L3Absolute_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL16APV/Summer19UL16APV_V7_DATA/Summer19UL16APV_Run@_V7_DATA_L2L3Residual_AK4PFchs.txt"}
//		},
//		{2017, {dataFilePath + "JEC/Summer19UL17/Summer19UL17_V5_DATA/Summer19UL17_Run@_V5_DATA_L1FastJet_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL17/Summer19UL17_V5_DATA/Summer19UL17_Run@_V5_DATA_L2Relative_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL17/Summer19UL17_V5_DATA/Summer19UL17_Run@_V5_DATA_L3Absolute_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL17/Summer19UL17_V5_DATA/Summer19UL17_Run@_V5_DATA_L2L3Residual_AK4PFchs.txt"}
//		},
//		{2018, {dataFilePath + "JEC/Summer19UL18/Summer19UL18_V5_DATA/Summer19UL18_Run@_V5_DATA_L1FastJet_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL18/Summer19UL18_V5_DATA/Summer19UL18_Run@_V5_DATA_L2Relative_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL18/Summer19UL18_V5_DATA/Summer19UL18_Run@_V5_DATA_L3Absolute_AK4PFchs.txt",
//			dataFilePath   + "JEC/Summer19UL18/Summer19UL18_V5_DATA/Summer19UL18_Run@_V5_DATA_L2L3Residual_AK4PFchs.txt"}
//		}
//	};
//
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	jecMCUnc = {
//		{2016, dataFilePath + "JEC/Summer19UL16/Summer19UL16_V7_MC/Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt"},
//		{20160, dataFilePath + "JEC/Summer19UL16APV/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt"},
//		{2017, dataFilePath + "JEC/Summer19UL17/Summer19UL17_V5_MC/Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt"},
//		{2018, dataFilePath + "JEC/Summer19UL18/Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt"}
//	};
//
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	jecDataUnc = {
//		{2016, dataFilePath + "JEC/Summer19UL16/Summer19UL16_V7_DATA/Summer19UL16_RunFGH_V7_DATA_UncertaintySources_AK4PFchs.txt"},
//		{20160, dataFilePath + "JEC/Summer19UL16APV/Summer19UL16APV_V7_DATA/Summer19UL16APV_Run@_V7_DATA_UncertaintySources_AK4PFchs.txt"},
//		{2017, dataFilePath + "JEC/Summer19UL17/Summer19UL17_V5_DATA/Summer19UL17_Run@_V5_DATA_UncertaintySources_AK4PFchs.txt"},
//		{2018, dataFilePath + "JEC/Summer19UL18/Summer19UL18_V5_DATA/Summer19UL18_Run@_V5_DATA_UncertaintySources_AK4PFchs.txt"}
//	};
//
//	// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
//	// https://github.com/cms-jet/JRDatabase
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	jmeSF = {
//		{2016, dataFilePath + "JME/Summer20UL16/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt"},
//		{20160, dataFilePath + "JME/Summer20UL16APV/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt"},
//		{2017, dataFilePath + "JME/Summer19UL17/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_SF_AK4PFchs.txt"},
//		{2018, dataFilePath + "JME/Summer19UL18/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt"}
//	};
//
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	jmePtReso = {
//		{2016, dataFilePath + "JME/Summer20UL16/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt"},
//		{20160, dataFilePath + "JME/Summer20UL16APV/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt"},
//		{2017, dataFilePath + "JME/Summer19UL17/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt"},
//		{2018, dataFilePath + "JME/Summer19UL18/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt"}
//	};
//
//	/* TODO UL BTAG SF
//	BTag Working Points
//	2016: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP
//	2016: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP
//+	2017: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
//	2018: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
//	DeepJet=DeepFlavour > DeepCSV, but keep both for comparison maybe
//	*/
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	deepFlavourBTag = {
//		{2016, {
//				{'l', 0.0480},
//				{'m', 0.2489},
//				{'t', 0.6377}
//			}
//		},
//		{20160, {
//				{'l', 0.0508},
//				{'m', 0.2598},
//				{'t', 0.6502}
//			}
//		},
//		{2017, {
//				{'l', 0.0532},
//				{'m', 0.3040},
//				{'t', 0.7476}
//			}
//		},
//		{2018, {
//				{'l', 0.0490},
//				{'m', 0.2783},
//				{'t', 0.7100}
//			}
//		}
//	};
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	deepCSVBTag = {
//		{2016, {
//				{'l', 0.1918},
//				{'m', 0.5847},
//				{'t', 0.8767}
//			}
//		},
//		{20160, {
//				{'l', 0.2027},
//				{'m', 0.6001},
//				{'t', 0.8819}
//			}
//		},
//		{2017, {
//				{'l', 0.1355},
//				{'m', 0.4506},
//				{'t', 0.7738}
//			}
//		},
//		{2018, {
//				{'l', 0.1208},
//				{'m', 0.4168},
//				{'t', 0.7665}
//			}
//		}
//	};
//
//	//BTag Scale Factor Files
//	// 2016 = postVFP = noHIPM = noAPV; 20160 = preVFP = HIPM = APV
//	//std::string btagSFFilePath = cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/btagSF/";
//	deepCSVTagSFFile = {
//		{2016, "BTagSF/wp_deepCSV_106XUL16postVFP_v3.csv"},
//		{20160, "BTagSF/wp_deepCSV_106XUL16preVFP_v2.csv"},
//		{2017, "BTagSF/wp_deepCSV_106XUL17_v3.csv"},
//		{2018, "BTagSF/wp_deepCSV_106XUL18_v2.csv"},
//	};
//
//	deepFlavourTagSFFile = {
//		{2016, "BTagSF/wp_deepJet_106XUL16postVFP_v3.csv"},
//		{20160, "BTagSF/wp_deepJet_106XUL16preVFP_v2.csv"},
//		{2017, "BTagSF/wp_deepJet_106XUL17_v3.csv"},
//		{2018, "BTagSF/wp_deepJet_106XUL18_v2.csv"},
//	};
//	if(!isData){
//		bTagReader = {
//			// FIXME deepjet for 2016 preVPF seems to make problems... Check other years as well!!
//			{"deepflavour", new BTagCSVReader(dataFilePath + deepFlavourTagSFFile[eraSelector])},
//			{"deepcsv",     new BTagCSVReader(dataFilePath + deepCSVTagSFFile[eraSelector])},
//		};
//	}
//
//	SetCorrector(runPeriod);
//
//	//jetNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nJet");
//	//jetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_pt");
//	//jetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_phi");
//	//jetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_eta");
//	//jetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_mass");
//	//jetArea = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_area");
//	//jetRawFactor = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_rawFactor");
//	//jetCSV = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_btagDeepB");
//	//jetDF = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_btagDeepFlavB");
//	//jetRho = std::make_unique<TTreeReaderValue<float>>(*reader, "fixedGridRhoFastjetAll");
//	//metPt = std::make_unique<TTreeReaderValue<float>>(*reader, "MET_pt");
//	//metPhi = std::make_unique<TTreeReaderValue<float>>(*reader, "MET_phi");
//
//	//fatJetNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nFatJet");
//	//fatJetDeepTagMDH4qvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_H4qvsQCD");
//	//fatJetDeepTagMDHbbvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_HbbvsQCD");
//	//fatJetDeepTagMDTvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_TvsQCD");
//	//fatJetDeepTagMDWvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_WvsQCD");
//	//fatJetDeepTagMDZHbbvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_ZHbbvsQCD");
//	//fatJetDeepTagMDZHccvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_ZHccvsQCD");
//	//fatJetDeepTagMDZbbvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_ZbbvsQCD");
//	//fatJetDeepTagMDZvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_ZvsQCD");
//	//fatJetDeepTagMDBbvsLight = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_bbvsLight");
//	//fatJetDeepTagMDCcvsLight = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTagMD_ccvsLight");
//	//fatJetDeepTagH = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTag_H");
//	//fatJetDeepTagQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTag_QCD");
//	//fatJetDeepTagQCDothers = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTag_QCDothers");
//	//fatJetDeepTagTvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTag_TvsQCD");
//	//fatJetDeepTagWvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTag_WvsQCD");
//	//fatJetDeepTagZvsQCD = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_deepTag_ZvsQCD");
//
//	//if (!this->isData) {
//	//	genJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_pt");
//	//	genJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_eta");
//	//	genJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_phi");
//	//	genJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_mass");
//
//	//	jetFlavour = std::make_unique<TTreeReaderArray<int>>(*reader, "Jet_partonFlavour");
//	//	jetGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Jet_genJetIdx");
//	//}
//
//	//Set configuration for JER tools
//	std::string fileName = jmePtReso[eraSelector];
//	resolution = JME::JetResolution(fileName);
//
//	fileName = jmeSF[eraSelector];
//	//fileName.replace(fileName.find("@"), 1, &runPeriod); // TODO Make sure the runperiod doesn't matter here
//	resolution_sf = JME::JetResolutionScaleFactor(fileName);
//
//	//Set object to get JEC uncertainty
//	fileName = isData? jecDataUnc[eraSelector] : jecMCUnc[eraSelector];
//	if (fileName.find("@") != std::string::npos) {
//		fileName.replace(fileName.find("@"), 1, &runPeriod);
//	}
//
//	if (isJECSystematic) {
//		jetCorrectionUncertainty = new JetCorrectionUncertainty(JetCorrectorParameters(fileName, "Total"));
//	}
//
//	//Set Branches of output tree
//	tree->Branch("JetPt", &JetPt);
//	tree->Branch("JetEta", &JetEta);
//	tree->Branch("JetPhi", &JetPhi);
//	tree->Branch("JetMass", &JetMass);
//	tree->Branch("JetRawFactor", &JetRawFactor);
//	tree->Branch("JetCSVBTag", &JetCSVBTag);
//	tree->Branch("JetDFBTag", &JetDFBTag);
//	tree->Branch("JetLooseDFBTag", &JetLooseDFBTag);
//	tree->Branch("JetMediumDFBTag", &JetMediumDFBTag);
//	tree->Branch("JetTightDFBTag", &JetTightDFBTag);
//	tree->Branch("JetLooseCSVBTag", &JetLooseCSVBTag);
//	tree->Branch("JetMediumCSVBTag", &JetMediumCSVBTag);
//	tree->Branch("JetTightCSVBTag", &JetTightCSVBTag);
//
//	if (!isData) {
//		tree->Branch("JetLooseDFBTagSF", &JetLooseDFBTagSF);
//		tree->Branch("JetMediumDFBTagSF", &JetMediumDFBTagSF);
//		tree->Branch("JetTightDFBTagSF", &JetTightDFBTagSF);
//		tree->Branch("JetLooseCSVBTagSF", &JetLooseCSVBTagSF);
//		tree->Branch("JetMediumCSVBTagSF", &JetMediumCSVBTagSF);
//		tree->Branch("JetTightCSVBTagSF", &JetTightCSVBTagSF);
//		if (!doSystematics) {
//			tree->Branch("JetLooseDFBTagSFUp", &JetLooseDFBTagSFUp);
//			tree->Branch("JetLooseDFBTagSFDown", &JetLooseDFBTagSFDown);
//			tree->Branch("JetMediumDFBTagSFUp", &JetMediumDFBTagSFUp);
//			tree->Branch("JetMediumDFBTagSFDown", &JetMediumDFBTagSFDown);
//			tree->Branch("JetTightDFBTagSFUp", &JetTightDFBTagSFUp);
//			tree->Branch("JetTightDFBTagSFDown", &JetTightDFBTagSFDown);
//			tree->Branch("JetLooseCSVBTagSFUp", &JetLooseCSVBTagSFUp);
//			tree->Branch("JetLooseCSVBTagSFDown", &JetLooseCSVBTagSFDown);
//			tree->Branch("JetMediumCSVBTagSFUp", &JetMediumCSVBTagSFUp);
//			tree->Branch("JetMediumCSVBTagSFDown", &JetMediumCSVBTagSFDown);
//			tree->Branch("JetTightCSVBTagSFUp", &JetTightCSVBTagSFUp);
//			tree->Branch("JetTightCSVBTagSFDown", &JetTightCSVBTagSFDown);
//		}
//	}
//
//	tree->Branch("METPt", &METPt);
//	tree->Branch("METPhi", &METPhi);
//
//	tree->Branch("nJet", &nJet);
//	tree->Branch("nLooseDFBTagJet", &nLooseDFBTagJet);
//	tree->Branch("nMediumDFBTagJet", &nMediumDFBTagJet);
//	tree->Branch("nTightDFBTagJet", &nTightDFBTagJet);
//	tree->Branch("nLooseCSVBTagJet", &nLooseCSVBTagJet);
//	tree->Branch("nMediumCSVBTagJet", &nMediumCSVBTagJet);
//	tree->Branch("nTightCSVBTagJet", &nTightCSVBTagJet);
//
//	tree->Branch("FatJetDeepTagMD_H4qvsQCD", &FatJetDeepTagMD_H4qvsQCD);
//	tree->Branch("FatJetDeepTagMD_HbbvsQCD", &FatJetDeepTagMD_HbbvsQCD);
//	tree->Branch("FatJetDeepTagMD_TvsQCD", &FatJetDeepTagMD_TvsQCD);
//	tree->Branch("FatJetDeepTagMD_WvsQCD", &FatJetDeepTagMD_WvsQCD);
//	tree->Branch("FatJetDeepTagMD_ZHbbvsQCD", &FatJetDeepTagMD_ZHbbvsQCD);
//	tree->Branch("FatJetDeepTagMD_ZHccvsQCD", &FatJetDeepTagMD_ZHccvsQCD);
//	tree->Branch("FatJetDeepTagMD_ZbbvsQCD", &FatJetDeepTagMD_ZbbvsQCD);
//	tree->Branch("FatJetDeepTagMD_ZvsQCD", &FatJetDeepTagMD_ZvsQCD);
//	tree->Branch("FatJetDeepTagMD_bbvsLight", &FatJetDeepTagMD_bbvsLight);
//	tree->Branch("FatJetDeepTagMD_ccvsLight", &FatJetDeepTagMD_ccvsLight);
//	tree->Branch("FatJetDeepTag_H", &FatJetDeepTag_H);
//	tree->Branch("FatJetDeepTag_QCD", &FatJetDeepTag_QCD);
//	tree->Branch("FatJetDeepTag_QCDothers", &FatJetDeepTag_QCDothers);
//	tree->Branch("FatJetDeepTag_TvsQCD", &FatJetDeepTag_TvsQCD);
//	tree->Branch("FatJetDeepTag_WvsQCD", &FatJetDeepTag_WvsQCD);
//	tree->Branch("FatJetDeepTag_ZvsQCD", &FatJetDeepTag_ZvsQCD);
//}

void JetProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	// clear vectors for each event, otherwise this creates a memory leak
	JetPt.clear();
	JetEta.clear();
	JetPhi.clear();
	JetMass.clear();
	JetRawPt.clear();
	JetRawMass.clear();
	JetRawFactor.clear();
	JetCSVBTag.clear();
	JetDFBTag.clear();

	JetLooseDFBTagSF.clear();
	JetMediumDFBTagSF.clear();
	JetTightDFBTagSF.clear();
	JetLooseCSVBTagSF.clear();
	JetMediumCSVBTagSF.clear();
	JetTightCSVBTagSF.clear();
	JetLooseDFBTagSFUp.clear();
	JetLooseDFBTagSFDown.clear();
	JetMediumDFBTagSFUp.clear();
	JetMediumDFBTagSFDown.clear();
	JetTightDFBTagSFUp.clear();
	JetTightDFBTagSFDown.clear();
	JetLooseCSVBTagSFUp.clear();
	JetLooseCSVBTagSFDown.clear();
	JetMediumCSVBTagSFUp.clear();
	JetMediumCSVBTagSFDown.clear();
	JetTightCSVBTagSFUp.clear();
	JetTightCSVBTagSFDown.clear();

	FatJetDeepTagMD_H4qvsQCD.clear();
	FatJetDeepTagMD_HbbvsQCD.clear();
	FatJetDeepTagMD_TvsQCD.clear();
	FatJetDeepTagMD_WvsQCD.clear();
	FatJetDeepTagMD_ZHbbvsQCD.clear();
	FatJetDeepTagMD_ZHccvsQCD.clear();
	FatJetDeepTagMD_ZbbvsQCD.clear();
	FatJetDeepTagMD_ZvsQCD.clear();
	FatJetDeepTagMD_bbvsLight.clear();
	FatJetDeepTagMD_ccvsLight.clear();
	FatJetDeepTag_H.clear();
	FatJetDeepTag_QCD.clear();
	FatJetDeepTag_QCDothers.clear();
	FatJetDeepTag_TvsQCD.clear();
	FatJetDeepTag_WvsQCD.clear();
	FatJetDeepTag_ZvsQCD.clear();

	JetLooseDFBTag.clear();
	JetMediumDFBTag.clear();
	JetTightDFBTag.clear();
	JetLooseCSVBTag.clear();
	JetMediumCSVBTag.clear();
	JetTightCSVBTag.clear();

	//Initialize all variables as -999
	METPt = -999;
	METPhi = -999;
	JetRho = -999;

	//unsigned int
	nJet = 0;
	nLooseCSVBTagJet = 0;
	nMediumCSVBTagJet = 0;
	nTightCSVBTagJet = 0;
	nLooseDFBTagJet = 0;
	nMediumDFBTagJet = 0;
	nTightDFBTagJet = 0;

	//float CSVBValue = 0, DeepBValue = 0;
	float metPx = 0, metPy = 0;

	nJet = -999;//*jetNumber->Get();
	nFatJet = -999;//*fatJetNumber->Get();

	const float &metpt = -999;//*metPt->Get();
	const float &metphi = -999;//*metPhi->Get();
	metPx = metpt * std::cos(metphi);
	metPy = metpt * std::sin(metphi);

	for (unsigned int i = 0; i < nJet; i++) {
		const float &pt = -999; //jetPt->At(i);
		const float &phi = -999; //jetPhi->At(i);
		const float &eta = -999; //jetEta->At(i);
		const float &rawFactor = -999; //jetRawFactor->At(i);
		const float &rawPt = pt * (1 - rawFactor);
		const int &flavour = isData ? -999 : -999; //jetFlavour->At(i); // only used for btagSF which is only applied to MC

		//float correctionFactor = CorrectEnergy(rawPt, eta, *jetRho->Get(), -999; jetArea->At(i));
		float correctionFactor = CorrectEnergy(rawPt, eta, -999, -999);

		if (isJECSystematic) {
			jetCorrectionUncertainty->setJetPt(correctionFactor * rawPt);
			jetCorrectionUncertainty->setJetEta(eta);
			const float uncertainty = jetCorrectionUncertainty->getUncertainty(isUp);
			correctionFactor *= isUp ? (1 + uncertainty) : (1 - uncertainty);
		}

		//const float &smearFactor = isData ? 1.0 : JetProducer::SmearEnergy(rawPt, eta, phi, *jetRho->Get(), -999; //jetArea->At(i));
		const float &smearFactor = isData ? 1.0 : JetProducer::SmearEnergy(rawPt, eta, phi, -999, -999); //jetArea->At(i));

		const float &correctedPt = smearFactor * correctionFactor * rawPt;

		if (correctedPt > ptCut && abs(eta) < etaCut) {
			const float &mass = -999; //jetMass->At(i);
			const float &rawMass = mass * (1 - rawFactor);
			const float &correctedMass = smearFactor * correctionFactor * rawMass;

			const float &csvBTagValue = -999; //jetCSV->At(i);
			JetCSVBTag.push_back(csvBTagValue);
			if (csvBTagValue > deepCSVBTag[eraSelector]['t']) {
				JetLooseCSVBTag.push_back(true);
				JetMediumCSVBTag.push_back(true);
				JetTightCSVBTag.push_back(true);

				if (!isData) {

					JetLooseCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 0, std::abs(flavour)));
					JetMediumCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 1, std::abs(flavour)));
					JetTightCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 2, std::abs(flavour)));

					if (!doSystematics) {
						JetLooseCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 0, std::abs(flavour)));
						JetMediumCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 1, std::abs(flavour)));
						JetTightCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 2, std::abs(flavour)));

						JetLooseCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 0, std::abs(flavour)));
						JetMediumCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 1, std::abs(flavour)));
						JetTightCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 2, std::abs(flavour)));
					}
				}
			} else if (csvBTagValue > deepCSVBTag[eraSelector]['m']) {
				JetLooseCSVBTag.push_back(true);
				JetMediumCSVBTag.push_back(true);
				JetTightCSVBTag.push_back(false);

				if (!isData) {
					JetLooseCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 0, std::abs(flavour)));
					JetMediumCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 1, std::abs(flavour)));
					JetTightCSVBTagSF.push_back(0);

					if(!doSystematics) {
						JetLooseCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 0, std::abs(flavour)));
						JetMediumCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 1, std::abs(flavour)));
						JetTightCSVBTagSFUp.push_back(0);

						JetLooseCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 0, std::abs(flavour)));
						JetMediumCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 1, std::abs(flavour)));
						JetTightCSVBTagSFDown.push_back(0);
					}
				}
			} else if (csvBTagValue > deepCSVBTag[eraSelector]['l']) {
				JetLooseCSVBTag.push_back(true);
				JetMediumCSVBTag.push_back(false);
				JetTightCSVBTag.push_back(false);

				if (!isData) {
					JetLooseCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 0, std::abs(flavour)));
					JetMediumCSVBTagSF.push_back(0);
					JetTightCSVBTagSF.push_back(0);

					if(!doSystematics) {
						JetLooseCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 0, std::abs(flavour)));
						JetMediumCSVBTagSFUp.push_back(0);
						JetTightCSVBTagSFUp.push_back(0);

						JetLooseCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 0, std::abs(flavour)));
						JetMediumCSVBTagSFDown.push_back(0);
						JetTightCSVBTagSFDown.push_back(0);
					}
				}
			} else {
				JetLooseCSVBTag.push_back(false);
				JetMediumCSVBTag.push_back(false);
				JetTightCSVBTag.push_back(false);

				if (!isData) {
					JetLooseCSVBTagSF.push_back(0);
					JetMediumCSVBTagSF.push_back(0);
					JetTightCSVBTagSF.push_back(0);

					if(!doSystematics) {
						JetLooseCSVBTagSFUp.push_back(0);
						JetMediumCSVBTagSFUp.push_back(0);
						JetTightCSVBTagSFUp.push_back(0);

						JetLooseCSVBTagSFDown.push_back(0);
						JetMediumCSVBTagSFDown.push_back(0);
						JetTightCSVBTagSFDown.push_back(0);
					}
				}
			}

			const float &dfBTagValue = -999; //jetDF->At(i);
			JetDFBTag.push_back(dfBTagValue);
			if (dfBTagValue > deepFlavourBTag[eraSelector]['t']) {
				JetLooseDFBTag.push_back(true);
				JetMediumDFBTag.push_back(true);
				JetTightDFBTag.push_back(true);

				if (!isData) {
					JetLooseDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 0, std::abs(flavour)));
					JetMediumDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 1, std::abs(flavour)));
					JetTightDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 2, std::abs(flavour)));

					if (!doSystematics) {
						JetLooseDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 0, std::abs(flavour)));
						JetMediumDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 1, std::abs(flavour)));
						JetTightDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 2, std::abs(flavour)));

						JetLooseDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 0, std::abs(flavour)));
						JetMediumDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 1, std::abs(flavour)));
						JetTightDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 2, std::abs(flavour)));
					}
				}
			} else if (dfBTagValue > deepFlavourBTag[eraSelector]['m']) {
				JetLooseDFBTag.push_back(true);
				JetMediumDFBTag.push_back(true);
				JetTightDFBTag.push_back(false);

				if (!isData) {
					JetLooseDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 0, std::abs(flavour)));
					JetMediumDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 1, std::abs(flavour)));
					JetTightDFBTagSF.push_back(0);

					if (!doSystematics) {
						JetLooseDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 0, std::abs(flavour)));
						JetMediumDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 1, std::abs(flavour)));
						JetTightDFBTagSFUp.push_back(0);

						JetLooseDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 0, std::abs(flavour)));
						JetMediumDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 1, std::abs(flavour)));
						JetTightDFBTagSFDown.push_back(0);
					}
				}
			} else if (dfBTagValue > deepFlavourBTag[eraSelector]['l']) {
				JetLooseDFBTag.push_back(true);
				JetMediumDFBTag.push_back(false);
				JetTightDFBTag.push_back(false);

				if (!isData) {
					JetLooseDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 0, std::abs(flavour)));
					JetMediumDFBTagSF.push_back(0);
					JetTightDFBTagSF.push_back(0);

					if (!doSystematics) {
						JetLooseDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 0, std::abs(flavour)));
						JetMediumDFBTagSFUp.push_back(0);
						JetTightDFBTagSFUp.push_back(0);

						JetLooseDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 0, std::abs(flavour)));
						JetMediumDFBTagSFDown.push_back(0);
						JetTightDFBTagSFDown.push_back(0);
					}
				}
			} else {
				JetLooseDFBTag.push_back(false);
				JetMediumDFBTag.push_back(false);
				JetTightDFBTag.push_back(false);

				if (!isData) {
					JetLooseDFBTagSF.push_back(0);
					JetMediumDFBTagSF.push_back(0);
					JetTightDFBTagSF.push_back(0);

					if (!doSystematics) {
						JetLooseDFBTagSFUp.push_back(0);
						JetMediumDFBTagSFUp.push_back(0);
						JetTightDFBTagSFUp.push_back(0);

						JetLooseDFBTagSFDown.push_back(0);
						JetMediumDFBTagSFDown.push_back(0);
						JetTightDFBTagSFDown.push_back(0);
					}
				}
			}

			//Jet four momentum components
			JetPt.push_back(correctedPt); JetEta.push_back(eta); JetPhi.push_back(phi); JetMass.push_back(correctedMass); JetRawPt.push_back(rawPt); JetRawMass.push_back(rawPt), JetRawFactor.push_back(rawFactor);

			//JECCorrection.push_back(correctionFactor);
			metPx += pt * std::cos(phi) - JetPt.back() * std::cos(phi);
			metPy += pt * std::sin(phi) - JetPt.back() * std::sin(phi);
		}
	}

	for (unsigned int i = 0; i < nFatJet; i++) {
		//FatJetDeepTagMD_H4qvsQCD.push_back(fatJetDeepTagMDH4qvsQCD->At(i));
		//FatJetDeepTagMD_HbbvsQCD.push_back(fatJetDeepTagMDHbbvsQCD->At(i));
		//FatJetDeepTagMD_TvsQCD.push_back(fatJetDeepTagMDTvsQCD->At(i));
		//FatJetDeepTagMD_WvsQCD.push_back(fatJetDeepTagMDWvsQCD->At(i));
		//FatJetDeepTagMD_ZHbbvsQCD.push_back(fatJetDeepTagMDZHbbvsQCD->At(i));
		//FatJetDeepTagMD_ZHccvsQCD.push_back(fatJetDeepTagMDZHccvsQCD->At(i));
		//FatJetDeepTagMD_ZbbvsQCD.push_back(fatJetDeepTagMDZbbvsQCD->At(i));
		//FatJetDeepTagMD_ZvsQCD.push_back(fatJetDeepTagMDZvsQCD->At(i));
		//FatJetDeepTagMD_bbvsLight.push_back(fatJetDeepTagMDBbvsLight->At(i));
		//FatJetDeepTagMD_ccvsLight.push_back(fatJetDeepTagMDCcvsLight->At(i));
		//FatJetDeepTag_H.push_back(fatJetDeepTagH->At(i));
		//FatJetDeepTag_QCD.push_back(fatJetDeepTagQCD->At(i));
		//FatJetDeepTag_QCDothers.push_back(fatJetDeepTagQCDothers->At(i));
		//FatJetDeepTag_TvsQCD.push_back(fatJetDeepTagTvsQCD->At(i));
		//FatJetDeepTag_WvsQCD.push_back(fatJetDeepTagWvsQCD->At(i));
		//FatJetDeepTag_ZvsQCD.push_back(fatJetDeepTagZvsQCD->At(i));
	}

	//Cleanup, remove at most 1 Jet (i.e. per lepton)
	float deltaRMin = 999;
	int nearestJetIndex = -1;
	for (unsigned int i = 0; i < JetPt.size(); i++) {
		float deltaR = Utility::DeltaR(JetEta.at(i), JetPhi.at(i), product.leptonEta, product.leptonPhi);
		if (deltaR < deltaRMin) {
			deltaRMin = deltaR;
			nearestJetIndex = i;
		}
	}
	if (deltaRMin < deltaRCut && nearestJetIndex > 0) { // TODO check if this is the correct cut
		JetPt.erase(JetPt.begin() + nearestJetIndex);
		JetEta.erase(JetEta.begin() + nearestJetIndex);
		JetPhi.erase(JetPhi.begin() + nearestJetIndex);
		JetMass.erase(JetMass.begin() + nearestJetIndex);
		JetRawFactor.erase(JetRawFactor.begin() + nearestJetIndex);

		JetLooseCSVBTag.erase(JetLooseCSVBTag.begin() + nearestJetIndex);
		JetMediumCSVBTag.erase(JetMediumCSVBTag.begin() + nearestJetIndex);
		JetTightCSVBTag.erase(JetTightCSVBTag.begin() + nearestJetIndex);
		JetLooseDFBTag.erase(JetLooseDFBTag.begin() + nearestJetIndex);
		JetMediumDFBTag.erase(JetMediumDFBTag.begin() + nearestJetIndex);
		JetTightDFBTag.erase(JetTightDFBTag.begin() + nearestJetIndex);

		if (!isData) {
			JetLooseCSVBTagSF.erase(JetLooseCSVBTagSF.begin() + nearestJetIndex);
			JetMediumCSVBTagSF.erase(JetMediumCSVBTagSF.begin() + nearestJetIndex);
			JetTightCSVBTagSF.erase(JetTightCSVBTagSF.begin() + nearestJetIndex);
			JetLooseDFBTagSF.erase(JetLooseDFBTagSF.begin() + nearestJetIndex);
			JetMediumDFBTagSF.erase(JetMediumDFBTagSF.begin() + nearestJetIndex);
			JetTightDFBTagSF.erase(JetTightDFBTagSF.begin() + nearestJetIndex);

			if (!doSystematics) {
				JetLooseCSVBTagSFUp.erase(JetLooseCSVBTagSFUp.begin() + nearestJetIndex);
				JetMediumCSVBTagSFUp.erase(JetMediumCSVBTagSFUp.begin() + nearestJetIndex);
				JetTightCSVBTagSFUp.erase(JetTightCSVBTagSFUp.begin() + nearestJetIndex);
				JetLooseDFBTagSFUp.erase(JetLooseDFBTagSFUp.begin() + nearestJetIndex);
				JetMediumDFBTagSFUp.erase(JetMediumDFBTagSFUp.begin() + nearestJetIndex);
				JetTightDFBTagSFUp.erase(JetTightDFBTagSFUp.begin() + nearestJetIndex);

				JetLooseCSVBTagSFDown.erase(JetLooseCSVBTagSFDown.begin() + nearestJetIndex);
				JetMediumCSVBTagSFDown.erase(JetMediumCSVBTagSFDown.begin() + nearestJetIndex);
				JetTightCSVBTagSFDown.erase(JetTightCSVBTagSFDown.begin() + nearestJetIndex);
				JetLooseDFBTagSFDown.erase(JetLooseDFBTagSFDown.begin() + nearestJetIndex);
				JetMediumDFBTagSFDown.erase(JetMediumDFBTagSFDown.begin() + nearestJetIndex);
				JetTightDFBTagSFDown.erase(JetTightDFBTagSFDown.begin() + nearestJetIndex);
			}
		}
	}

	nJet = JetPt.size();

	for (unsigned int i = 0; i < nJet; i++) {
		if (JetTightCSVBTag[i]) { nLooseCSVBTagJet++; nMediumCSVBTagJet++; nTightCSVBTagJet++;}
		else if (JetMediumCSVBTag[i]) { nLooseCSVBTagJet++; nMediumCSVBTagJet++;}
		else if (JetLooseCSVBTag[i]) { nLooseCSVBTagJet++;}

		if (JetTightDFBTag[i]) { nLooseDFBTagJet++; nMediumDFBTagJet++; nTightDFBTagJet++;}
		else if (JetMediumDFBTag[i]) { nLooseDFBTagJet++; nMediumDFBTagJet++;}
		else if (JetLooseDFBTag[i]) { nLooseDFBTagJet++;}
	}

	METPt = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
	METPhi = std::atan2(metPy, metPx);

	//Sort Vectors according to JetPt, since the correction could have changed the order
	std::vector<int> idx(nJet);
	std::iota(idx.begin(), idx.end(), 0);
	std::stable_sort(idx.begin(), idx.end(), [&](int i1, int i2) {return JetPt[i1] > JetPt[i2];});

	if (nJet != 0) {
		SortByIndex<std::vector<float>>(JetPt, idx, nJet);
		SortByIndex<std::vector<float>>(JetPhi, idx, nJet);
		SortByIndex<std::vector<float>>(JetEta, idx, nJet);
		SortByIndex<std::vector<float>>(JetMass, idx, nJet);
		SortByIndex<std::vector<float>>(JetRawFactor, idx, nJet);

		SortByIndex<std::vector<bool>>(JetLooseCSVBTag, idx, nJet);
		SortByIndex<std::vector<bool>>(JetMediumCSVBTag, idx, nJet);
		SortByIndex<std::vector<bool>>(JetTightCSVBTag, idx, nJet);
		SortByIndex<std::vector<bool>>(JetLooseDFBTag, idx, nJet);
		SortByIndex<std::vector<bool>>(JetMediumDFBTag, idx, nJet);
		SortByIndex<std::vector<bool>>(JetTightDFBTag, idx, nJet);

		if (!isData) {
			SortByIndex<std::vector<float>>(JetLooseCSVBTagSF, idx, nJet);
			SortByIndex<std::vector<float>>(JetMediumCSVBTagSF, idx, nJet);
			SortByIndex<std::vector<float>>(JetTightCSVBTagSF, idx, nJet);
			SortByIndex<std::vector<float>>(JetLooseDFBTagSF, idx, nJet);
			SortByIndex<std::vector<float>>(JetMediumDFBTagSF, idx, nJet);
			SortByIndex<std::vector<float>>(JetTightDFBTagSF, idx, nJet);

			if (!doSystematics) {
				SortByIndex<std::vector<float>>(JetLooseCSVBTagSFUp, idx, nJet);
				SortByIndex<std::vector<float>>(JetMediumCSVBTagSFUp, idx, nJet);
				SortByIndex<std::vector<float>>(JetTightCSVBTagSFUp, idx, nJet);
				SortByIndex<std::vector<float>>(JetLooseDFBTagSFUp, idx, nJet);
				SortByIndex<std::vector<float>>(JetMediumDFBTagSFUp, idx, nJet);
				SortByIndex<std::vector<float>>(JetTightDFBTagSFUp, idx, nJet);

				SortByIndex<std::vector<float>>(JetLooseCSVBTagSFDown, idx, nJet);
				SortByIndex<std::vector<float>>(JetMediumCSVBTagSFDown, idx, nJet);
				SortByIndex<std::vector<float>>(JetTightCSVBTagSFDown, idx, nJet);
				SortByIndex<std::vector<float>>(JetLooseDFBTagSFDown, idx, nJet);
				SortByIndex<std::vector<float>>(JetMediumDFBTagSFDown, idx, nJet);
				SortByIndex<std::vector<float>>(JetTightDFBTagSFDown, idx, nJet);
			}
		}
	}

	minMWjj = 999;
	minMWjjPt = -999;
	bestMWjj = -999;
	bestMWjjPt = -999;
	bestMTop = -999;
	bestMTopPt = -999;
	for (unsigned int i = 0; i < nJet; i++) {
		ROOT::Math::PtEtaPhiMVector jet1P4 = ROOT::Math::PtEtaPhiMVector(JetPt.at(i), JetEta.at(i), JetPhi.at(i), JetMass.at(i));
		if (!JetMediumCSVBTag.at(i)) {
			for (unsigned int j = i+1; j < nJet; j++) {
				ROOT::Math::PtEtaPhiMVector jet2P4 = ROOT::Math::PtEtaPhiMVector(JetPt.at(j), JetEta.at(j), JetPhi.at(j), JetMass.at(j));
				ROOT::Math::PtEtaPhiMVector diJetP4 = jet1P4 + jet2P4;
				float mjj = diJetP4.M();
				if (mjj > 30 && mjj < minMWjj) { minMWjj = mjj; minMWjjPt = diJetP4.Pt();}
				if (abs(mjj - 80.4) < abs(bestMWjj - 80.4)) {
					bestMWjj = mjj; bestMWjjPt = diJetP4.Pt();
				}
				for (unsigned int b = 0; b < nJet; b++) {
					if (JetMediumCSVBTag.at(b) && (Utility::DeltaR(JetEta.at(b), JetPhi.at(b), jet1P4.Eta(), jet1P4.Phi()) < 0.1 || Utility::DeltaR(JetEta.at(b), JetPhi.at(b), jet2P4.Eta(), jet2P4.Phi()) < 0.1)) {
						ROOT::Math::PtEtaPhiMVector bJetP4 = ROOT::Math::PtEtaPhiMVector(JetPt.at(j), JetEta.at(j), JetPhi.at(j), JetMass.at(j));
						ROOT::Math::PtEtaPhiMVector topP4 = diJetP4 + bJetP4;
						float topMass = topP4.M();
						if (abs(topMass - 172) < abs(bestMTop - 172)) { bestMTop = topMass; bestMTopPt = topP4.Pt();}
					}
				}
			}
		}
	}

	//Store values in product to calculate high level variables
	product.nJet = nJet;
	product.jetPt = JetPt;
	product.jetPhi = JetPhi;
	product.jetEta = JetEta;
	product.jetMass = JetMass;

	product.metPt = METPt;
	product.metPhi = METPhi;

	product.nLooseCSVBTagJet = nLooseCSVBTagJet;
	product.nMediumCSVBTagJet = nMediumCSVBTagJet;
	product.nTightCSVBTagJet = nTightCSVBTagJet;

	product.nLooseDFBTagJet = nLooseDFBTagJet;
	product.nMediumDFBTagJet = nMediumDFBTagJet;
	product.nTightDFBTagJet = nTightDFBTagJet;
	product.jetMediumCSVBTag = JetMediumCSVBTag;
	product.jetMediumDFBTag = JetMediumDFBTag;

	//if (nJet!=0) {
	//	cutflow.hist->Fill("Jets accepted", cutflow.weight);
	//} else {
	//	cutflow.passed = false;
	//}
}

void JetProducer::EndJob(TFile &file) {
	delete jetCorrector;
	delete bTagReader["deepcsv"];
	delete bTagReader["deepflavour"];
	if (doSystematics) {
		delete jetCorrectionUncertainty;
	}
}
