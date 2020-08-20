#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/JetProducer.h>

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>

#include <cmath>
#include <random>

JetProducer::JetProducer(const int& era, const float& ptCut, const float& etaCut, const float& deltaRCut, const char& runPeriod, TTreeReader& reader):
	BaseProducer(&reader),
	runPeriod(runPeriod),
	era(era),
	ptCut(ptCut),
	etaCut(etaCut),
	deltaRCut(deltaRCut)
	{}

template <typename T>
void JetProducer::SortByIndex(T& var, std::vector<int> idx) {
	//std::vector<float> tmp(nJet);
	T tmp(nJet);
	for (unsigned int i = 0; i < nJet; i++) {
		tmp.at(i) = var.at(idx[i]);
	}
	var = std::move(tmp);
}

void JetProducer::SetCorrector(const JetType& type, const char& runPeriod) {
	std::vector<JetCorrectorParameters> corrVec;

	for (std::string fileName : isData? jecData[era] : jecMC[era]) {
		if (fileName.find("@") != std::string::npos) {
			fileName.replace(fileName.find("@"), 1, runEras[era][runPeriod]);
		}

		fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4" : "AK8");

		corrVec.push_back(JetCorrectorParameters(fileName));
	}

	jetCorrector[type] = new FactorizedJetCorrector(corrVec);
}

//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
float JetProducer::CorrectEnergy(const float& pt, const float& eta, const float& rho, const float& area, const JetType &type) {
	jetCorrector[type]->setJetPt(pt);
	jetCorrector[type]->setJetEta(eta);
	jetCorrector[type]->setRho(rho);
	jetCorrector[type]->setJetA(area);
	return jetCorrector[type]->getCorrection();
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
//https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
std::map<char, float> JetProducer::SmearEnergy(const float& pt, const float& eta, const float& phi, const float& rho, const float& coneSize, const JetType& type) {
	jetParameter.setJetPt(pt).setJetEta(eta).setRho(rho);

	float reso = resolution[type].getResolution(jetParameter);
	float resoSF;
	float resoSFUp;
	float resoSFDown;
	//if (!isJERsyst) resoSF = resolution_sf[type].getScaleFactor(jetParameter);
	//else resoSF = resolution_sf[type].getScaleFactor(jetParameter, isUp ? Variation::UP : Variation::DOWN);
	resoSF     = resolution_sf[type].getScaleFactor(jetParameter);
	resoSFUp   = resolution_sf[type].getScaleFactor(jetParameter, Variation::UP);
	resoSFDown = resolution_sf[type].getScaleFactor(jetParameter, Variation::DOWN);

	float smearFactor  = 1.;
	float smearFactorUp = 1.;
	float smearFactorDown = 1.;
	float dR;
	float genPt, genPhi, genEta, genMass;
	//unsigned int size;
	bool isMatched = false;

	unsigned int size = (type == AK4) ? genJetPt->GetSize() : genFatJetPt->GetSize();

	//Loop over all gen jets and find match
	for (unsigned int i = 0; i < size; i++) {
		genPt = (type == AK4) ? genJetPt->At(i) : genFatJetPt->At(i);
		genPhi = (type == AK4) ? genJetPhi->At(i) : genFatJetPhi->At(i);
		genEta = (type == AK4) ? genJetEta->At(i) : genFatJetEta->At(i);
		genMass = (type == AK4) ? genJetMass->At(i) : genFatJetMass->At(i);

		dR = std::sqrt(std::pow(phi - genPhi, 2) + std::pow(eta - genEta, 2));

		//Check if jet and gen jet are matched
		if (dR < coneSize/2. and abs(pt - genPt) < 3. * reso * pt) {
			genJet[type] = ROOT::Math::PtEtaPhiMVector(genPt, genEta, genPhi, genMass);
			isMatched = true;
			break;
		}
	}

	//If you found gen matched
	if (isMatched) {
		smearFactor     = 1. + (resoSF-1) * (pt - genPt) / pt;
		smearFactorUp   = 1. + (resoSFUp-1) * (pt - genPt) / pt;
		smearFactorDown = 1. + (resoSFDown-1) * (pt - genPt) / pt;
	} else if (resoSF > 1.) {
		std::default_random_engine generator;
		std::normal_distribution<> gaus(0, reso * std::sqrt(resoSF * resoSF - 1));
		std::normal_distribution<> gausUp(0, reso * std::sqrt(resoSFUp * resoSFUp - 1));
		std::normal_distribution<> gausDown(0, reso * std::sqrt(resoSFDown * resoSFDown - 1));
		smearFactor     = 1. + gaus(generator);
		smearFactorUp   = 1. + gausUp(generator);
		smearFactorDown = 1. + gausDown(generator);
	}

	//Check if direction of jet not changed
	if (pt * smearFactor     < 1e-2) { smearFactor     = 1e-2 / pt;}
	if (pt * smearFactorUp   < 1e-2) { smearFactorUp   = 1e-2 / pt;}
	if (pt * smearFactorDown < 1e-2) { smearFactorDown = 1e-2 / pt;}

	return {{'c', smearFactor}, {'u', smearFactorUp}, {'d', smearFactorDown}};
}

void JetProducer::BeginJob(TTree* tree, bool &isData) {
	//Set data bool
	this->isData = isData;

	// Path to Correction files
	std::string cmsswBase = std::string(std::getenv("CMSSW_BASE"));
	std::string jmeFilePath = cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/jme/";
	jecMC = {
		{2016, {jmeFilePath + "Summer16/Summer16_07Aug2017_V11_MC_L1FastJet_&PFchs.txt",
			jmeFilePath + "Summer16/Summer16_07Aug2017_V11_MC_L2Relative_&PFchs.txt",
			jmeFilePath + "Summer16/Summer16_07Aug2017_V11_MC_L3Absolute_&PFchs.txt"}
		},
		{2017, {jmeFilePath + "Fall17/Fall17_17Nov2017_V32_MC_L1FastJet_&PFchs.txt",
			jmeFilePath + "Fall17/Fall17_17Nov2017_V32_MC_L2Relative_&PFchs.txt",
			jmeFilePath + "Fall17/Fall17_17Nov2017_V32_MC_L3Absolute_&PFchs.txt"}
		},
		{2018, {jmeFilePath + "Autumn18/Autumn18_V19_MC_L1FastJet_&PFchs.txt",
			jmeFilePath + "Autumn18/Autumn18_V19_MC_L2Relative_&PFchs.txt",
			jmeFilePath + "Autumn18/Autumn18_V19_MC_L3Absolute_&PFchs.txt"}
		},
	};

	jecData = {
		{2016, {jmeFilePath + "Summer16/Summer16_07Aug2017@_V11_DATA_L1FastJet_&PFchs.txt",
			jmeFilePath + "Summer16/Summer16_07Aug2017@_V11_DATA_L2Relative_&PFchs.txt",
			jmeFilePath + "Summer16/Summer16_07Aug2017@_V11_DATA_L3Absolute_&PFchs.txt",
			jmeFilePath + "Summer16/Summer16_07Aug2017@_V11_DATA_L2L3Residual_&PFchs.txt"}
		},
		{2017, {jmeFilePath + "Fall17/Fall17_17Nov2017@_V32_DATA_L1FastJet_&PFchs.txt",
			jmeFilePath + "Fall17/Fall17_17Nov2017@_V32_DATA_L2Relative_&PFchs.txt",
			jmeFilePath + "Fall17/Fall17_17Nov2017@_V32_DATA_L3Absolute_&PFchs.txt",
			jmeFilePath + "Fall17/Fall17_17Nov2017@_V32_DATA_L2L3Residual_&PFchs.txt"}
		},
		{2018, {jmeFilePath + "Autumn18/Autumn18_Run@_V9_DATA_L1FastJet_&PFchs.txt",
			jmeFilePath + "Autumn18/Autumn18_Run@_V9_DATA_L2Relative_&PFchs.txt",
			jmeFilePath + "Autumn18/Autumn18_Run@_V9_DATA_L3Absolute_&PFchs.txt",
			jmeFilePath + "Autumn18/Autumn18_Run@_V9_DATA_L2L3Residual_&PFchs.txt"}
		}
	};

	jecUnc = {
		{2016, jmeFilePath + "Summer16/Summer16_07Aug2017_V11_MC_UncertaintySources_&PFchs.txt"},
		{2017, jmeFilePath + "Fall17/Fall17_17Nov2017_V32_MC_UncertaintySources_&PFchs.txt"},
		{2018, jmeFilePath + "Autumn18/Autumn18_V19_MC_UncertaintySources_&PFchs.txt"},
	};

	jmeSF = {
		{2016, jmeFilePath + "Summer16/Summer16_25nsV1_MC_SF_&PFchs.txt"},
		{2017, jmeFilePath + "Fall17/Fall17_V3_MC_SF_&PFchs.txt"},
		{2018, jmeFilePath + "Autumn18/Autumn18_V7_MC_SF_&PFchs.txt"},
	};

	jmePtReso = {
		{2016, jmeFilePath + "Summer16/Summer16_25nsV1_MC_PtResolution_&PFchs.txt"},
		{2017, jmeFilePath + "Fall17/Fall17_V3_MC_PtResolution_&PFchs.txt"},
		{2018, jmeFilePath + "Autumn18/Autumn18_V7_MC_PtResolution_&PFchs.txt"},
	};

	jecFastSim = {
		{2016, {jmeFilePath + "Summer16/Summer16_FastSimV1_MC_L1FastJet_&PFchs.txt",
			jmeFilePath + "Summer16/Summer16_FastSimV1_MC_L2Relative_&PFchs.txt",
			jmeFilePath + "Summer16/Summer16_FastSimV1_MC_L3Absolute_&PFchs.txt"}
		},
		{2017, {jmeFilePath + "Fall17/Fall17_FastSimV1_MC_L1FastJet_&PFchs.txt",
			jmeFilePath + "Fall17/Fall17_FastSimV1_MC_L2Relative_&PFchs.txt",
			jmeFilePath + "Fall17/Fall17_FastSimV1_MC_L3Absolute_&PFchs.txt"}
		},
		{2018, {jmeFilePath + "Autumn18/Autumn18_FastSimV1_MC_L1FastJet_&PFchs.txt",
			jmeFilePath + "Autumn18/Autumn18_FastSimV1_MC_L2Relative_&PFchs.txt",
			jmeFilePath + "Autumn18/Autumn18_FastSimV1_MC_L3Absolute_&PFchs.txt"}
		}
	};// TODO also do uncertainties for FastSim

	/*
	BTag Working Points
	2016: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
+	2017: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
	2018: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
	*/
	deepFlavourBTag = {
		{2016, {
				{'l', 0.0614},
				{'m', 0.3093},
				{'t', 0.7221}
			}
		},
		{2017, {
				{'l', 0.0521},
				{'m', 0.3033},
				{'t', 0.7489}
			}
		},
		{2018, {
				{'l', 0.0494},
				{'m', 0.2770},
				{'t', 0.7264}
			}
		}
	};
	deepCSVBTag = {
		{2016, {
				{'l', 0.2219},
				{'m', 0.6324},
				{'t', 0.8958}
			}
		},
		{2017, {
				{'l', 0.1522},
				{'m', 0.4941},
				{'t', 0.8001}
			}
		},
		{2018, {
				{'l', 0.1241},
				{'m', 0.4184},
				{'t', 0.7527}
			}
		}
	};

	//BTag Scale Factor Files
	std::string btagSFFilePath = cmsswBase + "/src/PhysicsTools/NanoAODTools/data/btagSF/";
	deepCSVTagSFFile = {
		{2016, "DeepCSV_2016LegacySF_V1.csv"},
		{2017, "DeepCSV_94XSF_V4_B_F.csv"},
		{2018, "DeepCSV_102XSF_V1.csv"},
	};

	deepFlavourTagSFFile = {
		{2016, "DeepJet_2016LegacySF_V1.csv"},
		{2017, "DeepFlavour_94XSF_V3_B_F.csv"},
		{2018, "DeepJet_102XSF_V1.csv"},
	};

	if(!isData){
		bTagReader = {
			{"deepcsv", new BTagCSVReader(btagSFFilePath + deepCSVTagSFFile[era])},
			{"deepflavour", new BTagCSVReader(btagSFFilePath + deepFlavourTagSFFile[era])},
		};
	}

	//Set TTreeReader for genpart and trigger obj from BaseProducer
	SetCollection(this->isData);
	jetNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nJet");
	jetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_pt");
	jetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_phi");
	jetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_eta");
	jetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_mass");
	jetArea = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_area");
	jetRawFactor = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_rawFactor");
	jetCSV = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_btagDeepB");
	jetDF = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_btagDeepFlavB");
	jetRho = std::make_unique<TTreeReaderValue<float>>(*reader, "fixedGridRhoFastjetAll");
	metPt = std::make_unique<TTreeReaderValue<float>>(*reader, "MET_pt");
	metPhi = std::make_unique<TTreeReaderValue<float>>(*reader, "MET_phi");

	fatJetNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nFatJet");
	fatJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_pt");
	fatJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_eta");
	fatJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_phi");
	fatJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_mass");
	fatJetArea = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_area");
	fatJetCSV = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_btagDeepB");
	//fatJetDF = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_btagDeepFlavB");
	fatJetTau1 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau1");
	fatJetTau2 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau2");
	fatJetTau3 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau3");

	if (!this->isData) {
		genJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_pt");
		genJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_eta");
		genJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_phi");
		genJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_mass");

		genFatJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_pt");
		genFatJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_eta");
		genFatJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_phi");
		genFatJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_mass");

		jetFlavour = std::make_unique<TTreeReaderArray<int>>(*reader, "Jet_partonFlavour");
		jetGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Jet_genJetIdx");
	}

	for (JetType type : {AK4, AK8}) {
		//Set configuration for JER tools
		std::string fileName = jmePtReso[era];
		fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4" : "AK8");
		resolution[type] = JME::JetResolution(fileName);

		fileName = jmeSF[era];
		fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4" : "AK8");
		resolution_sf[type] = JME::JetResolutionScaleFactor(fileName);

		//Set object to get JEC uncertainty
		fileName = jecUnc[era];
		fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4" : "AK8");
		jetCorrectionUncertainty[type] = new JetCorrectionUncertainty(JetCorrectorParameters(fileName, "Total"));
	}

	//Set Branches of output tree
	tree->Branch("JetPt", &JetPt);
	tree->Branch("JetEta", &JetEta);
	tree->Branch("JetPhi", &JetPhi);
	tree->Branch("JetMass", &JetMass);
	tree->Branch("JetRawFactor", &JetRawFactor);
	tree->Branch("JetPt_jerUp", &JetPt_jerUp);
	tree->Branch("JetMass_jerUp", &JetMass_jerUp);
	tree->Branch("JetPt_jerDown", &JetPt_jerDown);
	tree->Branch("JetMass_jerDown", &JetMass_jerDown);
	tree->Branch("JetPt_jecUp", &JetPt_jecUp);
	tree->Branch("JetMass_jecUp", &JetMass_jecUp);
	tree->Branch("JetPt_jecDown", &JetPt_jecDown);
	tree->Branch("JetMass_jecDown", &JetMass_jecDown);
	tree->Branch("JetCSVBTag", &JetCSVBTag);
	tree->Branch("JetDFBTag", &JetDFBTag);
	tree->Branch("JetLooseDFBTag", &JetLooseDFBTag);
	tree->Branch("JetMediumDFBTag", &JetMediumDFBTag);
	tree->Branch("JetTightDFBTag", &JetTightDFBTag);
	tree->Branch("JetLooseCSVBTag", &JetLooseCSVBTag);
	tree->Branch("JetMediumCSVBTag", &JetMediumCSVBTag);
	tree->Branch("JetTightCSVBTag", &JetTightCSVBTag);
	tree->Branch("JetLooseDFBTagSF", &JetLooseDFBTagSF);
	tree->Branch("JetMediumDFBTagSF", &JetMediumDFBTagSF);
	tree->Branch("JetTightDFBTagSF", &JetTightDFBTagSF);
	tree->Branch("JetLooseCSVBTagSF", &JetLooseCSVBTagSF);
	tree->Branch("JetMediumCSVBTagSF", &JetMediumCSVBTagSF);
	tree->Branch("JetTightCSVBTagSF", &JetTightCSVBTagSF);
	tree->Branch("JetLooseDFBTagSFUp", &JetLooseDFBTagSFUp);
	tree->Branch("JetMediumDFBTagSFUp", &JetMediumDFBTagSFUp);
	tree->Branch("JetTightDFBTagSFUp", &JetTightDFBTagSFUp);
	tree->Branch("JetLooseCSVBTagSFUp", &JetLooseCSVBTagSFUp);
	tree->Branch("JetMediumCSVBTagSFUp", &JetMediumCSVBTagSFUp);
	tree->Branch("JetTightCSVBTagSFUp", &JetTightCSVBTagSFUp);
	tree->Branch("JetLooseDFBTagSFDown", &JetLooseDFBTagSFDown);
	tree->Branch("JetMediumDFBTagSFDown", &JetMediumDFBTagSFDown);
	tree->Branch("JetTightDFBTagSFDown", &JetTightDFBTagSFDown);
	tree->Branch("JetLooseCSVBTagSFDown", &JetLooseCSVBTagSFDown);
	tree->Branch("JetMediumCSVBTagSFDown", &JetMediumCSVBTagSFDown);
	tree->Branch("JetTightCSVBTagSFDown", &JetTightCSVBTagSFDown);

	tree->Branch("METPt", &METPt);
	tree->Branch("METPhi", &METPhi);

	tree->Branch("nJet", &nJet);
	tree->Branch("nFatJet", &nFatJet);
	tree->Branch("nLooseDFBTagJet", &nLooseDFBTagJet);
	tree->Branch("nMediumDFBTagJet", &nMediumDFBTagJet);
	tree->Branch("nTightDFBTagJet", &nTightDFBTagJet);
	tree->Branch("nLooseCSVBTagJet", &nLooseCSVBTagJet);
	tree->Branch("nMediumCSVBTagJet", &nMediumCSVBTagJet);
	tree->Branch("nTightCSVBTagJet ", &nTightCSVBTagJet);
}

void JetProducer::Produce(CutFlow& cutflow, Susy1LeptonProduct *product) {
	// clear vectors for each event, otherwise this creates a memory leak
	JetPt.clear();
	JetEta.clear();
	JetPhi.clear();
	JetMass.clear();
	JetRawFactor.clear();
	JetPt_jerUp.clear();
	JetMass_jerUp.clear();
	JetPt_jerDown.clear();
	JetMass_jerDown.clear();
	JetPt_jecUp.clear();
	JetMass_jecUp.clear();
	JetPt_jecDown.clear();
	JetMass_jecDown.clear();
	JetCSVBTag.clear();
	JetDFBTag.clear();
	JetLooseDFBTag.clear();
	JetMediumDFBTag.clear();
	JetTightDFBTag.clear();
	JetLooseCSVBTag.clear();
	JetMediumCSVBTag.clear();
	JetTightCSVBTag.clear();
	JetLooseDFBTagSF.clear();
	JetMediumDFBTagSF.clear();
	JetTightDFBTagSF.clear();
	JetLooseCSVBTagSF.clear();
	JetMediumCSVBTagSF.clear();
	JetTightCSVBTagSF.clear();
	JetLooseDFBTagSFUp.clear();
	JetMediumDFBTagSFUp.clear();
	JetTightDFBTagSFUp.clear();
	JetLooseCSVBTagSFUp.clear();
	JetMediumCSVBTagSFUp.clear();
	JetTightCSVBTagSFUp.clear();
	JetLooseDFBTagSFDown.clear();
	JetMediumDFBTagSFDown.clear();
	JetTightDFBTagSFDown.clear();
	JetLooseCSVBTagSFDown.clear();
	JetMediumCSVBTagSFDown.clear();
	JetTightCSVBTagSFDown.clear();

	//Initialize all variables as -999
	METPt = -999;
	METPhi = -999;
	JetRho = -999;

	//unsigned int
	nJet = 0;
	nFatJet = 0;
	nLooseCSVBTagJet = 0;
	nMediumCSVBTagJet = 0;
	nTightCSVBTagJet = 0;
	nLooseDFBTagJet = 0;
	nMediumDFBTagJet = 0;
	nTightDFBTagJet = 0;

	//float CSVBValue = 0, DeepBValue = 0;
	float metPx = 0, metPy = 0, metPx_jerUp = 0, metPy_jerUp = 0, metPx_jerDown = 0, metPy_jerDown = 0, metPx_jecUp = 0, metPy_jecUp = 0, metPx_jecDown = 0, metPy_jecDown = 0;

	nJet = *jetNumber->Get();
	nFatJet = *fatJetNumber->Get();

	for (const JetType& type : {AK4, AK8}) {
		SetCorrector(type, runPeriod);
	}

	const float& metpt = *metPt->Get();
	const float& metphi = *metPhi->Get();
	metPx = metpt * std::cos(metphi);
	metPy = metpt * std::sin(metphi);

	metPx_jerUp = metPx;
	metPy_jerUp = metPy;
	metPx_jerDown = metPx;
	metPy_jerDown = metPy;

	metPx_jecUp = metPx;
	metPy_jecUp = metPy;
	metPx_jecDown = metPx;
	metPy_jecDown = metPy;

	for (unsigned int i = 0; i < nJet; i++) {
		//do jet correction, smearing etc.
		const float& pt = jetPt->At(i);
		const float& phi = jetPhi->At(i);
		const float& eta = jetEta->At(i);
		const float& mass = jetMass->At(i);
		const float& rawFactor = jetRawFactor->At(i);
		const float& rawPt = pt * (1 - rawFactor);
		const float& rawMass = mass * (1 - rawFactor);

		std::map<char, float> smearFactor;
		float correctionFactor = CorrectEnergy(rawPt, eta, *jetRho->Get(), jetArea->At(i), AK4);
		float correctionFactorUp = 1, correctionFactorDown = 1;

		if (jecUnc[AK4] !=  nullptr) {
			jetCorrectionUncertainty[AK4]->setJetPt(correctionFactor * rawPt);
			jetCorrectionUncertainty[AK4]->setJetEta(eta);
			const float uncUp = jetCorrectionUncertainty[AK4]->getUncertainty(true);
			jetCorrectionUncertainty[AK4]->setJetPt(correctionFactor * rawPt);
			jetCorrectionUncertainty[AK4]->setJetEta(eta);
			const float uncDown = jetCorrectionUncertainty[AK4]->getUncertainty(false);
			correctionFactorUp   = (1 + uncUp) * correctionFactor ;
			correctionFactorDown = (1 - uncDown) * correctionFactor ;
		}

		if (isData) {
			smearFactor = {{'c', 1}, {'u', 1}, {'d', 1}};
		} else {
			smearFactor = JetProducer::SmearEnergy(rawPt, eta, phi, *jetRho->Get(), jetArea->At(i), AK4);
		}

		const float& correctedPt = smearFactor['c'] * correctionFactor * rawPt;
		const float& correctedMass = smearFactor['c'] * correctionFactor * rawMass;

		const float& correctedPt_jerUp = smearFactor['u'] * correctionFactor * rawPt;
		const float& correctedMass_jerUp = smearFactor['u'] * correctionFactor * rawMass;
		const float& correctedPt_jerDown = smearFactor['d'] * correctionFactor * rawPt;
		const float& correctedMass_jerDown = smearFactor['d'] * correctionFactor * rawMass;

		const float& correctedPt_jecUp = smearFactor['c'] * correctionFactorUp * rawPt;
		const float& correctedMass_jecUp = smearFactor['c'] * correctionFactorUp * rawMass;
		const float& correctedPt_jecDown = smearFactor['c'] * correctionFactorDown * rawPt;
		const float& correctedMass_jecDown = smearFactor['c'] * correctionFactorDown * rawMass;

		if (correctedPt > ptCut && abs(eta) < etaCut) {
			const float& csvBTagValue = jetCSV->At(i);
			JetCSVBTag.push_back(csvBTagValue);
			if (csvBTagValue > deepCSVBTag[era]['t']) {
				JetLooseCSVBTag.push_back(true);
				JetMediumCSVBTag.push_back(true);
				JetTightCSVBTag.push_back(true);

				if (!isData) {
					JetLooseCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 0));
					JetMediumCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 1));
					JetTightCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 2));

					JetLooseCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 0));
					JetMediumCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 1));
					JetTightCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 2));

					JetLooseCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 0));
					JetMediumCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 1));
					JetTightCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 2));
				}
			} else if (csvBTagValue > deepCSVBTag[era]['m']) {
				JetLooseCSVBTag.push_back(true);
				JetMediumCSVBTag.push_back(true);
				JetTightCSVBTag.push_back(false);

				if (!isData) {
					JetLooseCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 0));
					JetMediumCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 1));
					JetTightCSVBTagSF.push_back(0);

					JetLooseCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 0));
					JetMediumCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 1));
					JetTightCSVBTagSFUp.push_back(0);

					JetLooseCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 0));
					JetMediumCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 1));
					JetTightCSVBTagSFDown.push_back(0);
				}
			} else if (csvBTagValue > deepCSVBTag[era]['l']) {
				JetLooseCSVBTag.push_back(true);
				JetMediumCSVBTag.push_back(false);
				JetTightCSVBTag.push_back(false);

				if (!isData) {
					JetLooseCSVBTagSF.push_back(bTagReader["deepcsv"]->Get(correctedPt, 0));
					JetMediumCSVBTagSF.push_back(0);
					JetTightCSVBTagSF.push_back(0);

					JetLooseCSVBTagSFUp.push_back(bTagReader["deepcsv"]->GetUp(correctedPt, 0));
					JetMediumCSVBTagSFUp.push_back(0);
					JetTightCSVBTagSFUp.push_back(0);

					JetLooseCSVBTagSFDown.push_back(bTagReader["deepcsv"]->GetDown(correctedPt, 0));
					JetMediumCSVBTagSFDown.push_back(0);
					JetTightCSVBTagSFDown.push_back(0);
				}
			} else {
				JetLooseCSVBTag.push_back(false);
				JetMediumCSVBTag.push_back(false);
				JetTightCSVBTag.push_back(false);

				if (!isData) {
					JetLooseCSVBTagSF.push_back(0);
					JetMediumCSVBTagSF.push_back(0);
					JetTightCSVBTagSF.push_back(0);

					JetLooseCSVBTagSFUp.push_back(0);
					JetMediumCSVBTagSFUp.push_back(0);
					JetTightCSVBTagSFUp.push_back(0);

					JetLooseCSVBTagSFDown.push_back(0);
					JetMediumCSVBTagSFDown.push_back(0);
					JetTightCSVBTagSFDown.push_back(0);
				}
			}

			const float& dfBTagValue = jetDF->At(i);
			JetDFBTag.push_back(dfBTagValue);
			if (dfBTagValue > deepFlavourBTag[era]['t']) {
				JetLooseDFBTag.push_back(true);
				JetMediumDFBTag.push_back(true);
				JetTightDFBTag.push_back(true);

				if (!isData) {
					JetLooseDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 0));
					JetMediumDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 1));
					JetTightDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 2));

					JetLooseDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 0));
					JetMediumDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 1));
					JetTightDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 2));

					JetLooseDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 0));
					JetMediumDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 1));
					JetTightDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 2));
				}
			} else if (dfBTagValue > deepFlavourBTag[era]['m']) {
				JetLooseDFBTag.push_back(true);
				JetMediumDFBTag.push_back(true);
				JetTightDFBTag.push_back(false);

				if (!isData) {
					JetLooseDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 0));
					JetMediumDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 1));
					JetTightDFBTagSF.push_back(0);

					JetLooseDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 0));
					JetMediumDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 1));
					JetTightDFBTagSFUp.push_back(0);

					JetLooseDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 0));
					JetMediumDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 1));
					JetTightDFBTagSFDown.push_back(0);
				}
			} else if (dfBTagValue > deepFlavourBTag[era]['l']) {
				JetLooseDFBTag.push_back(true);
				JetMediumDFBTag.push_back(false);
				JetTightDFBTag.push_back(false);

				if (!isData) {
					JetLooseDFBTagSF.push_back(bTagReader["deepflavour"]->Get(correctedPt, 0));
					JetMediumDFBTagSF.push_back(0);
					JetTightDFBTagSF.push_back(0);

					JetLooseDFBTagSFUp.push_back(bTagReader["deepflavour"]->GetUp(correctedPt, 0));
					JetMediumDFBTagSFUp.push_back(0);
					JetTightDFBTagSFUp.push_back(0);

					JetLooseDFBTagSFDown.push_back(bTagReader["deepflavour"]->GetDown(correctedPt, 0));
					JetMediumDFBTagSFDown.push_back(0);
					JetTightDFBTagSFDown.push_back(0);
				}
			} else {
				JetLooseDFBTag.push_back(false);
				JetMediumDFBTag.push_back(false);
				JetTightDFBTag.push_back(false);

				if (!isData) {
					JetLooseDFBTagSF.push_back(0);
					JetMediumDFBTagSF.push_back(0);
					JetTightDFBTagSF.push_back(0);

					JetLooseDFBTagSFUp.push_back(0);
					JetMediumDFBTagSFUp.push_back(0);
					JetTightDFBTagSFUp.push_back(0);

					JetLooseDFBTagSFDown.push_back(0);
					JetMediumDFBTagSFDown.push_back(0);
					JetTightDFBTagSFDown.push_back(0);
				}
			}

			//Jet four momentum components
			JetPt.push_back(correctedPt); JetEta.push_back(eta); JetPhi.push_back(phi); JetMass.push_back(correctedMass); JetRawPt.push_back(rawPt); JetRawMass.push_back(rawPt), JetRawFactor.push_back(rawFactor);

			JetPt_jecUp.push_back(correctedPt_jecUp); JetMass_jecUp.push_back(correctedMass_jecUp); JetPt_jecDown.push_back(correctedPt_jecDown); JetMass_jecDown.push_back(correctedMass_jecDown);

			JetPt_jerUp.push_back(correctedPt_jerUp); JetMass_jerUp.push_back(correctedMass_jerUp); JetPt_jerDown.push_back(correctedPt_jerDown); JetMass_jerDown.push_back(correctedMass_jerDown);

			//JECCorrection.push_back(correctionFactor);
			metPx += pt * std::cos(phi) - JetPt.back() * std::cos(phi);
			metPy += pt * std::sin(phi) - JetPt.back() * std::sin(phi);

			metPx_jerUp += pt * std::cos(phi) - JetPt_jerUp.back() * std::cos(phi);
			metPy_jerUp += pt * std::sin(phi) - JetPt_jerUp.back() * std::sin(phi);
			metPx_jerDown += pt * std::cos(phi) - JetPt_jerDown.back() * std::cos(phi);
			metPy_jerDown += pt * std::sin(phi) - JetPt_jerDown.back() * std::sin(phi);

			metPx_jecUp += pt * std::cos(phi) - JetPt_jecUp.back() * std::cos(phi);
			metPy_jecUp += pt * std::sin(phi) - JetPt_jecUp.back() * std::sin(phi);
			metPx_jecDown += pt * std::cos(phi) - JetPt_jecDown.back() * std::cos(phi);
			metPy_jecDown += pt * std::sin(phi) - JetPt_jecDown.back() * std::sin(phi);
		}
	}

	//Cleanup, remove at most 1 Jet (i.e. per lepton)
	float minDeltaR = 999;
	int nearestJetIndex = -1;
	for (unsigned int i = 0; i < JetPt.size(); i++) {
		float deltaR = DeltaR(JetEta.at(i), JetPhi.at(i), product->leptonEta, product->leptonPhi);
		if (deltaR < minDeltaR) {
			minDeltaR = deltaR;
			nearestJetIndex = i;
		}
	}
	if (minDeltaR < deltaRCut && nearestJetIndex > 0) {
		JetPt.erase(JetPt.begin() + nearestJetIndex);
		JetEta.erase(JetEta.begin() + nearestJetIndex);
		JetPhi.erase(JetPhi.begin() + nearestJetIndex);
		JetMass.erase(JetMass.begin() + nearestJetIndex);
		JetRawFactor.erase(JetRawFactor.begin() + nearestJetIndex);

		JetPt_jecUp.erase(JetPt_jecUp.begin() + nearestJetIndex);
		JetPt_jecDown.erase(JetPt_jecDown.begin() + nearestJetIndex);
		JetPt_jerUp.erase(JetPt_jerUp.begin() + nearestJetIndex);
		JetPt_jerDown.erase(JetPt_jerDown.begin() + nearestJetIndex);
		JetMass_jecUp.erase(JetMass_jecUp.begin() + nearestJetIndex);
		JetMass_jecDown.erase(JetMass_jecDown.begin() + nearestJetIndex);
		JetMass_jerUp.erase(JetMass_jerUp.begin() + nearestJetIndex);
		JetMass_jerDown.erase(JetMass_jerDown.begin() + nearestJetIndex);

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

	for (unsigned int i = 0; i < nJet; i++) {
		if (JetTightCSVBTag[i]) { nLooseCSVBTagJet++; nMediumCSVBTagJet++; nTightCSVBTagJet++;}
		else if (JetMediumCSVBTag[i]) { nLooseCSVBTagJet++; nMediumCSVBTagJet++;}
		else if (JetLooseCSVBTag[i]) { nLooseCSVBTagJet++;}

		if (JetTightDFBTag[i]) { nLooseDFBTagJet++; nMediumDFBTagJet++; nTightDFBTagJet++;}
		else if (JetMediumDFBTag[i]) { nLooseDFBTagJet++; nMediumDFBTagJet++;}
		else if (JetLooseDFBTag[i]) { nLooseDFBTagJet++;}
	}

	nJet = JetPt.size();
	METPt = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
	METPhi = std::atan2(metPy, metPx);

	METPt_jerUp = std::sqrt(std::pow(metPx_jerUp, 2) + std::pow(metPy_jerUp, 2));
	METPhi_jerUp = std::atan2(metPy_jerUp, metPx_jerUp);
	METPt_jerDown = std::sqrt(std::pow(metPx_jerDown, 2) + std::pow(metPy_jerDown, 2));
	METPhi_jerDown = std::atan2(metPy_jerDown, metPx_jerDown);

	METPt_jecUp = std::sqrt(std::pow(metPx_jecUp, 2) + std::pow(metPy_jecUp, 2));
	METPhi_jecUp = std::atan2(metPy_jecUp, metPx_jecUp);
	METPt_jecDown = std::sqrt(std::pow(metPx_jecDown, 2) + std::pow(metPy_jecDown, 2));
	METPhi_jecDown = std::atan2(metPy_jecDown, metPx_jecDown);

	//Sort Vectors according to JetPt, since the correction could have changed the order
	std::vector<int> idx(nJet);
	std::iota(idx.begin(), idx.end(), 0);
	std::stable_sort(idx.begin(), idx.end(), [&](int i1, int i2) {return JetPt[i1] > JetPt[i2];});

	if (nJet != 0) {
		SortByIndex<std::vector<float>>(JetPt, idx);
		SortByIndex<std::vector<float>>(JetPhi, idx);
		SortByIndex<std::vector<float>>(JetEta, idx);
		SortByIndex<std::vector<float>>(JetMass, idx);
		SortByIndex<std::vector<float>>(JetRawFactor, idx);

		SortByIndex<std::vector<bool>>(JetLooseCSVBTag, idx);
		SortByIndex<std::vector<bool>>(JetMediumCSVBTag, idx);
		SortByIndex<std::vector<bool>>(JetTightCSVBTag, idx);
		SortByIndex<std::vector<bool>>(JetLooseDFBTag, idx);
		SortByIndex<std::vector<bool>>(JetMediumDFBTag, idx);
		SortByIndex<std::vector<bool>>(JetTightDFBTag, idx);

		if (!isData) {
			SortByIndex<std::vector<float>>(JetLooseCSVBTagSF, idx);
			SortByIndex<std::vector<float>>(JetMediumCSVBTagSF, idx);
			SortByIndex<std::vector<float>>(JetTightCSVBTagSF, idx);
			SortByIndex<std::vector<float>>(JetLooseDFBTagSF, idx);
			SortByIndex<std::vector<float>>(JetMediumDFBTagSF, idx);
			SortByIndex<std::vector<float>>(JetTightDFBTagSF, idx);

			SortByIndex<std::vector<float>>(JetLooseCSVBTagSFUp, idx);
			SortByIndex<std::vector<float>>(JetMediumCSVBTagSFUp, idx);
			SortByIndex<std::vector<float>>(JetTightCSVBTagSFUp, idx);
			SortByIndex<std::vector<float>>(JetLooseDFBTagSFUp, idx);
			SortByIndex<std::vector<float>>(JetMediumDFBTagSFUp, idx);
			SortByIndex<std::vector<float>>(JetTightDFBTagSFUp, idx);

			SortByIndex<std::vector<float>>(JetLooseCSVBTagSFDown, idx);
			SortByIndex<std::vector<float>>(JetMediumCSVBTagSFDown, idx);
			SortByIndex<std::vector<float>>(JetTightCSVBTagSFDown, idx);
			SortByIndex<std::vector<float>>(JetLooseDFBTagSFDown, idx);
			SortByIndex<std::vector<float>>(JetMediumDFBTagSFDown, idx);
			SortByIndex<std::vector<float>>(JetTightDFBTagSFDown, idx);
		}
	}

	//Store values in product to calculate high level variables
	product->nJet = nJet;
	product->jetPt = JetPt;
	product->jetPhi = JetPhi;
	product->jetEta = JetEta;
	product->jetMass = JetMass;

	product->metPt = METPt;
	product->metPhi = METPhi;

	product->nLooseCSVBTagJet = nLooseCSVBTagJet;
	product->nMediumCSVBTagJet = nMediumCSVBTagJet;
	product->nTightCSVBTagJet = nTightCSVBTagJet;

	product->nLooseDFBTagJet = nLooseDFBTagJet;
	product->nMediumDFBTagJet = nMediumDFBTagJet;
	product->nTightDFBTagJet = nTightDFBTagJet;


	/* TODO Loop over fat jets, will be done later
	for (unsigned int i = 0; i < nFatJet; i++) {
		//JER smearing

		float fatPt = fatJetPt->At(i);
		const float fatEta = fatJetEta->At(i);
		const float fatPhi = fatJetPhi->At(i);
		const float fatMass = fatJetMass->At(i);

		std::map<char, float> smearFactor;
		float correctionFactor = CorrectEnergy(fatPt, fatEta, *jetRho->Get(), fatJetArea->At(i), AK4);
		float correctionFactorUp, correctionFactorDown;

		if (isData) {
			smearFactor = {{'c', 1.0}, {'u', 1.0}, {'d', 1.0}};
		} else {
			smearFactor = JetProducer::SmearEnergy(fatPt, fatEta, fatPhi, *jetRho->Get(), fatJetArea->At(i), AK4);
		}

		//Get jet uncertainty
		if (jecUnc[AK8] !=  nullptr) {
			jetCorrectionUncertainty[AK8]->setJetPt(correctionFactor * fatPt);
			jetCorrectionUncertainty[AK8]->setJetEta(fatEta);
			const float uncUp = jetCorrectionUncertainty[AK4]->getUncertainty(true);
			const float uncDown = jetCorrectionUncertainty[AK4]->getUncertainty(false);
			correctionFactorUp   = (1 + uncUp) * correctionFactor ;
			correctionFactorDown = (1 - uncDown) * correctionFactor ;
		}
	}
	*/

	/* TODO
		->Proper runNumber of the datasets, maybe there is als a smarter way.. One could just parse the runPeriod as an argument to the wrapper or us a function that gets it from the filename
		->FatJet Loop
	*/

	if (nJet!=0) {
		std::string cutName("JetPt > 20, |JetEta| < 2.4");
		cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
	} else {
		cutflow.passed = false;
	}

	for (const JetType& type : {AK4, AK8}) {
		delete jetCorrector[type];
	}

}

void JetProducer::EndJob(TFile* file) {
	delete bTagReader["deepcsv"];
	delete bTagReader["deepflavour"];
	for (const JetType& type : {AK4, AK8}) {
		delete jetCorrectionUncertainty[type];
	}
}
