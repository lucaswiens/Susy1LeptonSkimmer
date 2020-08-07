#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/JetProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>

#include <cmath>
#include <random>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > RMFLV;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > CartesianRMFLV;

JetProducer::JetProducer(const int& era, const float& ptCut, const float& etaCut, const float& deltaRCut, TTreeReader& reader):
	BaseProducer(&reader),
	era(era),
	ptCut(ptCut),
	etaCut(etaCut),
	deltaRCut(deltaRCut)
	{}

template <typename T>
void JetProducer::SortByIndex(T& var, std::vector<int> idx){
	std::vector<float> tmp(nJet);
	for (unsigned int i = 0; i < nJet; i++){
		tmp.at(i) = var.at(idx[i]);
	}
	var = std::move(tmp);
}

void JetProducer::SetCorrector(const JetType& type, const int& runNumber){
	std::vector<JetCorrectorParameters> corrVec;

	for(std::string fileName: isData? jecData[era] : jecMC[era]){
		if(fileName.find("@") != std::string::npos){
			for(std::pair<std::string, std::pair<int, int>> eraNames: runEras[era]){
				if(eraNames.second.first <= runNumber and runNumber <= eraNames.second.second){
					fileName.replace(fileName.find("@"), 1, eraNames.first);
				}
			}
		}

		fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");

		corrVec.push_back(JetCorrectorParameters(fileName));
	}

	jetCorrector[type] = new FactorizedJetCorrector(corrVec);
}

//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
float JetProducer::CorrectEnergy(const float& pt, const float& eta, const float& rho, const float& area, const JetType &type){
	jetCorrector[type]->setJetPt(pt);
	jetCorrector[type]->setJetEta(eta);
	jetCorrector[type]->setRho(rho);
	jetCorrector[type]->setJetA(area);
	return jetCorrector[type]->getCorrection();
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
//https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
std::map<char, float> JetProducer::SmearEnergy(const float& pt, const float& eta, const float& phi, const float& rho, const float& coneSize, const JetType& type){
	jetParameter.setJetPt(pt).setJetEta(eta).setRho(rho);

	float reso = resolution[type].getResolution(jetParameter);
	float resoSF;
	float resoSFUp;
	float resoSFDown;
	//if(!isJERsyst) resoSF = resolution_sf[type].getScaleFactor(jetParameter);
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

	unsigned int size = (type == AK4) ? genJetPt->GetSize(): genFatJetPt->GetSize();

	//Loop over all gen jets and find match
	for(unsigned int i = 0; i < size; i++){
		genPt = (type == AK4) ? genJetPt->At(i): genFatJetPt->At(i);
		genPhi = (type == AK4) ? genJetPhi->At(i): genFatJetPhi->At(i);
		genEta = (type == AK4) ? genJetEta->At(i): genFatJetEta->At(i);
		genMass = (type == AK4) ? genJetMass->At(i): genFatJetMass->At(i);

		dR = std::sqrt(std::pow(phi - genPhi, 2) + std::pow(eta - genEta, 2));

		//Check if jet and gen jet are matched
		if(dR < coneSize/2. and abs(pt - genPt) < 3. * reso * pt){
			genJet[type] = ROOT::Math::PtEtaPhiMVector(genPt, genEta, genPhi, genMass);
			isMatched = true;
			break;
		}
	}

	//If you found gen matched
	if(isMatched){
		smearFactor     = 1. + (resoSF-1) * (pt - genPt) / pt;
		smearFactorUp   = 1. + (resoSFUp-1) * (pt - genPt) / pt;
		smearFactorDown = 1. + (resoSFDown-1) * (pt - genPt) / pt;
	} else if(resoSF > 1.) {
		std::default_random_engine generator;
		std::normal_distribution<> gaus(0, reso * std::sqrt(resoSF * resoSF - 1));
		std::normal_distribution<> gausUp(0, reso * std::sqrt(resoSFUp * resoSFUp - 1));
		std::normal_distribution<> gausDown(0, reso * std::sqrt(resoSFDown * resoSFDown - 1));
		smearFactor     = 1. + gaus(generator);
		smearFactorUp   = 1. + gausUp(generator);
		smearFactorDown = 1. + gausDown(generator);
	}

	//Check if direction of jet not changed
	if(pt * smearFactor     < 1e-2){ smearFactor     = 1e-2 / pt;}
	if(pt * smearFactorUp   < 1e-2){ smearFactorUp   = 1e-2 / pt;}
	if(pt * smearFactorDown < 1e-2){ smearFactorDown = 1e-2 / pt;}

	return {{'c', smearFactor}, {'u', smearFactorUp}, {'d', smearFactorDown}};
}

void JetProducer::SetGenParticles(const int& i, const float& pt, const float& eta, const float& phi, const std::vector<int>& pdgID, const JetType &type){
	float dR;

	//partID.at(type).push_back(-99.);
	//mothID.at(type).push_back(-99.);
	//grandID.at(type).push_back(-99.);

	partID[type].push_back(-99.);
	mothID[type].push_back(-99.);
	grandID[type].push_back(-99.);

	float bestDR = 30.; float bestPt = 100.;
	int index=-1;

	//Check if gen matched particle exist
	if(genJet[type].Pt() != 0){
		//Find Gen particle to gen Jet
		const char& size = genPt->GetSize();

		for(int i=0; i < size; i++){
			const int ID = abs(genID->At(i));

			if(std::find(pdgID.begin(), pdgID.end(), ID) != pdgID.end()){
				index = FirstCopy(i, ID);
			} else continue;

			const float& ptGen = genPt->At(index);
			const float& phiGen = genPhi->At(index);
			const float& etaGen = genEta->At(index);

			dR = BaseProducer::DeltaR(etaGen, phiGen, genJet[type].Eta(), genJet[type].Phi());
			const float& dPt = abs((genJet[type].Pt() - ptGen) / genJet[type].Pt());

			if(dR < bestDR and dPt < bestPt){
				if(std::find(alreadySeen.begin(), alreadySeen.end(), index) != alreadySeen.end()) continue;
				else alreadySeen.push_back(index);

				bestDR = dR;
				bestPt = dPt;

				const int motherID = abs(genID->At(genMotherIdx->At(index)));
				const int motherIdx = FirstCopy(genMotherIdx->At(index), motherID);
				const int grandMotherID = abs(genID->At(genMotherIdx->At(motherIdx)));

				partID[type].back() = ID;
				mothID[type].back() = motherID;
				grandID[type].back() = grandMotherID;
				break;
			}
		}
	}
	if(index != -1) alreadySeen.push_back(index);
}

void JetProducer::BeginJob(TTree* tree, bool &isData, Susy1LeptonProduct *product, const bool& isSyst){
	//Set data bool
	this->isData = isData;
	this->isSyst = isSyst;

	// Path to Correction files
	//std::string filePath = "$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/jme/";
	std::string filePath = "/nfs/dust/cms/user/wiens/CMSSW/CMSSW_10_6_8/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/jme/";
	jecMC = {
		{2016, {filePath + "Summer16/Summer16_07Aug2017_V11_MC_L1FastJet_&PFchs.txt",
			filePath + "Summer16/Summer16_07Aug2017_V11_MC_L2Relative_&PFchs.txt",
			filePath + "Summer16/Summer16_07Aug2017_V11_MC_L3Absolute_&PFchs.txt"}
		},
		{2017, {filePath + "Fall17/Fall17_17Nov2017_V32_MC_L1FastJet_&PFchs.txt",
			filePath + "Fall17/Fall17_17Nov2017_V32_MC_L2Relative_&PFchs.txt",
			filePath + "Fall17/Fall17_17Nov2017_V32_MC_L3Absolute_&PFchs.txt"}
		},
		{2018, {filePath + "Autumn18/Autumn18_V19_MC_L1FastJet_&PFchs.txt",
			filePath + "Autumn18/Autumn18_V19_MC_L2Relative_&PFchs.txt",
			filePath + "Autumn18/Autumn18_V19_MC_L3Absolute_&PFchs.txt"}
		},
	};

	jecData = {
		{2016, {filePath + "Summer16/Summer16_07Aug2017@_V11_DATA_L1FastJet_&PFchs.txt",
			filePath + "Summer16/Summer16_07Aug2017@_V11_DATA_L2Relative_&PFchs.txt",
			filePath + "Summer16/Summer16_07Aug2017@_V11_DATA_L3Absolute_&PFchs.txt",
			filePath + "Summer16/Summer16_07Aug2017@_V11_DATA_L2L3Residual_&PFchs.txt"}
		},
		{2017, {filePath + "Fall17/Fall17_17Nov2017@_V32_DATA_L1FastJet_&PFchs.txt",
			filePath + "Fall17/Fall17_17Nov2017@_V32_DATA_L2Relative_&PFchs.txt",
			filePath + "Fall17/Fall17_17Nov2017@_V32_DATA_L3Absolute_&PFchs.txt",
			filePath + "Fall17/Fall17_17Nov2017@_V32_DATA_L2L3Residual_&PFchs.txt"}
		},
		{2018, {filePath + "Autumn18/Autumn18_Run@_V9_DATA_L1FastJet_&PFchs.txt",
			filePath + "Autumn18/Autumn18_Run@_V9_DATA_L2Relative_&PFchs.txt",
			filePath + "Autumn18/Autumn18_Run@_V9_DATA_L3Absolute_&PFchs.txt",
			filePath + "Autumn18/Autumn18_Run@_V9_DATA_L2L3Residual_&PFchs.txt"}
		}
	};

	jecUnc = {
		{2016, filePath + "Summer16/Summer16_07Aug2017_V11_MC_UncertaintySources_&PFchs.txt"},
		{2017, filePath + "Fall17/Fall17_17Nov2017_V32_MC_UncertaintySources_&PFchs.txt"},
		{2018, filePath + "Autumn18/Autumn18_V19_MC_UncertaintySources_&PFchs.txt"},
	};

	jmeSF = {
		{2016, filePath + "Summer16/Summer16_25nsV1_MC_SF_&PFchs.txt"},
		{2017, filePath + "Fall17/Fall17_V3_MC_SF_&PFchs.txt"},
		{2018, filePath + "Autumn18/Autumn18_V7_MC_SF_&PFchs.txt"},
	};

	jmePtReso = {
		{2016, filePath + "Summer16/Summer16_25nsV1_MC_PtResolution_&PFchs.txt"},
		{2017, filePath + "Fall17/Fall17_V3_MC_PtResolution_&PFchs.txt"},
		{2018, filePath + "Autumn18/Autumn18_V7_MC_PtResolution_&PFchs.txt"},
	};

	jecFastSim = {
		{2016, {filePath + "Summer16/Summer16_FastSimV1_MC_L1FastJet_&PFchs.txt",
			filePath + "Summer16/Summer16_FastSimV1_MC_L2Relative_&PFchs.txt",
			filePath + "Summer16/Summer16_FastSimV1_MC_L3Absolute_&PFchs.txt"}
		},
		{2017, {filePath + "Fall17/Fall17_FastSimV1_MC_L1FastJet_&PFchs.txt",
			filePath + "Fall17/Fall17_FastSimV1_MC_L2Relative_&PFchs.txt",
			filePath + "Fall17/Fall17_FastSimV1_MC_L3Absolute_&PFchs.txt"}
		},
		{2018, {filePath + "Autumn18/Autumn18_FastSimV1_MC_L1FastJet_&PFchs.txt",
			filePath + "Autumn18/Autumn18_FastSimV1_MC_L2Relative_&PFchs.txt",
			filePath + "Autumn18/Autumn18_FastSimV1_MC_L3Absolute_&PFchs.txt"}
		}
	};//also do uncertainties for this

	//Set TTreeReader for genpart and trigger obj from BaseProducer
	SetCollection(this->isData);
	jetNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nJet");
	jetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_pt");
	jetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_phi");
	jetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_eta");
	jetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_mass");
	jetArea = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_area");
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
	fatJetTau1 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau1");
	fatJetTau2 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau2");
	fatJetTau3 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau3");

	if(!this->isData){
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

	for(JetType type: {AK4, AK8}){
		//Set configuration for JER tools
		std::string fileName = jmePtReso[era];
		fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
		resolution[type] = JME::JetResolution(fileName);

		fileName = jmeSF[era];
		fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
		resolution_sf[type] = JME::JetResolutionScaleFactor(fileName);

		//Set object to get JEC uncertainty
		fileName = jecUnc[era];
		fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
		jetCorrectionUncertainty[type] = new JetCorrectionUncertainty(JetCorrectorParameters(fileName, "Total"));
	}


	//Set Branches of output tree
	tree->Branch("nJet", &nJet);
	tree->Branch("JetPt", &JetPt);
	tree->Branch("JetPhi", &JetPhi);
	tree->Branch("JetEta", &JetEta);
	tree->Branch("JetRho", &JetRho);
	tree->Branch("JetMass", &JetMass);
	tree->Branch("METPt", &METPt);
	tree->Branch("METPhi", &METPhi);
}

void JetProducer::Produce(CutFlow& cutflow, Susy1LeptonProduct *product){
	//Initialize all variables as -999
	JetPt.clear();
	JetPhi.clear();
	JetEta.clear();
	METPt = -999;
	METPhi = -999;

	//float CSVBValue = 0, DeepBValue = 0;
	float metPx = 0, metPy;

	nJet = *jetNumber->Get();
	nFatJet = *fatJetNumber->Get();
	runNumber = *run->Get();

	for(const JetType& type: {AK4, AK8}){
			SetCorrector(type, runNumber);
	}

	//RMFLV metP4 = RMFLV(0, 0, 0, 0);
	//metP4.SetPt(1);

	const float& metpt = *metPt->Get();
	const float& metphi = *metPhi->Get();
	//const float& jetrho = *jetRho->Get();
	metPx = metpt * std::cos(metphi);
	metPy = metpt * std::sin(metphi);

	for (unsigned int i = 0; i < nJet; i++){
		//do jet correction, smearing etc.
		const float& pt = jetPt->At(i);
		const float& phi = jetPhi->At(i);
		const float& eta = jetEta->At(i);
		const float& mass = jetMass->At(i);

		std::map<char, float> smearFactor;
		float correctionFactor = CorrectEnergy(pt, eta, *jetRho->Get(), jetArea->At(i), AK4);
		float correctionFactorUp, correctionFactorDown;

		if(jecUnc[AK4] !=  nullptr){
			jetCorrectionUncertainty[AK4]->setJetPt(correctionFactor * pt);
			jetCorrectionUncertainty[AK4]->setJetEta(eta);
			const float uncUp = jetCorrectionUncertainty[AK4]->getUncertainty(true);
			jetCorrectionUncertainty[AK4]->setJetPt(correctionFactor * pt);
			jetCorrectionUncertainty[AK4]->setJetEta(eta);
			const float uncDown = jetCorrectionUncertainty[AK4]->getUncertainty(false);
			//correctionFactorUp   = (1 + uncUp) * correctionFactor ;
			//correctionFactorDown = (1 - uncDown) * correctionFactor ;
		}

		if(isData){
			smearFactor = {{'c', 1.0}, {'u', 1.0}, {'d', 1.0}};
		} else {
			smearFactor = JetProducer::SmearEnergy(pt, eta, phi, *jetRho->Get(), jetArea->At(i), AK4);
		}

		const float& correctedPt = smearFactor['c'] * correctionFactor * pt;
		const float& correctedMass = smearFactor['c'] * correctionFactor * mass;

		if( correctedPt > ptCut && abs(eta) < etaCut){
			//correctionFactor = CorrectEnergy(pt, eta, *jetRho->Get(), jetArea->At(i), AK4);

			//Jet four momentum components
			JetPt.push_back(correctedPt);
			JetEta.push_back(eta);
			JetPhi.push_back(phi);
			JetMass.push_back(correctedMass);
			//JECCorrection.push_back(correctionFactor);
			metPx += pt * std::cos(phi) - JetPt.back() * std::cos(phi);
			metPy += pt * std::sin(phi) - JetPt.back() * std::sin(phi);
			if(!isData){
				SetGenParticles(i, JetPt.back(), JetEta.back(), JetPhi.back(), {5}, AK4); //pdgid = 5 => b quark
			}
		}
	}

	//Cleanup, remove at most 1 Jet (i.e. per lepton)
	float minDeltaR = 999;
	int nearestJetIndex = -1;
	for (unsigned int i = 0; i < JetPt.size(); i++){
		float deltaR = DeltaR(JetEta.at(i), JetPhi.at(i), product->leptonEta, product->leptonPhi);
		if (deltaR < minDeltaR){
			minDeltaR = deltaR;
			nearestJetIndex = i;
		}
	}
	if (minDeltaR < deltaRCut && nearestJetIndex > 0){
		JetPt.erase(JetPt.begin() + nearestJetIndex);
		JetEta.erase(JetEta.begin() + nearestJetIndex);
		JetPhi.erase(JetPhi.begin() + nearestJetIndex);
		JetMass.erase(JetMass.begin() + nearestJetIndex);
	}

	/* TODO Loop over fat jets, will be done later
	for(unsigned int i = 0; i < nFatJet; i++){
		//JER smearing

		float fatPt = fatJetPt->At(i);
		const float fatEta = fatJetEta->At(i);
		const float fatPhi = fatJetPhi->At(i);
		const float fatMass = fatJetMass->At(i);

		std::map<char, float> smearFactor;
		float correctionFactor = CorrectEnergy(fatPt, fatEta, *jetRho->Get(), fatJetArea->At(i), AK4);
		float correctionFactorUp, correctionFactorDown;

		if(isData){
			smearFactor = {{'c', 1.0}, {'u', 1.0}, {'d', 1.0}};
		} else {
			smearFactor = JetProducer::SmearEnergy(fatPt, fatEta, fatPhi, *jetRho->Get(), fatJetArea->At(i), AK4);
		}

		//Get jet uncertainty
		if(jecUnc[AK8] !=  nullptr){
			jetCorrectionUncertainty[AK8]->setJetPt(correctionFactor * fatPt);
			jetCorrectionUncertainty[AK8]->setJetEta(fatEta);
			const float uncUp = jetCorrectionUncertainty[AK4]->getUncertainty(true);
			const float uncDown = jetCorrectionUncertainty[AK4]->getUncertainty(false);
			correctionFactorUp   = (1 + uncUp) * correctionFactor ;
			correctionFactorDown = (1 - uncDown) * correctionFactor ;
		}
	}
	*/

	nJet = JetPt.size();
	METPt = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
	METPhi = std::atan2(metPy, metPx);

	//Sort Vectors according to JetPt, since the correction could have changed the order
	std::vector<int> idx(nJet);
	std::iota(idx.begin(), idx.end(), 0);
	std::stable_sort(idx.begin(), idx.end(), [&](int i1, int i2) {return JetPt[i1] > JetPt[i2];});

	if (nJet != 0){
		SortByIndex<std::vector<float>>(JetPt, idx);
		SortByIndex<std::vector<float>>(JetPhi, idx);
		SortByIndex<std::vector<float>>(JetEta, idx);
		SortByIndex<std::vector<float>>(JetMass, idx);
	}

	//Store values in product to calculate high level variables
	product->nJet = nJet;
	product->jetPt = JetPt;
	product->jetPhi = JetPhi;
	product->jetEta = JetEta;
	product->jetMass = JetMass;

	product->metPt = METPt;
	product->metPhi = METPhi;

	if (nJet!=0){
		std::string cutName("JetPt > 20, |JetEta| < 2.4");
		cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
	} else {
		cutflow.passed = false;
	}

	for(const JetType& type: {AK4, AK8}){
		delete jetCorrector[type];
	}

}

void JetProducer::EndJob(TFile* file){
}