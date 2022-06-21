#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/DeltaPhiProducer.h>

//DeltaPhiProducer::DeltaPhiProducer(TTreeReader &reader): BaseProducer(&reader) {}


float DeltaPhiProducer::DeltaPhi(ROOT::Math::PtEtaPhiMVector v1, ROOT::Math::PtEtaPhiMVector v2) {
	float dPhi = v1.Phi() - v2.Phi();
	while (dPhi >= TMath::Pi()) dPhi -= TMath::TwoPi();
	while (dPhi < -TMath::Pi()) dPhi += TMath::TwoPi();
	return dPhi;
}

//void DeltaPhiProducer::BeginJob(std::shared_ptr<TTree> tree, bool &isData, bool &doSystematics) {
//	//Set data bool
//	this->isData = isData;
//	this->doSystematics = doSystematics;
//
//	//isoTrackNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nIsoTrack");
//	//isoTrackPt = std::make_unique<TTreeReaderArray<float>>(*reader, "IsoTrack_pt");
//	//isoTrackEta = std::make_unique<TTreeReaderArray<float>>(*reader, "IsoTrack_eta");
//	//isoTrackPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "IsoTrack_phi");
//	////isoTrackMass = std::make_unique<TTreeReaderArray<float>>(*reader, "IsoTrack_mass");
//	////isoTrackCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "IsoTrack_charge");
//	//isoTrackPdgId = std::make_unique<TTreeReaderArray<int>>(*reader, "IsoTrack_pdgId");
//
//	//Set Branches of output tree
//	tree->Branch("HT", &HT);
//	tree->Branch("LT", &LT);
//	tree->Branch("LP", &LP);
//	tree->Branch("DeltaPhi", &deltaPhi);
//	tree->Branch("dPhi", &dPhi);
//	tree->Branch("WBosonMt", &wBosonMt);
//	tree->Branch("signalRegion", &signalRegionCSV);
//	tree->Branch("signalRegionCSV", &signalRegionCSV);
//	tree->Branch("signalRegionDF", &signalRegionDF);
//
//	tree->Branch("IsoTrackMt2", &IsoTrackMt2);
//	tree->Branch("IsoTrackPt", &IsoTrackPt);
//	tree->Branch("IsoTrackPdgId", &IsoTrackPdgId);
//	tree->Branch("IsoTrackVeto", &IsoTrackVeto);
//	tree->Branch("IsoTrackHadronicDecay", &IsoTrackHadronicDecay);
//}

void DeltaPhiProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	//Initialize all variables as -999
	//IsoTrackMt2.clear();
	//IsoTrackPt.clear();
	//IsoTrackPdgId.clear();
	//IsoTrackHadronicDecay.clear();
	//IsoTrackVeto = false;

	//HT = -999;
	//LT = -999;
	//LP = -999;
	//deltaPhi = -999;
	//dPhi = -999;
	//wBosonMt = -999;
	//signalRegionCSV = -999;
	//signalRegionDF = -999;

	/*
	if (product.leptonPt != -999 && product.metPt != -999) {
		ROOT::Math::PtEtaPhiMVector leptonP4 = ROOT::Math::PtEtaPhiMVector(product.leptonPt, product.leptonEta, product.leptonPhi, product.leptonMass);
		ROOT::Math::PtEtaPhiMVector metP4 = ROOT::Math::PtEtaPhiMVector(product.metPt, 0, product.metPhi, 0);
		ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + metP4;

		LT = product.leptonPt + product.metPt;
		deltaPhi = DeltaPhi(leptonP4, wBosonP4);
		dPhi = std::abs(deltaPhi);
		LP = product.leptonPt / wBosonP4.Pt() * std::cos(deltaPhi);
		wBosonMt = wBosonP4.Mt();

		HT = 0;
		for (unsigned int i = 0; i < product.nJet; i++) {
			HT += product.jetPt.at(i);
		}

		if (product.nMediumCSVBTagJet == 0) {// 0-B SRs -- simplified dPhi
			if (LT < 250) {signalRegionCSV = 0;}
			else if (LT > 250) {signalRegionCSV = dPhi > 0.75;}
			// BLIND data
			if (isData && product.nJet >= 5) { signalRegionCSV = -signalRegionCSV;}
		} else if (product.nJet < 99) {// Multi-B SRs
			if (LT < 250) { signalRegionCSV = 0;}
			else if (LT < 350) { signalRegionCSV = dPhi > 1.0;}
			else if (LT < 600) { signalRegionCSV = dPhi > 0.75;}
			else if (LT > 600) { signalRegionCSV = dPhi > 0.5;}

			// BLIND data
			if (isData && product.nJet >= 6) { signalRegionCSV = -signalRegionCSV;}
		}

		if (product.nMediumDFBTagJet == 0) {//0-B SRs -- simplified dPhi
			if (LT < 250) {signalRegionDF = 0;}
			else if (LT > 250) {signalRegionDF = dPhi > 0.75;}
			// BLIND data
			if (isData && product.nJet >= 5) { signalRegionDF = -signalRegionDF;}
		} else if (product.nJet < 99) {// Multi-B SRs
			if (LT < 250) { signalRegionDF = 0;}
			else if (LT < 350) { signalRegionDF = dPhi > 1.0;}
			else if (LT < 600) { signalRegionDF = dPhi > 0.75;}
			else if (LT > 600) { signalRegionDF = dPhi > 0.5;}

			// BLIND data
			if (isData && product.nJet >= 6) { signalRegionDF = -signalRegionDF;}
		}

		heppy::Davismt2 mt2obj;
		nIsoTrack = -999;//*isoTrackNumber->Get();
		float deltaRMin = 0.1, isoMt2Cut = -999;

		for (unsigned int j = 0; j < nIsoTrack; j++) {
			const int &isotrackpdgid = -999;//isoTrackPdgId->At(j);
			// The isotrack veto is used to reject dilepton events, thus only calculate variables if its a lepton
			const float &isotrackpt = -999;//isoTrackPt->At(j);
			const float &isotracketa = -999;//isoTrackEta->At(j);
			const float &isotrackphi = -999;//isoTrackPhi->At(j);

			const int &isotrackcharge = (0 < isotrackpdgid) - (isotrackpdgid < 0); //only true for leptons but veto only matters for leptons

			if (isotrackcharge == product.leptonCharge) continue;

			float deltaR = Utility::DeltaR(product.leptonEta, product.leptonPhi, isotracketa, isotrackphi);

			if (deltaR < deltaRMin) continue;

			ROOT::Math::PtEtaPhiMVector isotrackP4 = ROOT::Math::PtEtaPhiMVector(isotrackpt, isotracketa, isotrackphi, 0);

			double a[3] = {leptonP4.M(), leptonP4.X(), leptonP4.Y()};
			double b[3] = {isotrackP4.M(), isotrackP4.X(), isotrackP4.Y()};
			double c[3] = {metP4.M(), metP4.X(), metP4.Y()};

			mt2obj.set_momenta(a, b, c);
			mt2obj.set_mn(0);

			IsoTrackMt2.push_back(mt2obj.get_mt2());
			IsoTrackPt.push_back(isotrackpt);
			IsoTrackPdgId.push_back(isotrackpdgid);

			if (10 < abs(isotrackpdgid) && abs(isotrackpdgid) < 14) {
				IsoTrackHadronicDecay.push_back(false); //leptonic
				isoMt2Cut = 80;
			} else {
				IsoTrackHadronicDecay.push_back(true); //hadronic track
				isoMt2Cut = 60;
			}

			if (IsoTrackMt2.back() <= isoMt2Cut) { IsoTrackVeto = true; break;}
		}
	}
	*/

	//if (product.leptonPt != -999) {
	//	std::string cutName("DeltaPhi Calculated!");
	//	cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
	//} else {
	//	cutflow.passed = false;
	//}
}

void DeltaPhiProducer::EndJob(TFile &file) {}
