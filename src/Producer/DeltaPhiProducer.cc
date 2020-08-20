#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/DeltaPhiProducer.h>

DeltaPhiProducer::DeltaPhiProducer(TTreeReader& reader):
	BaseProducer(&reader)
	{}


float DeltaPhiProducer::DeltaPhi(ROOT::Math::PtEtaPhiMVector v1, ROOT::Math::PtEtaPhiMVector v2) {
	float dPhi = v1.Phi() - v2.Phi();
	while (dPhi >= TMath::Pi()) dPhi -= TMath::TwoPi();
	while (dPhi < -TMath::Pi()) dPhi += TMath::TwoPi();
	return dPhi;
}

void DeltaPhiProducer::BeginJob(TTree* tree, bool &isData) {
	//Set data bool
	this->isData = isData;

	//Set TTreeReader for genpart and trigger obj from BaseProducer
	SetCollection(this->isData);

	//Set Branches of output tree
	tree->Branch("HT", &HT);
	tree->Branch("LT", &LT);
	tree->Branch("LP", &LP);
	tree->Branch("DeltaPhi", &deltaPhi);
	tree->Branch("dPhi", &dPhi);
	tree->Branch("WBosonMt", &wBosonMt);
	tree->Branch("signalRegion", &signalRegionCSV);
	tree->Branch("signalRegionCSV", &signalRegionCSV);
	tree->Branch("signalRegionDF", &signalRegionDF);
}

void DeltaPhiProducer::Produce(CutFlow& cutflow, Susy1LeptonProduct *product) {
	//Initialize all variables as -999
	HT = -999;
	LT = -999;
	LP = -999;
	deltaPhi = -999;
	dPhi = -999;
	wBosonMt = -999;
	signalRegionCSV = -999;
	signalRegionDF = -999;

	if (product->leptonPt != -999 && product->metPt != -999) {
		ROOT::Math::PtEtaPhiMVector leptonP4 = ROOT::Math::PtEtaPhiMVector(product->leptonPt, product->leptonEta, product->leptonPhi, product->leptonMass);
		ROOT::Math::PtEtaPhiMVector metP4 = ROOT::Math::PtEtaPhiMVector(product->metPt, 0, product->metPhi, 0);
		ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + metP4;

		LT = product->leptonPt + product->metPt;
		deltaPhi = DeltaPhi(leptonP4, wBosonP4);
		dPhi = std::abs(deltaPhi);
		LP = product->leptonPt / wBosonP4.Pt() * std::cos(deltaPhi);
		wBosonMt = wBosonP4.Mt();

		HT = 0;
		for (unsigned int i = 0; i < product->nJet; i++) {
			HT += product->jetPt.at(i);
		}

		if (product->nMediumCSVBTagJet == 0) {// 0-B SRs -- simplified dPhi
			if (LT < 250) {signalRegionCSV = 0;}
			else if (LT > 250) {signalRegionCSV = dPhi > 0.75;}
			// BLIND data
			if (isData && product->nJet >= 5) { signalRegionCSV = -signalRegionCSV;}
		} else if (product->nJet < 99) {// Multi-B SRs
			if (LT < 250) { signalRegionCSV = 0;}
			else if (LT < 350) { signalRegionCSV = dPhi > 1.0;}
			else if (LT < 600) { signalRegionCSV = dPhi > 0.75;}
			else if (LT > 600) { signalRegionCSV = dPhi > 0.5;}

			// BLIND data
			if (isData && product->nJet >= 6) { signalRegionCSV = -signalRegionCSV;}
		}

		if (product->nMediumDFBTagJet == 0) {//0-B SRs -- simplified dPhi
			if (LT < 250) {signalRegionDF = 0;}
			else if (LT > 250) {signalRegionDF = dPhi > 0.75;}
			// BLIND data
			if (isData && product->nJet >= 5) { signalRegionDF = -signalRegionDF;}
		} else if (product->nJet < 99) {// Multi-B SRs
			if (LT < 250) { signalRegionDF = 0;}
			else if (LT < 350) { signalRegionDF = dPhi > 1.0;}
			else if (LT < 600) { signalRegionDF = dPhi > 0.75;}
			else if (LT > 600) { signalRegionDF = dPhi > 0.5;}

			// BLIND data
			if (isData && product->nJet >= 6) { signalRegionDF = -signalRegionDF;}
		}
	}

	if (product->leptonPt != -999) {
		std::string cutName("DeltaPhi Calculated!");
		cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
	} else {
		cutflow.passed = false;
	} // This should probably check if both the lepton producer and the jet producer were successful
}

void DeltaPhiProducer::EndJob(TFile* file) {
}
