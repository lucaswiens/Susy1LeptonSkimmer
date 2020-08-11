#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/DeltaPhiProducer.h>

DeltaPhiProducer::DeltaPhiProducer(TTreeReader& reader):
	BaseProducer(&reader)
	{}


float DeltaPhiProducer::DeltaPhi(ROOT::Math::PtEtaPhiMVector v1, ROOT::Math::PtEtaPhiMVector v2){
	float dPhi = v1.Phi() - v2.Phi();
	while (dPhi >= TMath::Pi()) dPhi -= TMath::TwoPi();
	while (dPhi < -TMath::Pi()) dPhi += TMath::TwoPi();
	return dPhi;
}

void DeltaPhiProducer::BeginJob(TTree* tree, bool &isData){
	//Set data bool
	this->isData = isData;

	//Set TTreeReader for genpart and trigger obj from BaseProducer
	SetCollection(this->isData);

	//Set Branches of output tree
	tree->Branch("HT", &HT);
	tree->Branch("LT", &LT);
	tree->Branch("LP", &LP);
	tree->Branch("DeltaPhi", &deltaPhi);
	tree->Branch("deltaPhi", &deltaPhi);
	tree->Branch("dPhi", &dPhi);
	tree->Branch("WBosonMt", &wBosonMt);
	//tree->Branch("Lepton", &);
}

void DeltaPhiProducer::Produce(CutFlow& cutflow, Susy1LeptonProduct *product){
	//Initialize all variables as -999
	HT = -999;
	LT = -999;
	LP = -999;
	deltaPhi = -999;
	dPhi = -999;
	wBosonMt = -999;

	if(product->leptonPt != -999 && product->metPt != -999){
		ROOT::Math::PtEtaPhiMVector leptonP4 = ROOT::Math::PtEtaPhiMVector(product->leptonPt, product->leptonEta, product->leptonPhi, product->leptonMass);
		ROOT::Math::PtEtaPhiMVector metP4 = ROOT::Math::PtEtaPhiMVector(product->metPt, 0, product->metPhi, 0);
		ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + metP4;

		LT = product->leptonPt + product->metPt;
		deltaPhi = DeltaPhi(leptonP4, wBosonP4);
		dPhi = std::abs(deltaPhi);
		LP = product->leptonPt / wBosonP4.Pt() * std::cos(deltaPhi);
		wBosonMt = wBosonP4.Mt();

		HT = 0;
		for (unsigned int i = 0; i < product->nJet; i++){
			HT += product->jetPt.at(i);
		}

	}

	if (product->leptonPt != -999){
		std::string cutName("DeltaPhi Calculated!");
		cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
	} else {
		cutflow.passed = false;
	} // This should probably check if both the lepton producer and the jet producer were successful
}

void DeltaPhiProducer::EndJob(TFile* file){
}
