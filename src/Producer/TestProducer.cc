#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/TestProducer.h>

TestProducer::TestProducer(const int& era, const float& ptCut, const float& etaCut, TTreeReader& reader):
	BaseProducer(&reader),
	era(era),
	ptCut(ptCut),
	etaCut(etaCut)
	{}

void TestProducer::BeginJob(TTree* tree, bool &isData, const bool& isSyst){
	//Set data bool
	this->isData = isData;
	this->isSyst = isSyst;

	//Initiliaze TTreeReaderValues then using NANO AOD
	elePt = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_pt");
	eleEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_eta");
	elePhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_phi");

	//Set TTreeReader for genpart and trigger obj from BaseProducer
	SetCollection(this->isData);


	//Variable name mapping to branch name
	floatVar = {
			{"Pt", Pt},
			{"Eta", Eta},
			{"Phi", Phi},
	};

	//Set Branches of output tree
	tree->Branch("Electron_Size", &nElectrons);

	for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
		tree->Branch(("Electron_" + var.first).c_str(), &var.second);
	}
}

void TestProducer::Produce(CutFlow cutflow){
	//Clear variables vector
	for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
		var.second.clear();
	}

	const int& eleSize = elePt->GetSize();

	//Loop over all electrons
	for(int i = 0; i < eleSize; i++){
		const float& pt = elePt->At(i);
		const float& eta = eleEta->At(i);
		const float& phi = elePhi->At(i);

		if(pt > ptCut && abs(eta) < etaCut){
			//Electron four momentum components
			Pt.push_back(pt);
			Eta.push_back(eta);
			Phi.push_back(phi);
		}
	}

	nElectrons = Pt.size();

	if(cutflow.nMinElectron <= Pt.size()){
		if(cutflow.nMinElectron!=0 and cutflow.passed){
			std::string cutName("N_{e} >= " + std::to_string(cutflow.nMinElectron) + " (no iso/ID req)");
			cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
		}
	} else {
		cutflow.passed = false;
	}
}

void TestProducer::EndJob(TFile* file){
}
