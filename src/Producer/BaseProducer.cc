#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

//BaseProducer::BaseProducer(): isNANO(false){}
BaseProducer::BaseProducer(TTreeReader* reader): reader(reader){}

float BaseProducer::DeltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2){
	return std::sqrt(std::pow(phi1 - phi2, 2) + std::pow(eta1 - eta2, 2));
}

void BaseProducer::SetCollection(bool &isData){
	if(!isData){
		genPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_pt");
		genEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_eta");
		genPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_phi");
		genMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_mass");
		genID = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_pdgId");
		genMotherIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_genPartIdxMother");
		genStatus = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_statusFlags");
		eleGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_genPartIdx");
		muonGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_genPartIdx");
	}

	run = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "run");

	trigObjPt = std::make_unique<TTreeReaderArray<float>>(*reader, "TrigObj_pt");
	trigObjEta = std::make_unique<TTreeReaderArray<float>>(*reader, "TrigObj_eta");
	trigObjPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "TrigObj_phi");
	trigObjID = std::make_unique<TTreeReaderArray<int>>(*reader, "TrigObj_id");
	trigObjFilterBit = std::make_unique<TTreeReaderArray<int>>(*reader, "TrigObj_filterBits");
}

int BaseProducer::FirstCopy(const int& index, const int& pdgID){
	int daughterIdx = index;
	int motherIdx = genMotherIdx->At(index);

	while(abs(genID->At(motherIdx)) == pdgID){
		daughterIdx = motherIdx;
		motherIdx = genMotherIdx->At(daughterIdx);
	}

	return daughterIdx;
}

std::tuple<int, int, int> BaseProducer::SetGenParticles(const float& pt, const float& eta, const float& phi, const int &i, const int &pdgID){
	std::tuple<int, int, int> IDs = std::make_tuple(-99., -99., -99.);

	std::map<int, int> genIndex;
	genIndex = {{11, eleGenIdx->At(i)}, {13, muonGenIdx->At(i)}};

	bool isgenMatched = genIndex[pdgID] != -1;

	//Check if gen matched particle exist
	if(isgenMatched){
		int index=0, motherIdx=0;

		index = FirstCopy(genIndex[pdgID], pdgID);

		const int& motherID = abs(genID->At(genMotherIdx->At(index)));

		motherIdx = FirstCopy(genMotherIdx->At(index), motherID);

		const int& grandMotherID = abs(genID->At(genMotherIdx->At(motherIdx)));

		IDs = std::make_tuple(pdgID, motherID, grandMotherID);
	}

	return IDs;
}
