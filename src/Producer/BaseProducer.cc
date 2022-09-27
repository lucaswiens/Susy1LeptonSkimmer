#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

BaseProducer::BaseProducer() {}

int BaseProducer::FirstCopy(const int &index, const int &pdgID) {
	int daughterIdx = index;
	int motherIdx = 0;//genMotherIdx->At(index);

	//while(abs(genID->At(motherIdx)) == pdgID) {
	//	daughterIdx = motherIdx;
	//	motherIdx = genMotherIdx->At(daughterIdx);
	//}

	return daughterIdx;
}

std::tuple<int, int, int> BaseProducer::SetGenParticles(const float &pt, const float &eta, const float &phi, const int &i, const int &pdgID) {
	std::tuple<int, int, int> IDs = std::make_tuple(-99., -99., -99.);

	std::map<int, int> genIndex;
	//genIndex = {{11, eleGenIdx->At(i)}, {13, muonGenIdx->At(i)}};

	bool isgenMatched = genIndex[pdgID] != -1;

	//Check if gen matched particle exist
	if (isgenMatched) {
		int index=0, motherIdx=0;

		index = FirstCopy(genIndex[pdgID], pdgID);

		//const int &motherID = abs(genID->At(genMotherIdx->At(index)));

		//motherIdx = FirstCopy(genMotherIdx->At(index), motherID);

		//const int &grandMotherID = abs(genID->At(genMotherIdx->At(motherIdx)));

		//IDs = std::make_tuple(pdgID, motherID, grandMotherID);
	}

	return IDs;
}
