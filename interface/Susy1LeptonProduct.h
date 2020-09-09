#pragma once

class Susy1LeptonProduct {
	public:
		float leptonPt = -999;
		float leptonEta = -999;
		float leptonPhi = -999;
		float leptonMass = -999;
		int leptonPdgId = -999;
		int leptonCharge = -999;

		unsigned int nJet = 0;
		unsigned int nLooseCSVBTagJet = 0;
		unsigned int nMediumCSVBTagJet = 0;
		unsigned int nTightCSVBTagJet = 0;
		unsigned int nLooseDFBTagJet = 0;
		unsigned int nMediumDFBTagJet = 0;
		unsigned int nTightDFBTagJet = 0;
		std::vector<float> jetPt;
		std::vector<float> jetPhi;
		std::vector<float> jetEta;
		std::vector<float> jetMass;
		std::vector<bool> jetMediumCSVBTag;
		std::vector<bool> jetMediumDFBTag;
		//Required for Signal Region Flag

		float metPt = -999;
		float metPhi = -999;

		int nGenPart = -999;
		int nGenWBoson = -999;
		std::vector<float> genWBosonPt;
		std::vector<float> genWBosonEta;
		std::vector<float> genWBosonPhi;
		std::vector<float> genWBosonMass;

		int nGenTop = -999;
		std::vector<float> genTopPt;
		std::vector<float> genTopEta;
		std::vector<float> genTopPhi;
		std::vector<float> genTopMass;
};
