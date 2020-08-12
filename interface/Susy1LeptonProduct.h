#pragma once

class Susy1LeptonProduct {
	public:
		float leptonPt = -999;
		float leptonEta = -999;
		float leptonPhi = -999;
		float leptonMass = -999;

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
		//Required for Signal Region Flag

		float metPt = -999;
		float metPhi = -999;
};
