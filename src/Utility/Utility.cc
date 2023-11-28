#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Utility.h>

std::vector<std::string> Utility::GetKeys(const boost::property_tree::ptree &tree, const std::string path){
	std::vector<std::string> keys;
	boost::property_tree::ptree node = tree.get_child(path);

	for(const std::pair<const std::string, boost::property_tree::ptree> &p : node){
		keys.push_back(p.first);
	}

	return keys;
};

float Utility::DeltaPhi(const float &phi1, const float &phi2) {
	float deltaPhi = phi1 - phi2;
	while (deltaPhi >  M_PI) deltaPhi -= 2 * M_PI;
	while (deltaPhi < -M_PI) deltaPhi += 2 * M_PI;
	return deltaPhi;
}

float Utility::DeltaR(const float &eta1, const float &phi1, const float &eta2, const float &phi2) {
	const float &deltaEta = eta1 - eta2;
	const float &deltaPhi = Utility::DeltaPhi(phi1, phi2);
	return std::sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);
}

int Utility::GetChargeFromPdgId(const int &pdgId) {
	if (10 < std::abs(pdgId) && std::abs(pdgId) < 14) { // https://twiki.cern.ch/twiki/bin/view/Main/PdgId
		// leptonic tracks have negative charge for positive pdgId
		return (pdgId < 0) - (pdgId > 0);
	}
	// hadronic tracks have positive charge for positive pdgId
	return (pdgId > 0) - (pdgId < 0);
}

float Utility::GetWeight(const float &x, TH1F* histogram) {
	if(histogram==NULL) {
		std::cerr << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t bin = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	return histogram->GetBinContent(bin);
}

float Utility::GetWeightErr(const float &x, TH1F* histogram) {
	if(histogram==NULL) {
		std::cerr << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t bin = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	return histogram->GetBinError(bin);
}

float Utility::Get2DWeight(const float &x, const float &y, TH2F* histogram) {
	if(histogram==NULL) {
		std::cerr << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t binX = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	Int_t binY = std::max(1, std::min(histogram->GetNbinsY(), histogram->GetYaxis()->FindBin(y)));
	return histogram->GetBinContent(binX, binY);
}

float Utility::Get2DWeightErr(const float &x, const float &y, TH2F* histogram) {
	if(histogram==NULL) {
		std::cerr << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t binX = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	Int_t binY = std::max(1, std::min(histogram->GetNbinsY(), histogram->GetYaxis()->FindBin(y)));
	return histogram->GetBinError(binX, binY);
	return 1.;
}

std::vector<std::string> Utility::SplitString(const std::string &text, const std::string &delimiter, const bool removeSpaces) {
	int position = text.find(delimiter);
	int initialPosition = 0;
	std::vector<std::string> values;

	// Decompose statement
	while(position != std::string::npos) {
		std::string subString = text.substr(initialPosition, position - initialPosition);
		if (removeSpaces) {
			subString.erase(std::remove_if(subString.begin(), subString.end(), ::isspace), subString.end());
		}

		values.push_back(subString);
		initialPosition = position + 1;

		position = text.find(delimiter, initialPosition);
	}

	// Add the last one
	std::string subString = text.substr(initialPosition, std::min(position, (int)text.size()) - initialPosition + 1);
	if (removeSpaces) {
		subString.erase(std::remove_if(subString.begin(), subString.end(), ::isspace), subString.end());
	}
	values.push_back(subString);

	return values;
}
