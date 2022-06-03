#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Utility.h>

template <typename T>
std::vector<T> Utility::GetVector(const boost::property_tree::ptree &tree, const std::string &key){
	std::vector<T> vec;

	for(std::pair<std::string, boost::property_tree::ptree> values : tree.get_child(key)){
		vec.push_back(values.second.get_value<T>());
	}

	return vec;
};

std::vector<std::string> Utility::GetKeys(const boost::property_tree::ptree &tree, const std::string path){
	std::vector<std::string> keys;
	boost::property_tree::ptree node = tree.get_child(path);

	for(const std::pair<const std::string, boost::property_tree::ptree> &p : node){
		keys.push_back(p.first);
	}

	return keys;
};

double Utility::DeltaPhi(double phi1, double phi2) {
	double deltaPhi = phi1 - phi2;
	while (deltaPhi >  M_PI) deltaPhi -= 2 * M_PI;
	while (deltaPhi < -M_PI) deltaPhi += 2 * M_PI;
	return deltaPhi;
}

double Utility::DeltaR(const double &eta1, const double &phi1, const double &eta2, const double &phi2) {
	double deltaEta = eta1 - eta2;
	double deltaPhi = Utility::DeltaPhi(phi1, phi2);
	return std::sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);
}
