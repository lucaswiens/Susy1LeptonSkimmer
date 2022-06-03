#ifndef UTILITY_H
#define UTILITY_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace Utility{
	template <typename T>
	std::vector<T> GetVector(const boost::property_tree::ptree &tree, const std::string &key);

	std::vector<std::string> GetKeys(const boost::property_tree::ptree &tree, const std::string path);
	double DeltaPhi(double phi1, double phi2);
	double DeltaR(const double &eta1, const double &phi1, const double &eta2, const double &phi2);
}

#endif

