#ifndef UTILITY_H
#define UTILITY_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include<iostream>
#include <cmath>

#include<TH1F.h>
#include<TH2F.h>

namespace Utility{
	template <typename T>
	std::vector<T> GetVector(const boost::property_tree::ptree &tree, const std::string &key) {
		std::vector<T> vec;
		for(std::pair<std::string, boost::property_tree::ptree> values : tree.get_child(key)) {
			vec.push_back(values.second.get_value<T>());
		}
		return vec;
	}

	template <typename T>
	void RemoveByIndex(T &array, const int &removeIndex, const int &arraySize) {
		std::copy(array->begin() + removeIndex + 1, // copy everything starting here
			array->begin() + arraySize,         // and ending here, not including it,
			array->begin() + removeIndex        // to this destination
		);
	}

	template <typename T>
	void SortByIndex(T &array, const std::vector<int> &indices, const int &arraySize) {
		T tmp = array;
		for (int i = 0; i < arraySize; i++) {
			array[i] = tmp[indices[i]];
		}
	}

	std::vector<std::string> GetKeys(const boost::property_tree::ptree &tree, const std::string path);
	std::vector<std::string> SplitString(const std::string &text, const std::string &delimiter);
	float DeltaPhi(const float &phi1, const float &phi2);
	float DeltaR(const float &eta1, const float &phi1, const float &eta2, const float &phi2);
	// Function to read bin content or error by value (e.g. muonPt or muonEta)
	float GetWeight(const float &x, TH1F* histogram);
	float GetWeightErr(const float &x, TH1F* histogram);
	float Get2DWeight(const float &x, const float &y, TH2F* histogram);
	float Get2DWeightErr(const float &x, const float &y, TH2F* histogram);
}

#endif

