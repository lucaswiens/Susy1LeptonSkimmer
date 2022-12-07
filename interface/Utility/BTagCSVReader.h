#pragma once

#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Utility.h>

#include <TF1.h>

class BTagCSVReader{
	private:
		std::map<std::pair<int, int>, std::vector<std::shared_ptr<TF1>>> SF;
		std::map<std::pair<int, int>, std::vector<std::pair<float, float>>> ptRange;
		std::map<char, int> wpMap = { // To make it possible for mistakes to happen, the csv files encode the loose, medium and tight WP as int, instead of char
			{'L', 0}, // Loose
			{'M', 1}, // Medium
			{'T', 2}  // Tight
		};
	public:
		BTagCSVReader(){};
		BTagCSVReader(const std::string &fileName, const bool &isFastSim);

		float Get(const float &pt, const char &wp);
		float GetUp(const float &pt, const char &wp);
		float GetDown(const float &pt, const char &wp);
};

