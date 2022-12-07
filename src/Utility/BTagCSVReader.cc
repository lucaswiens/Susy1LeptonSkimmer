#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/BTagCSVReader.h>

BTagCSVReader::BTagCSVReader(const std::string &fileName, const bool &isFastSim) {
	//Open file and check if it is open
	std::ifstream CSV(fileName);
	if (!CSV.is_open()) throw std::runtime_error("File not exists: " + fileName);

	//Header, ignofire it
	int lineNr = 0;
	std::string line;
	std::getline(CSV, line);

	//Loop over file
	while (std::getline(CSV, line)) {
		//Split line and get respective parameter
		std::vector<std::string> lineVector = Utility::SplitString(line, ",", true);

		std::string type = lineVector.at(1),
			sysType = lineVector.at(2);
		int wp = std::stoi(lineVector.at(0)),
			jetFlavour = std::stoi(lineVector.at(3));
		std::string formula = lineVector.at(10).substr(1, lineVector.at(10).size()-2);
		float ptMin = std::stof(lineVector.at(6)),
			ptMax = std::stof(lineVector.at(7));

		if (isFastSim) {
			if (type != "fastsim" or jetFlavour != 0) continue;
		} else {
			if (type != "comb" or jetFlavour != 0) continue;
		}
		if (type == "iterativefit") break;

		std::shared_ptr<TF1> func = std::make_shared<TF1>(("Func" + std::to_string(lineNr)).c_str(), formula.c_str(), ptMin, ptMax);
		ptRange[{wp, sysType == "central" ? 0 : sysType == "up" ? 1 : 2}].push_back({ptMin, ptMax});
		SF[{wp, sysType == "central" ? 0 : sysType == "up" ? 1 : 2}].push_back(std::move(func));

		lineNr++;
	}
}

float BTagCSVReader::Get(const float &pt, const char &wp) {
	int i=0;
	int wpInt = wpMap.at(wp);
	for(const std::pair<float, float> range : ptRange[{wpInt, 0}]) {
		if (range.first < pt and pt < range.second) {
			return SF[{wpInt, 0}][i]->Eval(pt);
		}
		i++;
	}
	return 1.;
}

float BTagCSVReader::GetUp(const float &pt, const char &wp) {
	int i=0;
	int wpInt = wpMap.at(wp);
	for(const std::pair<float, float> range : ptRange[{wpInt, 1}]) {
		if (range.first < pt and pt < range.second) {
				return SF[{wpInt, 1}][i]->Eval(pt);
		}
		i++;
	}
	return 1.;
}

float BTagCSVReader::GetDown(const float &pt, const char &wp) {
	int i=0;
	int wpInt = wpMap.at(wp);
	for(const std::pair<float, float> range : ptRange[{wpInt, 2}]) {
		if (range.first < pt and pt < range.second) {
				return SF[{wpInt, 2}][i]->Eval(pt);
		}
		i++;
	}
	return 1.;
}
