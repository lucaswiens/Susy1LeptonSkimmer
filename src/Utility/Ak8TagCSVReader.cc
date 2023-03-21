#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Ak8TagCSVReader.h>

/*#################################################################################
#   TODO:                                                                         #
#   FastSim SF do not exist yet.                                                  #
#   We will have to calculate those ourselves and read them in during this step   #
#################################################################################*/

Ak8TagCSVReader::Ak8TagCSVReader(const std::string &fileName, const bool &isFastSim) {
	//Open file and check if it is open
	std::ifstream CSV(fileName);
	if (!CSV.is_open()) throw std::runtime_error("File not exists: " + fileName);

	//Header, ignofire it
	std::string line;
	std::getline(CSV, line);

	//Loop over file
	while (std::getline(CSV, line)) {
		std::vector<std::string> lineVector = Utility::SplitString(line, ",", true);
		std::string object         = lineVector.at(0);
		std::string year           = lineVector.at(1);
		std::string version        = lineVector.at(2);
		std::string mistaggingRate = lineVector.at(3);
		float ptMin                = std::stoi(lineVector.at(4));
		float ptMax                = std::stoi(lineVector.at(5));
		float scaleFactor          = std::stof(lineVector.at(6));
		float scaleFactorUp        = std::stof(lineVector.at(7));
		float scaleFactorDown      = std::stof(lineVector.at(8));

		// map sounds nice but makes the loop not so nice..
		std::string selectionString = object + "." + year + "." + version + "." + mistaggingRate;
		if (ptRanges.find(selectionString) == ptRanges.end() ) { // selectionString not found in ptRanges
			ptRanges.insert({selectionString, {{ptMin, ptMax}}});
		} else {
			ptRanges.at(selectionString).push_back({ptMin, ptMax});
		}

		if (sfMap.find(selectionString + ".n") == sfMap.end() ) { // selectionString not found in sfMap
			sfMap.insert({selectionString + ".n", {scaleFactor}});
			sfMap.insert({selectionString + ".u", {scaleFactorUp}});
			sfMap.insert({selectionString + ".d", {scaleFactorDown}});
		} else {
			sfMap.at(selectionString + ".n").push_back(scaleFactor);
			sfMap.at(selectionString + ".u").push_back(scaleFactorUp);
			sfMap.at(selectionString + ".d").push_back(scaleFactorDown);
		}
	}
}

float Ak8TagCSVReader::GetSf(const std::string &selectionString, const char &syst, const float &pt) {
	float sf = syst == 'n' ? 1 : 0;
	const std::vector<std::pair<float, float>> &ptRangeVector = ptRanges.at(selectionString);
	for(int iPt = 0; iPt < ptRangeVector.size(); iPt++) {
		if (ptRangeVector.at(iPt).first <= pt && pt < ptRangeVector.at(iPt).second) {
			sf = sfMap.at(selectionString + "." + syst).at(iPt);
			break;
		}
	}
	return sf;
}
