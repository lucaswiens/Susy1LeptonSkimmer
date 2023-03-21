#ifndef AK8TAGCSVREADER
#define AK8TAGCSVREADER

#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <functional>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Utility.h>

class Ak8TagCSVReader{
	private:
		std::map<std::string, std::vector<std::pair<float, float>>> ptRanges;
		std::map<std::string, std::vector<float>> sfMap;
	public:
		Ak8TagCSVReader(){};
		Ak8TagCSVReader(const std::string &fileName, const bool &isFastSim);
		float GetSf(const std::string &selectionString, const char &syst, const float &pt);
};

#endif

