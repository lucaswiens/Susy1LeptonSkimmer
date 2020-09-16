#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/NanoSkimmer.h>

#include <string>
#include <sstream>

float CrossSection(std::string);

int main(int argc, char* argv[]) {
	//Extract informations of command line
	std::string fileName = std::string(argv[1]);
	bool isData          = std::string(argv[2]) == "True" ? true : false;
	int era              = std::stoi(std::string(argv[3]));
	float xSec           = std::stof(std::string(argv[4]));
	std::string outName  = std::string(argv[5]);

	char runPeriod; //If it is MC, then runPeriod does not matter
	if (isData) {
		runPeriod = (char)*argv[6];
	}

	int nMaxEvents;
	if(argc == 8) {
		nMaxEvents = std::stoi(std::string(argv[7]));
	} else {
		nMaxEvents = -999;
	}

	NanoSkimmer skimmer(fileName, isData);
	skimmer.EventLoop(xSec, era, runPeriod, nMaxEvents);
	skimmer.WriteOutput(outName);
}
