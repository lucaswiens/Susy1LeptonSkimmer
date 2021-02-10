#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/NanoSkimmer.h>

#include <string>
#include <sstream>

float CrossSection(std::string);

int main(int argc, char* argv[]) {
	//Extract informations of command line
	std::string fileName = std::string(argv[1]);
	std::string outName  = std::string(argv[2]);
	bool isData          = std::string(argv[3]) == "True" ? true : false;
	bool doSystematics   = std::string(argv[4]) == "True" ? true : false;
	int era              = std::stoi(argv[5]);
	char runPeriod = (char)*argv[6]; //If it is MC, then runPeriod does not matter
	double xSection = std::stod(argv[7]);

	int nMaxEvents;
	if(argc == 9) {
		nMaxEvents = std::stoi(std::string(argv[8]));
	} else {
		nMaxEvents = -999;
	}

	NanoSkimmer skimmer(fileName, outName, isData, doSystematics, era, runPeriod, xSection);
	skimmer.EventLoop(nMaxEvents);
	skimmer.WriteOutput();

	return 0;
}
