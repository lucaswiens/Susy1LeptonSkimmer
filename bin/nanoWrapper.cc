#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/NanoSkimmer.h>

#include <string>

std::vector<std::string> SplitString(const std::string& splitString, const std::string& delimeter){
	//Function which handles splitting of string input
	std::vector<std::string> splittedString;
	std::string string;
	std::istringstream splittedStream(splitString);
	while (std::getline(splittedStream, string, delimeter.c_str()[0])){
		splittedString.push_back(string);
	}

int main(int argc, char* argv[]) {
	//Extract informations of command line
	std::string fileName = std::string(argv[1]);
	bool isData = std::string(argv[2]) == "True" ? true : false;
	int era = std::stoi(std::string(argv[3]));
	float xSec = std::stof(std::string(argv[4]));
	std::string outName = std::string(argv[5]);

	NanoSkimmer skimmer(fileName, isData);
	skimmer.EventLoop(xSec, era);
	skimmer.WriteOutput(outName);
}


