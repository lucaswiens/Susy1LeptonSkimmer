#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

#include <vector>
#include <string>
#include <chrono>
#include <memory>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

class NanoSkimmer{
	private:
		//Measure execution time
		std::chrono::steady_clock::time_point start;
		std::chrono::steady_clock::time_point end;

		//Input
		std::string inFile;
		bool isData;
		//Vector with wished producers
		std::vector<std::shared_ptr<BaseProducer>> producers;

		int finalNumberOfEvents;

		//Output Tree
		TTree* outputTree;

		//Vector of cutflow histograms for each analysis
		CutFlow cutflow;


		//Progress bar function
		void ProgressBar(const int &progress, const int &rate);

		//Configure analysis modules
		void Configure(const float &xSec, const int &era, const char &runPeriod, TTreeReader& reader);


	public:
		NanoSkimmer();
		NanoSkimmer(const std::string &inFile, const bool &isData);
		void EventLoop(const float &xSec = 1., const int &era = 2016, const char &runPeriod = 'm', const int &nMaxEvents = -999);
		void WriteOutput(const std::string &outFile);
};
