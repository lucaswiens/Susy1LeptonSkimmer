#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

#include <stdlib.h>
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
		std::string inFile, outFile;
		TFile *file;
		bool isData, doSystematics;
		//Vector with wished producers
		std::vector<std::vector<std::shared_ptr<BaseProducer>>> producers;
		//std::vector<std::shared_ptr<BaseProducer>> producers;

		int finalNumberOfEvents;

		//Output Trees
		std::vector<TTree*> outputTrees;

		//Vector of cutflow histograms for each systematic and each producer
		std::vector<CutFlow> cutflows;

		//Progress bar function
		void ProgressBar(const int &progress, const int &rate);

		//Configure analysis modules
		void Configure(const int &era, const char &runPeriod, TTreeReader &reader);


	public:
		NanoSkimmer();
		NanoSkimmer(const std::string &inFile, const std::string &outFile, const bool &isData, const bool &doSystematics);
		int EventLoop(const int &era = 2016, const char &runPeriod = 'm', const int &nMaxEvents = -999);
		void WriteOutput();
};
