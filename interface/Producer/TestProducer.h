#ifndef TESTPRODUCER_H
#define TESTPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

class TestProducer: public BaseProducer {
	/*
		Test Producer which can also be used as a template for writing a new Producer
	*/
	private:
		//Check if it is data or MC
		bool isData;
		//Vector with output varirables of the output tree
		std::map<std::string, std::vector<float>&> floatVar;

		// Cut Variables
		int era;
		float ptCut, etaCut;

		// Vector for the output variables
		std::vector<float> Pt, Eta, Phi;

		char nElectrons;

		//TTreeReader Values for NANO AOD analysis
		std::unique_ptr<TTreeReaderArray<float>> elePt, eleEta, elePhi, eleIso;
	public:
		TestProducer(const int& era, const float& ptCut, const float& etaCut, TTreeReader& reader);

		void BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst=false);
		void Produce(std::vector<CutFlow>& cutflows);
		void EndJob(TFile* file);
};

#endif
