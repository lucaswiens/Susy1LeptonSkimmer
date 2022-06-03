#ifndef PILEUPWEIGHTPRODUCER_H
#define PILEUPWEIGHTPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

class PileUpWeightProducer : public BaseProducer {
	private:
		//Check if it is data or MC
		bool isData, doSystematics;
		int era;

		//Histogram file names
		std::map<int, TString> pileupPathData, pileupPathMC;
		TH1D *pileupRatio, *pileupRatioPlus, *pileupRatioMinus;
		WeightCalculator *wc;

		//Variables to be stored in the output tree
		int nPV;
		double weight, weightPlus, weightMinus;

		//TTreeReader Values for NANO AOD analysis
		//std::unique_ptr<TTreeReaderValue<float>> pvNumber;
		//std::unique_ptr<TTreeReaderArray<float>> muonPt;

	public:
		std::string Name = "PileUpWeightProducer";
		PileUpWeightProducer(const int &era);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

