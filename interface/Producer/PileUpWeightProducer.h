#ifndef PILEUPWEIGHTPRODUCER_H
#define PILEUPWEIGHTPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <correction.h>

class PileUpWeightProducer : public BaseProducer {
	private:
		std::unique_ptr<correction::CorrectionSet> pileUpCorrectionSet;
		std::string goldenJsonString;
		std::shared_ptr<TH1D> pileUpMc;
	public:
		PileUpWeightProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

