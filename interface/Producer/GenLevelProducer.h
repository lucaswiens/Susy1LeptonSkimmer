#ifndef GENLEVELPRODUCER_H
#define GENLEVELPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <TMath.h>

class GenLevelProducer : public BaseProducer {
	private:
	public:
		GenLevelProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

