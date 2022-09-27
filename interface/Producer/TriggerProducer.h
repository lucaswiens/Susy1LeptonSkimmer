#ifndef TRIGGERPRODUCER_H
#define TRIGGERPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

class TriggerProducer : public BaseProducer {
	private:
		std::vector<std::string> triggerNames;
	public:
		TriggerProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product, const std::vector<std::string> &triggerNames);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

