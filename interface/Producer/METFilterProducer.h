#ifndef METFILTERPRODUCER_H
#define METFILTERPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

class METFilterProducer : public BaseProducer {
	private:
	public:
		std::string Name = "MetFilerProducer";
		METFilterProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

