#ifndef EVENTINFORMATIONPRODUCER_H
#define EVENTINFORMATIONPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <correction.h>

#include <boost/optional/optional.hpp>

class EventInformationProducer : public BaseProducer {
	private:
		//std::string goldenJsonString;
		pt::ptree goldenJsonTree;
	public:
		EventInformationProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

