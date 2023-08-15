#ifndef SIGNALPRODUCER_H
#define SIGNALPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/BTagCSVReader.h>

#include <TH2F.h>

#include <correction.h>

#include <boost/optional/optional.hpp>

class SignalProducer : public BaseProducer {
	private:
		int stopPdgId, gluinoPdgId, neutralinoPdgId, charginoPdgId, wBosonPdgId;

		// Cross Section Json File
		pt::ptree xSectionTree;
		std::shared_ptr<TH2F> numberOfGenEvents, numberOfGenEventsIsrWeighted, numberOfT5qqqqWWGenEvents, numberOfT5qqqqWWGenEventsIsrWeighted;
	public:
		SignalProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector, TFile &outputFile);
		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

