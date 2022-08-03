#ifndef SCALEFACTORPRODUCER_H
#define SCALEFACTORPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

#include <TH2D.h>

#include <correction.h>

class ScaleFactorProducer : public BaseProducer {
	private:
		//Era information
		std::string era, electronEraAlias, muonEraAlias;
		std::string muonTriggName;
		std::string run;

		bool isData;

		//Hist of MC btag efficiency
		/*
		std::shared_ptr<TH2F> bTagEffBLooseDeepJet, bTagEffBMediumDeepJet, bTagEffBTightDeepJet,
			bTagEffCLooseDeepJet, bTagEffCMediumDeepJet, bTagEffCTightDeepJet,
			bTagEffLightLooseDeepJet, bTagEffLightMediumDeepJet, bTagEffLightTightDeepJet,
			bTagEffBLooseDeepCSV, bTagEffBMediumDeepCSV, bTagEffBTightDeepCSV,
			bTagEffCLooseDeepCSV, bTagEffCMediumDeepCSV, bTagEffCTightDeepCSV,
			bTagEffLightLooseDeepCSV, bTagEffLightMediumDeepCSV, bTagEffLightTightDeepCSV,
			bTotal, cTotal, lightTotal;
		*/

		std::unique_ptr<correction::CorrectionSet> electronSf, muonSf, bTagSf;
		std::vector<std::string> bTagSyst;

	public:
		ScaleFactorProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

