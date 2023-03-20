#ifndef SCALEFACTORPRODUCER_H
#define SCALEFACTORPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Ak8TagCSVReader.h>

#include <TH2F.h>

#include <correction.h>

class ScaleFactorProducer : public BaseProducer {
	private:
		//Era information
		std::string deepAk8EraAlias, electronEraAlias, muonEraAlias;
		std::string muonTriggName;
		std::string run;

		bool isData;

		//Hist of MC btag efficiency
		std::shared_ptr<TH2F> bTagEffBLooseDeepJet, bTagEffBMediumDeepJet, bTagEffBTightDeepJet,
			bTagEffCLooseDeepJet, bTagEffCMediumDeepJet, bTagEffCTightDeepJet,
			bTagEffLightLooseDeepJet, bTagEffLightMediumDeepJet, bTagEffLightTightDeepJet,
			bTagEffBLooseDeepCsv, bTagEffBMediumDeepCsv, bTagEffBTightDeepCsv,
			bTagEffCLooseDeepCsv, bTagEffCMediumDeepCsv, bTagEffCTightDeepCsv,
			bTagEffLightLooseDeepCsv, bTagEffLightMediumDeepCsv, bTagEffLightTightDeepCsv,
			bTotal, cTotal, lightTotal;

		std::unique_ptr<correction::CorrectionSet> electronSf, muonSf, bTagSf;
		std::vector<std::string> bTagSyst, bTagSystLight;

		Ak8TagCSVReader ak8TagCSVReader;
	public:
		ScaleFactorProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, const std::string &eraSelector, const bool &isFastSim, TFile &outputFile);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

