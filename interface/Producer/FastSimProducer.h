#ifndef FASTSIMPRODUCER_H
#define FASTSIMPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

#include <TH2F.h>

#include <correction.h>

class FastSimProducer : public BaseProducer {
	private:
		//Era information
		std::string era, electronEraAlias, muonEraAlias;
		std::string muonTriggName;
		std::string run;
		int stopPdgId, gluinoPdgId, neutralinoPdgId, charginoPdgId;

	TFile *electronSfFile, *muonSfFile;

		// Histograms with SF information for all used working points
		std::shared_ptr<TH2F> electronVetoSf, electronTightSf,
			electronVetoMvaSf, electronTightMvaSf,
			muonLooseSf, muonMediumSf;

	public:
		FastSimProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector, TFile &outputFile);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

