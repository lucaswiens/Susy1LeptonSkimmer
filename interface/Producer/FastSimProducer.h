#ifndef FASTSIMPRODUCER_H
#define FASTSIMPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/BTagCSVReader.h>

#include <TH2F.h>

#include <correction.h>

#include <boost/optional/optional.hpp>

class FastSimProducer : public BaseProducer {
	private:
		//Era information
		std::string era, electronEraAlias, muonEraAlias;
		std::string muonTriggName;
		std::string run;
		int stopPdgId, gluinoPdgId, neutralinoPdgId, charginoPdgId;

		// Histograms with SF information for all used working points
		TFile *electronSfFile, *muonSfFile;
		std::shared_ptr<TH2F> electronVetoSf, electronTightSf,
			electronVetoMvaSf, electronTightMvaSf,
			muonLooseSf, muonMediumSf;

		// BTag ScaleFactors requrie the btag csv reader
		std::unique_ptr<BTagCSVReader> btagCsvReader;
	public:
		FastSimProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector, TFile &outputFile);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

