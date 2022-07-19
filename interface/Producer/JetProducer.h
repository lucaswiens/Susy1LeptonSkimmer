#ifndef JETPRODUCER_H
#define JETPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <JetMETCorrections/Modules/interface/JetCorrectionProducer.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>

//#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/BTagCSVReader.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Utility.h>

#include <correction.h>

class JetProducer : public BaseProducer {

	private:
		// Configuration Variables
		bool isData, doSystematics, isUp, isJERSystematic, isJECSystematic, preVFP;
		double jetPtCut, jetEtaCut;

		// Top quark and W boson mass
		double pdgTopMass     = 172.76,
			pdgWBosonMass = 80.379;

		// JEC and JER from json library
		std::unique_ptr<correction::CorrectionSet> ak4CorrectionSet, ak8CorrectionSet, jmeCorrectionSet;
		std::string era, dataType, runPeriod, jerVersion, ak4Algorithm, ak8Algorithm;


		//JEC/JER systematics
		std::shared_ptr<JetCorrectionUncertainty> jetCorrectionUncertainty;
		bool isJERsyst = false;

		// Btag Map
		std::map<char, double> deepCsvBTagMap, deepJetBTagMap;


		//Get jet energy correction
		//FactorizedJetCorrector *jetCorrector;
		//void SetCorrector(const char &runPeriod);

		double CorrectEnergy(const double &pt, const double &eta, const double &rho, const double &area);
		double SmearEnergy(DataReader &dataReader, const double &jetPtCorrected, const bool &isAk4);

	public:
		std::string Name = "JetProducer";
		JetProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

