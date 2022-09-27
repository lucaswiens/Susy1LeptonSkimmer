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
		float jetPtCut, jetEtaCut;

		// Top quark and W boson mass
		float pdgTopMass     = 172.76,
			pdgWBosonMass = 80.379;

		// JEC and JER from json library
		std::unique_ptr<correction::CorrectionSet> ak4CorrectionSet, ak8CorrectionSet, jmeCorrectionSet;
		std::string era, dataType, runPeriod, jerVersion, ak4Algorithm, ak8Algorithm;


		//JEC systematics
		std::vector<std::string> jecSystematics;
		std::map<std::string, std::shared_ptr<JetCorrectionUncertainty>> ak4CorrectionUncertainty, ak8CorrectionUncertainty;

		//Classes for reading jet energy SF (Only needed for AK8 jets as the ak4 can be done with the correction set)
		JME::JetParameters jetParameter;
		JME::JetResolution resolutionAk8;
		JME::JetResolutionScaleFactor resolutionSfAk8;

		// Btag Map
		std::map<char, float> deepCsvBTagMap, deepJetBTagMap;

		float CorrectEnergy(const float &pt, const float &eta, const float &rho, const float &area);
		float SmearEnergy(DataReader &dataReader, const float &jetPtCorrected, const bool &isAk4);
		float SmearFatEnergy(DataReader &dataReader, const float &jetPtCorrected);

	public:
		JetProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

