#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <JetMETCorrections/Modules/interface/JetCorrectionProducer.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/BTagCSVReader.h>

class JetProducer : public BaseProducer {

	private:
		//Check if it is data or MC
		bool isData, doSystematics, isUp, isJERSystematic, isJECSystematic ;
		char runPeriod;
		//bool isSignal // or maybe make a fastsim flag

		//Classes for reading jet energy SF
		JME::JetParameters jetParameter;
		JME::JetResolution resolution;
		JME::JetResolutionScaleFactor resolution_sf;

		//JEC/JER systematics
		JetCorrectionUncertainty *jetCorrectionUncertainty;
		bool isJERsyst = false;

		// JEC map
		std::map<int, std::vector<std::string>> jecMC, jecFastSim, jerMC, jecData;
		std::map<int, std::string> jmePtReso, jmeSF, jecUnc;

		// Btag Map
		std::map<int, std::map<char, float>> deepFlavourBTag, deepCSVBTag;
		std::map<char, std::string> deepCSVTagSFFile, deepFlavourTagSFFile;
		std::map<std::string, BTagCSVReader*> bTagReader;


		//Cut Variables
		int era;
		float ptCut, etaCut, deltaRCut;

		//Vector for the output variables
		std::vector<float> JetPt, JetEta, JetPhi, JetMass, JetRawPt, JetRawMass, JetRawFactor,
			JetCSVBTag, JetDFBTag,
			JetLooseDFBTagSF, JetMediumDFBTagSF, JetTightDFBTagSF, JetLooseCSVBTagSF, JetMediumCSVBTagSF, JetTightCSVBTagSF,
			JetLooseDFBTagSFUp, JetMediumDFBTagSFUp, JetTightDFBTagSFUp, JetLooseCSVBTagSFUp, JetMediumCSVBTagSFUp, JetTightCSVBTagSFUp,
			JetLooseDFBTagSFDown, JetMediumDFBTagSFDown, JetTightDFBTagSFDown, JetLooseCSVBTagSFDown, JetMediumCSVBTagSFDown, JetTightCSVBTagSFDown,
			FatJetDeepTagMD_H4qvsQCD, FatJetDeepTagMD_HbbvsQCD, FatJetDeepTagMD_TvsQCD, FatJetDeepTagMD_WvsQCD,
			FatJetDeepTagMD_ZHbbvsQCD, FatJetDeepTagMD_ZHccvsQCD, FatJetDeepTagMD_ZbbvsQCD, FatJetDeepTagMD_ZvsQCD,
			FatJetDeepTagMD_bbvsLight, FatJetDeepTagMD_ccvsLight, FatJetDeepTag_H, FatJetDeepTag_QCD,
			FatJetDeepTag_QCDothers, FatJetDeepTag_TvsQCD, FatJetDeepTag_WvsQCD, FatJetDeepTag_ZvsQCD;
		std::vector<bool> JetLooseDFBTag, JetMediumDFBTag, JetTightDFBTag, JetLooseCSVBTag, JetMediumCSVBTag, JetTightCSVBTag;
		float METPt, METPhi, JetRho, minMWjj, minMWjjPt, bestMWjj, bestMWjjPt, bestMTop, bestMTopPt;
		unsigned int nJet, nFatJet, nLooseDFBTagJet, nMediumDFBTagJet, nTightDFBTagJet, nLooseCSVBTagJet, nMediumCSVBTagJet, nTightCSVBTagJet;

		//TTreeReader Values
		std::unique_ptr<TTreeReaderValue<unsigned int>> jetNumber, fatJetNumber;
		std::unique_ptr<TTreeReaderArray<int>> jetGenIdx, jetFlavour;
		std::unique_ptr<TTreeReaderArray<float>> jetMass, jetPt, jetEta, jetPhi, jetArea, jetRawFactor, jetCSV, jetDF;

		std::unique_ptr<TTreeReaderArray<float>> genJetPt, genJetEta, genJetPhi, genJetMass;

		std::unique_ptr<TTreeReaderValue<float>> metPhi, metPt, jetRho;

		std::unique_ptr<TTreeReaderArray<float>> fatJetDeepTagMDH4qvsQCD, fatJetDeepTagMDHbbvsQCD, fatJetDeepTagMDTvsQCD, fatJetDeepTagMDWvsQCD,
			fatJetDeepTagMDZHbbvsQCD, fatJetDeepTagMDZHccvsQCD, fatJetDeepTagMDZbbvsQCD, fatJetDeepTagMDZvsQCD,
			fatJetDeepTagMDBbvsLight, fatJetDeepTagMDCcvsLight, fatJetDeepTagH, fatJetDeepTagQCD,
			fatJetDeepTagQCDothers, fatJetDeepTagTvsQCD, fatJetDeepTagWvsQCD, fatJetDeepTagZvsQCD;

		//Gen Level Information
		ROOT::Math::PtEtaPhiMVector genJet;
		//std::vector<int> alreadySeen; // is this needed?

		//Get jet energy correction
		FactorizedJetCorrector *jetCorrector;
		void SetCorrector(const char &runPeriod);

		float CorrectEnergy(const float &pt, const float &eta, const float &rho, const float &area);
		float SmearEnergy(const float &pt, const float &eta, const float &phi, const float &rho, const float &coneSize);

		template <typename T>
		void SortByIndex(T &var, std::vector<int> idx, unsigned int vectorSize);
	public:
		JetProducer(const int &era, const float &ptCut, const float &etaCut, const float &deltaRCut, const char &runPeriod, TTreeReader &reader);

		void BeginJob(TTree *tree, bool &isData, bool &doSystematics);
		void Produce(CutFlow &cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile *file);
};
