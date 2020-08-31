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
		bool isData;
		char runPeriod;
		//bool isSignal // or maybe make a fastsim flag

		//Classes for reading jet energy SF
		JME::JetParameters jetParameter;
		JME::JetResolution resolution;
		JME::JetResolutionScaleFactor resolution_sf;

		//JEC/JER systematics
		JetCorrectionUncertainty* jetCorrectionUncertainty;
		std::string jecSyst;
		bool isJERsyst = false;
		//std::int<systematic, float> smearFactor;

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
		std::vector<float> JetPt, JetEta, JetPhi, JetMass, JetRawPt, JetRawMass, JetRawFactor, JetCSVBTag, JetDFBTag, JetPt_jerUp, JetMass_jerUp, JetPt_jerDown, JetMass_jerDown, JetPt_jecUp, JetMass_jecUp, JetPt_jecDown, JetMass_jecDown, JetLooseDFBTagSF, JetMediumDFBTagSF, JetTightDFBTagSF, JetLooseCSVBTagSF, JetMediumCSVBTagSF, JetTightCSVBTagSF, JetLooseDFBTagSFUp, JetMediumDFBTagSFUp, JetTightDFBTagSFUp, JetLooseCSVBTagSFUp, JetMediumCSVBTagSFUp, JetTightCSVBTagSFUp, JetLooseDFBTagSFDown, JetMediumDFBTagSFDown, JetTightDFBTagSFDown, JetLooseCSVBTagSFDown, JetMediumCSVBTagSFDown, JetTightCSVBTagSFDown, FatJet_deepTagMD_H4qvsQCD, FatJet_deepTagMD_HbbvsQCD, FatJet_deepTagMD_TvsQCD, FatJet_deepTagMD_WvsQCD, FatJet_deepTagMD_ZHbbvsQCD, FatJet_deepTagMD_ZHccvsQCD, FatJet_deepTagMD_ZbbvsQCD, FatJet_deepTagMD_ZvsQCD, FatJet_deepTagMD_bbvsLight, FatJet_deepTagMD_ccvsLight, FatJet_deepTag_H, FatJet_deepTag_QCD, FatJet_deepTag_QCDothers, FatJet_deepTag_TvsQCD, FatJet_deepTag_WvsQCD, FatJet_deepTag_ZvsQCD;
		std::vector<bool> JetLooseDFBTag, JetMediumDFBTag, JetTightDFBTag, JetLooseCSVBTag, JetMediumCSVBTag, JetTightCSVBTag;
		float METPt, METPhi, JetRho, METPt_jerUp, METPhi_jerUp, METPt_jerDown, METPhi_jerDown, METPt_jecUp, METPhi_jecUp, METPt_jecDown, METPhi_jecDown, minMWjj, minMWjjPt, bestMWjj, bestMWjjPt, bestMTop, bestMTopPt;
		unsigned int nJet, nFatJet, nLooseDFBTagJet, nMediumDFBTagJet, nTightDFBTagJet, nLooseCSVBTagJet, nMediumCSVBTagJet, nTightCSVBTagJet;

		//TTreeReader Values
		std::unique_ptr<TTreeReaderValue<unsigned int>> jetNumber, fatJetNumber;
		std::unique_ptr<TTreeReaderArray<int>> jetGenIdx, jetFlavour;
		std::unique_ptr<TTreeReaderArray<float>> jetMass, jetPt, jetEta, jetPhi, jetArea, jetRawFactor, jetCSV, jetDF;

		std::unique_ptr<TTreeReaderArray<float>> genJetPt, genJetEta, genJetPhi, genJetMass;

		std::unique_ptr<TTreeReaderValue<float>> metPhi, metPt, jetRho;

		std::unique_ptr<TTreeReaderArray<float>> fatJetDeepTagMDH4qvsQCD, fatJetDeepTagMDHbbvsQCD, fatJetDeepTagMDTvsQCD, fatJetDeepTagMDWvsQCD, fatJetDeepTagMDZHbbvsQCD, fatJetDeepTagMDZHccvsQCD, fatJetDeepTagMDZbbvsQCD, fatJetDeepTagMDZvsQCD, fatJetDeepTagMDBbvsLight, fatJetDeepTagMDCcvsLight, fatJetDeepTagH, fatJetDeepTagQCD, fatJetDeepTagQCDothers, fatJetDeepTagTvsQCD, fatJetDeepTagWvsQCD, fatJetDeepTagZvsQCD;

		//Gen Level Information
		ROOT::Math::PtEtaPhiMVector genJet;
		std::vector<int> alreadySeen;

		//Get jet energy correction
		FactorizedJetCorrector* jetCorrector;
		void SetCorrector(const char& runPeriod);

		float CorrectEnergy(const float& pt, const float& eta, const float& rho, const float& area);
		std::map<char, float> SmearEnergy(const float& pt, const float& eta, const float& phi, const float& rho, const float& coneSize);

		template <typename T>
		void SortByIndex(T& var, std::vector<int> idx, unsigned int vectorSize);
	public:
		JetProducer(const int& era, const float& ptCut, const float& etaCut, const float& deltaRCut, const char& runPeriod, TTreeReader& reader);

		void BeginJob(TTree* tree, bool& isData);
		void Produce(CutFlow& cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile* file);
};
