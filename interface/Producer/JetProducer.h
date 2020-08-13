#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <JetMETCorrections/Modules/interface/JetCorrectionProducer.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>

class JetProducer : public BaseProducer {

	enum JetType {SUBAK4, AK4, AK8, PF, VTX};

	private:
		//Check if it is data or MC
		bool isData;
		//bool isSignal // or maybe make a fastsim flag

		//Classes for reading jet energy SF
		JME::JetParameters jetParameter;
		std::map<JetType, JME::JetResolution> resolution;
		std::map<JetType, JME::JetResolutionScaleFactor> resolution_sf;

		//JEC/JER systematics
		std::map<JetType, JetCorrectionUncertainty*> jetCorrectionUncertainty;
		std::string jecSyst;
		bool isJERsyst = false;
		//std::int<systematic, float> smearFactor;

		// JEC map
		std::map<int, std::vector<std::string>> jecMC, jecFastSim, jerMC, jecData;
		std::map<int, std::string> jmePtReso, jmeSF, jecUnc;

		// Btag Map
		std::map<int, std::map<char, float>> deepFlavourBTag, deepCSVBTag;

		//Cut Variables
		int era;
		float ptCut, etaCut, deltaRCut;

		//Vector for the output variables
		std::vector<float> JetPt, JetEta, JetPhi, JetMass, JetRawPt, JetRawMass, JetRawFactor, JetCSVBTag, JetDFBTag, JetPt_jerUp, JetMass_jerUp, JetPt_jerDown, JetMass_jerDown, JetPt_jecUp, JetMass_jecUp, JetPt_jecDown, JetMass_jecDown;
		std::vector<bool> JetLooseDFBTag, JetMediumDFBTag, JetTightDFBTag, JetLooseCSVBTag, JetMediumCSVBTag, JetTightCSVBTag;
		float METPt, METPhi, JetRho, METPt_jerUp, METPhi_jerUp, METPt_jerDown, METPhi_jerDown, METPt_jecUp, METPhi_jecUp, METPt_jecDown, METPhi_jecDown;
		int runNumber;
		unsigned int nJet, nFatJet, nLooseDFBTagJet, nMediumDFBTagJet, nTightDFBTagJet, nLooseCSVBTagJet, nMediumCSVBTagJet, nTightCSVBTagJet;

		//TTreeReader Values
		std::unique_ptr<TTreeReaderValue<unsigned int>> jetNumber, fatJetNumber;

		std::unique_ptr<TTreeReaderArray<float>> fatJetPt, fatJetEta, fatJetPhi, fatJetMass, fatJetArea, fatJetCSV, fatJetDF;
		std::unique_ptr<TTreeReaderArray<float>> fatJetTau1, fatJetTau2, fatJetTau3;

		std::unique_ptr<TTreeReaderArray<int>> jetGenIdx, jetFlavour;
		std::unique_ptr<TTreeReaderArray<float>> jetMass, jetPt, jetEta, jetPhi, jetArea, jetRawFactor, jetCSV, jetDF;

		std::unique_ptr<TTreeReaderArray<float>> genJetPt, genJetEta, genJetPhi, genJetMass;
		std::unique_ptr<TTreeReaderArray<float>> genFatJetPt, genFatJetEta, genFatJetPhi, genFatJetMass;

		std::unique_ptr<TTreeReaderValue<float>> metPhi, metPt, jetRho;

		//Gen Level Information
		std::map<JetType, ROOT::Math::PtEtaPhiMVector> genJet;
		std::vector<int> alreadySeen;

		//Get jet energy correction
		std::map<JetType, FactorizedJetCorrector*> jetCorrector;
		void SetCorrector(const JetType &type, const int& runNumber);

		float CorrectEnergy(const float& pt, const float& eta, const float& rho, const float& area, const JetType &type);
		std::map<char, float> SmearEnergy(const float& pt, const float& eta, const float& phi, const float& rho, const float& coneSize, const JetType& type);

		template <typename T>
		void SortByIndex(T& var, std::vector<int> idx);
	public:
		JetProducer(const int& era, const float& ptCut, const float& etaCut, const float& deltaRCut, TTreeReader& reader);

		void BeginJob(TTree* tree, bool& isData);
		void Produce(CutFlow& cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile* file);
};
