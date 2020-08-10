#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <JetMETCorrections/Modules/interface/JetCorrectionProducer.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>

//#include <DataFormats/PatCandidates/interface/Jet.h>
//#include <DataFormats/PatCandidates/interface/MET.h>
//#include <DataFormats/JetReco/interface/GenJet.h>

//typedef LorentzVector<PtEtaPhiM4D<double> > ROOT::Math::PtEtaPhiMVector

class JetProducer: public BaseProducer {

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

		// JEC dict
		//std::map<int, std::vector<std::string>> jecMC, jecFastSim, jerMC, jecData;
		//std::map<int, std::string> jmePtReso, jmeSF;
		std::map<int, std::vector<std::string>> jecMC, jecFastSim, jerMC, jecData;
		std::map<int, std::string> jmePtReso, jmeSF, jecUnc;

		//Cut Variables
		int era;
		float ptCut, etaCut, deltaRCut;

		//Vector for the output variables
		std::vector<float> JetPt, JetEta, JetPhi, JetMass, JetPtUp, JetEtaUp, JetPhiUp, JetMassUp, JetPtDown, JetEtaDown, JetPhiDown, JetMassDown;
		std::map<JetType, std::vector<int>> TrueFlavour, Charge, FatJetIdx, partID, mothID, grandID; //TODO Might be not needed
		float METPt, METEta, METPhi, METMass, JetRho;
		int runNumber;
		unsigned int nJet, nFatJet;

		//TTreeReader Values
		std::unique_ptr<TTreeReaderValue<unsigned int>> jetNumber, fatJetNumber;
		//std::unique_ptr<TTreeReaderArray<float>> jetPt, jetEta, jetPhi, jetMass;
		//std::unique_ptr<TTreeReaderValue<float>> metPt, metPhi;// jetRho;

		std::unique_ptr<TTreeReaderArray<float>> fatJetPt, fatJetEta, fatJetPhi, fatJetMass, fatJetArea, fatJetCSV;
		std::unique_ptr<TTreeReaderArray<float>> fatJetTau1, fatJetTau2, fatJetTau3;

		std::unique_ptr<TTreeReaderArray<int>> jetGenIdx, jetFlavour;
		std::unique_ptr<TTreeReaderArray<float>> jetMass, jetPt, jetEta, jetPhi, jetArea, jetDeepBValue;

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
		void SetGenParticles(const int& i, const float& pt, const float& eta, const float& phi, const std::vector<int>& pdgID, const JetType &type);

		template <typename T>
		void SortByIndex(T& var, std::vector<int> idx);
	public:
		JetProducer(const int& era, const float& ptCut, const float& etaCut, const float& deltaRCut, TTreeReader& reader);

		void BeginJob(TTree* tree, bool& isData);
		void Produce(CutFlow& cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile* file);
};
