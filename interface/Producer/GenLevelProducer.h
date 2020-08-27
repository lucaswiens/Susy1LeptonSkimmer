#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <TMath.h>

class GenLevelProducer : public BaseProducer {
	private:
		//Check if it is data or MC
		bool isData;
		int era;

		//Variables to be stored in the output tree
		std::vector<float> GenDeltaPhiLepWSum, GenDeltaPhiLepWDirect, GenWSumMass, GenWDirectMass, GenMTLepNu;
		std::vector<int> GrandMotherId, GenTauGrandMotherId, GenTauMotherId, GenLepGrandMotherId, GenLepMotherId;
		int NGenPart, NGenLepton, NGenTau, NGenLeptonFromTau, NGenMatchedW, LeptonDecayChannelFlag, nISR;
		bool isDiLeptonEvent, isHadTauEvent, LeptonsInAcceptance;
		float nISRWeight, nISRWeightUp, nISRWeightDown, nISRWeight_Mar17, nISRWeightUp_Mar17, nISRWeightDown_Mar17, nISRWeight_ICHEP16, nISRWeightUp_ICHEP16, nISRWeightDown_ICHEP16;

		std::map<int, float> ISRweights_Mar17, ISRweights_ICHEP16, ISRweightssyst_Mar17, ISRweightssyst_ICHEP16;

		//TTreeReader Values for NANO AOD analysis
		std::unique_ptr<TTreeReaderValue<unsigned int>> nGenPart;
		std::unique_ptr<TTreeReaderArray<float>> genPartPt, genPartEta, genPartMass, genPartPhi;
		std::unique_ptr<TTreeReaderArray<int>> genPartIdxMother, genPartPdgId, genPartStatus, genPartStatusFlag;

		float DeltaPhi(ROOT::Math::PtEtaPhiMVector v1, ROOT::Math::PtEtaPhiMVector v2);
	public:
		GenLevelProducer(const int& era, TTreeReader& reader);

		void BeginJob(TTree* tree, bool& isData);
		void Produce(CutFlow& cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile* file);
};
