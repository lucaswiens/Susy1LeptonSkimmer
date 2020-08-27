#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

class TriggerProducer : public BaseProducer {
	private:
		//Check if it is data or MC
		bool isData;

		//Cut Variables
		int era;

		//Vector for the output variables
		int nIsr, nISRweight, nISRttweightsystUp, nISRttweightsystDown;
		bool HLT_EleOr, HLT_MuOr, HLT_LepOr, HLT_MetOr, TriggerEfficiency;

		//TTreeReader Values for NANO AOD analysis
		std::unique_ptr<TTreeReaderValue<bool>> hlt_Ele105_CaloIdVT_GsfTrkIdT, hlt_Ele115_CaloIdVT_GsfTrkIdT, hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, hlt_Ele35_WPTight_Gsf, hlt_Ele27_WPTight_Gsf, hlt_Ele15_IsoVVVL_PFHT350, hlt_Ele15_IsoVVVL_PFHT400, hlt_Ele15_IsoVVVL_PFHT450, hlt_Mu15_IsoVVVL_PFHT450, hlt_Mu50, hlt_IsoMu24, hlt_IsoTkMu24, hlt_Mu15_IsoVVVL_PFHT350, hlt_Mu15_IsoVVVL_PFHT400, hlt_PFMET100_PFMHT100_IDTight, hlt_PFMETNoMu100_PFMHTNoMu100_IDTight, hlt_PFMET110_PFMHT110_IDTight, hlt_PFMETNoMu110_PFMHTNoMu110_IDTight, hlt_PFMET120_PFMHT120_IDTight, hlt_PFMETNoMu120_PFMHTNoMu120_IDTight;

	public:
		TriggerProducer(const int& era, TTreeReader& reader);

		void BeginJob(TTree* tree, bool& isData);
		void Produce(CutFlow& cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile* file);
};
