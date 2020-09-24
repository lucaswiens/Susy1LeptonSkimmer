#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

class METFilterProducer : public BaseProducer {
	private:
		//Check if it is data or MC
		bool isData, doSystematics;

		//Cut Variables
		int era;
		float ptCut, etaCut, dxyCut, dzCut, sip3dCut, isoCut;

		//Vector for the output variables
		unsigned int RunNumber, LuminosityNumber;
		unsigned long long EventNumber;
		bool PassFilters, PassFiltersMoriond2017Tight, PassCSCFilterList;

		//TTreeReader Values for NANO AOD analysis
		std::unique_ptr<TTreeReaderValue<unsigned int>> runNumber, luminosityNumber;
		std::unique_ptr<TTreeReaderValue<unsigned long long>> eventNumber;
		std::unique_ptr<TTreeReaderValue<bool>> flagEeBadScFilter, flagBadMuons, flagHBHENoiseFilter, flagHBHENoiseIsoFilter, flagEcalDeadCellTriggerPrimitiveFilter, flagGoodVertices, flagGlobalSuperTightHalo2016Filter;

	public:
		METFilterProducer(const int &era, TTreeReader &reader);

		void BeginJob(TTree *tree, bool &isData, bool &doSystematics);
		void Produce(CutFlow &cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile *file);
};
