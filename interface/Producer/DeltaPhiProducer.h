#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <TMath.h>

class DeltaPhiProducer: public BaseProducer {
	private:
		//Check if it is data or MC
		bool isData;

		//Vector for the output variables
		float HT, LT, LP, deltaPhi, wBosonMt;
		ROOT::Math::PtEtaPhiMVector leptonP4, metP4;

		float DeltaPhi(ROOT::Math::PtEtaPhiMVector v1, ROOT::Math::PtEtaPhiMVector v2);

	public:
		DeltaPhiProducer(TTreeReader& reader);//TODO reader needeed?

		void BeginJob(TTree* tree, bool& isData);
		void Produce(CutFlow& cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile* file);
};
