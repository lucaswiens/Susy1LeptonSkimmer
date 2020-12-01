#pragma once

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include "PhysicsTools/Heppy/interface/Davismt2.h"

#include <TMath.h>

class DeltaPhiProducer : public BaseProducer {
	private:
		//Check if it is data or MC
		bool isData, doSystematics;

		//Vector for the output variables
		float HT, LT, LP, deltaPhi, dPhi, wBosonMt;
		int signalRegionCSV, signalRegionDF;
		unsigned int nIsoTrack;
		bool IsoTrackVeto;
		std::vector<float> IsoTrackMt2, IsoTrackPt;
		std::vector<bool> IsoTrackHadronicDecay;
		std::vector<int> IsoTrackPdgId;

		std::unique_ptr<TTreeReaderValue<unsigned int>> isoTrackNumber;
		std::unique_ptr<TTreeReaderArray<float>> isoTrackPt, isoTrackEta, isoTrackPhi, isoTrackMass, isoTrackIso, isoTrackDxy, isoTrackDz, isoTrackSip3d, isoTrackMiniPFRelIsoAll;
		std::unique_ptr<TTreeReaderArray<int>> isoTrackPdgId, isoTrackCharge;

		float DeltaPhi(ROOT::Math::PtEtaPhiMVector v1, ROOT::Math::PtEtaPhiMVector v2);
	public:
		DeltaPhiProducer(TTreeReader &reader);//TODO reader needeed?

		void BeginJob(TTree *tree, bool &isData, bool &doSystematics);
		void Produce(CutFlow &cutflow, Susy1LeptonProduct *product);
		void EndJob(TFile *file);
};
