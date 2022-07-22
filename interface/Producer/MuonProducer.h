#ifndef MUONPRODUCER_H
#define MUONPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
//#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

#include <RoccoR/RoccoR.cc>

#include <Math/LorentzVector.h> // FIXME maybe not needed
#include <Math/PtEtaPhiM4D.h>

#include <cmath>
#include <random>


class MuonProducer : public BaseProducer {
	private:
		//Cut Variables
		double muonGoodPtCut, muonVetoPtCut, muonEtaCut, muonGoodIsoCut, muonVetoIsoCut, muonAntiIsoCut, muonDxyCut, muonDzCut, muonSip3dCut;
		char muonGoodCutBasedIdCut, muonVetoCutBasedIdCut, muonAntiCutBasedIdCut;

		//Muon scale corrector
		std::string eraSelector;
		RoccoR rc;

	public:
		std::string Name = "MuonProducer";
		MuonProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

