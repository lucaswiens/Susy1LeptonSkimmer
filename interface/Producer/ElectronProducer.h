#ifndef ELECTRONPRODUCER_H
#define ELECTRONPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>

#include <cmath>
#include <random>


class ElectronProducer : public BaseProducer {
	private:
		//Cut Variables
		int era;
		float electronGoodPtCut, electronVetoPtCut, electronEtaCut, electronGoodIsoCut, electronVetoIsoCut, electronAntiIsoCut;
		int electronGoodNumberOfLostHitsCut;
		char electronGoodCutBasedIdCut, electronVetoCutBasedIdCut, electronAntiIsCutBasedIdCut, electronAntiIsNotCutBasedIdCut;
	public:
		ElectronProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

