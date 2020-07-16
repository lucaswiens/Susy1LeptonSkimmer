#ifndef BASEPRODUCER_H
#define BASEPRODUCER_H

#include <memory>
#include <map>
#include <vector>
#include <cmath>
#include <bitset>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <Rtypes.h>
#include <Math/Vector4Dfwd.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <FWCore/Framework/interface/Event.h>

#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

//Struct for cutflow
struct CutFlow {
	TH1F* hist;
	Float_t weight = 1.;

	unsigned char nMinEle=0;
	unsigned char nMinMu=0;
	unsigned char nMinJet=0;
	unsigned char nMinFatjet=0;

	bool passed = true;
};

class BaseProducer {
	protected:
		//File path for SF etc.
		std::string filePath = std::string(std::getenv("CMSSW_BASE")) + "/src/Susy1LeptonAnalsis/Susy1LeptonSkimmer/data/";

		std::map<int, std::map<std::string, std::pair<int, int>>> runEras = {
			  {2017, {
						{"B", {297046, 299329}},
						{"C", {299368, 302029}},
						{"DE", {302030, 304797}},
						{"F", {305040, 306462}},
					 }
			  },
		};

		//Collection which are used in several producers if NANO AOD is produced
		TTreeReader* reader = NULL;

		std::unique_ptr<TTreeReaderValue<unsigned int>> run;

		std::unique_ptr<TTreeReaderArray<float>> trigObjPt, trigObjPhi, trigObjEta;
		std::unique_ptr<TTreeReaderArray<int>> trigObjID, trigObjFilterBit;

		std::unique_ptr<TTreeReaderArray<float>> genPhi, genEta, genPt, genMass;
		std::unique_ptr<TTreeReaderArray<int>> genID, genMotherIdx, genStatus, eleGenIdx, muonGenIdx;

		//Set trihObj and Gen particle collection
		void SetCollection(bool &isData);
		bool isSyst=false;

		//Check for gen particle if it is last copy
		int FirstCopy(const int& index, const int& pdgID); //NANOAOD

		//Match Reco to gen particles
		std::tuple<int, int, int> SetGenParticles(const float& Pt, const float& Eta, const float& Phi, const int &i, const int& pdgID);

	public:
		virtual ~BaseProducer(){};
		BaseProducer();
		BaseProducer(TTreeReader* reader);
		virtual void BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst=false) = 0;
		virtual void Produce(std::vector<CutFlow>& cutflows) = 0;
		virtual void EndJob(TFile* file) = 0;

		static float DeltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2);
};
#endif
