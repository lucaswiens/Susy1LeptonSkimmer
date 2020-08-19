#pragma once

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

#include <DataFormats/Candidate/interface/Candidate.h>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

//Struct for cutflow
struct CutFlow {
	TH1F* hist;
	Float_t weight = 1.;

	unsigned char nMinLepton = 1;
	unsigned char nMinElectron = 0;
	unsigned char nMinMuon = 0;
	unsigned char nMinJet = 0;
	unsigned char nMinFatjet = 0;

	bool passed = true;
};

class BaseProducer {
	protected:
		//File path for SF etc.
		std::string filePath = std::string(std::getenv("CMSSW_BASE")) + "/src/Susy1LeptonAnalsis/Susy1LeptonSkimmer/data/";

		std::map<int, std::map<char, std::string>> runEras = {
			{2016,	{
					{'B', "BCD"},
					{'C', "BCD"},
					{'D', "BCD"},
					{'E', "EF"},
					{'F', "GH"}
				}
			},
			{2017,	{
					{'B', "B"},
					{'C', "C"},
					{'D', "DE"},
					{'E', "DE"},
					{'F', "F"}
				}
			},
			{2018,	{
					{'A', "A"},
					{'B', "B"},
					{'C', "CD"},
					{'D', "CD"}
				}
			}
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
		virtual void BeginJob(TTree* tree, bool& isData) = 0;
		virtual void Produce(CutFlow& cutflow, Susy1LeptonProduct *product) = 0;
		virtual void EndJob(TFile* file) = 0;

		static float DeltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2);
};
