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
//#include <TTreeReader.h>
//#include <TTreeReaderValue.h>
//#include <TTreeReaderArray.h>

#include <DataFormats/Candidate/interface/Candidate.h>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/DataReader.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Utility.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

class BaseProducer {
	protected:
		//File path for SF etc.
		std::string filePath = std::string(std::getenv("CMSSW_BASE")) + "/src/Susy1LeptonAnalsis/Susy1LeptonSkimmer/data/";

		std::map<double, std::map<char, std::string>> runEras = { // TODO Remove and maybe put in config if it's still needed
			{2016.0, {
					{'B', "BCD"},
					{'C', "BCD"},
					{'D', "BCD"},
					{'E', "EF"},
					{'F', "EF"},
				 }
			},
			{2016.1, {
					{'F', "FGH"},
					{'G', "FGH"},
					{'H', "FGH"}
				}
			},
			{2017.0, {
					{'B', "B"},
					{'C', "C"},
					{'D', "D"},
					{'E', "E"},
					{'F', "F"}
				}
			},
			{2018.0, {
					{'A', "A"},
					{'B', "B"},
					{'C', "C"},
					{'D', "D"}
				}
			}
		};

		//Check for gen particle if it is last copy
		int FirstCopy(const int &index, const int &pdgID); //NANOAOD

		//Match Reco to gen particles
		std::tuple<int, int, int> SetGenParticles(const double &Pt, const double &Eta, const double &Phi, const int &i, const int &pdgID);

		static const int arrayMaxSize = 20;
	public:
		std::string Name;
		virtual ~BaseProducer(){};
		BaseProducer();
		//virtual void BeginJob(std::shared_ptr<TTree> tree, bool &isData, bool &doSystematics) = 0;
		virtual void Produce(DataReader &dataReader, Susy1LeptonProduct &product) = 0;
		virtual void EndJob(TFile &file) = 0;

};

#endif

