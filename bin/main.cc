//#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <memory>

#include <TFile.h>
#include <TTree.h>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/DataReader.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/CutFlow.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/Utility.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/MuonProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/ElectronProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/JetProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/DeltaPhiProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/PileUpWeightProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/METFilterProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/TriggerProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/GenLevelProducer.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

/*
	This is a skimmer for the 1 Lepton Delta Phi Analysis.
	The current design is designed to work on UltraLegacy samples:
	https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis
*/

int main(int argc, char *argv[]) {
	//Extract informations of command line
	if (argc < 7) {
		std::cerr << "It seems like you did not parse the arguments correctly..\nCheck the bin/main.cc to see which argument takes which position or use the python/createBatch.py that should handle everything automatically for you." << std::endl;
		std::exit(-1);
	}

	std::string inputFileName   = std::string(argv[1]);
	std::string outputFileName  = std::string(argv[2]);
	bool isData                 = std::string(argv[3]) == "True" ? true : false;
	bool doSystematics          = std::string(argv[4]) == "True" ? true : false;
	int era                     = std::stoi(argv[5]);
	char runPeriod              = (char)*argv[6]; //If it is MC, then runPeriod does not matter but still has to be given (e.g. just use "m")
	double xSection             = std::stod(argv[7]);

	int nMaxEvents;
	if(argc == 9) {
		nMaxEvents = std::stoi(std::string(argv[8]));
	} else {
		nMaxEvents = -999;
	}

	std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
	std::cout << "Input file for analysis: " << inputFileName << std::endl;

	std::string cmsswBase = std::getenv("CMSSW_BASE");
	pt::ptree scaleFactorTree, configTree;
	pt::read_json(std::string(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/config/config.json"), configTree);
	pt::read_json(std::string(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/config/scaleFactor.json"), scaleFactorTree);

	std::vector<CutFlow> cutflows;
	std::vector<std::shared_ptr<TTree>> outputTrees;
	std::vector<std::string> channels = {"Muon", "MuonIncl", "Electron", "ElectronIncl"};
	TFile outputFile(outputFileName.c_str(), "RECREATE");
	Susy1LeptonProduct product(era, isData, outputFileName, runPeriod, xSection, outputFile);

	std::cout << "\n" << outputFileName << "\npreVFP = " << product.GetIsPreVFP() << std::endl;

	for (const std::string &channel : channels) {
		outputFile.mkdir(channel.c_str());
		if (doSystematics) { // create a tree for each systematic variation
			for (std::string systematic : {"JEC", "JER"}) {
				for (std::string shift : {"Up", "Down"}) {
					outputFile.cd(channel.c_str());
					std::shared_ptr<TTree> tree = std::make_shared<TTree>((systematic + shift).c_str(), (systematic + shift).c_str());
					tree->SetDirectory(outputFile.GetDirectory(channel.c_str()));
					outputTrees.push_back(tree);

					cutflows.push_back(CutFlow(outputFile, channel, systematic, shift));
					std::string path = "Channel." + channel + ".Selection";
					for(const std::string part : Utility::GetKeys(configTree, path)) {
						cutflows.back().AddCut(part, product, configTree.get<std::string>(path + "." + part + ".operator"), configTree.get<short>(path + "." + part + ".threshold"));
					}
				}
			}
		} else { // create just the nominal tree
			outputFile.cd(channel.c_str());
			std::shared_ptr<TTree> tree = std::make_shared<TTree>("nominal", "nominal");
			tree->SetDirectory(outputFile.GetDirectory(channel.c_str()));
			outputTrees.push_back(tree);

			cutflows.push_back(CutFlow(outputFile, channel, "nominal", ""));
			std::string path = "Channel." + channel + ".Selection";
			for(const std::string part : Utility::GetKeys(configTree, path)) {
				std::cout << path << ": "<< part << configTree.get<std::string>(path + "." + part + ".operator")<<configTree.get<std::string>(path + "." + part + ".threshold") << std::endl;
				cutflows.back().AddCut(part, product, configTree.get<std::string>(path + "." + part + ".operator"), configTree.get<short>(path + "." + part + ".threshold"));
			}
		}
	}

	// Initialize Producers and register Product
	product.RegisterOutput(outputTrees);
	std::vector<std::shared_ptr<BaseProducer>> producers = {
		//std::shared_ptr<TriggerProducer>(new TriggerProducer(era, runPeriod)), // trigger names need update
		//std::shared_ptr<METFilterProducer>(new METFilterProducer(era)),
		//std::shared_ptr<LeptonProducer>(new LeptonProducer(era, vetoLeptonPtCut, etaCut, dxyCut, dzCut, sip3DCut, isoCut, preVFP)),
		std::shared_ptr<MuonProducer>(new MuonProducer(configTree, scaleFactorTree, product.GetEraSelector())),
		std::shared_ptr<ElectronProducer>(new ElectronProducer(configTree, scaleFactorTree, product.GetEraSelector())),
		//std::shared_ptr<JetProducer>(new JetProducer(era, jetPtCut, etaCut, deltaRCut, preVFP, runPeriod)), // exclude jet producer for now
		//std::shared_ptr<DeltaPhiProducer>(new DeltaPhiProducer()),
	};
	if (!isData) {
		//producers.push_back(std::shared_ptr<PileUpWeightProducer>(new PileUpWeightProducer(era))),
		//producers.push_back(std::shared_ptr<GenLevelProducer>(new GenLevelProducer(era)));
	}

	DataReader dataReader(inputFileName, "Events");
	std::size_t nEvents = dataReader.GetEntries();

	//ProgressBar(0, 0);
	std::cout << std::endl << "Starting Event loop with " << nEvents << " Events" << std::endl;
	for (std::size_t entry = 0; entry < dataReader.GetEntries(); ++entry) {
		// Stop Events Loop after nMaxEvents
		if (nMaxEvents > 0 && entry >= nMaxEvents) { break;}
		if (entry % 50000 == 0) {
			std::cout  << "Processed " << 100*(double)(entry + 1)/(nMaxEvents < 0? nEvents : nMaxEvents) << "% Events at a rate of " + std::to_string(nEvents / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count()) + " Hz." << std::endl;
		}

		dataReader.SetEntry(entry);
		for (std::shared_ptr<BaseProducer> producer : producers) {
			producer->Produce(dataReader, product);
		}

		for(std::size_t iTree = 0; iTree < outputTrees.size(); iTree++) {
			cutflows.at(iTree).Count();
			cutflows.at(iTree).FillCutflow();
			if(cutflows.at(iTree).Passed()) { outputTrees.at(iTree)->Fill();}
		}
	}

	std::size_t finalNumberOfEvents;
	std::cout  << "Processed Events at a rate of " + std::to_string(nEvents / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count()) + " Hz." << std::endl;
	for (std::shared_ptr<TTree> tree : outputTrees) {
		finalNumberOfEvents = tree->GetEntries();
		if (nMaxEvents > 0) {
			std::cout << std::setw(20) << tree->GetDirectory()->GetName() << " " << std::setw(10) << tree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nMaxEvents << " (" << 100*(double)finalNumberOfEvents/nMaxEvents << "%)" << std::endl;
		} else {
			std::cout << std::setw(20) << tree->GetDirectory()->GetName() << " " << std::setw(10) << tree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nEvents << " (" << 100*(double)finalNumberOfEvents/nEvents << "%)" << std::endl;
		}
	}

	for (std::size_t iTree = 0; iTree < outputTrees.size(); iTree++) {
		outputFile.cd(outputTrees.at(iTree)->GetDirectory()->GetName());
		outputTrees.at(iTree)->Write(0, TObject::kOverwrite);
		cutflows.at(iTree).WriteOutput();
	}

	outputFile.Write(0, TObject::kOverwrite);
	outputFile.Close();

	std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
	std::cout << "Finished event loop (in minutes): " << std::chrono::duration_cast<std::chrono::minutes>(end - start).count() << std::endl;
	std::cout << "Output file created: " + outputFileName << std::endl;

	return 0;
}
