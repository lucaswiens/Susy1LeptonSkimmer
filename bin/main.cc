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
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/ScaleFactorProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/DeltaPhiProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/PileUpWeightProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/METFilterProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/TriggerProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/GenLevelProducer.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

/*#####################################################################
#   This is a skimmer for the 1 Lepton Delta Phi Analysis.            #
#   The current design is designed to work on UltraLegacy samples:    #
#   https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis   #
#####################################################################*/

int main(int argc, char *argv[]) {
	//Extract informations of command line
	if (argc < 7) {
		std::cerr << "It seems like you did not parse the arguments correctly..\nCheck the bin/main.cc to see which argument takes which position or use the python/createBatch.py that should handle everything automatically for you." << std::endl;
		std::exit(-1);
	}

	std::string inputFileName  = std::string(argv[1]);
	std::string outputFileName = std::string(argv[2]);
	bool isData                = std::string(argv[3]) == "True" ? true : false;
	bool doSystematics         = std::string(argv[4]) == "True" ? true : false; // Maybe use a string instead that can be "Nominal", "Up" or "Down"
	int era                    = std::stoi(argv[5]);
	char runPeriod             = (char)*argv[6]; //If it is MC, then runPeriod does not matter but still has to be given (e.g. just use "m")
	double xSection            = std::stod(argv[7]);

	int nMaxEvents;
	if (argc == 9) {
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
	std::vector<std::string> triggerNames, metTriggerNames, metFilterNames;
	TFile outputFile(outputFileName.c_str(), "RECREATE");
	Susy1LeptonProduct product(era, isData, outputFileName, runPeriod, xSection, configTree, outputFile);
	for (const std::string &channel : channels) {
		outputFile.mkdir(channel.c_str());
		if (doSystematics) { // create a tree for each systematic variation #FIXME
			for (std::string systematic : {"JEC", "JER"}) {
				for (std::string shift : {"Up", "Down"}) {
					outputFile.cd(channel.c_str());
					std::shared_ptr<TTree> tree = std::make_shared<TTree>((systematic + shift).c_str(), (systematic + shift).c_str());
					tree->SetDirectory(outputFile.GetDirectory(channel.c_str()));
					tree->SetAutoFlush(10000);
					outputTrees.push_back(tree);

					cutflows.push_back(CutFlow(outputFile, channel, systematic, shift));
					std::string path = "Channel." + channel + ".Selection";
					for (const std::string part : Utility::GetKeys(configTree, path)) {
						cutflows.back().AddCut(part, product, configTree.get<std::string>(path + "." + part + ".Operator"), configTree.get<short>(path + "." + part + ".Threshold"));
					}
				}
			}
		} else { // create just the Nominal tree
			outputFile.cd(channel.c_str());
			std::shared_ptr<TTree> tree = std::make_shared<TTree>("Nominal", "Nominal");
			tree->SetDirectory(outputFile.GetDirectory(channel.c_str()));
			outputTrees.push_back(tree);

			for (const std::string &name : Utility::GetVector<std::string>(configTree, "Channel." + channel + ".Trigger." + product.GetEraSelector())) {
				if (std::find(triggerNames.begin(), triggerNames.end(), name) == triggerNames.end()) {
					triggerNames.push_back(name);
				}
			}

			for (const std::string &name : Utility::GetVector<std::string>(configTree, "Channel." + channel + ".METTrigger")) {
				if (std::find(metTriggerNames.begin(), metTriggerNames.end(), name) == metTriggerNames.end()) {
					metTriggerNames.push_back(name);
				}
			}

			for (const std::string &name : Utility::GetVector<std::string>(configTree, "METFilter." + product.GetEraSelector())) {
				if (std::find(metFilterNames.begin(), metFilterNames.end(), name) == metFilterNames.end()) {
					metFilterNames.push_back(name);
				}
			}

			cutflows.push_back(CutFlow(outputFile, channel, "Nominal", ""));
			std::string path = "Channel." + channel + ".Selection";
			for (const std::string part : Utility::GetKeys(configTree, path)) {
				cutflows.back().AddCut(part, product, configTree.get<std::string>(path + "." + part + ".Operator"), configTree.get<short>(path + "." + part + ".Threshold"));
			}
		}
	}

	for (int iChannel = 0; iChannel < channels.size(); iChannel++) {
		std::vector<int> triggerIndex;
		for(int iTrigger = 0; iTrigger < triggerNames.size(); iTrigger++){
			for(const std::string& triggerName : Utility::GetVector<std::string>(configTree, "Channel." + channels[iChannel] + ".Trigger." + product.GetEraSelector())) {
				if(triggerName == triggerNames[iTrigger]) triggerIndex.push_back(iTrigger);
			}
		}

		//Input class instead output class is used to register cut for trigger!
		cutflows[iChannel].AddTriggerOr(triggerIndex, product);
	}

	// Initialize Producers and register Product
	product.RegisterOutput(outputTrees, configTree);
	product.RegisterTrigger(triggerNames, metTriggerNames, outputTrees);
	product.RegisterMetFilter(metFilterNames, outputTrees);
	std::vector<std::shared_ptr<BaseProducer>> producers = {
		std::shared_ptr<TriggerProducer>(new TriggerProducer(configTree, scaleFactorTree, product)),
		std::shared_ptr<METFilterProducer>(new METFilterProducer(configTree, scaleFactorTree, product)),
		std::shared_ptr<MuonProducer>(new MuonProducer(configTree, scaleFactorTree, product.GetEraSelector())),
		std::shared_ptr<ElectronProducer>(new ElectronProducer(configTree, scaleFactorTree)),
		std::shared_ptr<JetProducer>(new JetProducer(configTree, scaleFactorTree, product)),
		std::shared_ptr<DeltaPhiProducer>(new DeltaPhiProducer(configTree, scaleFactorTree)),
	};
	if (!isData) {
		producers.push_back(std::shared_ptr<PileUpWeightProducer>(new PileUpWeightProducer(configTree, scaleFactorTree, product.GetEraSelector())));
		producers.push_back(std::shared_ptr<GenLevelProducer>(new GenLevelProducer(configTree, scaleFactorTree, product.GetEraSelector())));
		producers.push_back(std::shared_ptr<ScaleFactorProducer>(new ScaleFactorProducer(configTree, scaleFactorTree, product.GetEraSelector())));
	}

	DataReader dataReader(inputFileName, "Events", isData);
	int nEvents = dataReader.GetEntries();

	//Register trigger to data reader class
	dataReader.RegisterTrigger(triggerNames, metTriggerNames);
	dataReader.RegisterMetFilter(Utility::GetVector<std::string>(configTree, "METFilter." + product.GetEraSelector()));

	//ProgressBar(0, 0);
	if (nMaxEvents > 0) {
		std::cout << std::endl << "Starting Event loop with " << nMaxEvents << " Events" << std::endl;
	} else {
		std::cout << std::endl << "Starting Event loop with " << nEvents << " Events" << std::endl;
	}
	for (int entry = 0; entry < dataReader.GetEntries(); ++entry) {
		// Stop Events Loop after nMaxEvents
		if (nMaxEvents > 0 && entry >= nMaxEvents) { break;}
		if (entry % 20000 == 0) {
			std::cout  << "Processed " << int(100*(double)(entry + 1)/(nMaxEvents > 0? nMaxEvents : nEvents)) << "% Events at a rate of " + std::to_string(entry / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count()) + " Hz." << std::endl;
		}

		dataReader.SetEntry(entry);
		for (std::shared_ptr<BaseProducer> producer : producers) {
			producer->Produce(dataReader, product);
		}

		for (int iTree = 0; iTree < outputTrees.size(); iTree++) {
			cutflows.at(iTree).Count();
			cutflows.at(iTree).FillCutflow();
			if (cutflows.at(iTree).Passed()) { outputTrees.at(iTree)->Fill();}
		}
	}

	int finalNumberOfEvents;
	if (nMaxEvents > 0) {
		std::cout  << "Processed Events at a rate of " + std::to_string(nMaxEvents / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count()) + " Hz." << std::endl;
	} else {
		std::cout  << "Processed Events at a rate of " + std::to_string(nEvents / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count()) + " Hz." << std::endl;
	}
	for (std::shared_ptr<TTree> tree : outputTrees) {
		finalNumberOfEvents = tree->GetEntries();
		if (nMaxEvents > 0) {
			std::cout << std::setw(20) << tree->GetDirectory()->GetName() << tree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nMaxEvents << " (" << 100*(double)finalNumberOfEvents/nMaxEvents << "%)" << std::endl;
		} else {
			std::cout << std::setw(20) << tree->GetDirectory()->GetName() << tree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nEvents << " (" << 100*(double)finalNumberOfEvents/nEvents << "%)" << std::endl;
		}
	}

	for (int iTree = 0; iTree < outputTrees.size(); iTree++) {
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
