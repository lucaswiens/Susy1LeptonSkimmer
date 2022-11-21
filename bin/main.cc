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
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/FastSimProducer.h>
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
#   The current setup is designed to work on UltraLegacy samples:     #
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
	const int &era             = std::stoi(argv[3]);
	const char &runPeriod      = (char)*argv[4]; // Data : A-H; MC : M; FastSim : S
	const float &xSection      = std::stod(argv[5]);
	const bool &isData         = (runPeriod == 'M' || runPeriod == 'S') ? false : true;
	const bool &isFastSim      = runPeriod == 'S' ? true : false;

	std::cout << "Producing NTuples for a " << (isFastSim ? "fastsim " : "") << (isData ? "data" : "MC") << " sample." << std::endl <<
		"Year = " << era << std::endl <<
		"RunPeriod = " << runPeriod << std::endl <<
		"Cross Section  = " << xSection << std::endl <<
	std::endl;


	int nMaxEvents;
	if (argc >= 7) {
		nMaxEvents = std::stoi(std::string(argv[6]));
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
	std::vector<std::string> channels = {"Muon", "Electron", "LeptonIncl"};
	std::vector<std::string> triggerNames, metTriggerNames, metFilterNames;
	TFile outputFile(outputFileName.c_str(), "RECREATE");
	Susy1LeptonProduct product(era, isData, isFastSim, outputFileName, runPeriod, xSection, configTree, outputFile);

	// Create a TTree and CutFlow for each channel
	static const long autoFlush = - 50 * 1024 * 1024;
	for (const std::string &channel : channels) {
		std::shared_ptr<TTree> tree = std::make_shared<TTree>(channel.c_str(), channel.c_str());
		tree->SetAutoFlush(autoFlush);
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

		cutflows.push_back(CutFlow(outputFile, channel));
		std::string path = "Channel." + channel + ".Selection";
		for (const std::string part : Utility::GetKeys(configTree, path)) {
			std::cout << path << ": "<< part << configTree.get<std::string>(path + "." + part + ".Operator")<<configTree.get<std::string>(path + "." + part + ".Threshold") << std::endl;
			cutflows.back().AddCut(part, product, configTree.get<std::string>(path + "." + part + ".Operator"), configTree.get<short>(path + "." + part + ".Threshold"));
		}
		std::cout << std::endl;
	}

	// Register Trigger output
	for (int iChannel = 0; iChannel < channels.size(); iChannel++) {
		std::vector<int> triggerIndex;
		for(int iTrigger = 0; iTrigger < triggerNames.size(); iTrigger++) {
			for(const std::string& triggerName : Utility::GetVector<std::string>(configTree, "Channel." + channels.at(iChannel) + ".Trigger." + product.GetEraSelector())) {
				if(triggerName == triggerNames[iTrigger]) triggerIndex.push_back(iTrigger);
			}
		}

		cutflows[iChannel].AddTriggerOr(triggerIndex, product, channels.at(iChannel));
	}

	// Register branches that will be stored in the output
	product.RegisterTrigger(triggerNames, metTriggerNames, outputTrees);
	product.RegisterMetFilter(metFilterNames, outputTrees);
	product.RegisterOutput(outputTrees, configTree);

	// Initialize Producers
	std::vector<std::shared_ptr<BaseProducer>> producers = {
		std::shared_ptr<TriggerProducer>(new TriggerProducer(configTree, scaleFactorTree, product, triggerNames)),
		std::shared_ptr<METFilterProducer>(new METFilterProducer(configTree, scaleFactorTree, product)),
		std::shared_ptr<MuonProducer>(new MuonProducer(configTree, scaleFactorTree, product.GetEraSelector())),
		std::shared_ptr<ElectronProducer>(new ElectronProducer(configTree, scaleFactorTree)),
		std::shared_ptr<JetProducer>(new JetProducer(configTree, scaleFactorTree, product)),
		std::shared_ptr<DeltaPhiProducer>(new DeltaPhiProducer(configTree, scaleFactorTree)),
	};
	if (!isData) {
		producers.push_back(std::shared_ptr<PileUpWeightProducer>(new PileUpWeightProducer(configTree, scaleFactorTree, product.GetEraSelector())));
		producers.push_back(std::shared_ptr<GenLevelProducer>(new GenLevelProducer(configTree, scaleFactorTree, product.GetEraSelector())));
		producers.push_back(std::shared_ptr<ScaleFactorProducer>(new ScaleFactorProducer(configTree, scaleFactorTree, product.GetEraSelector(), outputFile)));
	}
	if (isFastSim) {
		producers.push_back(std::shared_ptr<FastSimProducer>(new FastSimProducer(configTree, scaleFactorTree, product.GetEraSelector(), outputFile)));
	}

	DataReader dataReader(inputFileName, "Events", isData, isFastSim);
	dataReader.RegisterTrigger(triggerNames, metTriggerNames);
	dataReader.RegisterMetFilter(Utility::GetVector<std::string>(configTree, "METFilter." + product.GetEraSelector()));
	int nEvents = nMaxEvents > 0 ? nMaxEvents : dataReader.GetEntries();

	std::cout << std::endl << "Starting Event loop over " << nEvents << " Events" << std::endl;
	for (int entry = 0; entry < nEvents; ++entry) {
		if (entry % 10000 == 0) {
			std::cout  << "Processed " << int(100*(float)(entry + 1)/nEvents) << "% Events at a rate of " + std::to_string(entry / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count()) + " Hz." << std::endl;
		}

		dataReader.SetEntry(entry);
		for (std::shared_ptr<BaseProducer> producer : producers) {
			//std::cout << producer->Name << std::endl; // Debug Information
			producer->Produce(dataReader, product);
		}

		for (int iTree = 0; iTree < outputTrees.size(); iTree++) {
			cutflows.at(iTree).Count();
			cutflows.at(iTree).FillCutflow();
			if (cutflows.at(iTree).Passed()) { outputTrees.at(iTree)->Fill();}
		}
	}

	std::cout  << "Processed Events at a rate of " + std::to_string(nEvents / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count()) + " Hz." << std::endl;

	for (std::shared_ptr<TTree> tree : outputTrees) {
		int finalNumberOfEvents = tree->GetEntries();
		std::cout << std::setw(20) << tree->GetDirectory()->GetName() << std::setw(10) << tree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nEvents << " (" << 100*(float)finalNumberOfEvents/nEvents << "%)" << std::endl;
	}

	for (int iTree = 0; iTree < outputTrees.size(); iTree++) {
		outputFile.cd();
		outputTrees.at(iTree)->Write(0, TObject::kOverwrite);
		cutflows.at(iTree).WriteOutput();
	}

	product.WriteMetaData(outputFile);
	outputFile.Write(0, TObject::kOverwrite);
	outputFile.Close();

	std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
	std::cout << "Finished event loop (in minutes): " << std::chrono::duration_cast<std::chrono::minutes>(end - start).count() << std::endl;
	std::cout << "Output file created: " + outputFileName << std::endl;

	return 0;
}
