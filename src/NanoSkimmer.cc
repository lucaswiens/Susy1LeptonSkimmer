#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/NanoSkimmer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/LeptonProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/JetProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/DeltaPhiProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/PileUpWeightProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/METFilterProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/TriggerProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/GenLevelProducer.h>

NanoSkimmer::NanoSkimmer() {}

NanoSkimmer::NanoSkimmer(const std::string &inFile, const std::string &outFile, const bool &isData, const bool &doSystematics):
	inFile(inFile),
	outFile(outFile),
	isData(isData),
	doSystematics(doSystematics)
	{
		start = std::chrono::steady_clock::now();
		std::cout << "Input file for analysis: " + inFile << std::endl;

		file = new TFile(this->outFile.c_str(), "RECREATE");

		if (doSystematics) { // create a tree for each systematic variation
			for (std::string systematic : {"JEC", "JER"}) {
				for (std::string shift : {"Up", "Down"}) {
					TTree* tree = new TTree((systematic + shift).c_str(), (systematic + shift).c_str());
					tree->SetDirectory(file);
					outputTrees.push_back(tree);
				}
			}
		} else { // create just the nominal tree
			TTree *tree = new TTree("nominal", "nominal");
			tree->SetDirectory(file);
			outputTrees.push_back(tree);
		}


	}

void NanoSkimmer::ProgressBar(const int &progress, const int &rate) {
	std::string progressBar = "[";

	for (int i = 0; i < progress - 1; i++) {
		progressBar += "·";
	}
	if (progress == 100) { progressBar += "·";}
	else if (progress % 2 == 0) { progressBar += "c";}
	else { progressBar += "C";}

	for (int i = 0; i < 100 - progress; i++) {
		progressBar += "•";
	}

	progressBar = progressBar + "] " + std::to_string(progress) + "% of Events processed at a rate of " + std::to_string(rate) + " Hz." ;
	std::cout << "\r" << progressBar << std::flush;

	if (progress == 100) std::cout << std::endl;

}

void NanoSkimmer::Configure(const int &era, const char &runPeriod, TTreeReader &reader) {
	for (unsigned int i = 0; i < outputTrees.size(); i++) {
		producers.push_back({
			std::shared_ptr<TriggerProducer>(new TriggerProducer(era, runPeriod, reader)),
			std::shared_ptr<METFilterProducer>(new METFilterProducer(era, reader)),
			std::shared_ptr<LeptonProducer>(new LeptonProducer(era, 10, 2.4, 0.5, 1, 4, 0.4, reader)),
			std::shared_ptr<JetProducer>(new JetProducer(era, 20, 2.4, 0.4, runPeriod, reader)),
			std::shared_ptr<DeltaPhiProducer>(new DeltaPhiProducer(reader)),
		});
		if (!isData) {
			producers.at(i).push_back(std::shared_ptr<PileUpWeightProducer>(new PileUpWeightProducer(era, reader))),
			producers.at(i).push_back(std::shared_ptr<GenLevelProducer>(new GenLevelProducer(era, reader)));
		}

		CutFlow cutflow;
		cutflow.hist = new TH1F();
		cutflow.hist->SetName("cutflow_" + *outputTrees.at(i)->GetName());
		cutflow.hist->GetYaxis()->SetName("N");
		cutflow.weight = 1;
		cutflows.push_back(cutflow);

		//Begin jobs for all producers
		for (std::shared_ptr<BaseProducer> producer: producers.at(i)) {
			producer->BeginJob(outputTrees.at(i), isData, doSystematics);
		}
	}
}

int NanoSkimmer::EventLoop(const int &era, const char &runPeriod, const int &nMaxEvents) {
	//TTreeReader preperation
	TFile *inputFile = TFile::Open(inFile.c_str(), "READ");
	if (inputFile == nullptr) {
		std::cerr << "Input File " << inFile << " is a zombie." << std::endl;
		return -1;
	}
	TTree *eventTree = (TTree*)inputFile->Get("Events");
	TTreeReader reader(eventTree);

	//Progress bar at 0%
	int processed = 0;
	ProgressBar(0, 0);
	int stepSize;
	if (nMaxEvents > 0) { stepSize = std::max<int>(1, (int)nMaxEvents/100);}
	else {stepSize = 10000;}
	int nEvents = eventTree->GetEntries();

	Configure(era, runPeriod, reader);

	while(reader.Next()) {
		//Stop Events Loop after nMaxEvents
		if (nMaxEvents > 0 && processed >= nMaxEvents) { break;}

		//Product for passing newly calculated variables to each producer
		Susy1LeptonProduct product;
		//Call each producer
		for (unsigned int i = 0; i < outputTrees.size(); i++) {
			for (unsigned int j = 0; j < producers.at(i).size(); j++) {
				producers.at(i).at(j)->Produce(cutflows.at(i), &product);
				//If the cutflow fails for one producer, reject the event
				if (!cutflows.at(i).passed) { break;}
			}

			//Check if event should be filled
			if (cutflows.at(i).passed) {
				outputTrees.at(i)->Fill();
			}
			cutflows.at(i).passed = true;
		}

		//progress bar
		processed++;
		if (processed % stepSize == 0) {
			int progress;
			if (nMaxEvents > 0) { progress = 100 * (float) processed / nMaxEvents;}
			else { progress = 100 * (float) processed / nEvents;}
			ProgressBar(progress, processed / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count());
		}

		product.clear();// probably not needed
	}

	ProgressBar(100, processed / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count());

	//Print stats
	for (TTree *tree : outputTrees) {
		finalNumberOfEvents = tree->GetEntries();
		if (nMaxEvents>0){
			std::cout << tree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nMaxEvents << " (" << 100*(float)finalNumberOfEvents/nMaxEvents << "%)" << std::endl;
		} else {
			std::cout << tree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nEvents << " (" << 100*(float)finalNumberOfEvents/nEvents << "%)" << std::endl;
		}
	}

	return 0;
}

void NanoSkimmer::WriteOutput() {
	//End jobs for all producers
	for (unsigned int i = 0; i < outputTrees.size(); i++) {
		for (unsigned int j = 0; j < producers.at(i).size(); j++) {
			producers.at(i).at(j)->EndJob(file);
		}
	}

	file->Write(0, TObject::kOverwrite);
	file->Close();
	delete file;

	end = std::chrono::steady_clock::now();
	std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;
	std::cout << "Output file created: " + outFile << std::endl;
}
