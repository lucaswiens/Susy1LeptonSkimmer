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

NanoSkimmer::NanoSkimmer(const std::string &inFile, const bool &isData):
	inFile(inFile),
	isData(isData)
	{
		start = std::chrono::steady_clock::now();

		this->outputTree = new TTree();
		this->outputTree->SetName("Events");
		outputTree->SetAutoFlush(400000); // Flush tree to Disk to avoid memory leaks

		std::cout << "Input file for analysis: " + inFile << std::endl;
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

void NanoSkimmer::Configure(const float &xSec, const int &era, const char &runPeriod, TTreeReader& reader) {
	producers = {
		std::shared_ptr<TriggerProducer>(new TriggerProducer(era, reader)),
		std::shared_ptr<METFilterProducer>(new METFilterProducer(era, reader)),
		std::shared_ptr<LeptonProducer>(new LeptonProducer(era, 10, 2.4, 0.5, 1, 4, 0.4, reader)),
		std::shared_ptr<JetProducer>(new JetProducer(era, 20, 2.4, 0.4, runPeriod, reader)),
		std::shared_ptr<DeltaPhiProducer>(new DeltaPhiProducer(reader)),
	};
	if (!isData) {
		producers.push_back(std::shared_ptr<PileUpWeightProducer>(new PileUpWeightProducer(era, reader))),
		producers.push_back(std::shared_ptr<GenLevelProducer>(new GenLevelProducer(era, reader)));
	}
}

void NanoSkimmer::EventLoop(const float &xSec, const int &era, const char &runPeriod, const int &nMaxEvents) {

	//TTreeReader preperation
	TFile* inputFile = TFile::Open(inFile.c_str(), "READ");
	TTree* eventTree = (TTree*)inputFile->Get("Events");
	TTreeReader reader(eventTree);

	Configure(xSec, era, runPeriod, reader);

	cutflow.hist = new TH1F();
	cutflow.hist->SetName("cutflow");
	cutflow.hist->SetName("cutflow");
	cutflow.hist->GetYaxis()->SetName("Events");
	cutflow.weight = 1;


	//Begin jobs for all producers
	for (std::shared_ptr<BaseProducer> producer: producers) {
		producer->BeginJob(outputTree, isData);
	}

	//Progress bar at 0%
	int processed = 0;
	ProgressBar(0., 0.);
	int nProcessedEvents = 0;
	int stepSize;
	if (nMaxEvents > 0) { stepSize = std::min((int)nMaxEvents/100, 0);}
	else {stepSize = 10000;}

	int nEvents = eventTree->GetEntries();
	while(reader.Next()) {
		//Stop Events Loop after nMaxEvents
		if (nMaxEvents > 0 && nProcessedEvents >= nMaxEvents){break;}
		nProcessedEvents++;
		//Product for passing newly calculated variables to each producer
		Susy1LeptonProduct product;
		//Call each producer
		for (unsigned int i = 0; i < producers.size(); i++) {
			producers[i]->Produce(cutflow, &product);

			//If the cutflow fails for one producer, reject the event
			if (!cutflow.passed) { break;}
		}

		//Check individual for each channel, if event should be filled
		if (cutflow.passed) {
			outputTree->Fill();
		}
		cutflow.passed = true;

		//progress bar
		processed++;
		if (processed % stepSize == 0) {
			int progress;
			if (nMaxEvents > 0) { progress = 100 * (float) processed / nMaxEvents;}
			else { progress = 100 * (float) processed / nEvents;}
			ProgressBar(progress, processed / std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count());
		}
	}

	ProgressBar(100, 0);

	finalNumberOfEvents = outputTree->GetEntries();

	//Print stats
	if (nMaxEvents>0){
		std::cout << outputTree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nMaxEvents << " (" << 100*(float)nEvents/nMaxEvents << "%)" << std::endl;
	} else {
		std::cout << outputTree->GetName() << " analysis: Selected " << finalNumberOfEvents << " events of " << nEvents << " (" << 100*(float)finalNumberOfEvents/nEvents << "%)" << std::endl;
	}
}

void NanoSkimmer::WriteOutput(const std::string &outFile) {
	TFile* file = TFile::Open(outFile.c_str(), "RECREATE");
	outputTree->Write();

	//End jobs for all producers
	for (unsigned int i = 0; i < producers.size(); i++) {
		producers[i]->EndJob(file);
	}

	cutflow.hist->Write();

	file->Write(0, TObject::kOverwrite);
	file->Close();

	end = std::chrono::steady_clock::now();
	std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;

	std::cout << "Output file created: " + outFile << std::endl;
}
