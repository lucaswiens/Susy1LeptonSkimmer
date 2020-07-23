#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/NanoSkimmer.h>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/MuonProducer.h>
//#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/TestProducer.h>

NanoSkimmer::NanoSkimmer(){}

NanoSkimmer::NanoSkimmer(const std::string &inFile, const bool &isData):
	inFile(inFile),
	isData(isData)
	{
		this->outputTree = new TTree();
		this->outputTree->SetName("Events");
		start = std::chrono::steady_clock::now();
		std::cout << "Input file for analysis: " + inFile << std::endl;
	}

void NanoSkimmer::ProgressBar(const int &progress){
	std::string progressBar = "[";

	for (int i = 0; i < progress; i++){
		if (i%2 == 0) progressBar += "#";
	}

	for (int i = 0; i < 100 - progress; i++){
		if (i%2 == 0) progressBar += " ";
	}

	progressBar = progressBar + "] " + std::to_string(progress) + "% of Events processed";
	std::cout << "\r" << progressBar << std::flush;

	if (progress == 100) std::cout << std::endl;

}

void NanoSkimmer::Configure(const float &xSec, const int &era, TTreeReader& reader){
	producers = {
		//std::shared_ptr<TestProducer>(new TestProducer(2017, 20., 2.4, reader)),
		std::shared_ptr<MuonProducer>(new MuonProducer(era, 10, 2.4, 0.5, 1, 4, 0.4, reader)),
	};
}

void NanoSkimmer::EventLoop(const float &xSec, const int &era){

	//TTreeReader preperation
	TFile* inputFile = TFile::Open(inFile.c_str(), "READ");
	TTree* eventTree = (TTree*)inputFile->Get("Events");
	TTreeReader reader(eventTree);

	Configure(xSec, era, reader);

	//Create cutflow histograms
	CutFlow cutflow;

	cutflow.hist = new TH1F();
	cutflow.hist->SetName("cutflow");
	cutflow.hist->SetName("cutflow");
	cutflow.hist->GetYaxis()->SetName("Events");
	cutflow.weight = 1;

	//Begin jobs for all producers
	for (std::shared_ptr<BaseProducer> producer: producers){
		producer->BeginJob(outputTree, isData);
	}

	//Progress bar at 0%
	int processed = 0;
	ProgressBar(0.);

	int nEvents = eventTree->GetEntries();
	while(reader.Next()){
		//Call each producer
		for (unsigned int i = 0; i < producers.size(); i++){
			producers[i]->Produce(cutflow);

			//If the cutflow fails for one producer, reject the event
			if (!cutflow.passed){ break;}
		}

		//Check individual for each channel, if event should be filled
		if (cutflow.passed){
			outputTree->Fill();
		}
		cutflow.passed = true;

		//progress bar
		processed++;
		if (processed % 10000 == 0){
			//int progress = 100*(float)processed/eventTree->GetEntries();
			int progress = 100 * (float) processed / nEvents;
			ProgressBar(progress);
		}
	}

	ProgressBar(100);

	//Print stats
	std::cout << outputTree->GetName() << " analysis: Selected " << outputTree->GetEntries() << " events of " << eventTree->GetEntries() << " (" << 100*(float)outputTree->GetEntries()/eventTree->GetEntries() << "%)" << std::endl;
}

void NanoSkimmer::WriteOutput(const std::string &outFile){
	TFile* file = TFile::Open(outFile.c_str(), "RECREATE");
	outputTree->Write();

	//End jobs for all producers
	for (unsigned int i = 0; i < producers.size(); i++){
		producers[i]->EndJob(file);
	}

	for (CutFlow& cutflow: cutflows){
		cutflow.hist->Write();
		delete cutflow.hist;
	}

	file->Write();
	file->Close();

	end = std::chrono::steady_clock::now();
	std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;

	std::cout << "Output file created: " + outFile << std::endl;
}
