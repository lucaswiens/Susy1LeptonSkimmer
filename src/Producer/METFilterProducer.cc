#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/METFilterProducer.h>

METFilterProducer::METFilterProducer(const int &era): era(era) {}

//void METFilterProducer::BeginJob(std::shared_ptr<TTree> tree, bool &isData, bool &doSystematics) {
//	//Set data bool
//	this->isData = isData;
//	this->doSystematics = doSystematics;
//
//	//Initiliaze TTreeReaderValues
//	//runNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "run");
//	//luminosityNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "luminosityBlock");
//	//eventNumber = std::make_unique<TTreeReaderValue<unsigned long long>>(*reader, "event");
//
//	//flagEeBadScFilter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_eeBadScFilter");
//	//flagHBHENoiseFilter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_HBHENoiseFilter");
//	//flagHBHENoiseIsoFilter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_HBHENoiseIsoFilter");
//	//flagEcalDeadCellTriggerPrimitiveFilter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_EcalDeadCellTriggerPrimitiveFilter");
//	//flagGoodVertices = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_goodVertices");
//	//flagGlobalSuperTightHalo2016Filter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_globalSuperTightHalo2016Filter");
//	//flagBadMuons = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_BadPFMuonFilter");
//
//	//Set Branches of output tree
//	tree->Branch("PassFilters", &PassFilters);
//	tree->Branch("PassFiltersMoriond2017Tight", &PassFiltersMoriond2017Tight);
//}

void METFilterProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	if (isData) {
		const bool &eeBadScFilter = -999;//*flagEeBadScFilter->Get();
		const bool &hbheNoiseFilter = -999;//*flagHBHENoiseFilter->Get();
		const bool &hbheNoiseIsoFilter = -999;//*flagHBHENoiseIsoFilter->Get();
		const bool &ecalDeadCellTriggerPrimitiveFilter = -999;//*flagEcalDeadCellTriggerPrimitiveFilter->Get();
		const bool &goodVertices = -999;//*flagGoodVertices->Get();
		const bool &globalSuperTightHalo2016Filter = -999;//*flagGlobalSuperTightHalo2016Filter->Get();
		const bool &badMuons = -999;//*flagBadMuons->Get();
		PassFilters = eeBadScFilter && hbheNoiseFilter && hbheNoiseIsoFilter && ecalDeadCellTriggerPrimitiveFilter && goodVertices && globalSuperTightHalo2016Filter;
		PassFiltersMoriond2017Tight = PassFilters && !badMuons;
	} else {
		PassFilters = true;
		PassFiltersMoriond2017Tight = true;
	}

	//if (PassFilters) {
	//	cutflow.hist->Fill("Passes MET Filters", cutflow.weight);
	//} else {
	//	cutflow.passed = false;
	//}
}

void METFilterProducer::EndJob(TFile &file) { }
