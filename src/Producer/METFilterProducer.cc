#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/METFilterProducer.h>

METFilterProducer::METFilterProducer(const int& era, TTreeReader& reader):
	BaseProducer(&reader),
	era(era)
	{}

void METFilterProducer::BeginJob(TTree* tree, bool &isData) {
	//Set data bool
	this->isData = isData;

	//Initiliaze TTreeReaderValues
	runNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "run");
	luminosityNumber = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "luminosityBlock");
	eventNumber = std::make_unique<TTreeReaderValue<unsigned long long>>(*reader, "event");

	flagEeBadScFilter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_eeBadScFilter");
	flagHBHENoiseFilter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_HBHENoiseFilter");
	flagHBHENoiseIsoFilter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_HBHENoiseIsoFilter");
	flagEcalDeadCellTriggerPrimitiveFilter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_EcalDeadCellTriggerPrimitiveFilter");
	flagGoodVertices = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_goodVertices");
	flagGlobalSuperTightHalo2016Filter = std::make_unique<TTreeReaderValue<bool>>(*reader, "Flag_globalSuperTightHalo2016Filter");
	/* badMuon flag does not exist anymore?
	if hasattr(event, "Flag_badMuons"):
		flag_badMuons = getattr(event, "Flag_badMuons")
		flag_duplicateMuons = getattr(event, "Flag_duplicateMuons")
	else :
		flag_duplicateMuons = True
		flag_badMuons = getattr(event,"Flag_BadPFMuonFilter")
	*/

	SetCollection(this->isData);

	//Set Branches of output tree
	tree->Branch("PassFilters", &PassFilters);
	tree->Branch("PassFiltersMoriond2017Tight", &PassFiltersMoriond2017Tight);
	tree->Branch("PassCSCFilterList", &PassCSCFilterList);
}

void METFilterProducer::Produce(CutFlow& cutflow, Susy1LeptonProduct *product) {
	if (isData) {
		const bool& eeBadScFilter = *flagEeBadScFilter->Get();
		const bool& hBHENoiseFilter = *flagHBHENoiseFilter->Get();
		const bool& hBHENoiseIsoFilter = *flagHBHENoiseIsoFilter->Get();
		const bool& ecalDeadCellTriggerPrimitiveFilter = *flagEcalDeadCellTriggerPrimitiveFilter->Get();
		const bool& goodVertices = *flagGoodVertices->Get();
		const bool& globalSuperTightHalo2016Filter = *flagGlobalSuperTightHalo2016Filter->Get();
		PassFilters = eeBadScFilter && hBHENoiseFilter && hBHENoiseIsoFilter && ecalDeadCellTriggerPrimitiveFilter && goodVertices && globalSuperTightHalo2016Filter;
		PassCSCFilterList = true; //What is the point?
		//PassFiltersMoriond2017Tight = PassFilters && !badMuons && !duplicateMuons; //not in NanoAOD?
		/*
		# check MET text filter files
		if (runNr,lumiNr,eventNr) in filterList:
			#print "yes", runNr,lumiNr,eventNr
			PassCSCFilterList = False
		else:
			#print "no", runNr,lumiNr,eventNr
			PassCSCFilterList = True

		# check filters present in event
		*/
	} else {
		PassFilters = true;
		PassFiltersMoriond2017Tight = true;
		PassCSCFilterList = true;
	}

	if (PassFilters) {
		std::string cutName("Passes MET Filters");
		cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
	} else {
		cutflow.passed = false;
	}
}

void METFilterProducer::EndJob(TFile* file) {
}
