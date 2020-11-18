#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/TriggerProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

TriggerProducer::TriggerProducer(const int &era, const char &runPeriod, TTreeReader &reader):
	BaseProducer(&reader),
	era(era),
	runPeriod(runPeriod)
	{}

void TriggerProducer::BeginJob(TTree *tree, bool &isData, bool &doSystematics) {
	//Set data bool
	this->isData = isData;
	this->doSystematics = doSystematics;

	//Initiliaze TTreeReaderValues
	if(era == 2016){
		hlt_Ele105_CaloIdVT_GsfTrkIdT = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele105_CaloIdVT_GsfTrkIdT");
		hlt_Ele115_CaloIdVT_GsfTrkIdT = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele115_CaloIdVT_GsfTrkIdT");
		hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165");
		hlt_Ele27_WPTight_Gsf = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele27_WPTight_Gsf");
		if (isData) {
			hlt_Ele15_IsoVVVL_PFHT400 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400");
		} else { //Only exists in MC
			hlt_Ele15_IsoVVVL_PFHT400 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele15_IsoVVVL_PFHT400");
		}

		hlt_Mu50 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Mu50");
		hlt_IsoMu24 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_IsoMu24");
		hlt_IsoTkMu24 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_IsoTkMu24");
		if (isData) {
			hlt_Mu15_IsoVVVL_PFHT400 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400");
		} else { //Only exists in MC
			hlt_Mu15_IsoVVVL_PFHT400 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Mu15_IsoVVVL_PFHT400");
		}

		hlt_PFMET110_PFMHT110_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET110_PFMHT110_IDTight");
		hlt_PFMETNoMu110_PFMHTNoMu110_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight");
		hlt_PFMET120_PFMHT120_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET120_PFMHT120_IDTight");
		hlt_PFMETNoMu120_PFMHTNoMu120_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
		if (!(runPeriod == 'G' || runPeriod == 'H')) {
			hlt_Ele15_IsoVVVL_PFHT350 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele15_IsoVVVL_PFHT350");
			hlt_Mu15_IsoVVVL_PFHT350 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Mu15_IsoVVVL_PFHT350");
			hlt_PFMET100_PFMHT100_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET100_PFMHT100_IDTight");
			hlt_PFMETNoMu100_PFMHTNoMu100_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight");
		}

	} else if (era == 2017){
		hlt_Ele115_CaloIdVT_GsfTrkIdT = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele115_CaloIdVT_GsfTrkIdT");
		hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165");
		hlt_Ele35_WPTight_Gsf = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele35_WPTight_Gsf");
		hlt_Ele15_IsoVVVL_PFHT450 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele15_IsoVVVL_PFHT450");

		hlt_Mu50 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Mu50");
		hlt_IsoMu24 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_IsoMu24");
		hlt_IsoTkMu24 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_IsoTkMu24");
		hlt_Mu15_IsoVVVL_PFHT450 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Mu15_IsoVVVL_PFHT450");

		hlt_PFMET100_PFMHT100_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET100_PFMHT100_IDTight");
		hlt_PFMETNoMu100_PFMHTNoMu100_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight");
		hlt_PFMET110_PFMHT110_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET110_PFMHT110_IDTight");
		hlt_PFMETNoMu110_PFMHTNoMu110_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight");
		hlt_PFMET120_PFMHT120_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET120_PFMHT120_IDTight");
		hlt_PFMETNoMu120_PFMHTNoMu120_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");

	} else if (era ==2018) {
		hlt_Ele115_CaloIdVT_GsfTrkIdT = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele115_CaloIdVT_GsfTrkIdT");
		hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165");
		hlt_Ele35_WPTight_Gsf = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele35_WPTight_Gsf");
		hlt_Ele15_IsoVVVL_PFHT450 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Ele15_IsoVVVL_PFHT450");

		hlt_Mu50 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Mu50");
		hlt_IsoMu24 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_IsoMu24");
		hlt_IsoTkMu24 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_IsoTkMu24");
		hlt_Mu15_IsoVVVL_PFHT450 = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_Mu15_IsoVVVL_PFHT450");

		hlt_PFMET100_PFMHT100_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET100_PFMHT100_IDTight");
		hlt_PFMETNoMu100_PFMHTNoMu100_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight");
		hlt_PFMET110_PFMHT110_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET110_PFMHT110_IDTight");
		hlt_PFMETNoMu110_PFMHTNoMu110_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight");
		hlt_PFMET120_PFMHT120_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMET120_PFMHT120_IDTight");
		hlt_PFMETNoMu120_PFMHTNoMu120_IDTight = std::make_unique<TTreeReaderValue<bool>>(*reader, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
	}

	//Set Branches of output tree
	tree->Branch("HLTElectronOr", &HLT_EleOr);
	tree->Branch("HLTMuonOr", &HLT_MuOr);
	tree->Branch("HLTLeptonOr", &HLT_LepOr);
	tree->Branch("HLTMETOr", &HLT_MetOr);
}

void TriggerProducer::Produce(CutFlow &cutflow, Susy1LeptonProduct *product) {
	//Initialize all variables as -999
	HLT_EleOr = false;
	HLT_MuOr = false;
	HLT_LepOr = false;
	HLT_MetOr = false;

	if (era ==2016) {
		HLT_EleOr = *hlt_Ele105_CaloIdVT_GsfTrkIdT->Get() || *hlt_Ele115_CaloIdVT_GsfTrkIdT->Get() || *hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165->Get() || *hlt_Ele27_WPTight_Gsf->Get() || *hlt_Ele15_IsoVVVL_PFHT400->Get();
		HLT_MuOr = *hlt_Mu50->Get() || *hlt_IsoMu24->Get() || *hlt_IsoTkMu24->Get() || *hlt_Mu15_IsoVVVL_PFHT400->Get();
		HLT_MetOr = *hlt_PFMET110_PFMHT110_IDTight->Get() || *hlt_PFMETNoMu110_PFMHTNoMu110_IDTight->Get() || *hlt_PFMET120_PFMHT120_IDTight->Get() || *hlt_PFMETNoMu120_PFMHTNoMu120_IDTight->Get();
		if (!(runPeriod == 'G' || runPeriod == 'H')) { // These Triggers only exist for run before runPeriod G or H
			HLT_EleOr = HLT_EleOr || *hlt_Ele15_IsoVVVL_PFHT350->Get();
			HLT_MuOr = HLT_MuOr || *hlt_Mu15_IsoVVVL_PFHT350->Get();
			HLT_MetOr = HLT_MetOr || *hlt_PFMET100_PFMHT100_IDTight->Get() || *hlt_PFMETNoMu100_PFMHTNoMu100_IDTight->Get();;
		}
	} else if (era ==2017) {
		HLT_EleOr = *hlt_Ele115_CaloIdVT_GsfTrkIdT->Get() || *hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165->Get() || *hlt_Ele35_WPTight_Gsf->Get() || *hlt_Ele15_IsoVVVL_PFHT450->Get();
		HLT_MuOr = *hlt_Mu50->Get() || *hlt_IsoMu24->Get() || *hlt_IsoTkMu24->Get() || *hlt_Mu15_IsoVVVL_PFHT450->Get();
		HLT_MetOr = *hlt_PFMET100_PFMHT100_IDTight->Get() || *hlt_PFMETNoMu100_PFMHTNoMu100_IDTight->Get() || *hlt_PFMET110_PFMHT110_IDTight->Get() || *hlt_PFMETNoMu110_PFMHTNoMu110_IDTight->Get() || *hlt_PFMET120_PFMHT120_IDTight->Get() || *hlt_PFMETNoMu120_PFMHTNoMu120_IDTight->Get();
	} else if (era ==2018) {
		HLT_EleOr = *hlt_Ele115_CaloIdVT_GsfTrkIdT->Get() || *hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165->Get() || *hlt_Ele35_WPTight_Gsf->Get() || *hlt_Ele15_IsoVVVL_PFHT450->Get();
		HLT_MuOr = *hlt_Mu50->Get() || *hlt_IsoMu24->Get() || *hlt_IsoTkMu24->Get() || *hlt_Mu15_IsoVVVL_PFHT450->Get();
		HLT_MetOr = *hlt_PFMET100_PFMHT100_IDTight->Get() || *hlt_PFMETNoMu100_PFMHTNoMu100_IDTight->Get() || *hlt_PFMET110_PFMHT110_IDTight->Get() || *hlt_PFMETNoMu110_PFMHTNoMu110_IDTight->Get() || *hlt_PFMET120_PFMHT120_IDTight->Get() || *hlt_PFMETNoMu120_PFMHTNoMu120_IDTight->Get();
	}

	HLT_LepOr = HLT_EleOr || HLT_MuOr;

	std::string cutName("Lepton Passed Trigger");
	cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
}

void TriggerProducer::EndJob(TFile *file) {
}
