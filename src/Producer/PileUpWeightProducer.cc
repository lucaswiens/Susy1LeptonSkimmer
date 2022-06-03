#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/PileUpWeightProducer.h>

PileUpWeightProducer::PileUpWeightProducer(const int &era): era(era) {}

//void PileUpWeightProducer::BeginJob(std::shared_ptr<TTree> tree, bool &isData, bool &doSystematics) {
//	//Set data bool
//	this->isData = isData;
//	this->doSystematics = doSystematics;
//
//	if (!isData) {
//		TString puFilePath = "$CMSSW_BASE/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/pileup/";
//
//		pileupPathData = {
//			{2016, puFilePath + "PileupData_GoldenJSON_Full2016.root"},
//			{2017, puFilePath + "PileupHistogram-goldenJSON-13tev-2017-99bins_withVar.root"},
//			{2018, puFilePath + "PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root"}
//		};
//		pileupPathMC = {
//			{2016, puFilePath + "pileup_profile_Summer16.root"},
//			{2017, puFilePath + "mcPileup2017.root"},
//			{2018, puFilePath + "mcPileup2018.root"}
//		};
//
//		TFile *pileupFileData = TFile::Open(pileupPathData[era], "READ");
//		TH1D *pileupData = static_cast<TH1D*>(pileupFileData->Get("pileup"));
//		TH1D *pileupDataPlus = static_cast<TH1D*>(pileupFileData->Get("pileup_plus"));
//		TH1D *pileupDataMinus = static_cast<TH1D*>(pileupFileData->Get("pileup_minus"));
//
//		TFile *pileupFileMC = TFile::Open(pileupPathMC[era], "READ");
//		TH1D *pileupMC = static_cast<TH1D*>(pileupFileMC->Get("pu_mc"));
//
//		wc = new WeightCalculator(true, false);
//
//		pileupRatio = wc->Ratio(pileupMC, pileupData, false);
//		pileupRatioPlus = wc->Ratio(pileupMC, pileupDataPlus, false);
//		pileupRatioMinus = wc->Ratio(pileupMC, pileupDataMinus, false);
//
//		pileupFileData->Close();
//		pileupFileMC->Close();
//
//		//Initiliaze TTreeReaderValues
//		//pvNumber = std::make_unique<TTreeReaderValue<float>>(*reader, "Pileup_nTrueInt");
//
//		//Set Branches of output tree
//		tree->Branch("PileUpWeight", &weight);
//		tree->Branch("PileUpWeightPlus", &weightPlus);
//		tree->Branch("PileUpWeightMinus", &weightMinus);
//	}
//
//}

void PileUpWeightProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	//Initialize all variables as -999
	nPV = -999;
	weight = -999;
	weightPlus = -999;
	weightMinus = -999;

	if (!isData) {
		nPV = -999;//*pvNumber->Get();
		weight = 1; weightPlus = 1; weightMinus = 1;
		weight = wc->GetWeight(nPV, pileupRatio);
		weightPlus = wc->GetWeight(nPV, pileupRatioPlus);
		weightMinus = wc->GetWeight(nPV, pileupRatioMinus);
	}

	//cutflow.hist->Fill("PileUpWeight", cutflow.weight);
}

void PileUpWeightProducer::EndJob(TFile &file) {
	if (!isData) {
		delete wc;
	}
}
