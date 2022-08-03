#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/PileUpWeightProducer.h>

PileUpWeightProducer::PileUpWeightProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {
	Name = "PileUpWeightProducer";
	std::string cmsswBase = std::getenv("CMSSW_BASE");
	pileUpCorrectionSet = correction::CorrectionSet::from_file(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("PileUp." + eraSelector + ".JSON"));
	goldenJsonString = scaleFactorTree.get<std::string>("PileUp." + eraSelector + ".GoldenJSON");
}

void PileUpWeightProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	dataReader.ReadPileUpEntry();
	dataReader.GetPileUpValues(); // No indices required, all weights are read at once

	product.nPdfWeight = dataReader.nPdfWeight;
	std::copy(std::begin(dataReader.pdfWeight), std::end(dataReader.pdfWeight), std::begin(product.pdfWeight));

	product.nScaleWeight = dataReader.nScaleWeight;
	std::copy(std::begin(dataReader.scaleWeight), std::end(dataReader.scaleWeight), std::begin(product.scaleWeight));

	product.nTrueInt = dataReader.nTrueInt;
	product.pileUpWeight     = pileUpCorrectionSet->at(goldenJsonString)->evaluate({dataReader.nTrueInt, "nominal"});
	product.pileUpWeightUp   = pileUpCorrectionSet->at(goldenJsonString)->evaluate({dataReader.nTrueInt, "up"});
	product.pileUpWeightDown = pileUpCorrectionSet->at(goldenJsonString)->evaluate({dataReader.nTrueInt, "down"});

	product.preFire     = dataReader.preFire;
	product.preFireUp   = dataReader.preFireUp;
	product.preFireDown = dataReader.preFireDown;
}

void PileUpWeightProducer::EndJob(TFile &file) {}
