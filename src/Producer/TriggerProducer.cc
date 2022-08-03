#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/TriggerProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

TriggerProducer::TriggerProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product) {
	Name = "TriggerProducer";
}

//TODO Think about missing triggers
void TriggerProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	dataReader.ReadTrigger();
	dataReader.GetTrigger();

	product.hltEleOr = false;
	product.hltMuOr = false;
	product.hltMetOr = false;

	for(int i = 0; i < dataReader.triggerValues.size(); ++i){
		product.triggerValues[i] = dataReader.triggerValues[i];
		product.hltEleOr = product.hltEleOr || dataReader.triggerValues[i];
		product.hltMuOr = product.hltMuOr || dataReader.triggerValues[i];
	}

	for(int i = 0; i < dataReader.metTriggerValues.size(); ++i){
		product.metTriggerValues[i] = dataReader.metTriggerValues[i];
		product.hltMetOr = product.hltMetOr || dataReader.metTriggerValues[i];
	}
}

void TriggerProducer::EndJob(TFile &file) {}
