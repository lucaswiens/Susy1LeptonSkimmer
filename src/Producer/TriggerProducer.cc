#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/TriggerProducer.h>

TriggerProducer::TriggerProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product) {
	Name = "TriggerProducer";
}

void TriggerProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	dataReader.ReadTrigger();
	dataReader.GetTrigger();

	product.hltEleOr = false;
	product.hltMuOr = false;
	product.hltMetOr = false;

	for(int i = 0; i < dataReader.triggerValues.size(); ++i){
		product.triggerValues[i] = dataReader.triggerValues.at(i);

		std::string triggerName = dataReader.triggerNames[i];
		if (triggerName.find("Ele") != std::string::npos) {
			product.hltEleOr = product.hltEleOr || dataReader.triggerValues[i];
		} else if (triggerName.find("Mu") != std::string::npos) {
			product.hltMuOr = product.hltMuOr || dataReader.triggerValues[i];
		}
	}

	for(int i = 0; i < dataReader.metTriggerValues.size(); ++i){
		product.metTriggerValues[i] = dataReader.metTriggerValues[i];
		product.hltMetOr = product.hltMetOr || dataReader.metTriggerValues[i];
	}
}

void TriggerProducer::EndJob(TFile &file) {}
