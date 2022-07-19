#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/METFilterProducer.h>

METFilterProducer::METFilterProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, Susy1LeptonProduct &product) {}

void METFilterProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	dataReader.ReadMetFilter();
	dataReader.GetMetFilter();

	for(int i = 0; i < dataReader.metFilterValues.size(); ++i){
		product.metFilterValues[i] = dataReader.metFilterValues[i];
	}
}

void METFilterProducer::EndJob(TFile &file) { }
