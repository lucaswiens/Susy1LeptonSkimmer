#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/EventInformationProducer.h>

EventInformationProducer::EventInformationProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {
	Name = "EventInformationProducer";
	std::string cmsswBase = std::getenv("CMSSW_BASE");
	pt::read_json((cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/" + configTree.get<std::string>("GoldenJSON." + eraSelector)).c_str(), goldenJsonTree);
}

void EventInformationProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	dataReader.ReadEventInfo();

	product.eventNumber        = dataReader.eventNumber;
	product.runNumber          = dataReader.runNumber;
	product.isInLumiBlockRange = false;

	// check if run exist in golden json
	boost::optional<boost::property_tree::ptree&> goldenLumiRanges = goldenJsonTree.get_child_optional(std::to_string(product.runNumber));
	if (goldenLumiRanges  == boost::none) { return;}

	// iterate over good lumi sections
	for (const std::pair<std::string, boost::property_tree::ptree> &luminosityRanges : goldenLumiRanges.get()) {
		std::vector<int> luminosityRange;
		for(const std::pair<std::string, boost::property_tree::ptree> &lumiRange : luminosityRanges.second.get_child("")) {
			luminosityRange.push_back(lumiRange.second.get<int>(""));
		}
		assert(luminosityRange.size() == 2); // there should be only the min and max values per entry

		if (luminosityRange.at(0) <= dataReader.luminosityBlock && dataReader.luminosityBlock <= luminosityRange.at(1)) {
			product.isInLumiBlockRange = true;
			break;
		}
	}
}

void EventInformationProducer::EndJob(TFile &file) {}
