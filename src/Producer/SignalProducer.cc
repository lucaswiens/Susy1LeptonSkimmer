#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/SignalProducer.h>

SignalProducer::SignalProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector, TFile &outputFile) {
	Name = "SignalProducer";
	std::string cmsswBase = std::getenv("CMSSW_BASE");

	stopPdgId       = configTree.get<int>("Producer.Signal.StopPdgId");
	gluinoPdgId     = configTree.get<int>("Producer.Signal.GluinoPdgId");
	neutralinoPdgId = configTree.get<int>("Producer.Signal.NeutralinoPdgId");
	charginoPdgId   = configTree.get<int>("Producer.Signal.CharginoPdgId");

	// Open json file with cross section
	pt::read_json(std::string(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/fastsim/" + configTree.get<std::string>("Producer.Signal.CrossSection")).c_str(), xSectionTree);
}

void SignalProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	dataReader.ReadGenEntry();
	product.stopMass       = -999;
	product.gluinoMass     = -999;
	product.neutralinoMass = -999;
	product.charginoMass   = -999;

	product.susyXSectionNLO      = -999;
	product.susyXSectionNLLO     = -999;
	product.susyXSectionNLOUp    = -999;
	product.susyXSectionNLLOUp   = -999;
	product.susyXSectionNLODown  = -999;
	product.susyXSectionNLLODown = -999;

	/*******************************************************************************************************************************************************************
	*   Method taken from:                                                                                                                                             *
	*   https://github.com/CERN-PH-CMG/cmgtools-lite/blob/ad5a23dd087e6e0ba826c8e03b39050a46a18c38/TTHAnalysis/python/analyzers/susyParameterScanAnalyzer.py#L63-L77   *
	*******************************************************************************************************************************************************************/
	int stopCounter, gluinoCounter, neutralinoCounter, charginoCounter;
	float stopMass, gluinoMass, neutralinoMass, charginoMass;
	for (int iGen = 0; iGen < dataReader.nGenPart; iGen++) {
		dataReader.GetGenValues(iGen);
		if (std::abs(dataReader.genPdgId) == stopPdgId) {
			stopCounter++;
			stopMass += dataReader.genMass;
		} else if (std::abs(dataReader.genPdgId) == gluinoPdgId) {
			gluinoCounter++;
			gluinoMass += dataReader.genMass;
		} else if (std::abs(dataReader.genPdgId) == neutralinoPdgId) {
			neutralinoCounter++;
			neutralinoMass += dataReader.genMass;
		} else if (std::abs(dataReader.genPdgId) == charginoPdgId) {
			charginoCounter++;
			charginoMass += dataReader.genMass;
		}
	}

	product.stopMass       = floor(stopMass       / stopCounter + 0.5);
	product.gluinoMass     = floor(gluinoMass     / gluinoCounter + 0.5);
	product.neutralinoMass = floor(neutralinoMass / neutralinoCounter + 0.5);
	product.charginoMass   = floor(charginoMass   / charginoCounter + 0.5);

	// Check if the mass exists in the cross section json file... Sometimes the masses are set to strange unexpected values, one might needs to correct this later or discard those events
	boost::optional< pt::ptree& > child = xSectionTree.get_child_optional("NLO+NLL Cross section [pb]." + std::to_string((int)product.gluinoMass));
	if (!child) { return;}

	product.susyXSectionNLO                  = xSectionTree.get<float>("NLO+NLL Cross section [pb]." + std::to_string(product.gluinoMass));
	const float &susyXSectionNLOUncertainty  = xSectionTree.get<float>("Uncertainty (NLO + NLL) [%]." + std::to_string(product.gluinoMass));
	product.susyXSectionNLLO                 = xSectionTree.get<float>("NNLO+NNLL Cross section [pb]." + std::to_string(product.gluinoMass));
	const float &susyXSectionNLLOUncertainty = xSectionTree.get<float>("Uncertainty (NNLOapprox + NNLL) [%]." + std::to_string(product.gluinoMass));

	product.susyXSectionNLOUp  = product.susyXSectionNLO  + susyXSectionNLOUncertainty;
	product.susyXSectionNLLOUp = product.susyXSectionNLLO + susyXSectionNLLOUncertainty;

	product.susyXSectionNLODown  = product.susyXSectionNLO  - susyXSectionNLOUncertainty;
	product.susyXSectionNLLODown = product.susyXSectionNLLO - susyXSectionNLLOUncertainty;
}

void SignalProducer::EndJob(TFile &file) {}
