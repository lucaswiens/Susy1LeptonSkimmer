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

	// Number of Generated events per mass points (bin size is 25 GeV, bin edges are shifted by 12.5 GeV to center the mass values)
	numberOfGenEvents = std::make_shared<TH2F>("numberOfGenEvents", "numberOfGenEvents", 120, -12.5, 3000-12.5, 120, -12.5, 3000-12.5);
	numberOfGenEventsIsrWeighted = std::make_shared<TH2F>("numberOfGenEventsIsrWeighted", "numberOfGenEventsIsrWeighted", 120, -12.5, 3000-12.5, 120, -12.5, 3000-12.5);
}

void SignalProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	dataReader.ReadGenEntry();
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
	int stopCounter = 0, gluinoCounter = 0, neutralinoCounter = 0, charginoCounter = 0;
	float stopMass = 0., gluinoMass = 0., neutralinoMass = 0., charginoMass = 0.;
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

	product.stopMass       = stopCounter       == 0 ? -999 : floor(stopMass       / stopCounter       + 0.5);
	product.gluinoMass     = gluinoCounter     == 0 ? -999 : floor(gluinoMass     / gluinoCounter     + 0.5);
	product.neutralinoMass = neutralinoCounter == 0 ? -999 : floor(neutralinoMass / neutralinoCounter + 0.5);
	product.charginoMass   = charginoCounter   == 0 ? -999 : floor(charginoMass   / charginoCounter   + 0.5);

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

	numberOfGenEvents->Fill(product.gluinoMass, product.neutralinoMass);
	numberOfGenEventsIsrWeighted->Fill(product.gluinoMass, product.neutralinoMass, product.nIsrWeight);
}

void SignalProducer::EndJob(TFile &file) {
	numberOfGenEvents->Write(0, TObject::kOverwrite);
	numberOfGenEventsIsrWeighted->Write(0, TObject::kOverwrite);
}
