#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/ElectronProducer.h>

ElectronProducer::ElectronProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree) {
	Name = "ElectronProducer";
	electronGoodPtCut               = configTree.get<float>("Producer.Electron.Pt.Good");
	electronVetoPtCut               = configTree.get<float>("Producer.Electron.Pt.Veto");
	electronEtaCut                  = configTree.get<float>("Producer.Electron.Eta");
	electronGoodIsoCut              = configTree.get<float>("Producer.Electron.Iso.Good");
	electronVetoIsoCut              = configTree.get<float>("Producer.Electron.Iso.Veto");
	electronAntiIsoCut              = configTree.get<float>("Producer.Electron.Iso.Anti");
	electronGoodCutBasedIdCut       = configTree.get<char>("Producer.Electron.CutBasedId.Good");
	electronVetoCutBasedIdCut       = configTree.get<char>("Producer.Electron.CutBasedId.Veto");
	electronAntiIsCutBasedIdCut     = configTree.get<char>("Producer.Electron.CutBasedId.Anti.Is");
	electronAntiIsNotCutBasedIdCut  = configTree.get<char>("Producer.Electron.CutBasedId.Anti.Not");
	electronGoodNumberOfLostHitsCut = configTree.get<int>("Producer.Electron.NumberOfLostHits.Good");

	std::cout << std::endl <<
		"The following cuts are applied to Electrons:" << std::endl <<
		"|Eta|                < " << electronEtaCut                  << std::endl <<
		"GoodPt               > " << electronGoodPtCut               << std::endl <<
		"GoodIso              < " << electronGoodIsoCut              << std::endl <<
		"GoodCutBasedId       = " << electronGoodCutBasedIdCut       << std::endl <<
		"VetoPt               > " << electronVetoPtCut               << std::endl <<
		"VetoIso              < " << electronVetoIsoCut              << std::endl <<
		"VetoCutBasedId       = " << electronVetoCutBasedIdCut       << std::endl <<
		"AntiIso              < " << electronAntiIsoCut              << std::endl <<
		"AntiCutBasedId       = " << electronAntiIsCutBasedIdCut     << std::endl <<
		"AntiCutBasedId      != " << electronAntiIsNotCutBasedIdCut  << std::endl <<
		"GoodNumberOfLostHits = " << electronGoodNumberOfLostHitsCut << std::endl;
}

void ElectronProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	dataReader.ReadElectronEntry();
	assert(dataReader.nElectron < product.nMax);

	int electronCounter = 0, goodElectronCounter = 0, vetoElectronCounter = 0, antiSelectedElectronCounter = 0;
	for (int iElectron = 0; iElectron < dataReader.nElectron; iElectron++) {
		dataReader.GetElectronValues(iElectron);

		if (dataReader.electronPt < electronVetoPtCut || abs(dataReader.electronEta) > electronEtaCut) { continue;}

		product.electronPt[electronCounter] = dataReader.electronPt;
		product.electronEta[electronCounter] = dataReader.electronEta;
		product.electronPhi[electronCounter] = dataReader.electronPhi;
		product.electronMass[electronCounter] = dataReader.electronMass;
		product.electronDxy[electronCounter] = dataReader.electronDxy;
		product.electronDz[electronCounter] = dataReader.electronDz;
		product.electronECorr[electronCounter] = dataReader.electronECorr;
		product.electronMiniIso[electronCounter] = dataReader.electronMiniIso;
		//product.electronIso03[electronCounter] = dataReader.electronIso03;
		//product.electronIso04[electronCounter] = dataReader.electronIso04;
		product.electronRelJetIso[electronCounter] = dataReader.electronRelJetIso;
		if (product.GetIsFastSim()) {
			product.electronEnergyScaleUp[electronCounter] = dataReader.electronEnergyScaleUp;
			product.electronEnergyScaleDown[electronCounter] = dataReader.electronEnergyScaleDown;
			product.electronEnergySigmaUp[electronCounter] = dataReader.electronEnergySigmaUp;
			product.electronEnergySigmaDown[electronCounter] = dataReader.electronEnergySigmaDown;
		}

		product.electronCharge[electronCounter] = dataReader.electronCharge;
		product.electronCutBasedId[electronCounter] = dataReader.electronCutBasedId;
		product.electronConvVeto[electronCounter] = dataReader.electronConvVeto;

		product.electronTightMvaId[electronCounter] = dataReader.electronTightMvaId;
		product.electronMediumMvaId[electronCounter] = dataReader.electronMediumMvaId;
		product.electronLooseMvaId[electronCounter] = dataReader.electronLooseMvaId;

		product.electronTightId[electronCounter] = dataReader.electronCutBasedId >= 4;
		product.electronMediumId[electronCounter] = dataReader.electronCutBasedId >= 3;
		product.electronLooseId[electronCounter] = dataReader.electronCutBasedId >= 2;
		product.electronVetoId[electronCounter] = dataReader.electronCutBasedId >= 1;

		product.electronIsGood[electronCounter] = dataReader.electronPt > electronGoodPtCut &&
								dataReader.electronMiniIso < electronGoodIsoCut &&
								dataReader.electronConvVeto &&
								dataReader.electronNLostHits == electronGoodNumberOfLostHitsCut &&
								dataReader.electronIdMap.at(electronGoodCutBasedIdCut);

		product.electronIsVeto[electronCounter] = dataReader.electronPt <= electronGoodPtCut &&
								dataReader.electronMiniIso < electronVetoIsoCut &&
								dataReader.electronIdMap.at(electronVetoCutBasedIdCut);

		product.electronIsAntiSelected[electronCounter] = dataReader.electronMiniIso < electronAntiIsoCut && // FIXME Check if is Medium but not Tight is correct
								!dataReader.electronIdMap.at(electronAntiIsNotCutBasedIdCut);

		if (product.electronIsGood[electronCounter]) { goodElectronCounter++;}
		if (product.electronIsVeto[electronCounter]) { vetoElectronCounter++;}
		if (product.electronIsAntiSelected[electronCounter]) { antiSelectedElectronCounter++;}
		electronCounter++;
	}

	//Overwrite number of leptons by excluding leptons that do not survive the cuts
	product.nElectron = electronCounter;
	product.nGoodElectron = goodElectronCounter;
	product.nVetoElectron = vetoElectronCounter;
	product.nAntiSelectedElectron = antiSelectedElectronCounter;

	product.nLepton = product.nMuon + product.nElectron;
	product.nGoodLepton = product.nGoodMuon + product.nGoodElectron;
	product.nVetoLepton = product.nVetoMuon + product.nVetoElectron;
}

void ElectronProducer::EndJob(TFile &file) {}
