#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/ElectronProducer.h>

ElectronProducer::ElectronProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {
	electronGoodPtCut               = configTree.get<double>("Producer.Electron.Pt.Good");
	electronVetoPtCut               = configTree.get<double>("Producer.Electron.Pt.Veto");
	electronEtaCut                  = configTree.get<double>("Producer.Electron.Eta");
	electronGoodIsoCut              = configTree.get<double>("Producer.Electron.Iso.Good");
	electronVetoIsoCut              = configTree.get<double>("Producer.Electron.Iso.Veto");
	electronAntiIsoCut              = configTree.get<double>("Producer.Electron.Iso.Anti");
	electronGoodCutBasedIdCut       = configTree.get<char>("Producer.Electron.CutBasedId.Good");
	electronVetoCutBasedIdCut       = configTree.get<char>("Producer.Electron.CutBasedId.Veto");
	electronAntiIsCutBasedIdCut     = configTree.get<char>("Producer.Electron.CutBasedId.Anti.Is");
	electronAntiIsNotCutBasedIdCut  = configTree.get<char>("Producer.Electron.CutBasedId.Anti.Not");
	electronGoodNumberOfLostHitsCut = configTree.get<int>("Producer.Electron.NumberOfLostHits.Good");

	std::cout << std::endl <<
		"The following cuts are applied to Electrons:" << std::endl <<
		"|Eta|                < " << electronEtaCut                 << std::endl <<
		"GoodPt               > " << electronGoodPtCut              << std::endl <<
		"GoodIso              < " << electronGoodIsoCut             << std::endl <<
		"GoodCutBasedId       = " << electronGoodCutBasedIdCut      << std::endl <<
		"VetoPt               > " << electronVetoPtCut              << std::endl <<
		"VetoIso              < " << electronVetoIsoCut             << std::endl <<
		"VetoCutBasedId       = " << electronVetoCutBasedIdCut      << std::endl <<
		"AntiIso              < " << electronAntiIsoCut             << std::endl <<
		"AntiCutBasedId       = " << electronAntiIsCutBasedIdCut    << std::endl <<
		"AntiCutBasedId      != " << electronAntiIsNotCutBasedIdCut << std::endl <<
		"GoodNumberOfLostHits = " << electronGoodNumberOfLostHitsCut   << std::endl;
}

void ElectronProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	//product.nElectron = 0;
	//product.nGoodElectron = 0;
	//product.nVetoElectron = 0;

	dataReader.ReadElectronEntry();
	product.nElectron = dataReader.nElectron;
	assert(product.nElectron < product.nMax);

	std::size_t electronCounter = 0, goodElectronCounter = 0, vetoElectronCounter = 0, antiSelectedElectronCounter = 0;
	if (product.nElectron > 0) {
		for (std::size_t iElectron = 0; iElectron < dataReader.nElectron; iElectron++) {
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
			product.electronEnergyScaleUp[electronCounter] = dataReader.electronEnergyScaleUp;
			product.electronEnergyScaleDown[electronCounter] = dataReader.electronEnergyScaleDown;
			product.electronEnergySigmaUp[electronCounter] = dataReader.electronEnergySigmaUp;
			product.electronEnergySigmaDown[electronCounter] = dataReader.electronEnergySigmaDown;

			product.electronCharge[electronCounter] = dataReader.electronCharge;
			product.electronCutBasedId[electronCounter] = dataReader.electronCutBasedId;
			product.electronConvVeto[electronCounter] = dataReader.electronConvVeto;

			product.electronLooseMvaId[electronCounter] = dataReader.electronLooseMvaId;
			product.electronMediumMvaId[electronCounter] = dataReader.electronMediumMvaId;
			product.electronTightMvaId[electronCounter] = dataReader.electronTightMvaId;

			product.electronTightId[electronCounter] = dataReader.electronCutBasedId >= 4;
			product.electronMediumId[electronCounter] = dataReader.electronCutBasedId >= 3;
			product.electronLooseId[electronCounter] = dataReader.electronCutBasedId >= 2;
			product.electronVetoId[electronCounter] = dataReader.electronCutBasedId >= 1;

			product.electronIsGood[electronCounter] = dataReader.electronPt > electronGoodPtCut &&
									dataReader.electronMiniIso < electronGoodIsoCut &&
									dataReader.electronConvVeto &&
									dataReader.electronNLostHits == electronGoodNumberOfLostHitsCut &&
									(electronGoodCutBasedIdCut == 'T'? product.electronTightId[electronCounter] :
									electronGoodCutBasedIdCut  == 'M'? product.electronMediumId[electronCounter] :
									electronGoodCutBasedIdCut  == 'L'? product.electronLooseId[electronCounter] :
									electronGoodCutBasedIdCut  == 'V'? product.electronVetoId[electronCounter] : false);

			product.electronIsVeto[electronCounter] = dataReader.electronPt <= electronGoodPtCut &&
									dataReader.electronMiniIso < electronVetoIsoCut &&
									(electronVetoCutBasedIdCut == 'T'? product.electronTightId[electronCounter] :
									electronVetoCutBasedIdCut  == 'M'? product.electronMediumId[electronCounter] :
									electronVetoCutBasedIdCut  == 'L'? product.electronLooseId[electronCounter] :
									electronVetoCutBasedIdCut  == 'V'? product.electronVetoId[electronCounter] : false);

			product.electronIsAntiSelected[electronCounter] = dataReader.electronMiniIso < electronAntiIsoCut && // FIXME Check if is Medium but not Tight is correct
									(electronAntiIsCutBasedIdCut == 'T'? product.electronTightId[electronCounter] :
									electronAntiIsCutBasedIdCut  == 'M'? product.electronMediumId[electronCounter] :
									electronAntiIsCutBasedIdCut  == 'L'? product.electronLooseId[electronCounter] :
									electronAntiIsCutBasedIdCut  == 'V'? product.electronVetoId[electronCounter] : false) &&
									(electronAntiIsNotCutBasedIdCut == 'T'? !product.electronTightId[electronCounter] :
									electronAntiIsNotCutBasedIdCut  == 'M'? !product.electronMediumId[electronCounter] :
									electronAntiIsNotCutBasedIdCut  == 'L'? !product.electronLooseId[electronCounter] :
									electronAntiIsNotCutBasedIdCut  == 'V'? !product.electronVetoId[electronCounter] : false);

			if (product.electronIsGood[electronCounter]) { goodElectronCounter++;}
			if (product.electronIsVeto[electronCounter]) { vetoElectronCounter++;}
			if (product.electronIsAntiSelected[electronCounter]) { antiSelectedElectronCounter++;}
			electronCounter++;
		}

		//Overwrite number of leptons by excluding leptons that do not survive the cuts
		product.nElectron = electronCounter;
		product.nGoodElectron = goodElectronCounter;
		product.nVetoElectron = vetoElectronCounter;

		/*
		if (product.nLepton!=0 && false) { //nLepton can be 0 since unselected leptons are not counted
			ROOT::Math::PtEtaPhiMVector leadingLeptonP4 = ROOT::Math::PtEtaPhiMVector(Pt.at(0), Eta.at(0), Phi.at(0), Mass.at(0));
			for (std::size_t i = 1; i < product.nLepton; i++){
				ROOT::Math::PtEtaPhiMVector otherLeptonP4 = ROOT::Math::PtEtaPhiMVector(Pt.at(i), Eta.at(i), Phi.at(i), Mass.at(i));
				ROOT::Math::PtEtaPhiMVector diLeptonP4 = leadingLeptonP4 + otherLeptonP4;
				dileptonMass[muonCounter] = diLeptonP4.M();
			}
		}
		*/
	}
}

void ElectronProducer::EndJob(TFile &file) {
	//muonIdSFFile->Close();
	//muonIsolationSFFile->Close();
	//muonTriggerSFFile->Close();

	//electronGSFSFFile->Close();
	//electronMVASFFile->Close();
}
