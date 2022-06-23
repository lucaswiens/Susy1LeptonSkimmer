#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/MuonProducer.h>

MuonProducer::MuonProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {
	std::string cmsswBase = std::getenv("CMSSW_BASE");
	//Name = "MuonProducer";
	rc.init(std::string(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Muon.Scale." + eraSelector)));

	muonGoodPtCut         = configTree.get<double>("Producer.Muon.Pt.Good");
	muonVetoPtCut         = configTree.get<double>("Producer.Muon.Pt.Veto");
	muonEtaCut            = configTree.get<double>("Producer.Muon.Eta");
	muonGoodIsoCut        = configTree.get<double>("Producer.Muon.Iso.Good");
	muonVetoIsoCut        = configTree.get<double>("Producer.Muon.Iso.Veto");
	muonAntiIsoCut        = configTree.get<double>("Producer.Muon.Iso.Anti");
	muonGoodCutBasedIdCut = configTree.get<char>("Producer.Muon.CutBasedId.Good");
	muonVetoCutBasedIdCut = configTree.get<char>("Producer.Muon.CutBasedId.Veto");
	muonAntiCutBasedIdCut = configTree.get<char>("Producer.Muon.CutBasedId.Anti");
	//muonDxyCut          = configTree.get<double>("Producer.Muon.Dxy");
	//muonDzCut           = configTree.get<double>("Producer.Muon.Dz");
	//muonSip3dCut        = configTree.get<double>("Producer.Muon.Sip3DCut");

	std::cout << std::endl <<
		"The following cuts are applied to Muons:"   << std::endl <<
		"|Eta|          < " << muonEtaCut          << std::endl <<
		"GoodPt         > " << muonGoodPtCut         << std::endl <<
		"GoodIso        < " << muonGoodIsoCut        << std::endl <<
		"GoodCutBasedId = " << muonGoodCutBasedIdCut << std::endl <<
		"VetoPt         > " << muonVetoPtCut         << std::endl <<
		"VetoIso        < " << muonVetoIsoCut        << std::endl <<
		"VetoCutBasedId = " << muonVetoCutBasedIdCut << std::endl <<
		"AntiIso       >= " << muonAntiIsoCut        << std::endl <<
		"AntiCutBasedId = " << muonAntiCutBasedIdCut << std::endl << std::endl;
}

void MuonProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	//Initialize all variables as -999
	product.muonPtVector.clear(); // Just for testing, this is required for vectors
	//product.nMuon = 0;
	//product.nGoodMuon = 0;
	//product.nVetoMuon = 0;

	dataReader.ReadMuonEntry();
	product.nMuon = dataReader.nMuon;
	assert(product.nMuon < product.nMax);

	int muonCounter = 0, goodMuonCounter = 0, vetoMuonCounter = 0, antiSelectedMuonCounter = 0;
	for (int iMuon = 0; iMuon < dataReader.nMuon; iMuon++) {
		dataReader.GetMuonValues(iMuon);

		// Rochester Correction
		int muonMatchedGenIndex = -999;
		double dataScaleFactor = 1., mcScaleFactor = 1., scaleFactorUnc = 0.; // TODO think about how to do this better
		if (!product.GetIsData()) {
			muonMatchedGenIndex = dataReader.GetGenMatchedIndex(dataReader.muonPt, dataReader.muonPhi, dataReader.muonEta, 13, 0.4, 0.4);
			if(muonMatchedGenIndex > 0){
				dataReader.alreadyMatchedIndex.push_back(muonMatchedGenIndex);
				dataReader.GetGenValues(muonMatchedGenIndex);

				mcScaleFactor = rc.kSpreadMC(dataReader.muonCharge, dataReader.muonPt, dataReader.muonEta, dataReader.muonPhi, dataReader.genPt, 0, 0);
				scaleFactorUnc = rc.kSpreadMCerror(dataReader.muonCharge, dataReader.muonPt, dataReader.muonEta, dataReader.muonPhi, dataReader.genPt);
			} else {

				mcScaleFactor = rc.kSmearMC(dataReader.muonCharge, dataReader.muonPt, dataReader.muonEta, dataReader.muonPhi, dataReader.muonNTrackerLayers, dataReader.muonRandomNumber, 0, 0);
				scaleFactorUnc = rc.kSmearMCerror(dataReader.muonCharge, dataReader.muonPt, dataReader.muonEta, dataReader.muonPhi, dataReader.muonNTrackerLayers, dataReader.muonRandomNumber);
			}
		} else {
			dataScaleFactor = rc.kScaleDT(dataReader.muonCharge, dataReader.muonPt, dataReader.muonEta, dataReader.muonPhi, 0, 0);
		}

		double muonPt      = dataReader.muonPt * (product.GetIsData() ? dataScaleFactor : mcScaleFactor),
			muonPtUp   = dataReader.muonPt * (product.GetIsData() ? (dataScaleFactor + scaleFactorUnc) : (mcScaleFactor + scaleFactorUnc)),
			muonPtDown = dataReader.muonPt * (product.GetIsData() ? (dataScaleFactor - scaleFactorUnc) : (mcScaleFactor - scaleFactorUnc));

		if (muonPt < muonVetoPtCut ||
				std::abs(dataReader.muonEta) > muonEtaCut ||
				//dataReader.muonDxy > muonDxyCut ||
				//dataReader.muonDz > muonDzCut ||
				//dataReader.muonSip3d > muonSip3dCut ||
				//dataReader.muonMiniIso > muonVetoIsoCut ||
				!dataReader.muonIsPfCand)
		{continue; }


		product.muonPtVector.push_back(muonPt);
		product.muonPt[muonCounter] = muonPt;
		product.muonEta[muonCounter] = dataReader.muonEta;
		product.muonPhi[muonCounter] = dataReader.muonPhi;
		product.muonMass[muonCounter] = dataReader.muonMass;
		product.muonCharge[muonCounter] = dataReader.muonCharge;
		product.muonMiniIso[muonCounter] = dataReader.muonMiniIso;
		product.muonPdgId[muonCounter] = dataReader.muonPdgId;
		product.muonGenMatchedIndex[muonCounter] = muonMatchedGenIndex;

		product.muonLooseId[muonCounter]    = dataReader.muonLooseId;
		product.muonMediumId[muonCounter]   = dataReader.muonMediumId;
		product.muonLooseId[muonCounter]    = dataReader.muonLooseId;
		product.muonCutBasedId[muonCounter] = dataReader.muonTightId ? 4 : dataReader.muonMediumId ? 3 : dataReader.muonLooseId ? 2 : 1;

		product.muonIsGood[muonCounter] = muonPt > muonGoodPtCut &&
							dataReader.muonMiniIso < muonGoodIsoCut &&
							dataReader.muonIdMap.at(muonGoodCutBasedIdCut);
							//(muonGoodCutBasedIdCut == 'T'?  dataReader.muonTightId :
							//muonGoodCutBasedIdCut  == 'M'? dataReader.muonMediumId :
							//muonGoodCutBasedIdCut  == 'L'?  dataReader.muonLooseId : false);

		product.muonIsVeto[muonCounter] = muonPt <= muonGoodPtCut &&
							dataReader.muonMiniIso < muonVetoIsoCut &&
							dataReader.muonIdMap.at(muonVetoCutBasedIdCut);
							//(muonVetoCutBasedIdCut == 'T'?  dataReader.muonTightId :
							//muonVetoCutBasedIdCut  == 'M'? dataReader.muonMediumId :
							//muonVetoCutBasedIdCut  == 'L'?  dataReader.muonLooseId : false);

		product.muonIsAntiSelected[muonCounter] = dataReader.muonMiniIso >= muonAntiIsoCut &&
							dataReader.muonIdMap.at(muonAntiCutBasedIdCut);
							//(muonAntiCutBasedIdCut == 'T'?  dataReader.muonTightId :
							//muonAntiCutBasedIdCut  == 'M'? dataReader.muonMediumId :
							//muonAntiCutBasedIdCut  == 'L'?  dataReader.muonLooseId : false);

		if (product.muonIsGood[muonCounter]) { goodMuonCounter++;}
		if (product.muonIsVeto[muonCounter]) { vetoMuonCounter++;}
		if (product.muonIsAntiSelected[muonCounter]) { antiSelectedMuonCounter++;}
		muonCounter++;
	}

	//Overwrite number of leptons by excluding leptons that do not survive the cuts
	product.nMuon = muonCounter;
	product.nGoodMuon = goodMuonCounter;
	product.nVetoMuon = vetoMuonCounter;
	product.nAntiSelectedMuon = vetoMuonCounter;

	/*
	if (product.nLepton!=0 && false) { //nLepton can be 0 since unselected leptons are not counted
		ROOT::Math::PtEtaPhiMVector leadingLeptonP4 = ROOT::Math::PtEtaPhiMVector(Pt.at(0), Eta.at(0), Phi.at(0), Mass.at(0));
		for (int i = 1; i < product.nLepton; i++){
			ROOT::Math::PtEtaPhiMVector otherLeptonP4 = ROOT::Math::PtEtaPhiMVector(Pt.at(i), Eta.at(i), Phi.at(i), Mass.at(i));
			ROOT::Math::PtEtaPhiMVector diLeptonP4 = leadingLeptonP4 + otherLeptonP4;
			dileptonMass[muonCounter] = diLeptonP4.M();
		}
	}
	*/

}

void MuonProducer::EndJob(TFile &file) {
	//muonIdSFFile->Close();
	//muonIsolationSFFile->Close();
	//muonTriggerSFFile->Close();

	//electronGSFSFFile->Close();
	//electronMVASFFile->Close();
}
