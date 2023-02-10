#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/MuonProducer.h>

MuonProducer::MuonProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {
	Name = "MuonProducer";
	std::string cmsswBase = std::getenv("CMSSW_BASE");
	rc.init(std::string(cmsswBase + "/src/" + scaleFactorTree.get<std::string>("Muon." + eraSelector + ".RoccoR")));

	muonGoodPtCut         = configTree.get<float>("Producer.Muon.Pt.Good");
	muonVetoPtCut         = configTree.get<float>("Producer.Muon.Pt.Veto");
	muonEtaCut            = configTree.get<float>("Producer.Muon.Eta");
	muonGoodIsoCut        = configTree.get<float>("Producer.Muon.Iso.Good");
	muonVetoIsoCut        = configTree.get<float>("Producer.Muon.Iso.Veto");
	muonAntiIsoCut        = configTree.get<float>("Producer.Muon.Iso.Anti");
	muonGoodCutBasedIdCut = configTree.get<char>("Producer.Muon.CutBasedId.Good");
	muonVetoCutBasedIdCut = configTree.get<char>("Producer.Muon.CutBasedId.Veto");
	muonAntiCutBasedIdCut = configTree.get<char>("Producer.Muon.CutBasedId.Anti");
	//muonDxyCut          = configTree.get<float>("Producer.Muon.Dxy");
	//muonDzCut           = configTree.get<float>("Producer.Muon.Dz");
	//muonSip3dCut        = configTree.get<float>("Producer.Muon.Sip3DCut");

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
	dataReader.ReadMuonEntry();
	product.nMuon = dataReader.nMuon;
	assert(product.nMuon < product.nMax);

	int muonCounter = 0, goodMuonCounter = 0, vetoMuonCounter = 0, antiSelectedMuonCounter = 0;
	for (int iMuon = 0; iMuon < dataReader.nMuon; iMuon++) {
		dataReader.GetMuonValues(iMuon);

		// Rochester Correction
		int muonMatchedGenIndex = -999;
		float dataScaleFactor = 1., mcScaleFactor = 1., scaleFactorUnc = 0.;
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

		float muonPt      = dataReader.muonPt * (product.GetIsData() ? dataScaleFactor : mcScaleFactor),
			muonPtUp   = dataReader.muonPt * (product.GetIsData() ? (dataScaleFactor + scaleFactorUnc) : (mcScaleFactor + scaleFactorUnc)),
			muonPtDown = dataReader.muonPt * (product.GetIsData() ? (dataScaleFactor - scaleFactorUnc) : (mcScaleFactor - scaleFactorUnc));

		if (muonPt < muonVetoPtCut ||
				std::abs(dataReader.muonEta) > muonEtaCut ||
				!dataReader.muonIsPfCand)
		{ continue;}

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

		product.muonIsVeto[muonCounter] = muonPt <= muonGoodPtCut &&
							muonPt >= muonVetoPtCut &&
							dataReader.muonMiniIso < muonVetoIsoCut &&
							dataReader.muonIdMap.at(muonVetoCutBasedIdCut);

		product.muonIsAntiSelected[muonCounter] = dataReader.muonMiniIso >= muonAntiIsoCut &&
							dataReader.muonIdMap.at(muonAntiCutBasedIdCut);

		if (product.muonIsGood[muonCounter]) { goodMuonCounter++;}
		if (product.muonIsVeto[muonCounter]) { vetoMuonCounter++;}
		if (product.muonIsAntiSelected[muonCounter]) { antiSelectedMuonCounter++;}
		muonCounter++;
	}

	product.nMuon = muonCounter;
	product.nGoodMuon = goodMuonCounter;
	product.nVetoMuon = vetoMuonCounter;
	product.nAntiSelectedMuon = vetoMuonCounter;
}

void MuonProducer::EndJob(TFile &file) { }
