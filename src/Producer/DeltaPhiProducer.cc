#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/DeltaPhiProducer.h>

DeltaPhiProducer::DeltaPhiProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree) {
	Name = "DeltaPhiProducer";
	deltaRCut = configTree.get<float>("Producer.IsoTrack.DeltaR");
	isoTrackPtCut = configTree.get<float>("Producer.IsoTrack.Pt");
	hadronicMt2Cut = configTree.get<float>("Producer.IsoTrack.Mt2.Hadronic");
	leptonicMt2Cut = configTree.get<float>("Producer.IsoTrack.Mt2.Leptonic");
}

void DeltaPhiProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	int isoTrackCounter = 0;
	bool isGoodMuon = false;
	int goodMuonIndex = 0, goodElectronIndex = 0;
	/*###########################################################################################
	#   Most of these variables are only well defined if the event has exactly one good lepton. #
	#   For inclusive channels they are just calculated using the leading good lepton.          #
	###########################################################################################*/
	if (product.nMuon > 0 && product.nElectron > 0) {
		for (int iMuon = 0; iMuon < product.nMuon; iMuon++) {
			if (product.muonIsGood[iMuon]){
				goodMuonIndex = iMuon;
				break;
			}
		}

		for (int iElectron = 0; iElectron < product.nElectron; iElectron++) {
			if (product.electronIsGood[iElectron]){
				goodElectronIndex = iElectron;
				break;
			}
		}
		isGoodMuon  = product.muonPt[goodMuonIndex] > product.electronPt[goodElectronIndex];
	} else if (product.nMuon > 0) {
		for (int iMuon = 0; iMuon < product.nMuon; iMuon++) {
			if (product.muonIsGood[iMuon]){
				goodMuonIndex = iMuon;
				break;
			}
		}
		isGoodMuon  = true;
	} else {
		for (int iElectron = 0; iElectron < product.nElectron; iElectron++) {
			if (product.electronIsGood[iElectron]){
				goodElectronIndex = iElectron;
				break;
			}
		}
		isGoodMuon  = false;
	}

	product.LT = -999;
	product.deltaPhi = -999;
	product.LP = -999;
	product.wBosonMt = -999;

	product.leptonPt    = -999;
	product.leptonEta   = -999;
	product.leptonPhi   = -999;
	product.leptonMass  = -999;
	product.leptonPdgId = -999;
	if (product.nLepton > 0) {
		product.leptonPt    = isGoodMuon ? product.muonPt[goodMuonIndex]    : product.electronPt[goodElectronIndex];
		product.leptonEta   = isGoodMuon ? product.muonEta[goodMuonIndex]   : product.electronEta[goodElectronIndex];
		product.leptonPhi   = isGoodMuon ? product.muonPhi[goodMuonIndex]   : product.electronPhi[goodElectronIndex];
		product.leptonMass  = isGoodMuon ? product.muonMass[goodMuonIndex]  : product.electronMass[goodElectronIndex];
		product.leptonPdgId = isGoodMuon ? product.muonPdgId[goodMuonIndex] : product.electronPdgId[goodElectronIndex];

		int leptonCharge  = isGoodMuon ? product.muonCharge[goodMuonIndex] : product.electronCharge[goodElectronIndex];
		ROOT::Math::PtEtaPhiMVector leptonP4 = ROOT::Math::PtEtaPhiMVector(product.leptonPt, product.leptonEta, product.leptonPhi, product.leptonMass);
		ROOT::Math::PtEtaPhiMVector metP4 = ROOT::Math::PtEtaPhiMVector(product.metPt, 0, product.metPhi, 0);
		ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + metP4;

		product.LT = product.leptonPt + product.metPt;
		product.deltaPhi = Utility::DeltaPhi(leptonP4.Phi(), wBosonP4.Phi());
		product.LP = product.leptonPt / wBosonP4.Pt() * std::cos(product.deltaPhi);
		product.wBosonMt = wBosonP4.Mt();

		dataReader.ReadIsoTrackEntry();
		assert(dataReader.nIsoTrack < product.nMax);
		heppy::Davismt2 mt2obj;
		product.isoTrackVeto = false;
		isoMt2Cut = -999;
		for (int iTrack = 0; iTrack < dataReader.nIsoTrack; iTrack++) {
			dataReader.GetIsoTrackValues(iTrack);
			if (10 < std::abs(dataReader.isoTrackPdgId) && std::abs(dataReader.isoTrackPdgId) < 14 && (dataReader.isoTrackCharge == leptonCharge)) continue;

			float deltaR = Utility::DeltaR(product.leptonEta, product.leptonPhi, dataReader.isoTrackEta, dataReader.isoTrackPhi);
			if (deltaR < deltaRCut && dataReader.isoTrackPt > isoTrackPtCut) continue;

			ROOT::Math::PtEtaPhiMVector isotrackP4 = ROOT::Math::PtEtaPhiMVector(dataReader.isoTrackPt, dataReader.isoTrackEta, dataReader.isoTrackPhi, 0);

			double a[3] = {leptonP4.M(), leptonP4.X(), leptonP4.Y()};
			double b[3] = {isotrackP4.M(), isotrackP4.X(), isotrackP4.Y()};
			double c[3] = {metP4.M(), metP4.X(), metP4.Y()};

			mt2obj.set_momenta(a, b, c);
			mt2obj.set_mn(0);

			product.isoTrackMt2[isoTrackCounter] = mt2obj.get_mt2();
			product.isoTrackPt[isoTrackCounter] = dataReader.isoTrackPt;
			product.isoTrackPdgId[isoTrackCounter] = dataReader.isoTrackPdgId;
			product.isoTrackCharge[isoTrackCounter] = dataReader.isoTrackCharge;

			if (10 < std::abs(dataReader.isoTrackPdgId) && std::abs(dataReader.isoTrackPdgId) < 14) { // https://twiki.cern.ch/twiki/bin/view/Main/PdgId
				product.isoTrackIsHadronicDecay[isoTrackCounter] = false; //leptonic track
				isoMt2Cut = leptonicMt2Cut;
			} else {
				product.isoTrackIsHadronicDecay[isoTrackCounter] = true; //hadronic track
				isoMt2Cut = hadronicMt2Cut;
			}

			if (product.isoTrackMt2[isoTrackCounter] <= isoMt2Cut) { product.isoTrackVeto = true;}
			isoTrackCounter++;
		}
	}
	product.nIsoTrack = isoTrackCounter;

	// count the number of btags
	for (int *nBtag : {&product.nDeepJetLooseBTag, &product.nDeepJetMediumBTag, &product.nDeepJetTightBTag}) {
		*nBtag = 0;
	}

	product.HT = 0;
	for (int iJet = 0; iJet < product.nJet; iJet++) {
		product.HT += product.jetPt.at(iJet);

		if (product.jetDeepJetLooseId[iJet])  { product.nDeepJetLooseBTag++;}
		if (product.jetDeepJetMediumId[iJet]) { product.nDeepJetMediumBTag++;}
		if (product.jetDeepJetTightId[iJet])  { product.nDeepJetTightBTag++;}
	}

	// count the number of t and W tags
	for (int *nDeepAk8Tag : {&product.nDeepAk8TopLooseId, &product.nDeepAk8TopMediumId, &product.nDeepAk8TopTightId, &product.nDeepAk8TopVeryTightId, &product.nDeepAk8TopMDLooseId, &product.nDeepAk8TopMDMediumId, &product.nDeepAk8TopMDTightId, &product.nDeepAk8TopMDVeryTightId, &product.nDeepAk8WVeryLooseId, &product.nDeepAk8WLooseId, &product.nDeepAk8WMediumId, &product.nDeepAk8WTightId, &product.nDeepAk8WMDVeryLooseId, &product.nDeepAk8WMDLooseId, &product.nDeepAk8WMDMediumId, &product.nDeepAk8WMDTightId}) {
		*nDeepAk8Tag = 0;
	}

	for (int iFatJet = 0; iFatJet < product.nFatJet; iFatJet++) {
		if (product.fatJetDeepAk8TopLooseId[iFatJet]) {product.nDeepAk8TopLooseId++;}
		if (product.fatJetDeepAk8TopMediumId[iFatJet]) {product.nDeepAk8TopMediumId++;}
		if (product.fatJetDeepAk8TopTightId[iFatJet]) {product.nDeepAk8TopTightId++;}
		if (product.fatJetDeepAk8TopVeryTightId[iFatJet]) {product.nDeepAk8TopVeryTightId++;}

		if (product.fatJetDeepAk8TopMDLooseId[iFatJet]) {product.nDeepAk8TopMDLooseId++;}
		if (product.fatJetDeepAk8TopMDMediumId[iFatJet]) {product.nDeepAk8TopMDMediumId++;}
		if (product.fatJetDeepAk8TopMDTightId[iFatJet]) {product.nDeepAk8TopMDTightId++;}
		if (product.fatJetDeepAk8TopMDVeryTightId[iFatJet]) {product.nDeepAk8TopMDVeryTightId++;}

		if (product.fatJetDeepAk8WVeryLooseId[iFatJet]) {product.nDeepAk8WVeryLooseId++;}
		if (product.fatJetDeepAk8WLooseId[iFatJet]) {product.nDeepAk8WLooseId++;}
		if (product.fatJetDeepAk8WMediumId[iFatJet]) {product.nDeepAk8WMediumId++;}
		if (product.fatJetDeepAk8WTightId[iFatJet]) {product.nDeepAk8WTightId++;}

		if (product.fatJetDeepAk8WMDVeryLooseId[iFatJet]) {product.nDeepAk8WMDVeryLooseId++;}
		if (product.fatJetDeepAk8WMDLooseId[iFatJet]) {product.nDeepAk8WMDLooseId++;}
		if (product.fatJetDeepAk8WMDMediumId[iFatJet]) {product.nDeepAk8WMDMediumId++;}
		if (product.fatJetDeepAk8WMDTightId[iFatJet]) {product.nDeepAk8WMDTightId++;}
	}
}

void DeltaPhiProducer::EndJob(TFile &file) {}
