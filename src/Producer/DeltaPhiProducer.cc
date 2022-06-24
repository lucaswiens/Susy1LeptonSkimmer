#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/DeltaPhiProducer.h>

DeltaPhiProducer::DeltaPhiProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree) {
	hadronicMt2Cut = configTree.get<double>("Producer.IsoTrack.Mt2.Hadronic");
	leptonicMt2Cut = configTree.get<double>("Producer.IsoTrack.Mt2.Leptonic");
}

void DeltaPhiProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	int isoTrackCounter = 0;
	bool isGoodMuon = false;
	int goodMuonIndex = 0, goodElectronIndex = 0;
	/*##########################################################
	#   Calculate the DeltaPhi using the leading good lepton   #
	#   In most cases this should be the leading lepton        #
	##########################################################*/
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
	product.absoluteDeltaPhi = -999;
	product.LP = -999;
	product.wBosonMt = -999;
	if (product.nMuon > 0 || product.nElectron > 0) {
		double leptonPt   = isGoodMuon ? product.muonPt[goodMuonIndex]     : product.electronPt[goodElectronIndex];
		double leptonEta  = isGoodMuon ? product.muonEta[goodMuonIndex]    : product.electronEta[goodElectronIndex];
		double leptonPhi  = isGoodMuon ? product.muonPhi[goodMuonIndex]    : product.electronPhi[goodElectronIndex];
		double leptonMass = isGoodMuon ? product.muonMass[goodMuonIndex]   : product.electronMass[goodElectronIndex];
		int leptonCharge  = isGoodMuon ? product.muonCharge[goodMuonIndex] : product.electronCharge[goodElectronIndex];
		ROOT::Math::PtEtaPhiMVector leptonP4 = ROOT::Math::PtEtaPhiMVector(leptonPt, leptonEta, leptonPhi, leptonMass);
		ROOT::Math::PtEtaPhiMVector metP4 = ROOT::Math::PtEtaPhiMVector(product.metPt, 0, product.metPhi, 0);
		ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + metP4;

		product.LT = leptonPt + product.metPt;
		product.deltaPhi = Utility::DeltaPhi(leptonP4.Phi(), wBosonP4.Phi());
		product.absoluteDeltaPhi = std::abs(product.deltaPhi);
		product.LP = leptonPt / wBosonP4.Pt() * std::cos(product.deltaPhi);
		product.wBosonMt = wBosonP4.Mt();

		dataReader.ReadIsoTrackEntry();
		assert(dataReader.nIsoTrack < product.nMax);
		heppy::Davismt2 mt2obj;
		product.isoTrackVeto = false;
		double deltaRMin = 0.1, isoMt2Cut = -999;
		for (int iTrack = 0; iTrack < dataReader.nIsoTrack; iTrack++) {
			dataReader.GetIsoTrackValues(iTrack);

			const int &isotrackcharge = (0 < dataReader.isoTrackPdgId) - (dataReader.isoTrackPdgId < 0); //only true for leptons but veto only matters for leptons

			if (10 < std::abs(dataReader.isoTrackPdgId) && std::abs(dataReader.isoTrackPdgId) < 14 && (isotrackcharge == leptonCharge)) continue;

			double deltaR = Utility::DeltaR(leptonEta, leptonPhi, dataReader.isoTrackEta, dataReader.isoTrackPhi);

			if (deltaR < deltaRMin) continue;

			ROOT::Math::PtEtaPhiMVector isotrackP4 = ROOT::Math::PtEtaPhiMVector(dataReader.isoTrackPt, dataReader.isoTrackEta, dataReader.isoTrackPhi, 0);

			double a[3] = {leptonP4.M(), leptonP4.X(), leptonP4.Y()};
			double b[3] = {isotrackP4.M(), isotrackP4.X(), isotrackP4.Y()};
			double c[3] = {metP4.M(), metP4.X(), metP4.Y()};

			mt2obj.set_momenta(a, b, c);
			mt2obj.set_mn(0);

			product.isoTrackMt2[isoTrackCounter] = mt2obj.get_mt2();
			product.isoTrackPt[isoTrackCounter] = dataReader.isoTrackPt;
			product.isoTrackPdgId[isoTrackCounter] = dataReader.isoTrackPdgId;

			if (10 < std::abs(dataReader.isoTrackPdgId) && std::abs(dataReader.isoTrackPdgId) < 14) { // TODO Add these to config
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

	product.HT = 0;
	for (int *nBtag : {&product.nDeepCsvLooseBTag, &product.nDeepCsvMediumBTag, &product.nDeepCsvTightBTag, &product.nDeepJetLooseBTag, &product.nDeepJetMediumBTag, &product.nDeepJetTightBTag}) {

		*nBtag = 0;
	}

	for (int iJet = 0; iJet < product.nJet; iJet++) {
		product.HT += product.jetPt.at(iJet);

		if (product.jetDeepCsvLooseId[iJet])  { product.nDeepCsvLooseBTag++;}
		if (product.jetDeepCsvMediumId[iJet]) { product.nDeepCsvMediumBTag++;}
		if (product.jetDeepCsvTightId[iJet])  { product.nDeepCsvTightBTag++;}
		if (product.jetDeepJetLooseId[iJet])  { product.nDeepJetLooseBTag++;}
		if (product.jetDeepJetMediumId[iJet]) { product.nDeepJetMediumBTag++;}
		if (product.jetDeepJetTightId[iJet])  { product.nDeepJetTightBTag++;}
	}

	// TODO this is kind of stupid.. Probably best to do this during the anaysis step
	if (product.nDeepJetMediumBTag == 0) {// 0-B SRs
		if (product.LT < 250) { product.isSignalRegion = false;}
		else if (product.LT > 250) { product.isSignalRegion = product.absoluteDeltaPhi > 0.75;}
	} else {// Multi-B SRs
		if (product.LT < 250) { product.isSignalRegion = false;}
		else if (product.LT < 350) { product.isSignalRegion = product.absoluteDeltaPhi > 1.0;}
		else if (product.LT < 600) { product.isSignalRegion = product.absoluteDeltaPhi > 0.75;}
		else if (product.LT > 600) { product.isSignalRegion = product.absoluteDeltaPhi > 0.5;}
	}
}

void DeltaPhiProducer::EndJob(TFile &file) {}
