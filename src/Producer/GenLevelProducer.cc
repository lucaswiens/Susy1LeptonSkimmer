#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/GenLevelProducer.h>

GenLevelProducer::GenLevelProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {

}

void GenLevelProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	ROOT::Math::PtEtaPhiMVector leptonP4;
	ROOT::Math::PtEtaPhiMVector neutrinoP4;
	double maxLeptonPt;
	int genPartCounter = 0, genLeptonCounter = 0, genTauCounter = 0, genLeptonFromTauCounter = 0, genWCounter = 0, genMatchedWCounter = 0, genNeutrinoCounter = 0, genTopCounter = 0;
	std::vector<int> tauLeptonIndices;
	dataReader.ReadGenEntry();
	for (int iGen = 0; iGen < dataReader.nGenPart; iGen++) {
		dataReader.GetGenValues(iGen);
		const double leptonPt = dataReader.genPt;


		/*#################################################################################################
		#   Meaning of genStatus                                                                          #
		#   https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc102X_doc.html#GenPart   #
		#################################################################################################*/
		if (std::abs(dataReader.genPdgId) == 14 || std::abs(dataReader.genPdgId) == 12) {
			product.genNeutrinoPt[genNeutrinoCounter] = dataReader.genPt;
			genNeutrinoCounter++;
		}

		if ((std::abs(dataReader.genPdgId) == 13 || std::abs(dataReader.genPdgId) == 11 ) && dataReader.genMotherIndex >= 0) { // 13:Muon, 11:Electron and mother exists
			if (leptonPt > maxLeptonPt) { maxLeptonPt = leptonPt;}
			if (dataReader.genStatus == 2) { // 2 : isTauDecayProduct
				tauLeptonIndices.push_back(iGen);
				const int motherIndex = dataReader.genMotherIndex;
				dataReader.GetGenValues(motherIndex); // Read Mother Info
				const int motherPdgId = dataReader.genPdgId;
				const int grandMotherIndex = dataReader.genMotherIndex;
				dataReader.GetGenValues(motherIndex); // Read GrandMother Info
				const int grandMotherPdgId = dataReader.genPdgId;
				product.genTauMotherPdgId[genTauCounter] = motherPdgId;
				product.genTauGrandMotherPdgId[genTauCounter] = grandMotherPdgId;
				genLeptonFromTauCounter++;
				genTauCounter++;
			} else {
				// lepton Pt already read
				const double leptonPhi  = dataReader.genPhi;
				const double leptonEta  = dataReader.genEta;
				const double leptonMass = dataReader.genMass;
				const int leptonStatus = dataReader.genStatus;

				const int motherIndex = dataReader.genMotherIndex;
				dataReader.GetGenValues(motherIndex); // Read Mother Info
				const double motherPt   = dataReader.genPt;
				const double motherPhi  = dataReader.genPhi;
				const double motherEta  = dataReader.genEta;
				const double motherMass = dataReader.genMass;
				const int motherStatus  = dataReader.genStatus;
				const int motherPdgId   = dataReader.genPdgId;
				product.genLepMotherPdgId[genLeptonCounter] = motherPdgId;

				const int grandMotherIndex = dataReader.genMotherIndex;
				dataReader.GetGenValues(grandMotherIndex); // Read GrandMother Info
				const double grandMotherPt   = dataReader.genPt;
				const double grandMotherPhi  = dataReader.genPhi;
				const double grandMotherEta  = dataReader.genEta;
				const double grandMotherMass = dataReader.genMass;
				const int grandMotherStatus  = dataReader.genStatus;
				const int grandMotherPdgId   = dataReader.genPdgId;
				product.genLepGrandMotherPdgId[genLeptonCounter] = grandMotherPdgId;

				if (std::abs(motherPdgId) == 24 && leptonStatus == 23) { // genLep is outgoing and has W as mother
					ROOT::Math::PtEtaPhiMVector motherP4(motherPt, motherEta, motherPhi, motherMass);
					product.genWDirectMass[genWCounter] = motherMass;

					if (grandMotherIndex > 0 && std::abs(grandMotherPdgId) == 6) {
						genTopCounter++;
					}

					for (int iGen2 = iGen + 1; iGen2 < dataReader.nGenPart; iGen2++) {
						dataReader.GetGenValues(iGen2);
						const int neutrinoMotherIndex = dataReader.genMotherIndex;
						const int statusNeutrino = dataReader.genStatus;
						if (neutrinoMotherIndex == motherIndex && statusNeutrino == 23) { // Has the same mother as the lepton and the neutrino status
							const double neutrinoPt   = dataReader.genPt;
							const double neutrinoPhi  = dataReader.genPhi;
							const double neutrinoEta  = dataReader.genEta;
							const double neutrinoMass = dataReader.genMass;

							leptonP4 = ROOT::Math::PtEtaPhiMVector(leptonPt, leptonEta, leptonPhi, leptonMass);
							neutrinoP4 = ROOT::Math::PtEtaPhiMVector(neutrinoPt, neutrinoEta, neutrinoPhi, neutrinoMass);
							ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + neutrinoP4;

							product.genDeltaPhiLepWSum[genMatchedWCounter] = Utility::DeltaPhi(leptonP4.Phi(), wBosonP4.Phi());
							product.genWSumMass[genMatchedWCounter] = wBosonP4.M();
							product.genDeltaPhiLepWDirect[genMatchedWCounter] = Utility::DeltaPhi(leptonP4.Phi(), motherP4.Phi());
							product.genMTLepNu[genMatchedWCounter] = wBosonP4.Mt();

							genMatchedWCounter++;
						} // if (idxNeutrinoMother == dataReader.genMotherIndex && statusNeutrino == 23)
					} // for (int iGen2 = i+1; iGen2 < dataReader.nGenPart; iGen2++)
					genWCounter++;
				} // if (std::abs(motherPdgId) == 24 && dataReader.genStatus == 23)
				genLeptonCounter++;
			} // if (dataReader.genStatus != 2), i.e. is not a tau decay product
		} // if is a muon or electron
	}

	// TODO figure out if this is used at all
	if (genLeptonCounter + genLeptonFromTauCounter == 2) {
		product.leptonsInAcceptance = true;
		product.isDiLeptonEvent = true;

		product.isHadTauEvent = genTauCounter > genLeptonFromTauCounter;
		if (product.isHadTauEvent) { product.leptonDecayChannelFlag =1;}

		for (int iTau : tauLeptonIndices) {
			dataReader.GetGenValues(iTau);
			const double &tauPt   = dataReader.genPt;
			const double &tauPhi  = dataReader.genPhi;
			const double &tauEta  = dataReader.genEta;
			const double &tauMass = dataReader.genMass;
			const double &tauPdgId = dataReader.genPdgId;

			if (maxLeptonPt >= 25 && tauPt < 10) { product.leptonsInAcceptance = false;}
			if (maxLeptonPt < 25 && tauPt < 5) { product.leptonsInAcceptance = false;}
			if (std::abs(tauEta) > 2.4) { product.leptonsInAcceptance = false;}
			if (std::abs(tauEta) > 2.4) { product.leptonsInAcceptance = false;}
			if (tauPdgId == 11 &&
					1.44 <= std::abs(tauEta) &&
					std::abs(tauEta) < 1.57
			) {
				product.leptonsInAcceptance = false;
			}

			if (product.isHadTauEvent && !product.leptonsInAcceptance) { product.leptonDecayChannelFlag = 0;}
			else if (product.isHadTauEvent) { product.leptonDecayChannelFlag = 1;}
			else if (!product.leptonsInAcceptance) { product.leptonDecayChannelFlag = 2;}
			else { product.leptonDecayChannelFlag = 3;}
		}
	} else {
		product.leptonDecayChannelFlag = -999;
	}

	product.nGenLepton = genLeptonCounter;
	product.nGenTau = genTauCounter;
	product.nGenLeptonFromTau = genLeptonFromTauCounter;
	product.nGenMatchedW = genMatchedWCounter;
	product.nGenNeutrino = genNeutrinoCounter;
}

void GenLevelProducer::EndJob(TFile &file) {}
