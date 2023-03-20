#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/GenLevelProducer.h>

GenLevelProducer::GenLevelProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {
	Name = "GenLevelProducer";

	jetPtCut  = configTree.get<float>("Producer.Jet.Pt");
	jetEtaCut = configTree.get<float>("Producer.Jet.Eta");
}

void GenLevelProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	product.genDeltaPhi = -999;
	product.genLP = -999;
	product.genMuonPt = -999;
	product.genMuonEta = -999;
	product.genMuonPhi = -999;
	product.genMuonMass = -999;
	product.genElectronPt = -999;
	product.genElectronEta =-999;
	product.genElectronPhi = -999;
	product.genElectronMass = -999;

	ROOT::Math::PtEtaPhiMVector leptonP4;
	ROOT::Math::PtEtaPhiMVector neutrinoP4;
	float maxLeptonPt;
	int genPartCounter = 0, genLeptonCounter = 0, genTauCounter = 0, genLeptonFromTauCounter = 0, genWCounter = 0, genMatchedWCounter = 0, genNeutrinoCounter = 0, genTopCounter = 0;
	std::vector<int> tauLeptonIndices;

	dataReader.ReadGenEntry();
	product.genMetPt = dataReader.genMetPt;
	product.genMetPhi = dataReader.genMetPhi;
	for (int iGen = 0; iGen < dataReader.nGenPart; iGen++) {
		dataReader.GetGenValues(iGen);
		const float leptonPt = dataReader.genPt;


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
				const float leptonPhi  = dataReader.genPhi;
				const float leptonEta  = dataReader.genEta;
				const float leptonMass = dataReader.genMass;
				const int leptonStatus = dataReader.genStatus;
				if (std::abs(dataReader.genPdgId) == 13) {
					product.genMuonPt   = leptonPt;
					product.genMuonEta  = leptonEta;
					product.genMuonPhi  = leptonPhi;
					product.genMuonMass = leptonMass;
				} else {
					product.genElectronPt   = leptonPt;
					product.genElectronEta  = leptonEta;
					product.genElectronPhi  = leptonPhi;
					product.genElectronMass = leptonMass;
				}

				const int motherIndex = dataReader.genMotherIndex;
				dataReader.GetGenValues(motherIndex); // Read Mother Info
				const float motherPt   = dataReader.genPt;
				const float motherPhi  = dataReader.genPhi;
				const float motherEta  = dataReader.genEta;
				const float motherMass = dataReader.genMass;
				const int motherStatus = dataReader.genStatus;
				const int motherPdgId  = dataReader.genPdgId;
				product.genLepMotherPdgId[genLeptonCounter] = motherPdgId;

				const int grandMotherIndex = dataReader.genMotherIndex;
				dataReader.GetGenValues(grandMotherIndex); // Read GrandMother Info
				const float grandMotherPt   = dataReader.genPt;
				const float grandMotherPhi  = dataReader.genPhi;
				const float grandMotherEta  = dataReader.genEta;
				const float grandMotherMass = dataReader.genMass;
				const int grandMotherStatus = dataReader.genStatus;
				const int grandMotherPdgId  = dataReader.genPdgId;
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
							const float neutrinoPt   = dataReader.genPt;
							const float neutrinoPhi  = dataReader.genPhi;
							const float neutrinoEta  = dataReader.genEta;
							const float neutrinoMass = dataReader.genMass;

							leptonP4 = ROOT::Math::PtEtaPhiMVector(leptonPt, leptonEta, leptonPhi, leptonMass);
							neutrinoP4 = ROOT::Math::PtEtaPhiMVector(neutrinoPt, neutrinoEta, neutrinoPhi, neutrinoMass);
							ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + neutrinoP4;

							product.genDeltaPhiLepWSum[genMatchedWCounter] = Utility::DeltaPhi(leptonP4.Phi(), wBosonP4.Phi());
							product.genWSumMass[genMatchedWCounter] = wBosonP4.M();
							product.genDeltaPhiLepWDirect[genMatchedWCounter] = Utility::DeltaPhi(leptonP4.Phi(), motherP4.Phi());
							product.genMTLepNu[genMatchedWCounter] = wBosonP4.Mt();

							product.genLT = leptonPt + product.genMetPt;
							product.genDeltaPhi = Utility::DeltaPhi(leptonP4.Phi(), wBosonP4.Phi());
							product.genLP = leptonPt / wBosonP4.Pt() * std::cos(product.deltaPhi);
							product.genWBosonMt = wBosonP4.Mt();

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
			const float &tauPt   = dataReader.genPt;
			const float &tauPhi  = dataReader.genPhi;
			const float &tauEta  = dataReader.genEta;
			const float &tauMass = dataReader.genMass;
			const float &tauPdgId = dataReader.genPdgId;

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
	product.genWeight = dataReader.genWeight;

	dataReader.ReadGenJetEntry();
	product.genHT = 0;
	for (int iGen = 0; iGen < dataReader.nGenJet; iGen++) {
		dataReader.GetGenJetValues(iGen);
		if (dataReader.genJetPt < jetPtCut || std::abs(dataReader.genJetEta) > jetEtaCut) { continue;}
		product.genHT += dataReader.genJetPt;
	}
}

void GenLevelProducer::EndJob(TFile &file) {}
