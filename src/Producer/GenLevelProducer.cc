#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/GenLevelProducer.h>

GenLevelProducer::GenLevelProducer(const pt::ptree &configTree, const pt::ptree &scaleFactorTree, std::string eraSelector) {

}

void GenLevelProducer::Produce(DataReader &dataReader, Susy1LeptonProduct &product) {
	/*
	ROOT::Math::PtEtaPhiMVector leptonP4;
	ROOT::Math::PtEtaPhiMVector neutrinoP4;
	for (int i = 0; i < NGenPart; i++){
		const double &leptonPt = -999;//genPartPt->At(i);
		if (leptonPt > maxLeptonPt) {maxLeptonPt = leptonPt;}

		if (dataReader.genMotherIndex >= 0) { //store daughter indices for nISR
			nDaughters++;
			idxDaughter.push_back(i);
			//GenTauMotherId.push_back(motherId);
		}

		if (abs(dataReader.genPdgId) == 14 || abs(dataReader.genPdgId) == 12) {
			NGenNeutrino++;
			GenNeutrinoPt.push_back(-999);
		}//genPartPt->At(i));
		if ((abs(dataReader.genPdgId) == 13 || abs(dataReader.genPdgId) == 11 ) && dataReader.genMotherIndex >= 0) { // 13:Muon, 11:Electron
			if (dataReader.StatusFlag == 2) { // 2 : isTauDecayProduct
				NGenLeptonFromTau++;
				NGenTau++;
				idxTauLepton.push_back(i);
			} else {
				NGenLepton++;
				//const int &motherStatus = genPartStatus->At(dataReader.genMotherIndex);
				const int &motherId = -999;//genPartPdgId->At(dataReader.genMotherIndex);
				const int &idxGrandMother = -999;//genPartIdxMother->At(dataReader.genMotherIndex);
				GenLepMotherId.push_back(motherId);
				if(idxGrandMother > 0){
					const int &grandMotherId = -999;//genPartPdgId->At(idxGrandMother);
					GenLepGrandMotherId.push_back(grandMotherId);
				}

				if (abs(motherId) == 24 && dataReader.genStatus == 23) { // genLep is outgoing and has W as mother
					const double &motherPt = -999;//genPartPt->At(dataReader.genMotherIndex);
					const double &motherPhi = -999;//genPartPhi->At(dataReader.genMotherIndex);
					const double &motherEta = -999;//genPartEta->At(dataReader.genMotherIndex);
					const double &motherMass = -999;//genPartMass->At(dataReader.genMotherIndex);

					ROOT::Math::PtEtaPhiMVector motherP4(motherPt, motherEta, motherPhi, motherMass);
					product.genWBosonPt.push_back(motherPt);
					product.genWBosonEta.push_back(motherEta);
					product.genWBosonPhi.push_back(motherPhi);
					product.genWBosonMass.push_back(motherMass);
					GenWDirectMass.push_back(motherMass);

					if(idxGrandMother > 0 && abs(GenLepGrandMotherId.back()) == 6){
						//Used for sytematic variations, store in product for now
						const double &grandMotherPt = -999;//genPartPt->At(idxGrandMother);
						const double &grandMotherPhi = -999;//genPartPhi->At(idxGrandMother);
						const double &grandMotherEta = -999;//genPartEta->At(idxGrandMother);
						const double &grandMotherMass = -999;//genPartMass->At(idxGrandMother);
						product.genTopPt.push_back(grandMotherPt);
						product.genTopEta.push_back(grandMotherEta);
						product.genTopPhi.push_back(grandMotherPhi);
						product.genTopMass.push_back(grandMotherMass);
					}

					for (int j = i+1; j < NGenPart; j++){
						const int &idxNeutrinoMother = -999;//genPartIdxMother->At(j);
						const int &statusNeutrino = -999;//genPartStatus->At(j);
						if (idxNeutrinoMother == dataReader.genMotherIndex && statusNeutrino == 23) {
							NGenMatchedW++;
							//leptonPt already read
							const double &leptonPhi = -999;//genPartPhi->At(i);
							const double &leptonEta = -999;//genPartEta->At(i);
							const double &leptonMass = -999;//genPartMass->At(i);

							const double &neutrinoPt = -999;//genPartPt->At(j);
							const double &neutrinoPhi = -999;//genPartPhi->At(j);
							const double &neutrinoEta = -999;//genPartEta->At(j);
							const double &neutrinoMass = -999;//genPartMass->At(j);

							leptonP4 = ROOT::Math::PtEtaPhiMVector(leptonPt, leptonEta, leptonPhi, leptonMass);
							neutrinoP4 = ROOT::Math::PtEtaPhiMVector(neutrinoPt, neutrinoEta, neutrinoPhi, neutrinoMass);
							ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + neutrinoP4;

							GenDeltaPhiLepWSum.push_back(Utility::Utility::DeltaPhi(leptonP4, wBosonP4));
							GenWSumMass.push_back(wBosonP4.M());
							GenDeltaPhiLepWDirect.push_back(Utility::Utility::DeltaPhi(leptonP4, motherP4));
							GenMTLepNu.push_back(wBosonP4.Mt());
						}
					}
				}
			}
		}
		product.nGenWBoson = product.genWBosonPt.size();
		product.nGenTop = product.genTopPt.size();

		if (NGenLepton + NGenLeptonFromTau == 2) {
			LeptonsInAcceptance = true;
			isDiLeptonEvent = true;
			isHadTauEvent = NGenTau > NGenLeptonFromTau;
			if (isHadTauEvent) { LeptonDecayChannelFlag =1;}
			for (int idx : idxTauLepton) {
				const double &pt = -999;//genPartPt->At(idx);
				//const double &phi = -999;//genPartPhi->At(idx);
				const double &eta = -999;//genPartEta->At(idx);
				//const double &mass = -999;//genPartMass->At(idx);
				const double &dataReader.genPdgId = -999;//genPartPdgId->At(idx);

				if (maxLeptonPt >= 25 && pt < 10) { LeptonsInAcceptance = false;}
				if (maxLeptonPt < 25 && pt < 5) { LeptonsInAcceptance = false;}
				if (abs(eta) > 2.5) { LeptonsInAcceptance = false;}
				if (abs(eta) > 2.5) { LeptonsInAcceptance = false;}
				if (dataReader.genPdgId == 11 && abs(eta) >= 1.44 && abs(eta) < 1.57) { LeptonsInAcceptance = false;}

				if (isHadTauEvent && !LeptonsInAcceptance) { LeptonDecayChannelFlag = 0;}
				else if (isHadTauEvent) { LeptonDecayChannelFlag = 1;}
				else if (!LeptonsInAcceptance) { LeptonDecayChannelFlag = 2;}
				else { LeptonDecayChannelFlag = 3;}
			}
		}

	}
	//product.nGenPart = NGenPart;
	*/
}

void GenLevelProducer::EndJob(TFile &file) {}
