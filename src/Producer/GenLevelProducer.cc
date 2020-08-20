#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/GenLevelProducer.h>

GenLevelProducer::GenLevelProducer(const int& era, TTreeReader& reader):
	BaseProducer(&reader),
	era(era)
	{}

float GenLevelProducer::DeltaPhi(ROOT::Math::PtEtaPhiMVector v1, ROOT::Math::PtEtaPhiMVector v2) {
	float dPhi = v1.Phi() - v2.Phi();
	while (dPhi >= TMath::Pi()) dPhi -= TMath::TwoPi();
	while (dPhi < -TMath::Pi()) dPhi += TMath::TwoPi();
	return dPhi;
}
void GenLevelProducer::BeginJob(TTree* tree, bool &isData) {
	//Set data bool
	this->isData = isData;

	nGenPart = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "nGenPart");
	genPartPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_pt");
	genPartPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_phi");
	genPartEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_eta");
	genPartMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_mass");
	genPartPdgId = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_pdgId");
	genPartIdxMother = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_genPartIdxMother");
	genPartStatus = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_status");
	genPartStatusFlag = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_statusFlags");

	//Set Branches of output tree
	//tree->Branch("", &);
	tree->Branch("GenDeltaPhiLepWSum", &GenDeltaPhiLepWSum);
	tree->Branch("GenDeltaPhiLepWDirect", &GenDeltaPhiLepWDirect);
	tree->Branch("GenWSumMass", &GenWSumMass);
	tree->Branch("GenWDirectMass", &GenWDirectMass);
	tree->Branch("nGenMatchedW", &NGenMatchedW);
	tree->Branch("GenMTLepNu", &GenMTLepNu);
	tree->Branch("LeptonDecayChannelFlag", &LeptonDecayChannelFlag);
	tree->Branch("GenTauGrandMotherId", &GenTauGrandMotherId);
	tree->Branch("GenTauMotherId", &GenTauMotherId);
	tree->Branch("GenLepGrandMotherId", &GenLepGrandMotherId);
	tree->Branch("GenLepMotherId", &GenLepMotherId);
	tree->Branch("GenLeptonsInAcceptance", &LeptonsInAcceptance);


	//Set TTreeReader for genpart and trigger obj from BaseProducer
	SetCollection(this->isData);
}

void GenLevelProducer::Produce(CutFlow& cutflow, Susy1LeptonProduct *product) {
	//Initialize all variables as -999
	GenDeltaPhiLepWSum.clear();
	GenDeltaPhiLepWDirect.clear();
	GenWSumMass.clear();
	GenWDirectMass.clear();
	GenMTLepNu.clear();
	GenTauGrandMotherId.clear();
	GenTauMotherId.clear();
	GenLepGrandMotherId.clear();
	GenLepMotherId.clear();

	LeptonDecayChannelFlag = -999;
	NGenPart = -999;
	NGenLepton = 0;
	NGenTau = 0;
	NGenLeptonFromTau = 0;
	NGenMatchedW = 0;

	isHadTauEvent = false;
	LeptonsInAcceptance = false;
	//std::unique_ptr<TTreeReaderArray<float>> nGenPart, GenPartEta, GenPartMass, GenPartPhi, GenPartPt;
	//std::unique_ptr<TTreeReaderArray<int>> genPartIdxMother, GenPartPdgId, GenPartStatus, GenPartStatusFlags;

	ROOT::Math::PtEtaPhiMVector leptonP4;
	ROOT::Math::PtEtaPhiMVector neutrinoP4;
	float maxLeptonPt = -999;
	std::vector<int> idxTauLepton;

	NGenPart = *nGenPart->Get();
	for (int i = 0; i < NGenPart; i++){
		const int& pdgId = genPartPdgId->At(i);
		const int& status = genPartStatus->At(i);
		const int& statusFlag = genPartStatusFlag->At(i);
		const int& idxMother = genPartIdxMother->At(i);

		const float& leptonPt = genPartPt->At(i);
		if (leptonPt > maxLeptonPt) {maxLeptonPt = leptonPt;}

		if ((abs(pdgId) == 13 || abs(pdgId) == 11 ) && idxMother >= 0) { // 13:Muon, 11:Electron
			if (statusFlag == 2) { // 2 : isTauDecayProduct
				NGenLeptonFromTau++;
				NGenTau++;
				idxTauLepton.push_back(i);
			} else {
				NGenLepton++;
				//const int& motherStatus = genPartStatus->At(idxMother);
				const int& motherId = genPartPdgId->At(idxMother);
				const int& idxGrandMother = genPartIdxMother->At(idxMother);
				if(idxGrandMother > 0){
					const int& grandMotherId = genPartPdgId->At(idxGrandMother);
					GenLepGrandMotherId.push_back(grandMotherId);
				}

				if (abs(motherId) == 24 && status == 23) { // genLep is outgoing and has W as mother
					const float& motherPt = genPartPt->At(idxMother);
					const float& motherPhi = genPartPhi->At(idxMother);
					const float& motherEta = genPartEta->At(idxMother);
					const float& motherMass = genPartMass->At(idxMother);

					ROOT::Math::PtEtaPhiMVector motherP4(motherPt, motherEta, motherPhi, motherMass);
					GenWDirectMass.push_back(motherMass);

					for (int j = 0; j < NGenPart; j++){
						const int& idxNeutrinoMother = genPartIdxMother->At(j);
						const int& statusNeutrino = genPartStatus->At(j);
						if (idxNeutrinoMother == motherId && statusNeutrino == 23) {
							NGenMatchedW++;
							//leptonPt already read
							const float& leptonPhi = genPartPhi->At(i);
							const float& leptonEta = genPartEta->At(i);
							const float& leptonMass = genPartMass->At(i);

							const float& neutrinoPt = genPartPt->At(j);
							const float& neutrinoPhi = genPartPhi->At(j);
							const float& neutrinoEta = genPartEta->At(j);
							const float& neutrinoMass = genPartMass->At(j);

							leptonP4 = ROOT::Math::PtEtaPhiMVector(leptonPt, leptonEta, leptonPhi, leptonMass);
							neutrinoP4 = ROOT::Math::PtEtaPhiMVector(neutrinoPt, neutrinoEta, neutrinoPhi, neutrinoMass);
							ROOT::Math::PtEtaPhiMVector wBosonP4 = leptonP4 + neutrinoP4;

							GenDeltaPhiLepWSum.push_back(DeltaPhi(leptonP4, wBosonP4));
							GenWSumMass.push_back(wBosonP4.M());
							GenDeltaPhiLepWDirect.push_back(DeltaPhi(leptonP4, motherP4));
							GenMTLepNu.push_back(wBosonP4.Mt());
						}
					}
				}
				GenLepMotherId.push_back(motherId);
			}
		} else if (abs(pdgId) == 15 && idxMother > 0) { //15:Tau
			NGenTau++;
			const int& motherId = genPartPdgId->At(idxMother);
			const int& idxGrandMother = genPartIdxMother->At(idxMother);
			if(idxGrandMother > 0){
				const int& grandMotherId = genPartPdgId->At(idxGrandMother);
				GenTauGrandMotherId.push_back(grandMotherId);
			}

			GenTauMotherId.push_back(motherId);
		}
	}

	if (NGenLepton + NGenLeptonFromTau== 2) {
		isDiLeptonEvent = true;
		isHadTauEvent = NGenTau > NGenLeptonFromTau;
		if (isHadTauEvent) { LeptonDecayChannelFlag =1;}
		for (int idx : idxTauLepton) {
			const float& pt = genPartPt->At(idx);
			//const float& phi = genPartPhi->At(idx);
			const float& eta = genPartEta->At(idx);
			//const float& mass = genPartMass->At(idx);
			const float& pdgId = genPartPdgId->At(idx);

			if (maxLeptonPt >= 25 && pt < 10) { LeptonsInAcceptance = false;}
			if (maxLeptonPt < 25 && pt < 5) { LeptonsInAcceptance = false;}
			if (abs(eta) > 2.5) { LeptonsInAcceptance = false;}
			if (abs(eta) > 2.5) { LeptonsInAcceptance = false;}
			if (pdgId == 11 && abs(eta) >= 1.44 && abs(eta) < 1.57) { LeptonsInAcceptance = false;}

			if (isHadTauEvent && !LeptonsInAcceptance) { LeptonDecayChannelFlag = 0;}
			else if (isHadTauEvent) { LeptonDecayChannelFlag = 1;}
			else if (!LeptonsInAcceptance) { LeptonDecayChannelFlag = 2;}
			else { LeptonDecayChannelFlag = 3;}
		}
	}

	std::string cutName("No Cuts");
	cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
}

void GenLevelProducer::EndJob(TFile* file) {
}
