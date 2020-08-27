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

	//https://indico.cern.ch/event/592621/contributions/2398559/attachments/1383909/2105089/16-12-05_ana_manuelf_isr.pdf
	ISRweights_Mar17       = {{0, 1}, {1, 0.920},  {2, 0.821}, {3, 0.715}, {4, 0.662}, {5, 0.561}, {6, 0.511}};
	ISRweights_ICHEP16     = {{0, 1}, {1, 0.882}, {2, 0.792}, {3, 0.702}, {4, 0.648}, {5, 0.601}, {6, 0.515}};
	ISRweightssyst_Mar17   = {{0, 0}, {1, 0.040},  {2, 0.090}, {3, 0.143}, {4, 0.169}, {5, 0.219}, {6, 0.244}};
	ISRweightssyst_ICHEP16 = {{0, 0}, {1, 0.059},  {2, 0.104}, {3, 0.149}, {4, 0.176}, {5, 0.199}, {6, 0.242}};

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

	//nISR Reweighting
	tree->Branch("nISR", &nISR);
	tree->Branch("nISRWeight", &nISRWeight);
	tree->Branch("nISRttweightsystUp", &nISRWeightUp);
	tree->Branch("nISRttweightsystDown", &nISRWeightDown);

	tree->Branch("nISRWeight_Mar17", &nISRWeight_Mar17);
	tree->Branch("nISRWeightUp_Mar17", &nISRWeightUp_Mar17);
	tree->Branch("nISRWeightDown_Mar17", &nISRWeightDown_Mar17);

	tree->Branch("nISRWeight_ICHEP16", &nISRWeight_ICHEP16);
	tree->Branch("nISRWeightUp_ICHEP16", &nISRWeightUp_ICHEP16);
	tree->Branch("nISRWeightDown_ICHEP16", &nISRWeightDown_ICHEP16);


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

	unsigned int nDaughters = 0;
	std::vector<int> idxDaughter;
	NGenPart = *nGenPart->Get();
	for (int i = 0; i < NGenPart; i++){
		const int& pdgId = genPartPdgId->At(i);
		const int& status = genPartStatus->At(i);
		const int& statusFlag = genPartStatusFlag->At(i);
		const int& idxMother = genPartIdxMother->At(i);

		const float& leptonPt = genPartPt->At(i);
		if (leptonPt > maxLeptonPt) {maxLeptonPt = leptonPt;}

		if (idxMother >= 0) { //store daughter indices for nISR
			nDaughters++;
			idxDaughter.push_back(i);
		}

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
				GenLepMotherId.push_back(idxMother);
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
					product->genWBosonPt.push_back(motherPt);
					product->genWBosonEta.push_back(motherEta);
					product->genWBosonPhi.push_back(motherPhi);
					product->genWBosonMass.push_back(motherMass);
					GenWDirectMass.push_back(motherMass);

					if(idxGrandMother > 0 && abs(GenLepGrandMotherId.back()) == 6){
						//Used for sytematic variations, store in product for now
						const float& grandMotherPt = genPartPt->At(idxGrandMother);
						const float& grandMotherPhi = genPartPhi->At(idxGrandMother);
						const float& grandMotherEta = genPartEta->At(idxGrandMother);
						const float& grandMotherMass = genPartMass->At(idxGrandMother);
						product->genTopPt.push_back(grandMotherPt);
						product->genTopEta.push_back(grandMotherEta);
						product->genTopPhi.push_back(grandMotherPhi);
						product->genTopMass.push_back(grandMotherMass);
					}

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
			GenTauMotherId.push_back(motherId);
			}
		}
		product->nGenWBoson = product->genWBosonPt.size();
		product->nGenTop = product->genTopPt.size();

		if (NGenLepton + NGenLeptonFromTau == 2) {
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

	}

	//https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295
	nISR = 0;
	for (unsigned int ijet = 0; ijet < product->nJet; ijet++) {
		bool matched=false;
		for (int imc = 0; imc < NGenPart; imc++) {
			if (matched) {
				break;
			} else {
				const int& pdgId = genPartPdgId->At(imc);
				const int& status = genPartStatus->At(imc);
				if (status !=23 || abs(pdgId) > 5) continue;
				const int& idxMother = genPartIdxMother->At(imc);
				int motherId = genPartPdgId->At(idxMother);
				if (!(motherId == 6 || motherId == 23 || motherId == 24 || motherId == 25 || motherId>1e6)) continue;
				//check against daughter in case of hard initial splitting
				for (unsigned int idau = 0; idau < nDaughters; idau++) {
					const float& genPhi = genPartPhi->At(idxDaughter.at(idau));
					const float& genEta = genPartEta->At(idxDaughter.at(idau));
					float dR = DeltaR(product->jetEta.at(ijet), product->jetPhi.at(ijet), genEta, genPhi);
					if (dR < 0.3 ) {
						matched = true;
						break;
					}
				}
			}
		} // Loop over MC particles
		if (!matched) { nISR++;}
	} // Loop over jets

	//if (nISR > 6) {nISR = 6;}
	float C_ISR = 1.090;
	float C_ISRUp   = 1.043;
	float C_ISRDown = 1.141;

	nISRWeight_Mar17 = C_ISR * ISRweights_Mar17[nISR];
	nISRWeightUp_Mar17   =  C_ISRUp   * (ISRweights_Mar17[nISR] + ISRweightssyst_Mar17[nISR]);
	nISRWeightDown_Mar17 =  C_ISRDown * (ISRweights_Mar17[nISR] - ISRweightssyst_Mar17[nISR]);

	nISRWeight_ICHEP16 = C_ISR * ISRweights_ICHEP16[nISR];
	nISRWeightUp_ICHEP16   =  C_ISRUp   * (ISRweights_ICHEP16[nISR] + ISRweightssyst_ICHEP16[nISR]);
	nISRWeightDown_ICHEP16 =  C_ISRDown * (ISRweights_ICHEP16[nISR] - ISRweightssyst_ICHEP16[nISR]);

	nISRWeight = nISRWeight_Mar17;
	nISRWeightUp   = nISRWeightUp_Mar17;
	nISRWeightDown = nISRWeightDown_Mar17;

	product->nGenPart = NGenPart;
	std::string cutName("No Cuts");
	cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
}

void GenLevelProducer::EndJob(TFile* file) {
}
