#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/CutFlow.h>
#include <iostream>

CutFlow::CutFlow(TFile &outputFile, const std::string &channel, const std::string &systematic, const std::string &shift) {
	outputFile.cd(channel.c_str());
	hist = std::make_shared<TH1D>();
	hist->SetName(("Cutflow_" + systematic + shift).c_str());
	hist->SetTitle(("Cutflow_" + systematic + shift).c_str());
	hist->SetDirectory(outputFile.GetDirectory(channel.c_str()));
}

std::function<bool()> CutFlow::ConstructCut(int &value, const std::string &op, const int &threshold) {
	if(op == "==") {
		return [&value, threshold]() {return value == threshold;};
	} else if(op == ">=") {
		return [&value, threshold]() {return value >= threshold;};
	} else if(op == "<=") {
		return [&value, threshold]() {return value <= threshold;};
	} else {
		throw std::runtime_error(("Unknow cut operator: '" + op + "'").c_str());
	}
}

void CutFlow::AddCut(const std::string &part, Susy1LeptonProduct &product, const std::string &op, const int &threshold) {
	if(part == "Electron") {
		cuts.push_back(ConstructCut(product.nElectron, op, threshold));
		cutNames.push_back("N_{e} " + op + std::to_string(threshold) + " (No ID.)");
	} else if(part == "GoodElectron") {
		cuts.push_back(ConstructCut(product.nGoodElectron, op, threshold));
		cutNames.push_back("N_{e}^{good} " + op + std::to_string(threshold) + " (No ID.)");
	} else if(part == "VetoElectron") {
		cuts.push_back(ConstructCut(product.nVetoElectron, op, threshold));
		cutNames.push_back("N_{e}^{veto} " + op + std::to_string(threshold) + " (No ID.)");
	} else if(part == "Muon") {
		cuts.push_back(ConstructCut(product.nMuon, op, threshold));
		cutNames.push_back("N_{#mu} " + op + std::to_string(threshold) + " (No ID.)");
	} else if(part == "GoodMuon") {
		cuts.push_back(ConstructCut(product.nGoodMuon, op, threshold));
		cutNames.push_back("N_{#mu}^{good} " + op + std::to_string(threshold));
	} else if(part == "VetoMuon") {
		cuts.push_back(ConstructCut(product.nVetoMuon, op, threshold));
		cutNames.push_back("N_{#mu}^{veto} " + op + std::to_string(threshold));
	} else if(part == "Jet") {
		cuts.push_back(ConstructCut(product.nJet, op, threshold));
		//cutNames.push_back("N_{j} " + op + std::to_string(threshold) + " (Not clean)");
		cutNames.push_back("N_{j} " + op + std::to_string(threshold) + " (Cleaned)");
	} else {
		throw std::runtime_error(("Unknow part: '" + op + "'").c_str());
	}
}

bool CutFlow::Passed() {
	for(std::size_t i = 0; i < cuts.size(); ++i) {
		if(!cuts[i]()) return false;
	}

	return true;
}

void CutFlow::FillCutflow() {
	for(std::size_t i = 0; i < cuts.size(); ++i) {
		if(!cuts[i]()) return;
		hist->Fill(cutNames[i].c_str(), 1);
	}
}
