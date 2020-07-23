#include <iostream>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

float WeightCalculator::GetWeight(float x, TH1F* histogram) {
	if(histogram==NULL) {
		std::cout << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t bin = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	return histogram->GetBinContent(bin);
}

float WeightCalculator::GetWeightErr(float x, TH1F* histogram) {
	if(histogram==NULL) {
		std::cout << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t bin = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	return histogram->GetBinError(bin);
}

float WeightCalculator::Get2DWeight(float x, float y, TH2F* histogram) {
	if(histogram==NULL) {
		std::cout << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t binX = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	Int_t binY = std::max(1, std::min(histogram->GetNbinsY(), histogram->GetYaxis()->FindBin(y)));
	return histogram->GetBinContent(binX, binY);
}

float WeightCalculator::Get2DWeightErr(float x, float y, TH2F* histogram) {
	if(histogram==NULL) {
		std::cout << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t binX = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	Int_t binY = std::max(1, std::min(histogram->GetNbinsY(), histogram->GetYaxis()->FindBin(y)));
	return histogram->GetBinError(binX, binY);
	return 1.;
}

