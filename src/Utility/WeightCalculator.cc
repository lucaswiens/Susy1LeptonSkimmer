#include <iostream>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Utility/WeightCalculator.h>

WeightCalculator::WeightCalculator(){};
WeightCalculator::WeightCalculator(const bool& norm, const bool& verbose):
	norm(norm),
	verbose(verbose)
	{};

double WeightCalculator::GetWeight(float x, TH1D* histogram) {
	if(histogram==NULL) {
		std::cout << "ERROR! The weights input* histogram is not loaded. Returning weight 0!" << std::endl;
		return 0.;
	}
	Int_t bin = std::max(1, std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(x)));
	return histogram->GetBinContent(bin);
}

double WeightCalculator::GetWeightErr(float x, TH1D* histogram) {
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

TH1D* WeightCalculator::Ratio(TH1D* hist, TH1D* targethist, bool fixLargeWgts) {
	TH1D* ret = (TH1D*)hist->Clone("hweights");
	ret->SetDirectory(0);

	std::vector<float> vals = loadVals(hist,norm);
	std::vector<float> targetvals = loadVals(targethist,norm);
	std::vector<float> weights;
	int nbins = vals.size();
	if(verbose) std::cout << "Weights for variable " << hist->GetName() << " with a number of bins equal to " << nbins << ":" << std::endl;
	for(int i=0; i<nbins; ++i) {
		float weight = vals[i] !=0 ? targetvals[i]/vals[i] : 1.;
		if(verbose) std::cout << std::setprecision(3) << weight << " ";
		weights.push_back(weight);
	}
	if(verbose) std::cout << "." << std::endl;
	if(fixLargeWgts) FixLargeWeights(weights);
	if(verbose) std::cout << "Final weights: " << std::endl;
	for(int i=0; i<(int)weights.size(); ++i) {
		ret->SetBinContent(i,weights[i]);
		if(verbose) std::cout << std::setprecision(3) << weights[i] << " ";
	}
	if(verbose) std::cout << "." << std::endl;
	return ret;
}

std::vector<float> WeightCalculator::loadVals(TH1D* hist, bool norm) {
	int nbins=hist->GetNcells();
	std::vector<float> vals;
	for(int i=0; i<nbins; ++i) {
		double bc=hist->GetBinContent(i);
		//double val = (i>0 && bc==0 && hist->GetBinContent(i-1)>0 && hist->GetBinContent(i+1)>0) ? 0.5*(hist->GetBinContent(i-1)+hist->GetBinContent(i+1)) : bc;
		vals.push_back(std::max(bc,0.));
	}
	if(verbose) std::cout << "Normalization of " << hist->GetName() << ": " << hist->Integral() << std::endl;
	if(norm) {
		float scale = 1.0/hist->Integral();
		for(int i=0; i<nbins; ++i) vals[i] *= scale;
	}
	return vals;
}

void WeightCalculator::FixLargeWeights(std::vector<float> &weights, float maxshift,float hardmax) {
	float maxw = std::min(*(std::max_element(weights.begin(),weights.end())),float(5.));
	std::vector<float> cropped;
	while (maxw > hardmax) {
		cropped.clear();
		for(int i=0; i<(int)weights.size(); ++i) cropped.push_back(std::min(maxw,weights[i]));
		float shift = CheckIntegral(cropped,weights);
		if(verbose) std::cout << "For maximum weight " << maxw << ": integral relative change: " << shift << std::endl;
		if(fabs(shift) > maxshift) break;
		maxw *= 0.95;
	}
	maxw /= 0.95;
	if (cropped.size()>0) {
			for(int i=0; i<(int)weights.size(); ++i) cropped[i] = std::min(maxw,weights[i]);
			float normshift = CheckIntegral(cropped,weights);
			for(int i=0; i<(int)weights.size(); ++i) weights[i] = cropped[i]*(1-normshift);
	}
}

float WeightCalculator::CheckIntegral(std::vector<float> wgt1, std::vector<float> wgt2) {
	float myint=0;
	float refint=0;
	for(int i=0; i<(int)wgt1.size(); ++i) {
		myint += wgt1[i]*refvals_[i];
		refint += wgt2[i]*refvals_[i];
	}
	return (myint-refint)/refint;
}
