#pragma once


#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TH2D.h>

#include <iomanip>

class WeightCalculator {
private:
	bool norm, verbose;
	std::vector<float> refvals_, targetvals_;
public:
	WeightCalculator();
	WeightCalculator(const bool& norm, const bool& verbose=false);

	double GetWeight(float x, TH1D* histogram);
	double GetWeightErr(float x, TH1D* histogram);
	float Get2DWeight(float x, float y, TH2F* histogram);
	float Get2DWeightErr(float x, float y, TH2F* histogram);
	TH1D* Ratio(TH1D* hist, TH1D* targethist, bool fixLargeWgts = false);
	std::vector<float> loadVals(TH1D* hist, bool norm);
	void FixLargeWeights(std::vector<float> &weights, float maxshift=0.0025, float hardmax=3);
	float CheckIntegral(std::vector<float> wgt1, std::vector<float> wgt2);
};
