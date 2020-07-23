#pragma once

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

class WeightCalculator {
public:
	float GetWeight(float x, TH1F* histogram);
	float GetWeightErr(float x, TH1F* histogram);
	float Get2DWeight(float x, float y, TH2F* histogram);
	float Get2DWeightErr(float x, float y, TH2F* histogram);
};
