#ifndef TRIGGERPRODUCER_H
#define TRIGGERPRODUCER_H

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Producer/BaseProducer.h>
#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

class TriggerProducer : public BaseProducer {
	private:
		//Check if it is data or MC
		//bool isData, doSystematics;

		//Cut Variables
		//int era;
		//char runPeriod;
		//unsigned int run;

		//Vector for the output variables
		//int nIsr, nISRweight, nISRttweightsystUp, nISRttweightsystDown;
		//bool HLT_EleOr, HLT_MuOr, HLT_LepOr, HLT_MetOr;


	public:
		std::string Name = "TriggerProducer";
		TriggerProducer(const int &era, const char &runPeriod);

		void Produce(DataReader &dataReader, Susy1LeptonProduct &product);
		void EndJob(TFile &file);
};

#endif

