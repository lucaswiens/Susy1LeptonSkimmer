#ifndef CUTFLOW_H
#define CUTFLOW_H

#include <vector>
#include <string>
#include <functional>
#include <memory>

#include <TFile.h>
#include <TH1F.h>

#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/Susy1LeptonProduct.h>

class CutFlow{
	private:
		std::shared_ptr<TH1D> hist;
		std::vector<std::function<bool()>> cuts;
		std::vector<std::string> cutNames;

		std::function<bool()> ConstructCut(int &value, const std::string &op, const int &threshold);

	public:
		CutFlow(){}
		CutFlow(TFile &outputFile, const std::string &channel);
		void AddCut(const std::string &part, Susy1LeptonProduct& product, const std::string &op, const int &threshold);

		// The Lambda function checks if any of the triggers are true and only returns false if all triggers fail
		void AddTrigger(const std::vector<int>& triggerIndex, Susy1LeptonProduct &product){
			cuts.insert(cuts.begin(), [&product, triggerIndex](){for(const int& index : triggerIndex){if(product.triggerValues[index]) return true;} return false;});
			cutNames.insert(cutNames.begin(), "Trigger");
		}

		// The Lambda function checks if any of the MET triggers are true and only returns false if all triggers fail
		void AddMetFilter(Susy1LeptonProduct &product){
			cuts.insert(cuts.begin(), [&product](){for(const bool& passed : product.metFilterValues){if(!passed) return false;} return true;});
			cutNames.insert(cutNames.begin(), "MET Filter");
		}

		// The Lambda function checks if either the MET triggers or Lepton triggers are true and only return false if all of them are false
		void AddTriggerOr(const std::vector<int>& triggerIndex, Susy1LeptonProduct &product, const std::string &channel){
			cuts.insert(cuts.begin(), [&product, triggerIndex](){
					bool hltLepOr = false, hltMetOr = false;
					for (const int& index : triggerIndex) {
						hltLepOr = hltLepOr || product.triggerValues[index];
					}
					for (const bool& passed : product.metFilterValues) {
						hltMetOr = hltMetOr || passed;
					}
					return hltLepOr || hltMetOr;
				}
			);
			cutNames.insert(cutNames.begin(), "Trigger(HLT_" + channel + "Or || HLT_METOr)");
		}

		bool Passed();
		void Count(){hist->Fill("No cuts", 1);}; // TODO maybe just use GetEntries()?
		void FillCutflow();
		void WriteOutput(){ hist->Write(0, TObject::kOverwrite);};
};

#endif

