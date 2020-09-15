#include <Susy1LeptonAnalysis/Susy1LeptonSkimmer/interface/NanoSkimmer.h>

#include <string>
#include <sstream>

float CrossSection(std::string);

int main(int argc, char* argv[]) {
	//Extract informations of command line
	std::string fileName = std::string(argv[1]);
	bool isData = std::string(argv[2]) == "True" ? true : false;
	int era = std::stoi(std::string(argv[3]));
	std::string outName = std::string(argv[4]);
	char runPeriod; //If it is MC, then runPeriod does not matter

	if (isData) {
		runPeriod = (char)*argv[5];
	}

	int nMaxEvents;
	if(argc == 7) {
		nMaxEvents = std::stoi(std::string(argv[6]));
	} else {
		nMaxEvents = -999;
	}

	std::string sampleName = fileName;
	float xSec = CrossSection(sampleName);

	NanoSkimmer skimmer(fileName, isData);
	skimmer.EventLoop(xSec, era, runPeriod, nMaxEvents);
	skimmer.WriteOutput(outName);
}

float CrossSection(std::string sample) {
	if      (sample.find("WJetsToQQ_HT"                           ) != std::string::npos) { return 95.14;}
	else if (sample.find("ZJetsToQQ_HT600toInf"                   ) != std::string::npos) { return 41.34;}
	else if (sample.find("QCD_Pt_170to300_"                       ) != std::string::npos) { return 117276.;}
	else if (sample.find("QCD_Pt_300to470_"                       ) != std::string::npos) { return 7823.;}
	else if (sample.find("QCD_Pt_470to600_"                       ) != std::string::npos) { return 648.2;}
	else if (sample.find("QCD_Pt_600to800_"                       ) != std::string::npos) { return 186.9;}
	else if (sample.find("QCD_Pt_800to1000_"                      ) != std::string::npos) { return 32.293;}
	else if (sample.find("QCD_Pt_1000to1400_"                     ) != std::string::npos) { return 9.4183;}
	else if (sample.find("QCD_Pt_1400to1800_"                     ) != std::string::npos) { return 0.84265;}
	else if (sample.find("QCD_Pt_1800to2400_"                     ) != std::string::npos) { return 0.114943;}
	else if (sample.find("QCD_Pt_2400to3200_"                     ) != std::string::npos) { return 0.006830;}
	else if (sample.find("QCD_Pt_3200toInf_"                      ) != std::string::npos) { return 0.000165445;}
	else if (sample.find("QCD_HT100to200"                         ) != std::string::npos) { return 27990000;}
	else if (sample.find("QCD_HT200to300"                         ) != std::string::npos) { return 1712000.;}
	else if (sample.find("QCD_HT300to500"                         ) != std::string::npos) { return 347700.;}
	else if (sample.find("QCD_HT500to700"                         ) != std::string::npos) { return 32100.;}
	else if (sample.find("QCD_HT700to1000"                        ) != std::string::npos) { return 6831.;}
	else if (sample.find("QCD_HT1000to1500"                       ) != std::string::npos) { return 1207.;}
	else if (sample.find("QCD_HT1500to2000"                       ) != std::string::npos) { return 119.9;}
	else if (sample.find("QCD_HT2000toInf"                        ) != std::string::npos) { return 25.24;}
	else if (sample.find("QCD_Pt-15to7000" ) != std::string::npos || sample.find( "QCD_Pt_15to7000" ) != std::string::npos) { return 2.022100000e+09*60.5387252324;}
	else if (sample.find("TT_TuneCUETP8M2T4"                      ) != std::string::npos) { return 831.76;}
	else if (sample.find("TT_TuneEE5C_13TeV-powheg-herwigpp"      ) != std::string::npos) { return 831.76;}
	else if (sample.find("TT_"                                    ) != std::string::npos) { return 831.76;}
	else if (sample.find("TTJets_"                                ) != std::string::npos) { return 831.76;}
	else if (sample.find("WJetsToLNu_HT-100To200"                 ) != std::string::npos) { return 1347*1.21;}
	else if (sample.find("WJetsToLNu_HT-200To400"                 ) != std::string::npos) { return 360*1.21;}
	else if (sample.find("WJetsToLNu_HT-400To600"                 ) != std::string::npos) { return 48.9*1.21;}
	else if (sample.find("WJetsToLNu_HT-600To800"                 ) != std::string::npos) { return 12.08*1.21;}
	else if (sample.find("WJetsToLNu_HT-800To1200"                ) != std::string::npos) { return 5.26*1.21;}
	else if (sample.find("WJetsToLNu_HT-1200To2500"               ) != std::string::npos) { return 1.33*1.21;}
	else if (sample.find("WJetsToLNu_HT-2500ToInf"                ) != std::string::npos) { return 0.03089*1.21 ;}
	else if (sample.find("WJetsToLNu_TuneCUETP8M1"                ) != std::string::npos) { return 50380.0*1.22 ;}
	else if (sample.find("W1JetsToLNu_TuneCUETP8M1"               ) != std::string::npos) { return 9644.5*1.22 ;}
	else if (sample.find("W2JetsToLNu_TuneCUETP8M1"               ) != std::string::npos) { return 3144.5*1.22 ;}
	else if (sample.find("W3JetsToLNu_TuneCUETP8M1"               ) != std::string::npos) { return 954.8*1.22 ;}
	else if (sample.find("W4JetsToLNu_TuneCUETP8M1"               ) != std::string::npos) { return 485.6*1.22 ;}
	else if (sample.find("WW_TuneCUETP8M1"                        ) != std::string::npos) { return 118.7;}
	else if (sample.find("WZ_TuneCUETP8M1"                        ) != std::string::npos) { return 47.13;}
	else if (sample.find("ZZ_TuneCUETP8M1"                        ) != std::string::npos) { return 16.5;}
	else if (sample.find("ST_s-channel_4f_leptonDecays"           ) != std::string::npos) { return 11.36*0.3272;}
	else if (sample.find("ST_t-channel_top_4f_leptonDecays"       ) != std::string::npos) { return 136.02*0.322;}
	else if (sample.find("ST_t-channel_antitop_4f_leptonDecays"   ) != std::string::npos) { return 80.95*0.322;}
	else if (sample.find("ST_t-channel_antitop_4f_inclusiveDecays") != std::string::npos) { return 136.02;}
	else if (sample.find("ST_t-channel_top_4f_inclusiveDecays"    ) != std::string::npos) { return 80.95;}
	else if (sample.find("ST_tW_antitop_5f_inclusiveDecays"       ) != std::string::npos) { return 35.6;}
	else if (sample.find("ST_tW_top_5f_inclusiveDecays_"          ) != std::string::npos) { return 35.6;}
	else if (sample.find("SMS-T1tttt_mGluino-1200_mLSP-800"       ) != std::string::npos) { return 0.04129;}
	else if (sample.find("SMS-T1tttt_mGluino-1500_mLSP-100"       ) != std::string::npos) { return 0.006889;}
	else if (sample.find("SMS-T1tttt_mGluino-2000_mLSP-100"       ) != std::string::npos) { return 0.0004488;}
	else if (sample.find("SingleMuon")!= std::string::npos  || sample.find("SingleElectron") != std::string::npos || sample.find("JetHT")  != std::string::npos || sample.find("MET") != std::string::npos || sample.find("MTHT") != std::string::npos) { return 1;}
	else {
		std::cout <<  "Cross section not defined for this sample!!" << std::endl;
		return 0.;
	}
}
