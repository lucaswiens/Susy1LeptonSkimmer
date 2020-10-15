#!/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_8/external/slc7_amd64_gcc700/bin/python
import sys, os
import argparse
import subprocess
import shutil
from re import findall

def createFileList(sampleFileName):
	fileList = []
	sampleFile = open(sampleFileName, "r")
	return findFileLocation(sample)

def findFileLocation(sample):
	return subprocess.check_output("dasgoclient -query=\"file dataset=" + str(sample) + "\"", shell=True).split()

def prepareArguments(sample):
	isData = True
	isSignal = False
	isUSER = False
	year = "unknown"
	runPeriod = "m"
	isFastSim = False
	if "RunIISummer16" in str(sample) or "Run2016" in str(sample):
		year = 2016
		if "PUSummer16v3Fast" in str(sample):
			isFastSim == True
		else:
			isFastSim == False
	elif "RunIIFall17" in str(sample) or "Run2017" in str(sample):
		year = 2017
		if "TuneCP2" in str(sample):
			isFastSim == True
		else:
			isFastSim == False
	if "/NANOAODSIM" in sample:
		isData = False
	if "/SMS-T1tttt" in sample:
		isSignal = True
	if isData:
		runPeriod = findall("Run201..", sample)[0][-1]
	return isData, isSignal, year, runPeriod, isFastSim

def getXSec(sample):
  if sample.find( "WJetsToQQ_HT"                           ) !=-1 : return 95.14;
  elif sample.find( "ZJetsToQQ_HT600toInf"                 ) !=-1 : return 41.34;
  elif sample.find( "QCD_Pt_170to300_"                     ) !=-1 : return 117276.;
  elif sample.find( "QCD_Pt_300to470_"                     ) !=-1 : return 7823.;
  elif sample.find( "QCD_Pt_470to600_"                     ) !=-1 : return 648.2;
  elif sample.find( "QCD_Pt_600to800_"                     ) !=-1 : return 186.9;
  elif sample.find( "QCD_Pt_800to1000_"                    ) !=-1 : return 32.293;
  elif sample.find( "QCD_Pt_1000to1400_"                   ) !=-1 : return 9.4183;
  elif sample.find( "QCD_Pt_1400to1800_"                   ) !=-1 : return 0.84265;
  elif sample.find( "QCD_Pt_1800to2400_"                   ) !=-1 : return 0.114943;
  elif sample.find( "QCD_Pt_2400to3200_"                   ) !=-1 : return 0.006830;
  elif sample.find( "QCD_Pt_3200toInf_"                    ) !=-1 : return 0.000165445;
  elif sample.find( "QCD_HT100to200"                       ) !=-1 : return 27990000;
  elif sample.find( "QCD_HT200to300"                       ) !=-1 : return 1712000.;
  elif sample.find( "QCD_HT300to500"                       ) !=-1 : return 347700.;
  elif sample.find( "QCD_HT500to700"                       ) !=-1 : return 32100.;
  elif sample.find( "QCD_HT700to1000"                      ) !=-1 : return 6831.;
  elif sample.find( "QCD_HT1000to1500"                     ) !=-1 : return 1207.;
  elif sample.find( "QCD_HT1500to2000"                     ) !=-1 : return 119.9;
  elif sample.find( "QCD_HT2000toInf"                      ) !=-1 : return 25.24;
  elif sample.find("QCD_Pt-15to7000" ) !=-1 or sample.find( "QCD_Pt_15to7000" ) !=-1: return  2.022100000e+09*60.5387252324;
  elif sample.find("TT_TuneCUETP8M2T4"                    ) !=-1 : return  831.76      ;
  elif sample.find("TT_TuneEE5C_13TeV-powheg-herwigpp"    ) !=-1 : return  831.76      ;
  elif sample.find("TT_"    							   ) !=-1 : return  831.76      ; #for testing
  elif sample.find("TTJets_"                              ) !=-1 : return  831.76      ;
  elif sample.find("WJetsToLNu_HT-100To200"                ) !=-1 : return 1347*1.21    ;
  elif sample.find("WJetsToLNu_HT-200To400"                ) !=-1 : return 360*1.21     ;
  elif sample.find("WJetsToLNu_HT-400To600"                ) !=-1 : return 48.9*1.21    ;
  elif sample.find("WJetsToLNu_HT-600To800"                ) !=-1 : return 12.08*1.21   ;
  elif sample.find("WJetsToLNu_HT-800To1200"               ) !=-1 : return 5.26*1.21    ;
  elif sample.find("WJetsToLNu_HT-1200To2500"              ) !=-1 : return 1.33*1.21    ;
  elif sample.find("WJetsToLNu_HT-2500ToInf"               ) !=-1 : return 0.03089*1.21 ;
  elif sample.find("WJetsToLNu_TuneCUETP8M1"              ) !=-1 : return 50380.0*1.22 ;
  elif sample.find("W1JetsToLNu_TuneCUETP8M1"             ) !=-1 : return 9644.5*1.22 ;
  elif sample.find("W2JetsToLNu_TuneCUETP8M1"             ) !=-1 : return 3144.5*1.22 ;
  elif sample.find("W3JetsToLNu_TuneCUETP8M1"             ) !=-1 : return  954.8*1.22 ;
  elif sample.find("W4JetsToLNu_TuneCUETP8M1"             ) !=-1 : return  485.6*1.22 ;
  elif sample.find("WW_TuneCUETP8M1"                       ) !=-1 : return 118.7        ;
  elif sample.find("WZ_TuneCUETP8M1"                       ) !=-1 : return 47.13        ;
  elif sample.find("ZZ_TuneCUETP8M1"                       ) !=-1 : return 16.5         ;
  elif sample.find("ST_s-channel_4f_leptonDecays"          ) !=-1 : return 11.36*0.3272 ;
  elif sample.find("ST_t-channel_top_4f_leptonDecays"      ) !=-1 : return 136.02*0.322 ;
  elif sample.find("ST_t-channel_antitop_4f_leptonDecays"  ) !=-1 : return 80.95*0.322  ;
  elif sample.find("ST_t-channel_antitop_4f_inclusiveDecays") !=-1 : return 136.02      ;
  elif sample.find("ST_t-channel_top_4f_inclusiveDecays"   ) !=-1 : return 80.95        ;
  elif sample.find("ST_tW_antitop_5f_inclusiveDecays"      ) !=-1 : return 35.6         ;
  elif sample.find("ST_tW_top_5f_inclusiveDecays_"         ) !=-1 : return 35.6         ;
  elif sample.find("SMS-T1tttt_mGluino-1200_mLSP-800"         ) !=-1 : return 0.04129   ;
  elif sample.find("SMS-T1tttt_mGluino-1500_mLSP-100"         ) !=-1 : return 0.006889  ;
  elif sample.find("SMS-T1tttt_mGluino-2000_mLSP-100"         ) !=-1 : return 0.0004488	;
  elif sample.find("SingleMuon")!=-1  or sample.find("SingleElectron") !=-1 or sample.find("JetHT")  !=-1 or sample.find("MET") !=-1 or sample.find("MTHT") !=-1: return 1.
  else:
	print "Cross section not defined for this sample!!"
	return 0.

def getOSVariable(Var):
	try:
		variable = os.environ[Var]
	except KeyError:
		print "Please set the environment variable " + Var
		sys.exit(1)
	return variable

if __name__=="__main__":
	date = subprocess.check_output("date +\"%Y_%m_%d\"", shell=True).replace("\n", "")
	cmsswBase = getOSVariable("CMSSW_BASE")

	parser = argparse.ArgumentParser(description="Runs a NAF batch system for nanoAOD", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", "--input-file", required=True, help="Path to the file containing a list of samples.")
	parser.add_argument("-o", "--output", help="Path to the output directory", default = cmsswBase + "/" "Batch/" + date)
	parser.add_argument("--do-systematics", help="Perform systematic variations", default = "False") #vs. localSkim for faster skimming

	args = parser.parse_args()

	#Make specific subdirectory depending on input file (expects input file to be path/to/input/inputFile.md)
	outputDirName = args.input_file.split("/")[-1][:-3]
	args.output = args.output + "/" + outputDirName
	executable = cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/scripts/produceSkim"

	if  os.path.exists(args.output):
		keepDirectory = raw_input("Output directory already exists: " + str(args.output) + " Do you want to remove it [y/n]: ")
		if ( "y" in keepDirectory or "Y" in keepDirectory or "Yes" in keepDirectory or "yes" in keepDirectory):
			shutil.rmtree(str(args.output))
			os.makedirs(str(args.output))
			os.makedirs(str(args.output) + "/logs")
			os.makedirs(str(args.output) + "/error")
			os.makedirs(str(args.output) + "/output")
		elif ( "N" in keepDirectory or  "n" in keepDirectory or  "No" in keepDirectory or "no" in keepDirectory): print str(args.output) , "will be ovewritten by the job output -- take care"
		else:
			raise ValueError( "invalid input, answer with \"Yes\" or \"No\"")
	else:
		os.makedirs(str(args.output))
		os.makedirs(str(args.output) + "/logs")
		os.makedirs(str(args.output) + "/error")
		os.makedirs(str(args.output) + "/output")


	submitFileContent = open(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/scripts/condor.submit", "r").read()
	submitFileContent = submitFileContent.replace("@EXECUTABLE", executable)
	submitFileContent = submitFileContent.replace("@OUT", args.output)

	submitFile = open(args.output + "/condor.submit", "w")
	submitFile.write(submitFileContent)
	submitFile.close()

	sampleFile = open(args.input_file, "r")
	argumentFile = open(args.output + "/arguments.md", "w")
	for sample in sampleFile:
		sample = sample.strip()
		if not (sample.startswith("#") or sample in ["", "\n", "\r\n"]):
			fileList = findFileLocation(sample)
			isData, isSignal, year, runPeriod, isFastSim = prepareArguments(sample)
			sampleName = sample.replace("/", "_")[1:]

			for filename in fileList:
				argumentFile.write("root://cms-xrd-global.cern.ch/" + filename + " " + str(sampleName) + " " + str(isData) + " " + str(args.do_systematics) + " " + str(year) + " " + str(runPeriod) + " " + cmsswBase + "/src\n")
	sampleFile.close()
	argumentFile.close()

	print "Condor submission file created. Can be submitted via:"
	print "cd " + args.output
	print "condor_submit condor.submit"
