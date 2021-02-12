import json
import os

"""
	Short python script to creat a GetXSection function used for plotting.
	Input json files are taken from:
	https://cms-gen-dev.cern.ch/xsdb/
"""

def GetOSVariable(Var):
	try:
		variable = os.environ[Var]
	except KeyError:
		print "Please set the environment variable " + Var
		sys.exit(1)
	return variable

if __name__ == "__main__":
	cmsswBase = GetOSVariable("CMSSW_BASE")
	xsDict = {}
	#List of json files download from xsdb, edit manually
	jsonList = ["DY.json",  "QCD.json",  "ST.json",  "TTJets.json",  "TTW.json",  "TTZ.json",  "WJets.json",  "WW.json",  "WZ.json",  "ZZ.json",  "tZq.json", "TTTo.json"]

	xsFile = open("python/getXSection.py", "w")
	xsFile.write("def GetXSection(fileName): #[pb]" + "\n")
	xsFile.write("\t#Cross Section derived from sample name using https://cms-gen-dev.cern.ch/xsdb/\n")
	xsFile.write("\t#TODO UL QCD files have PSWeights in their name, but xsdb does not include it in its name\n")
	xsFile.write("\tfileName = fileName.replace(\"PSWeights_\", \"\")" + "\n")
	xsFile.write("\t#TODO: UL Single Top xSection only defined for filename without inclusive decays specified\n")
	xsFile.write("\tfileName = fileName.replace(\"5f_InclusiveDecays_\", \"5f_\")" + "\n")
	xsFile.write("\t#TODO: UL tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8 only defined with PSWeights...\n")
	xsFile.write("\tfileName = fileName.replace(\"tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8\", \"tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8\")" + "\n")

	first = True
	for jsonFile in jsonList:
		with open(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/xSection/" + jsonFile) as xsecFile:
			xSecDict = json.load(xsecFile)

		for dict in xSecDict:
			if dict["process_name"] not in xsDict.keys():
				xsDict[dict["process_name"]] = dict["cross_section"]
	xsKeys = xsDict.keys()
	for key in sorted(xsKeys, key = lambda process : -len(process)):
		if first:
			xsFile.write("\tif   fileName.find(\"" + key + "\") !=-1 : return " + str(xsDict[key]) + "\n")
			first = False
		else:
			xsFile.write("\telif fileName.find(\"" + key + "\") !=-1 : return " + str(xsDict[key]) + "\n")

	xsFile.write("\telif fileName.find(\"SingleMuon\")!=-1  or fileName.find(\"SingleElectron\") !=-1 or fileName.find(\"JetHT\")  !=-1 or fileName.find(\"MET\") !=-1 or fileName.find(\"MTHT\") !=-1: return 1." + "\n")
	xsFile.write("\telse:" + "\n")
	xsFile.write("\t\tprint(\"Cross section not defined! Returning 0 and skipping sample:\\n{}\\n\".format(fileName))" + "\n")
	xsFile.write("\t\treturn 0" + "\n")
