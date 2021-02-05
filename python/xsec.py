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
	jsonList = ["DY.json",  "qcd.json",  "st.json",  "ttjets.json",  "ttw.json",  "ttz.json",  "vv.json", "wjets.json"]

	xsFile = open("python/getXSection.py", "w")
	xsFile.write("def GetXSection(fileName): #[pb]" + "\n")
	first = True

	for jsonFile in jsonList:
		with open(cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/data/xSection/" + jsonFile) as xsecFile:
			xSecDict = json.load(xsecFile)

		for dir in xSecDict:
			if dir["process_name"] not in xsDict.keys():
				xsDict[dir["process_name"]] = dir["cross_section"]
				#print("\telif fileName.find(\"" + dir["process_name"] + "\") !=-1 : return " + str(dir["cross_section"]))
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
