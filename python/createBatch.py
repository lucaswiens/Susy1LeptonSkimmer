#!/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_8/external/slc7_amd64_gcc700/bin/python
import sys, os
import argparse
import subprocess
import shutil
from re import findall

import getXSection

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
	xSection = getXSection.GetXSection(sample)
	return isData, isSignal, year, runPeriod, isFastSim, xSection

def GetOSVariable(Var):
	try:
		variable = os.environ[Var]
	except KeyError:
		print "Please set the environment variable " + Var
		sys.exit(1)
	return variable

if __name__=="__main__":
	date = subprocess.check_output("date +\"%Y_%m_%d\"", shell=True).replace("\n", "")
	cmsswBase = GetOSVariable("CMSSW_BASE")
	#redirector = "root://cms-xrd-global.cern.ch/"
	redirector = "root://xrootd-cms.infn.it/"

	parser = argparse.ArgumentParser(description="Runs a NAF batch system for nanoAOD", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", "--input-file", required=True, help="Path to the file containing a list of samples.")
	parser.add_argument("-o", "--output", help="Path to the output directory", default = cmsswBase + "/" "Batch/" + date)
	parser.add_argument("--do-systematics", help="Perform systematic variations", default = "False") #vs. localSkim for faster skimming

	args = parser.parse_args()

	#Make specific subdirectory depending on input file (expects input file to be path/to/input/inputFile.md)
	outputDirName = args.input_file.split("/")[-2] + "_" + args.input_file.split("/")[-1][:-3]
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
			isData, isSignal, year, runPeriod, isFastSim, xSection = prepareArguments(sample)
			sampleName = sample.replace("/", "_")[1:]

			for filename in fileList:
				argumentFile.write(redirector + filename + " " + str(sampleName) + " " + str(isData) + " " + str(args.do_systematics) + " " + str(year) + " " + str(runPeriod) + " " + str(xSection) + " " + cmsswBase + "/src\n")
	sampleFile.close()
	argumentFile.close()

	print "Condor submission file created. Can be submitted via:"
	print "cd " + args.output
	print "condor_submit condor.submit"
