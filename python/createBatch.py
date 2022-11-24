#!/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_8/external/slc7_amd64_gcc700/bin/python3
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
	runPeriod = "M"
	if "RunIISummer20UL16" in str(sample) or "Run2016" in str(sample):
		year = 2016
		if "PUSummer16v3Fast" in str(sample): #FastSim
			runPeriod = "S"
	elif "RunIISummer20UL17" in str(sample) or "Run2017" in str(sample):
		year = 2017
		if "TuneCP2" in str(sample): #FastSim
			runPeriod = "S"
	elif "RunIISummer20UL18" in str(sample) or "Run2018" in str(sample):
		year = 2018
		if "TuneCP2" in str(sample): #FastSim
			runPeriod = "S"
	else:
		print("Could not determine era of the sample:\n" + str(sample) + "\nCheck the prepareArguments function to configure it properly")
		sys.exit(1)

	if "/NANOAODSIM" in sample:
		isData = False
	if "/SMS-T1tttt" in sample or "SMS-T5qqqqWW" in sample:
		isSignal = True
	if isData:
		runPeriod = findall("Run201..", sample)[0][-1]
	xSection = getXSection.GetXSection(sample)
	return year, runPeriod, xSection

def GetOSVariable(Var):
	try:
		variable = os.environ[Var]
	except KeyError:
		print("Please set the environment variable " + Var)
		sys.exit(1)
	return variable

if __name__=="__main__":
	date = subprocess.check_output("date +\"%Y_%m_%d\"", shell=True).decode().replace("\n", "")#GetOSVariable("DATE")
	cmsswBase = GetOSVariable("CMSSW_BASE")
	#redirector = "root://cms-xrd-global.cern.ch/"
	redirector = "root://xrootd-cms.infn.it/"

	parser = argparse.ArgumentParser(description="Runs a NAF batch system for nanoAOD", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", "--input-file", required=True, help="Path to the file containing a list of samples.")
	parser.add_argument("-o", "--output", help="Path to the output directory", default = cmsswBase + "/" "Batch/" + date)

	args = parser.parse_args()

	#Make specific subdirectory depending on input file (expects input file to be path/to/input/inputFile.md)
	splitInputFile = args.input_file.split("/")
	outputDirName = splitInputFile[-2] + "/" + splitInputFile[-1].replace(".txt", "")
	args.output = args.output + "/" + outputDirName
	executable = cmsswBase + "/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/scripts/produceSkim"

	try:
		os.makedirs(str(args.output))
		os.makedirs(str(args.output) + "/logs")
		os.makedirs(str(args.output) + "/error")
		os.makedirs(str(args.output) + "/output")
		os.makedirs(str(args.output) + "/root")
	except:
		print("Output Dir already exists")

	submitFile = open(args.output + "/condor.submit", "w")
	submitFile.write(f"executable = {cmsswBase}/src/Susy1LeptonAnalysis/Susy1LeptonSkimmer/scripts/produceSkim\n")
	submitFile.write(f"arguments  = $(inputFile) $(outputFile) $(year) $(runPeriod) $(Process) $(cmsswBase) $(outputDir)\n")
	submitFile.write(f"\n")
	submitFile.write(f"universe       = vanilla\n")
	submitFile.write(f"request_memory = 500 MB\n")
	submitFile.write(f"\n")
	submitFile.write(f"output = {args.output}/output/job$(Cluster)_$(Process).stdout\n")
	submitFile.write(f"error  = {args.output}/error/job$(Cluster)_$(Process).stderr\n")
	submitFile.write(f"log    = {args.output}/logs/job$(Cluster)_$(Process).log\n")
	submitFile.write(f"\n")
	submitFile.write(f"queue inputFile outputFile year runPeriod cmsswBase outputDir from arguments.md\n")
	submitFile.close()

	sampleFile = open(args.input_file, "r")
	argumentFile = open(args.output + "/arguments.md", "w")
	for sample in sampleFile:
		sample = sample.strip()
		if not (sample.startswith("#") or sample in ["", "\n", "\r\n"]):
			fileList = findFileLocation(sample)
			year, runPeriod, xSection = prepareArguments(sample)
			sampleName = sample.replace("/", "_")[1:]

			for filename in fileList:
				argumentFile.write(redirector + filename.decode('utf-8') + " " + str(sampleName) + " " + str(year) + " " + str(runPeriod) + " " + str(xSection) + " " + cmsswBase + "/src " + args.output + "\n")
	sampleFile.close()
	argumentFile.close()

	print("Condor submission file created. Can be submitted via:")
	print("cd " + args.output)
	print("condor_submit condor.submit")
