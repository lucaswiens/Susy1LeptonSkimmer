#!/bin/sh
for sample in DY QCD ST TTJets TTW TTZ WW WZ ZZ tZq W1Jets W2Jets W3Jets W4Jets
do
	echo \#${sample} preVFP > ${sample}.md
	dasgoclient -query="dataset=/${sample}*/RunIISummer19UL16NanoAODAPVv2-*/*" >> ${sample}.md
	echo \#${sample} postVFP >> ${sample}.md
	dasgoclient -query="dataset=/${sample}*/RunIISummer19UL16NanoAODv2-*/*" >> ${sample}.md
done

for sample in SingleMuon SingleElectron
	#echo \#${sample} preVFP > ${sample}.md
	dasgoclient -query="dataset=/${sample}*/Run2016*UL2016_MiniAODv1_NanoAODv2*/NANOAOD" >> ${sample}.md
	#echo \#${sample} postVFP >> ${sample}.md
	#dasgoclient -query="dataset=/${sample}*/Run2016*UL2016_MiniAODv1_NanoAODv2*/NANOAOD" >> ${sample}.md
do
