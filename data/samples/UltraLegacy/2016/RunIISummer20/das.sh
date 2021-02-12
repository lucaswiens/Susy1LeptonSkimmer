#!/bin/sh
for sample in DY QCD ST TTJets TTTo2L2Nu TTToSemiLeptonic TTToW TTToZ WW WZ ZZ tZq W1Jets W2Jets W3Jets W4Jets
do
	echo \#${sample} preVFP > ${sample}.md
	dasgoclient -query="dataset=/${sample}*/RunIISummer20UL16NanoAODAPVv2-*/*" >> ${sample}.md
	echo \#${sample} postVFP >> ${sample}.md
	dasgoclient -query="dataset=/${sample}*/RunIISummer20UL16NanoAODv2-*/*" >> ${sample}.md
done
cat W[0-9]*Jets.md > WJets.md
rm W[0-9]*Jets.md
rm MC.md

cat *md > MC.md

for sample in SingleMuon SingleElectron
do
	#echo \#${sample} preVFP > ${sample}.md
	dasgoclient -query="dataset=/${sample}*/Run2016*UL2016_MiniAODv1_NanoAODv2*/NANOAOD" >> ${sample}.md
	#echo \#${sample} postVFP >> ${sample}.md
	#dasgoclient -query="dataset=/${sample}*/Run2016*UL2016_MiniAODv1_NanoAODv2*/NANOAOD" >> ${sample}.md
done

cat SingleMuon.md SingleElectron.md > Data.md
