#!/usr/bin/env python

import ROOT, getopt, sys, math
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetCombinedTH1FsWeightRek, DrawTDRStyle


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:'
keywords = ['path']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

path_tag=""
xaxis=""
secondaries = ""

for o,p in opts: 
    if o in ['-p','--path']:
        path_tag = p

print "" 

#_____________________
#
# READ INPUT FILES
#_____________________

file_help = GetListOfFilesWithEnding(path_tag ,".root")

h3 = ()

f = ()

for file_ite,file_name in enumerate(file_help):
	
	f+= ROOT.TFile.Open(file_name),
	print "File " + file_name + " ...",
	
	f[file_ite].cd("demo")

	liste = ROOT.gDirectory.GetListOfKeys()

	for i in range(len(liste)):
		
		h = f[file_ite].Get("demo/" + liste[i].GetName())
		
		if "h3_" in str(liste[i].GetName()):
			h3+= h,
	
	print " done"
	
print ""
	
h_eff_num = ()
h_eff_den = ()
h_pur_num = ()
h_pur_den = ()

for file_ite in range(0,len(f)):
	
	h_eff_num+= h3[file_ite].ProjectionZ(h3[file_ite].GetName() + "en_"+str(file_ite), 1, 1, 1, 1,"e"),
	h_eff_den+= h3[file_ite].ProjectionZ(h3[file_ite].GetName() + "ed_"+str(file_ite), 1, 1, 0,-1,"e"),
	h_pur_num+= h3[file_ite].ProjectionZ(h3[file_ite].GetName() + "pn_"+str(file_ite), 1, 1, 1, 1,"e"),
	h_pur_den+= h3[file_ite].ProjectionZ(h3[file_ite].GetName() + "pd_"+str(file_ite), 0,-1, 1, 1,"e"),
	
h_eff = GetCombinedTH1FsWeightRek(h_eff_num)
h_eff_den = GetCombinedTH1FsWeightRek(h_eff_den)
h_pur = GetCombinedTH1FsWeightRek(h_pur_num)
h_pur_den = GetCombinedTH1FsWeightRek(h_pur_den)
	
h_eff.Sumw2()
h_eff_den.Sumw2()	
h_pur.Sumw2()
h_pur_den.Sumw2()
	
h_eff.Divide(h_eff_den)
h_pur.Divide(h_pur_den)

numBin = h_eff.GetXaxis().GetNbins()

for i in range(1,numBin+1):
	val = h_eff.GetBinContent(i)
	den = h_eff_den.GetBinContent(i)
	err = 0.
	if not den==0.:
		err = math.sqrt( val*(1.-val)/den )
	h_eff.SetBinError(i, err )
	
	val = h_pur.GetBinContent(i)
	den = h_pur_den.GetBinContent(i)
	err = 0.
	if not den==0.:
		err = math.sqrt( val*(1.-val)/den )
	h_pur.SetBinError(i, err )
	
			
outputname = "TrackWeightAnalysis_"

## 3d histos

pur_name = outputname + 'Purity'
eff_name = outputname + 'Efficiency'
	
h_pur.GetYaxis().SetTitle("purity")
h_eff.GetYaxis().SetTitle("efficiency")
	
h_pur.GetXaxis().SetTitle("track weight")
h_eff.GetXaxis().SetTitle("track weight")

DrawTDRStyle(histos=(h_pur,), profiles=(), additional=(), logs=(0,0,0), legend=0, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=pur_name, moreTickz=False, usePalette='default')

DrawTDRStyle(histos=(h_eff,), profiles=(), additional=(), logs=(0,0,0), legend=0, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=eff_name, moreTickz=False, usePalette='default')

print ""

