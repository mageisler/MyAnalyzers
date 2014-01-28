#!/usr/bin/env python

import ROOT, getopt, sys
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, Create1DProf, GetCummulativeProf, GetMergedProf1DRek, DrawTDRStyle


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:s:'
keywords = ['path','secondaries']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

path_tag = ""
secondaries = ""
xaxis = ""

for o,p in opts: 
    if o in ['-p','--path']:
        path_tag = p 
    if o in ['-s','--secondaries']:
        secondaries = p

print "" 

file_help = GetListOfFilesWithEnding(path_tag, ".root")

px = ()
py = ()
pz = ()

ex = ()
ey = ()
ez = ()

f = ()

for file_ite,file_name in enumerate(file_help):
	
	f+= ROOT.TFile.Open(file_name),
	print "File " + file_name + " ...",
	
	f[file_ite].cd("demo")

	liste = ROOT.gDirectory.GetListOfKeys()

	for i in range(len(liste)):
		
		if secondaries in str(liste[i].GetName()) and "p3_" in str(liste[i].GetName()):
		
			help = f[file_ite].Get("demo/" +liste[i].GetName())
			
			xb = help.GetNbinsX()
			yb = help.GetNbinsY()
			zb = help.GetNbinsZ()
			
			if "pur" in str(liste[i].GetName()):
			
				px+= Create1DProf(help,0,yb,0,zb,"x"),
				py+= Create1DProf(help,0,xb,0,zb,"y"),
				pz+= Create1DProf(help,0,xb,0,yb,"z"),
			
			if "eff" in str(liste[i].GetName()):
			
				ex+= Create1DProf(help,0,yb,0,zb,"x"),
				ey+= Create1DProf(help,0,xb,0,zb,"y"),
				ez+= Create1DProf(help,0,xb,0,yb,"z"),
			
	print " done"

print ""

px_merged = GetCummulativeProf(GetMergedProf1DRek(px),reverse=False)
py_merged = GetCummulativeProf(GetMergedProf1DRek(py),reverse=False)
pz_merged = GetCummulativeProf(GetMergedProf1DRek(pz),reverse=True)

ex_merged = GetCummulativeProf(GetMergedProf1DRek(ex),reverse=False)
ey_merged = GetCummulativeProf(GetMergedProf1DRek(ey),reverse=False)
ez_merged = GetCummulativeProf(GetMergedProf1DRek(ez),reverse=True)

ex_merged.SetLineColor(1)
ey_merged.SetLineColor(1)
ez_merged.SetLineColor(1)

DrawTDRStyle(histos=(), profiles=(px_merged,ex_merged), additional=(), logs=(0,0,0), legend=0, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name = secondaries + "_x", moreTickz=False, usePalette='default')

DrawTDRStyle(histos=(), profiles=(py_merged,ey_merged), additional=(), logs=(0,0,0), legend=0, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name = secondaries + "_y", moreTickz=False, usePalette='default')

DrawTDRStyle(histos=(), profiles=(pz_merged,ez_merged), additional=(), logs=(0,0,0), legend=0, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name = secondaries + "_z", moreTickz=False, usePalette='default')


