#!/usr/bin/env python

import ROOT, getopt, sys
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetMergedProf1DRek, DrawTDRStyle


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

hs_c = ()
hr_c = ()
pur_c = ()
eff_c = ()

hs_u = ()
hr_u = ()
pur_u = ()
eff_u = ()

f = ()

for file_ite,file_name in enumerate(file_help):
	
	f+= ROOT.TFile.Open(file_name),
	print "File " + file_name + " ...",
	
	f[file_ite].cd("clean")

	liste = ROOT.gDirectory.GetListOfKeys()

	for i in range(len(liste)):
		
		if secondaries in str(liste[i].GetName()):
		
			help = f[file_ite].Get("clean/" +liste[i].GetName())
		
			if "_numSV_" in str(liste[i].GetName()):
				hr_c+= help,
		
			if "_numTV_" in str(liste[i].GetName()):
				hs_c+= help,
		
			if "_pur_" in str(liste[i].GetName()):
				pur_c+= help,
				
			if "_eff_" in str(liste[i].GetName()):
				eff_c+= help,

	for i in range(len(liste)):
		
		if secondaries in str(liste[i].GetName()):
		
			help = f[file_ite].Get("unclean/" +liste[i].GetName())
		
			if "_numSV_" in str(liste[i].GetName()):
				hr_u+= help,
		
			if "_numTV_" in str(liste[i].GetName()):
				hs_u+= help,
		
			if "_pur_" in str(liste[i].GetName()):
				pur_u+= help,
				
			if "_eff_" in str(liste[i].GetName()):
				eff_u+= help,
	
	print " done"
	
print ""

hs_c_add = hs_c[0]
hr_c_add = hr_c[0]	

hs_u_add = hs_u[0]
hr_u_add = hr_u[0]	

for fi in range(len(f)):

	hs_c_add.Add(hs_c[fi])
	hr_c_add.Add(hr_c[fi])

	hs_u_add.Add(hs_u[fi])
	hr_u_add.Add(hr_u[fi])
	
h1_c = (hs_c_add,hr_c_add)
h1_u = (hs_u_add,hr_u_add)

pur_c_add = GetMergedProf1DRek(pur_c)
eff_c_add = GetMergedProf1DRek(eff_c)

pur_u_add = GetMergedProf1DRek(pur_u)
eff_u_add = GetMergedProf1DRek(eff_u)

##########
			
outputname = "SecondaryValidator_" + secondaries + "_"

## 1d histos = TrackNum

h1_name = outputname + 'TrackNum'

integ = 1

for i in range(len(h1_c)):
	
	bn = h1_c[i].GetXaxis().GetNbins()
	
	for bi in range(1,bn+1):
		bw = h1_c[i].GetBinWidth(bi)
		bc = h1_c[i].GetBinContent(bi)
		bc = h1_c[i].SetBinContent(bi, bc * 1./bw)
		
		bw = h1_u[i].GetBinWidth(bi)
		bc = h1_u[i].GetBinContent(bi)
		bc = h1_u[i].SetBinContent(bi, bc * 1./bw)
	
	#h1_c[i].GetXaxis().SetRangeUser(0., 25.)
	#h1_u[i].GetXaxis().SetRangeUser(0., 25.)
	
	h1_c[i].Scale( 1./h1_u[i].Integral() )
	h1_u[i].Scale( 1./h1_u[i].Integral() )

h1_c[1].SetName("Filtered")
h1_u[1].SetName("Unfiltered")
h1_u[0].SetName("Simulated")
	
h1_c[1].GetYaxis().SetTitle("number of secondary vertices / cm")
h1_u[1].GetYaxis().SetTitle("number of secondary vertices / cm")

h1_c[1].GetXaxis().SetTitle("rho / cm")
h1_u[1].GetXaxis().SetTitle("rho / cm")

#DrawTDRStyle(histos=(h1_u[1],h1_u[0],h1_c[1]), profiles=(), additional=(), logs=(1,0,0), legend=2, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=h1_name, moreTickz=False, usePalette='default')

#sys.exit()

pur_name = outputname + 'Purity'
eff_name = outputname + 'Efficiency'

pur_c_add.SetName("Filtered")
pur_u_add.SetName("Unfiltered")

eff_c_add.SetName("Filtered")
eff_u_add.SetName("Unfiltered")
	
pur_c_add.GetYaxis().SetTitle("purity")
pur_u_add.GetYaxis().SetTitle("purity")

pur_c_add.GetXaxis().SetTitle("rho / cm")
pur_u_add.GetXaxis().SetTitle("rho / cm")
	
eff_c_add.GetYaxis().SetTitle("efficiency")
eff_u_add.GetYaxis().SetTitle("efficiency")

eff_c_add.GetXaxis().SetTitle("rho / cm")
eff_u_add.GetXaxis().SetTitle("rho / cm")
	
pur_c_add.GetXaxis().SetRangeUser(0.05, 150.)
pur_u_add.GetXaxis().SetRangeUser(0.05, 150.)
	
eff_c_add.GetXaxis().SetRangeUser(0.05, 150.)
eff_u_add.GetXaxis().SetRangeUser(0.05, 150.)

DrawTDRStyle(histos=(pur_c_add,pur_u_add), profiles=(), additional=(), logs=(1,0,0), legend=2, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=pur_name, moreTickz=False, usePalette='default')

#DrawTDRStyle(histos=(eff_c_add,eff_u_add), profiles=(), additional=(), logs=(1,0,0), legend=2, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=eff_name, moreTickz=False, usePalette='default')
