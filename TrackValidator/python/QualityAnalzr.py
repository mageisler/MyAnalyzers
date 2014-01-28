#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetCummulativeHisto, GetCummulativeYaxisHisto, GetCombinedTH2FsWeightRek, CreateFraction2DHistosRek, DrawTDRStyle, GetXYrelativeDistance, GetHistoFromProf, GetCummulativeProf


#_____________________
#
# READ INPUT PARAMETERS
#_____________________


letters = 'p:s:a:'
keywords = ['path','step','association']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)


path_tag =""
step = ""
asstag = "_0"

steps = ["tw","sc","fa"]
analyzers = ["qualanlzrFV","qualanlzr3D","qualanlzrZ"]

for o,p in opts: 
    if o in ['-p','--path']:
        path_tag = p 
    if o in ['-a','--association']:
        asstag = "_" + p 
    if o in ['-s','--step']:
        steps.append(p)

print "" 

file_help = GetListOfFilesWithEnding(path_tag ,".root")
    

#_____________________
#
# READ INPUT FILES
#_____________________


f = ()

h1_num = ()
h1_den = ()

h2_num = ()
h2_den = ()

for file_ite,file_name in enumerate(file_help):

	h1_num_f = ()
	h1_den_f = ()

	h2_num_f = ()
	h2_den_f = ()
	
	f+= ROOT.TFile.Open(file_name),
	print "File " + file_name + " ...",
	
	f[file_ite].cd("qualanlzrFV")

	liste = ROOT.gDirectory.GetListOfKeys()
	
	for analzr in analyzers:

		h1_num_f_a = ()
		h1_den_f_a = ()

		h2_num_f_a = ()
		h2_den_f_a = ()
	
		for step in steps:

			for i in range(len(liste)):
		
				if asstag in str(liste[i].GetName()) and step in str(liste[i].GetName()):
		
					h = f[file_ite].Get(analzr + "/" +liste[i].GetName())
		
					if "h_" in str(liste[i].GetName()) and "matched" in str(liste[i].GetName()):
						h.SetName("h_matched_"+str(file_ite)+"_"+analzr+"_"+step)
						h1_num_f_a+= h.Rebin(4),		
					if "h_" in str(liste[i].GetName()) and "all" in str(liste[i].GetName()):
						h.SetName("h_all_"+str(file_ite)+"_"+analzr+"_"+step)
						h1_den_f_a+= h.Rebin(4),
		
					if "h2_" in str(liste[i].GetName()) and "matched_a1distance_a2distance" in str(liste[i].GetName()):
						h.SetName("h2_matched_"+str(file_ite)+"_"+analzr+"_"+step)
						h2_num_f_a+= h,		
					if "h2_" in str(liste[i].GetName()) and "all_a1distance_a2distance" in str(liste[i].GetName()):	
						h.SetName("h2_all_"+str(file_ite)+"_"+analzr+"_"+step)
						h2_den_f_a+= h,
						
		h1_num_f+= h1_num_f_a,
		h1_den_f+= h1_den_f_a,
	
		h2_num_f+= h2_num_f_a,
		h2_den_f+= h2_den_f_a,	

	h1_num+= h1_num_f,
	h1_den+= h1_den_f,
	
	h2_num+= h2_num_f,
	h2_den+= h2_den_f,
	
	print " done"
	
print ""

## old: file - ana - step
## new: ana - step - file

h2_num_c = ()
h2_den_c = ()

for ana_ite in range(len(h2_num[0])):

	h2_num_c_a = ()
	h2_den_c_a = ()

	for step_ite in range(len(h2_num[0][0])): 

		h2_num_c_a_s = ()
		h2_den_c_a_s = ()

		for file_ite in range(len(f)):

			h2_num_c_a_s+= h2_num[file_ite][ana_ite][step_ite],
			h2_den_c_a_s+= h2_den[file_ite][ana_ite][step_ite], 

		h2_num_c_a+= h2_num_c_a_s,
		h2_den_c_a+= h2_den_c_a_s,

	h2_num_c+= h2_num_c_a,
	h2_den_c+= h2_den_c_a,
	
	
## new: ana - step - file

h1_pur = ()
h2_pur = ()

for ana_ite in range(len(h1_num[0])):

	h1_pur_a = ()
	h2_pur_a = ()

	for step_ite in range(len(h1_num[0][0])):

		h1_pur_i = h1_num[0][ana_ite][step_ite]
		h1_den_s = h1_den[0][ana_ite][step_ite]

		for file_ite in range(1,len(f)):

			h1_pur_i.Add(h1_num[file_ite][ana_ite][step_ite])
			h1_den_s.Add(h1_den[file_ite][ana_ite][step_ite])
	
		h1_pur_cum = GetCummulativeHisto(h1_pur_i)
		h1_den_cum = GetCummulativeHisto(h1_den_s)
	
		h1_pur_cum.Divide(h1_den_cum)

		numBin = h1_pur_cum.GetXaxis().GetNbins()
		
		for i in range(1,numBin+1):		
			val = h1_pur_cum.GetBinContent(i)
			den = h1_den_cum.GetBinContent(i)
			err = 0.
			if not den==0.:
				err = math.sqrt( val*(1.-val)/den )
			h1_pur_cum.SetBinError(i, err )
	
		h1_pur_a+= h1_pur_cum,

		h2_num_s = GetCombinedTH2FsWeightRek(h2_num_c[ana_ite][step_ite])
		h2_den_s = GetCombinedTH2FsWeightRek(h2_den_c[ana_ite][step_ite])

		h2_num_cum = GetCummulativeYaxisHisto(h2_num_s)
		h2_den_cum = GetCummulativeYaxisHisto(h2_den_s)

		h2_pur_cum = CreateFraction2DHistosRek( (h2_den_cum,h2_num_cum) )
	
		h2_pur_a+= h2_pur_cum[0],

	h1_pur+= h1_pur_a,
	h2_pur+= h2_pur_a,
	
	
## new: ana - step - file
	
p1_pur = ()	
h1_sum = ()

for i1 in range(len(h2_pur)):

	p1_pur_a = ()
	h1_sum_a = ()

	for i2 in range(len(h2_pur[0])):
	
		p1_pur_a+= GetCummulativeProf(GetXYrelativeDistance(h2_pur[i1][i2])),
		h1_help = GetHistoFromProf(p1_pur_a[i2])
		if not  h1_help.Integral()==0:
			h1_help.Scale(1./ h1_help.Integral())
		h1_sum_a+= h1_help,

	p1_pur+= p1_pur_a,
	h1_sum+= h1_sum_a,
    

#________________________
#
# DRAW THEM ALL	
#________________________

line90 = ROOT.TLine(0.0015,0.9,45.,0.9)
line90.SetLineColor(921)
line90.SetLineWidth(2)
line70 = ROOT.TLine(0.0015,0.7,45.,0.7)
line70.SetLineColor(921)
line70.SetLineWidth(2)
line50 = ROOT.TLine(0.0015,0.5,45.,0.5)
line50.SetLineColor(921)
line50.SetLineWidth(2)
line30 = ROOT.TLine(-0.98,0.3,1.98,0.3)
line30.SetLineColor(921)
line30.SetLineWidth(2)

for i1 in range(len(h1_pur)):

	h1_pur[i1][0].GetYaxis().SetTitle("cumulative purity")

	h2_pur[i1][0].GetXaxis().SetTitle("distance^{first} / cm")
	h2_pur[i1][0].GetYaxis().SetTitle("distance^{second} / cm")
	h2_pur[i1][0].GetZaxis().SetTitle("purity")

	p1_pur[i1][0].GetXaxis().SetTitle("relative distance")
	p1_pur[i1][0].GetYaxis().SetTitle("cumulative purity")

	for i in range(len(h1_pur[i1])):
		
		if "tw" in h1_pur[i1][i].GetName():
			h1_pur[i1][i].SetName("Step 1")
		
		if "sc" in h1_pur[i1][i].GetName():
			h1_pur[i1][i].SetName("Step 2")
		
		if "fa" in h1_pur[i1][i].GetName():
			h1_pur[i1][i].SetName("Step 3 - 1st")
		
		if "QA3" in h1_pur[i1][i].GetName():
			h1_pur[i1][i].SetName("Step 3 - 3D")
		
		if "QAZ" in h1_pur[i1][i].GetName():
			h1_pur[i1][i].SetName("Step 3 - z")
	
	for i in range(len(h2_pur[i1])):
		
		if "tw" in h2_pur[i1][i].GetName():
			h2_pur[i1][i].SetName("Step 1")
		
		if "sc" in h2_pur[i1][i].GetName():
			h2_pur[i1][i].SetName("Step 2")
		
		if "fa" in h2_pur[i1][i].GetName():
			h2_pur[i1][i].SetName("Step 3 - 1st")
		
		if "QA3" in h2_pur[i1][i].GetName():
			h2_pur[i1][i].SetName("Step 3 - 3D")
		
		if "QAZ" in h2_pur[i1][i].GetName():
			h2_pur[i1][i].SetName("Step 3 - z")
	
	for i in range(len(p1_pur[i1])):
		
		if "tw" in p1_pur[i1][i].GetName():
			p1_pur[i1][i].SetName("Step 1")
			h1_sum[i1][i].SetName("Step 1")
		
		if "sc" in p1_pur[i1][i].GetName():
			p1_pur[i1][i].SetName("Step 2")
			h1_sum[i1][i].SetName("Step 2")
		
		if "fa" in p1_pur[i1][i].GetName():
			p1_pur[i1][i].SetName("Step 3 - 1st")
			h1_sum[i1][i].SetName("Step 3 - 1st")
		
		if "QA3" in p1_pur[i1][i].GetName():
			p1_pur[i1][i].SetName("Step 3 - 3D")
			h1_sum[i1][i].SetName("Step 3 - 3D")
		
		if "QAZ" in p1_pur[i1][i].GetName():
			p1_pur[i1][i].SetName("Step 3 - z")
			h1_sum[i1][i].SetName("Step 3 - z")
			
outputname = "QualityAnalyzer" + asstag + "_"
h1_name = outputname + "1D"
h2_name = outputname + "2D"

#DrawTDRStyle(histos=h1_pur, profiles=(), additional=(line90,line70,line50), logs=(1,0,0), legend=3, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=h1_name, moreTickz=False, usePalette='default')

#sys.exit()

for i1 in range(len(h1_pur)):
	DrawTDRStyle(histos=(), profiles=(p1_pur[i1]), additional=(line30,), logs=(0,0,0), legend=3, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=h2_name+"_"+ analyzers[i1].split("qualanlzr")[1]+"_Profiles", moreTickz=False, usePalette='default')

	DrawTDRStyle(histos=(h1_sum[i1]), profiles=(), additional=(), logs=(0,0,0), legend=3, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=h2_name+"_"+ analyzers[i1].split("qualanlzr")[1]+"_Histos", moreTickz=False, usePalette='default')

