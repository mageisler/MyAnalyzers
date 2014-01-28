#!/usr/bin/env python

import ROOT, getopt, sys, math
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, DrawTDRStyle

steps = ["Secondary","Final"]

#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:x:s:'
keywords = ['path','xaxis','secondaries']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

path_tag=""
xaxis=""
secondaries = ""

for o,p in opts: 
    if o in ['-p','--path']:
        path_tag = p
    if o in ['-x','--xaxis']:
        xaxis = p
    if o in ['-s','--secondaries']:
        secondaries = p

print "" 

#_____________________
#
# READ INPUT FILES
#_____________________

file_help = GetListOfFilesWithEnding(path_tag ,".root")

h1 = ()
h3 = ()

f = ()

for file_ite,file_name in enumerate(file_help):

	h1_f = ()
	h3_f = ()
	
	f+= ROOT.TFile.Open(file_name),
	print "File " + file_name + " ...",
	
	f[file_ite].cd("demo")

	liste = ROOT.gDirectory.GetListOfKeys()

	for i in range(len(liste)):
		
		if secondaries in str(liste[i].GetName()) and xaxis in str(liste[i].GetName()):
		
			h = f[file_ite].Get("demo/" +liste[i].GetName())
		
			if "h_" in str(liste[i].GetName()):
				h1_f+= h,
		
			if "h3_" in str(liste[i].GetName()):
				h3_f+= h,

	h1+= h1_f,
	h3+= h3_f,
	
	print " done"
	
print ""
	
h_eff_num = ()
h_eff_den = ()
h_pur_num = ()
h_pur_den = ()

for file_ite in range(0,len(f)):
	
	h_eff_num_f = ()
	h_eff_den_f = ()
	h_pur_num_f = ()
	h_pur_den_f = ()
	
	for sec_ite in range(len(h3[0])):
	
		h_eff_num_f+= h3[file_ite][sec_ite].ProjectionZ(h3[file_ite][sec_ite].GetName() + "en_"+str(file_ite)+str(sec_ite), 1, 1, 1, 1,"e"),
		h_eff_den_f+= h3[file_ite][sec_ite].ProjectionZ(h3[file_ite][sec_ite].GetName() + "ed_"+str(file_ite)+str(sec_ite), 1, 1, 0,-1,"e"),
		h_pur_num_f+= h3[file_ite][sec_ite].ProjectionZ(h3[file_ite][sec_ite].GetName() + "pn_"+str(file_ite)+str(sec_ite), 1, 1, 1, 1,"e"),
		h_pur_den_f+= h3[file_ite][sec_ite].ProjectionZ(h3[file_ite][sec_ite].GetName() + "pd_"+str(file_ite)+str(sec_ite), 0,-1, 1, 1,"e"),
	
	h_eff_num+= h_eff_num_f,
	h_eff_den+= h_eff_den_f,
	h_pur_num+= h_pur_num_f,
	h_pur_den+= h_pur_den_f,



h1_sum = ()

h_eff = ()
h_pur = ()
	
for sec_ite in range(len(h1[0])):

	h1_sum_s = h1[0][sec_ite]

	for file_ite in range(1,len(f)):

		h1_sum_s.Add(h1[file_ite][sec_ite])
		
	#h1_sum_s.Rebin(2)
	
	h1_sum+= h1_sum_s,
	
	
for sec_ite in range(len(h3[0])):
	
	h_eff_num_s = h_eff_num[0][sec_ite]
	h_eff_den_s = h_eff_den[0][sec_ite]
	h_pur_num_s = h_pur_num[0][sec_ite]
	h_pur_den_s = h_pur_den[0][sec_ite]

	for file_ite in range(1,len(f)):
		
		h_eff_num_s.Add(h_eff_num[file_ite][sec_ite])
		h_eff_den_s.Add(h_eff_den[file_ite][sec_ite])
		h_pur_num_s.Add(h_pur_num[file_ite][sec_ite])
		h_pur_den_s.Add(h_pur_den[file_ite][sec_ite])

	effsum = 0.
	pursum = 0.
	esum = 0.
	psum = 0.

	numBin = h_eff_num_s.GetXaxis().GetNbins()

	for bi in range(0,numBin+2):
		
		e_val = h_eff_num_s.GetBinContent(bi)
		e_num = h_eff_den_s.GetBinContent(bi)
	
		effsum+= e_val
		esum+= e_num
		
		p_val = h_pur_num_s.GetBinContent(bi)
		p_num = h_pur_den_s.GetBinContent(bi)
	
		pursum+= p_val
		psum+= p_num

	pur = pursum*1./psum
	eff = effsum*1./esum
	
	print steps[sec_ite] , ": Purity:", pur, "  - Efficiency:", eff, "  - Product:" , pur*eff
		
	#h_eff_num_s.Rebin(2)		
	#h_eff_den_s.Rebin(2)		
	#h_pur_num_s.Rebin(2)		
	#h_pur_den_s.Rebin(2)
		
	h_eff_num_s.Divide(h_eff_den_s)
	h_pur_num_s.Divide(h_pur_den_s)

	numBin = h_eff_num_s.GetXaxis().GetNbins()
	
	for i in range(1,numBin+1):
		val = h_eff_num_s.GetBinContent(i)
		den = h_eff_den_s.GetBinContent(i)
		err = 0.
		if not den==0.:
			err = math.sqrt( val*(1.-val)/den )
		h_eff_num_s.SetBinError(i, err )
		
		val = h_pur_num_s.GetBinContent(i)
		den = h_pur_den_s.GetBinContent(i)
		err = 0.
		if not den==0.:
			err = math.sqrt( val*(1.-val)/den )
		h_pur_num_s.SetBinError(i, err )
	
	h_eff+= h_eff_num_s,
	h_pur+= h_pur_num_s,
	
print ""
			
outputname = "SecondaryComparison_" + secondaries + "_"

## 1d histos = TrackNum

#h1_name = outputname + 'TrackNum_Vs_' + xaxis

#integ = 1

#for i in range(len(h1_sum)):
	
	#integ = h1_sum[i].Integral()
	
	#bn = h1_sum[i].GetXaxis().GetNbins()
	
	#for bi in range(1,bn+1):
		#bw = h1_sum[i].GetBinWidth(bi)
		#bc = h1_sum[i].GetBinContent(bi)
		#bc = h1_sum[i].SetBinContent(bi, bc * 1./bw)
	
	##h1_sum[i].GetXaxis().SetRangeUser(0., 25.)

#for i in range(len(h1_sum)):
	#h1_sum[i].Scale(1./integ)	
	
	#if "Sec" in h1_sum[i].GetName():
		#h1_sum[i].SetName("Secondary")

	#if "Fin" in h1_sum[i].GetName():
		#h1_sum[i].SetName("Final")
	
#if xaxis=="Val":
	#h1_sum[0].GetYaxis().SetTitle("number of tracks / cm")
#elif xaxis=="nMissHits":
	#h1_sum[0].GetYaxis().SetTitle("number of tracks")

l = (1,0,0)
#if "Val" in xaxis:
	#h1_sum[0].GetXaxis().SetTitle("rho / cm")
#elif "nMissHits" in xaxis:
	#h1_sum[0].GetXaxis().SetTitle("number of missed inner hits")
	#l = (0,0,0)
#else:
	#h1_sum[0].GetXaxis().SetTitle("#sigma_{rho}")

#DrawTDRStyle(histos=h1_sum, profiles=(), additional=(), logs=l, legend=0, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=h1_name, moreTickz=False, usePalette='default')

## 3d histos = Purity

pur_name = outputname + 'Purity_Vs_' + xaxis
eff_name = outputname + 'Efficiency_Vs_' + xaxis

for i in range(len(h_pur)):
	
	#h_pur[i].GetXaxis().SetRangeUser(0., 25.)
	#h_eff[i].GetXaxis().SetRangeUser(0., 25.)		
	
	if "Sec" in h_pur[i].GetName():
		h_pur[i].SetName("Secondary")
	if "Sec" in h_eff[i].GetName():
		h_eff[i].SetName("Secondary")
	
	if "Fin" in h_pur[i].GetName():
		h_pur[i].SetName("Final")
	if "Fin" in h_eff[i].GetName():
		h_eff[i].SetName("Final")
	
h_pur[0].GetYaxis().SetTitle("purity")
h_eff[0].GetYaxis().SetTitle("efficiency")

if "Val" in xaxis:
	h_pur[0].GetXaxis().SetTitle("rho / cm")
	h_eff[0].GetXaxis().SetTitle("rho / cm")
elif "nMissHits" in xaxis:
	h_pur[0].GetXaxis().SetTitle("number of missed inner hits")
	h_eff[0].GetXaxis().SetTitle("number of missed inner hits")
	l = (0,0,0)
else:
	h_pur[0].GetXaxis().SetTitle("#sigma_{rho}")
	h_eff[0].GetXaxis().SetTitle("#sigma_{rho}")

#DrawTDRStyle(histos=h_pur, profiles=(), additional=(), logs=l, legend=2, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=pur_name, moreTickz=False, usePalette='default')

DrawTDRStyle(histos=h_eff, profiles=(), additional=(), logs=l, legend=2, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name=eff_name, moreTickz=False, usePalette='default')
