#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetListOfSubdirectories, DrawTDRStyle

	
def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1", "True", "Yes")


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:s:t:'
keywords = ['path','searchOption','tracksweight']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

#weights = ["0","01","001","0001","00001"]
#weightsS = ["0","01S","001S","0001S","00001S"]

weights = ["1_0","3_001S","z_0001"]
weightsName = ["1st","3D","z"]

path=""
searchOption=""
sqrt = False

for o,p in opts: 
    if o in ['-p','--path']:
        path = p
    if o in ['-s','--searchOption']:
        searchOption = p
    if o in ['-t','--tracksweight']:
        sqrt = str2bool(p) 
		
#if sqrt:
	#weights = weightsS

print "" 

file_help = GetListOfFilesWithEnding(path, ".root")

		
output_name = 'FinalAssociationAnalysis_'
    

#_____________________
#
# READ INPUT FILES
#_____________________


print ""
print "Read input file"
print ""

h1_v = ()
 
h_eff_r_num_v = ()
h_eff_r_den_v = ()
 
h_pur_r_num_v = ()
h_pur_r_den_v = ()

f = ()

for file_ite,file_name in enumerate(file_help):
	
	h_eff_r_num_v_f = ()
	h_eff_r_den_v_f = ()
	
	h_pur_r_num_v_f = ()
	h_pur_r_den_v_f = ()
	
	f+= ROOT.TFile.Open(file_name),
	print "File " + file_name + " ...",
	
	h1_v+= f[file_ite].Get("demo/h_NumTracks_rhoVal"),
	
	for w in weights:
		
		#h = f[file_ite].Get("demo/h3_Fin_pur_" + searchOption + "_" + w)
		h = f[file_ite].Get("demo/h3_Fin_pur_" + w)
	
		h_eff_r_num_v_f+= h.ProjectionZ("en"+w+str(file_ite), 1, 1, 1, 1,"e"),
		h_eff_r_den_v_f+= h.ProjectionZ("ed"+w+str(file_ite), 1, 1, 0,-1,"e"),
				
		h_pur_r_num_v_f+= h.ProjectionZ("pn"+w+str(file_ite), 1, 1, 1, 1,"e"),
		h_pur_r_den_v_f+= h.ProjectionZ("pd"+w+str(file_ite), 0,-1, 1, 1,"e"),
		
 
	h_eff_r_num_v+= h_eff_r_num_v_f,
	h_eff_r_den_v+= h_eff_r_den_v_f,
	
	h_pur_r_num_v+= h_pur_r_num_v_f,
	h_pur_r_den_v+= h_pur_r_den_v_f,
	
	print " done"
	
h1 = h1_v[0]

for fi in range(1,len(f)):
		
	h1.Add(h1_v[fi])
	

bn = h1.GetXaxis().GetNbins()
	
for bi in range(1,bn+1):
	bw = h1.GetBinWidth(bi)
	bc = h1.GetBinContent(bi)
	bc = h1.SetBinContent(bi, bc * 1./bw)
	
h1.GetXaxis().SetTitle("d_{xy} / cm")
h1.GetYaxis().SetTitle("number of tracks / cm")

#DrawTDRStyle(histos=(h1,), profiles=(), additional=(), logs=(1,0,0), legend=0, scale=True, zeroSuppression=False, changeLineStyle=False, drawOption='p', name="FinalAssociation_NumTracks_Vs_Rho", moreTickz=False, usePalette='default')

#sys.exit()
	
h_eff_r = ()
h_pur_r = ()

print ""
	
	
for wi in range(len(weights)):
		
	h_eff_r_num_w = h_eff_r_num_v[0][wi]
	h_eff_r_den_w = h_eff_r_den_v[0][wi]
	
	h_pur_r_num_w = h_pur_r_num_v[0][wi]
	h_pur_r_den_w = h_pur_r_den_v[0][wi]

	for fi in range(1,len(f)):
		
		h_eff_r_num_w.Add(h_eff_r_num_v[fi][wi])
		h_eff_r_den_w.Add(h_eff_r_den_v[fi][wi])
		
		h_pur_r_num_w.Add(h_pur_r_num_v[fi][wi])
		h_pur_r_den_w.Add(h_pur_r_den_v[fi][wi])

	effsum = 0.
	pursum = 0.
	esum = 0.
	psum = 0.

	numBin = h_eff_r_num_w.GetXaxis().GetNbins()

	for bi in range(0,numBin+2):
		
		e_val = h_eff_r_num_w.GetBinContent(bi)
		e_num = h_eff_r_den_w.GetBinContent(bi)
	
		effsum+= e_val
		esum+= e_num
		
		p_val = h_pur_r_num_w.GetBinContent(bi)
		p_num = h_pur_r_den_w.GetBinContent(bi)
	
		pursum+= p_val
		psum+= p_num

	pur = pursum*1./psum
	eff = effsum*1./esum
	
	print weights[wi] , ": Purity:", pur, "  - Efficiency:", eff, "  - Product:" , pur*eff
	
	h_eff_r_num_w.Divide(h_eff_r_den_w)
	h_pur_r_num_w.Divide(h_pur_r_den_w)

	numBin = h_pur_r_num_w.GetXaxis().GetNbins()
	
	for i in range(1,numBin+1):
		val = h_eff_r_num_w.GetBinContent(i)
		den = h_eff_r_den_w.GetBinContent(i)
		err = 0.
		if not den==0.:
			err = math.sqrt( val*(1.-val)/den ) + 0.0001
		h_eff_r_num_w.SetBinError(i, err )
		
		val = h_pur_r_num_w.GetBinContent(i)
		den = h_pur_r_den_w.GetBinContent(i)
		err = 0.
		if not den==0.:
			err = math.sqrt( val*(1.-val)/den ) + 0.0001
		h_pur_r_num_w.SetBinError(i, err )
	
	h_eff_r_num_w.GetXaxis().SetRangeUser(0., 25.)	
	h_pur_r_num_w.GetXaxis().SetRangeUser(0., 25.)
		
	h_eff_r+= h_eff_r_num_w,	
	h_pur_r+= h_pur_r_num_w,

h_eff_r[0].GetXaxis().SetTitle("d_{xy} / cm")
h_eff_r[0].GetYaxis().SetTitle("efficiency")

h_pur_r[0].GetXaxis().SetTitle("d_{xy} / cm")
h_pur_r[0].GetYaxis().SetTitle("purity")
	
	
#for wi in range(len(weights)):
	#h_eff_r[wi].SetName(weights[wi].split("S")[0].replace("0","0.",1))
	#h_pur_r[wi].SetName(weights[wi].split("S")[0].replace("0","0.",1))
	
for wi in range(len(weightsName)):
	h_eff_r[wi].SetName(weightsName[wi])
	h_pur_r[wi].SetName(weightsName[wi])

DrawTDRStyle(histos=h_pur_r, profiles=(), additional=(), logs=(1,0,0), legend=3, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name="FinalAssociation_Pur_Third", moreTickz=False, usePalette='default')

DrawTDRStyle(histos=h_eff_r, profiles=(), additional=(), logs=(1,0,0), legend=3, scale=False, zeroSuppression=False, changeLineStyle=False, drawOption='p', name="FinalAssociation_Eff_Third", moreTickz=False, usePalette='default')

	