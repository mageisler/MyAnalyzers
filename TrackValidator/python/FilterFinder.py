#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os

sys.path.append('/.automount/home/home__home2/institut_3b/geisler/Phd-Study/CMSSW/Helpers/')
from AnalysisTools import *


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:'
keywords = ['path']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)



path=""
name="nuclearInteraction"

for o,p in opts: 
    if o in ['-p','--path']:
        path = p


file_list = GetListOfFilesWithEnding(path,".root")
file_number = len(file_list)
#file_number = 1

if file_number<1:
    print "Too less files"
    sys.exit()
    
for dataset in possibleDatasets:
    if dataset in file_list[0]:
	name = dataset

#_____________________
#
# READ INPUT FILES
#_____________________

dir_list = GetListOfSubdirectories(file_list[0],"","demo")
dir_number = len(dir_list)
#dir_number = 1

zaxis_list = GetListOfKeyWords(file_list[0],dir_list[0],name)
zaxis_number = len(zaxis_list)
#zaxis_number = 1
	
print ""

File_ref = ()
histo_queue_2D = ()

for file_ite in range(file_number):

    File_ref+= ROOT.TFile.Open(file_list[file_ite]),
    
    print "File " + str(file_ite+1) + ": " + str(file_list[file_ite]) + " ...",

    histo_2D_dir = ()
    
    for dir_ite in range(dir_number):
 
        histo_2D_dir_zaxis = ()       
	
	for zaxis_ite in range(zaxis_number):
		
	    help = int(random.random()*1000)
	 
	    histo_help = File_ref[file_ite].Get(dir_list[dir_ite] + zaxis_list[zaxis_ite])
	    
	    zbin_number = histo_help.GetNbinsZ()

            histo_2D_dir_zaxis_zaxisrange = () 	        
	
	    for zbin_ite in range(1,zbin_number+1):
	    
	        histo_help.GetZaxis().SetRange(zbin_ite,zbin_ite)
	        histo_help.SetName(dir_list[dir_ite] + zaxis_list[zaxis_ite] + str(zbin_ite) + str(help))
	                
	        histo_2D_dir_zaxis_zaxisrange+= histo_help.Project3D("yxo"),
	    
	    histo_2D_dir_zaxis+= histo_2D_dir_zaxis_zaxisrange,
	    
	histo_2D_dir+= histo_2D_dir_zaxis,
	
    histo_queue_2D+= histo_2D_dir,

    print "done" 

print ""

#_____________________________________________
#
# REARRANGE TUPLE SO THAT FILE IS AT THE BACK
#_____________________________________________   
    
##the former order is file - dir - zaxis (- zaxisbin)
##the new order is dir - zaxis (- zaxisbin) - file

print "Rearranging"
     
histo_queue_2D_rearranged = ()  
	
for dir_ite in range(dir_number):
	    
    histo_queue_2D_rearranged_zaxis = ()
    
    for zaxis_ite in range(zaxis_number):
    
        histo_queue_2D_rearranged_zaxis_zaxisrange = () 
	    
	zbin_number = len(histo_queue_2D[0][dir_ite][zaxis_ite])
	
	for zbin_ite in range(zbin_number):
    
            histo_queue_2D_rearranged_zaxis_zaxisrange_file = ()
     
            for file_ite in range(file_number):  
			
                histo_queue_2D_rearranged_zaxis_zaxisrange_file+= histo_queue_2D[file_ite][dir_ite][zaxis_ite][zbin_ite],
			
            histo_queue_2D_rearranged_zaxis_zaxisrange+= histo_queue_2D_rearranged_zaxis_zaxisrange_file,
			
        histo_queue_2D_rearranged_zaxis+= histo_queue_2D_rearranged_zaxis_zaxisrange,
	
    histo_queue_2D_rearranged+= histo_queue_2D_rearranged_zaxis,

print ""
 

#________________________________
#
# COMBINE THE HISTOS WITH WEIGHT
#_______________________________   
    
##the former order is dir - zaxis - zaxisbin - file
##the new order is dir - zaxis - zaxisbin
	
print "Combining TH2Fs with weight..."

histo_com_2D = GetCombinedTH2FsWeightRek(histo_queue_2D_rearranged,file_list)

print ""
 

#________________________________________
#
# CALCULATE OVERALL, PURITY, EFFICIENCY
#________________________________________
    
    
##the former order is dir - zaxis - zaxisbin
##the new order is dir - zaxis - zaxisbin

histos_fracs = () 
histos_refs = ()  
	
for dir_ite in range(dir_number):
	
    print "Directory ", dir_list[dir_ite], ":"
	    
    histos_fracs_zaxis = ()
    histos_refs_zaxis = ()
    
    for zaxis_ite in range(zaxis_number):
	    
        histos_fracs_zaxis_zaxisrange = ()
     
        print "Z-Axis ", zaxis_list[zaxis_ite]
	    
	all_right = 0.
	all_entries = 0.
	    
	all_signal_right = 0.
	all_signal_entries = 0.
	    
	all_effic_num = 0.
	all_effic_denum = 0.
	 
	histo_help = File_ref[0].Get(dir_list[0] + zaxis_list[zaxis_ite])
	    
	zbin_number = len(histo_com_2D[dir_ite][zaxis_ite])	
	zaxis_title = histo_help.GetZaxis().GetTitle()
	
	new_min = 0
	
	for zbin_ite in range(zbin_number):
		
	    if ( (histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetEntries()==0) and (zbin_ite==new_min) ):
		new_min+= 1
		
	new_max = zbin_number-1
	
	for zbin_ite in range(zbin_number-1,new_min,-1):
		
	    if ( (histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetEntries()==0) and (zbin_ite==new_max) and (zbin_ite>new_min+1)):
		new_max= zbin_ite-1
		
	zaxis_min = histo_help.GetZaxis().GetBinUpEdge(new_min)
	zaxis_max = histo_help.GetZaxis().GetBinUpEdge(new_max)
	zbin_number_new = new_max-new_min
	
	histo_help = ROOT.TH1F(dir_list[dir_ite]+zaxis_list[zaxis_ite]+"1D_Over","FF Overall; " + zaxis_title + "; fraction of correctly ass. tracks", zbin_number_new, zaxis_min, zaxis_max)
	
	histo_ref_Over = ROOT.TH1F(dir_list[dir_ite]+zaxis_list[zaxis_ite]+"1D_Ref_Over","", zbin_number_new, zaxis_min, zaxis_max)
	
	histo_help_Sig = ROOT.TH1F(dir_list[dir_ite]+zaxis_list[zaxis_ite]+"1D_Sig","FF Purity; " + zaxis_title + "; purity", zbin_number_new, zaxis_min, zaxis_max)
	
	histo_ref_Sig = ROOT.TH1F(dir_list[dir_ite]+zaxis_list[zaxis_ite]+"1D_Ref_Sig","", zbin_number_new, zaxis_min, zaxis_max)
	
	histo_effic_sig = ROOT.TH1F(dir_list[dir_ite]+zaxis_list[zaxis_ite]+"1D_Effic","FF Efficiency; " + zaxis_title + "; efficiency", zbin_number_new, zaxis_min, zaxis_max)
	
	histo_ref_eff = ROOT.TH1F(dir_list[dir_ite]+zaxis_list[zaxis_ite]+"1D_Ref_Eff","", zbin_number_new, zaxis_min, zaxis_max)
	
	for zbin_ite in range(new_min,new_max):
	    
	    right_ones = histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(1,1) + histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(2,2)
	    
	    just_signal_right = histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(1,1)
	    just_signal = histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(1,1) + histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(2,1) + histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(3,1)
	    
	    effic_num = histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(1,1)
	    effic_denum = histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(1,1) + histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetBinContent(1,2)
	    
	    all_right+= right_ones
	    all_entries+= histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetEntries()
	    
	    all_signal_right+= just_signal_right
	    all_signal_entries+= just_signal
	    
	    all_effic_num+= effic_num
	    all_effic_denum+= effic_denum
	
	    if histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetEntries()==0:
		frac = 1.1
		frac_err = 100.
	    else:
		frac = right_ones * 1./ histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetEntries()
		frac_err = math.sqrt( frac*(1.-frac)/histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetEntries() )
	
	    if just_signal==0:
		frac_sig = 1.1
		frac_sig_err = 100.
	    else:
		frac_sig = just_signal_right * 1./ just_signal
		frac_sig_err = math.sqrt( frac_sig*(1.-frac_sig)/just_signal )
	
	    if effic_denum==0:
		effic = 1.1
		effic_err = 100.
	    else:
		effic = effic_num * 1./ effic_denum
		effic_err = math.sqrt( effic*(1.-effic)/effic_denum )
		
	    histo_help.SetBinContent(zbin_ite+1-new_min,frac)
	    histo_help.SetBinError(zbin_ite+1-new_min,frac_err)
		
	    histo_ref_Over.SetBinContent(zbin_ite+1-new_min,histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetEntries())
	    #histo_ref_Over.SetBinError(zbin_ite+1-new_min,math.sqrt.histo_com_2D[dir_ite][zaxis_ite][zbin_ite].GetEntries()))
	    
	    
	    histo_help_Sig.SetBinContent(zbin_ite+1-new_min,frac_sig)
	    histo_help_Sig.SetBinError(zbin_ite+1-new_min,frac_sig_err)
	    
	    histo_ref_Sig.SetBinContent(zbin_ite+1-new_min,just_signal)
	    #histo_ref_Sig.SetBinError(zbin_ite+1-new_min,math.sqrt(just_signal))
	    
	    
	    histo_effic_sig.SetBinContent(zbin_ite+1-new_min,effic)
	    histo_effic_sig.SetBinError(zbin_ite+1-new_min,effic_err) 
	    
	    histo_ref_eff.SetBinContent(zbin_ite+1-new_min,effic_denum)
	    #histo_ref_eff.SetBinError(zbin_ite+1-new_min,math.sqrt(effic_denum)) 
	
	histo_all = ()
	histo_all+= histo_help,
	histo_all+= histo_ref_Over,
	histo_all+= histo_help_Sig,
	histo_all+= histo_ref_Sig,
	histo_all+= histo_effic_sig,
	histo_all+= histo_ref_eff,
	
	RebinTH1Fs(histo_all,5,30)
	    
	histos_fracs_zaxis+= histo_all[0],
	histos_fracs_zaxis+= histo_all[2],
	histos_fracs_zaxis+= histo_all[4],
	    
	histos_refs_zaxis+= histo_all[1],
	histos_refs_zaxis+= histo_all[3],
	histos_refs_zaxis+= histo_all[5],
	
    histos_fracs+= histos_fracs_zaxis,
    histos_refs+= histos_refs_zaxis,

print ""


#_______________
#
# DRAW IT ALL	
#_______________
    
##the order is dir - zaxis
    
    
canv = ROOT.TCanvas("canv","Guck mal",10,10,1200,850)
ROOT.gStyle.SetNumberContours(255)
ROOT.gStyle.SetTitleX(0.58)
ROOT.gStyle.SetTitleY(0.93)
ROOT.gStyle.SetTitleW(1.)
ROOT.gStyle.SetOptStat(0)

canv.Divide(3,zaxis_number)

ROOT.gPad.SetGrid()
ROOT.gPad.SetBorderMode(0)
ROOT.gPad.SetFillColor(0)

for zaxis_ite in range(zaxis_number):
	
    for ite in range(3):
	    
	pad = (zaxis_ite*3)+ite
		    
        canv.cd(pad+1)
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetLeftMargin(0.2)
        ROOT.gPad.SetTopMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetFillColor(0)
	
        histos_fracs[0][pad].GetXaxis().CenterTitle()
        histos_fracs[0][pad].GetYaxis().CenterTitle()
   
        histos_fracs[0][pad].GetXaxis().SetTitleOffset(1.1)
        histos_fracs[0][pad].GetYaxis().SetTitleOffset(1.1)
        histos_fracs[0][pad].GetYaxis().SetRangeUser(0.,1.5)
	
        color = my_colors[0]
        histos_fracs[0][pad].SetLineColor(color)
        histos_fracs[0][pad].SetLineWidth(2)

        histos_fracs[0][pad].Draw()
	
	histos_refs[0][pad].Scale(1./histos_refs[0][pad].GetMaximum())
	
        histos_refs[0][pad].SetLineColor(922)
        histos_refs[0][pad].SetLineWidth(2)

        histos_refs[0][pad].Draw("same")
	
	
canv.SaveAs("../root/FilterFinder_"+name+".root")
canv.SaveAs("../pdfs/FilterFinder_"+name+".pdf")
canv.SaveAs("../pictures/FilterFinder_"+name+".png")