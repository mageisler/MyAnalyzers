#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os
sys.path.append('/.automount/home/home__home2/institut_3b/geisler/Phd-Study/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetListOfSubdirectories, GetMergedProf3DRek, npu_ranges, pt_ranges, eta_ranges,  Create1DProf, npu_namesP, pt_namesP, eta_namesP, DrawProfs, DrawProfs2D


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:k:'
keywords = ['path','keyword']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)


    
output_name = 'JetValidator_efficiency'

path=""
histo_name=""

for o,p in opts: 
    if o in ['-p','--path']:
        path = p 
    if o in ['-k','--keyword']:
        histo_name = p
	
histo_name="p_effic"
    

#_____________________
#
# READ INPUT FILES
#_____________________

# order is file - directory - association

print "Read input file"


file_names = GetListOfFilesWithEnding(path,".root")
file_number = len(file_names)

if file_number<1:
    print "Too less files"
    sys.exit()
dir_names = GetListOfSubdirectories(file_names[0],"","jetvalidator")
dir_number = len(dir_names)

ass_names = ()

for i in range(dir_number):
	
    ass_names+=GetListOfSubdirectories(file_names[0],dir_names[i],"pfJets"),
    
ass_number = len(ass_names[0])
	
print ""

File_ref = ()
prof_effic_queue_in = ()

for file_ite in range(file_number):
	
    File_ref+= ROOT.TFile.Open(file_names[file_ite]),    
    print "File " + str(file_ite+1)+ ": " + str(file_names[file_ite]) + " ...",
    
    prof_effic_file = () 
       
    for dir_ite in range(dir_number):

        prof_effic_file_dir = ()
       
        for ass_ite in range(ass_number):
	
	    help = int(random.random() *1000)    

            prof_help = File_ref[file_ite].Get(ass_names[dir_ite][ass_ite]+histo_name)
            prof_help.SetName(histo_name+"_"+str(help))
            prof_effic_file_dir+= prof_help,

        prof_effic_file+= prof_effic_file_dir,

    prof_effic_queue_in+= prof_effic_file,
    
    print " done"
    
print ""
    

#_____________________
#
# CHANGE ORDER	
#_____________________

# old order is file - directory - association
# new order is directory - association - file

print "Change order ...",

prof_effic_queue = ()

for dir_ite in range(dir_number):
    
    prof_effic_dir = () 
       
    for ass_ite in range(ass_number):

        prof_effic_dir_ass = ()
       
        for file_ite in range(file_number):

            prof_effic_dir_ass+= prof_effic_queue_in[file_ite][dir_ite][ass_ite],

        prof_effic_dir+= prof_effic_dir_ass,

    prof_effic_queue+= prof_effic_dir,
    
print " done"
print ""
    

#_____________________
#
# Merge files	
#_____________________

# old order is directory - association - file
# new order is directory - association

print "Merge files ...",

prof_effic_merged = GetMergedProf3DRek(prof_effic_queue)
    
print " done"
print ""
    

#________________________
#
# CREATE SEPERATED PLOTS	
#________________________

# old order is directory - association
# new order is directory - association - axis1 - axis2

print "Create seperated plots ...",

prof_effic_eta = ()
prof_effic_pt  = ()
prof_effic_npu = ()

for dir_ite in range(dir_number):
    
    prof_effic_eta_dir = () 
    prof_effic_pt_dir = ()
    prof_effic_npu_dir = ()  
       
    for ass_ite in range(ass_number):
    
        prof_effic_eta_dir_ass = ()
        prof_effic_pt_dir_ass = ()
        prof_effic_npu_dir_ass = ()

        ##eta && pt
        for npu_ite in range(len(npu_ranges)):
            if npu_ite==(len(npu_ranges)-1):
                npu_min = npu_ranges[0]
                npu_max = npu_ranges[len(npu_ranges)-1]
            else:
                npu_min = npu_ranges[npu_ite]
                npu_max = npu_ranges[npu_ite+1]
			
            ##eta    
            prof_effic_eta_dir_ass_npu = ()	 
  
            for pt_ite in range(len(pt_ranges)):
                if pt_ite==(len(pt_ranges)-1):
                    pt_min = pt_ranges[0]
                    pt_max = pt_ranges[len(pt_ranges)-1]
                else:
                    pt_min = pt_ranges[pt_ite]
                    pt_max = pt_ranges[pt_ite+1]
	    	    
	        prof_effic_eta_dir_ass_npu+= Create1DProf(prof_effic_merged[dir_ite][ass_ite], pt_min,pt_max, npu_min,npu_max,"x"),
	    	    
	    prof_effic_eta_dir_ass+= prof_effic_eta_dir_ass_npu,
    
            ##pt
            prof_effic_pt_dir_ass_npu = ()	 
  
            for eta_ite in range(len(eta_ranges)):
                if eta_ite==(len(eta_ranges)-1):
                    eta_min = eta_ranges[0]
                    eta_max = eta_ranges[len(eta_ranges)-1]
                else:
                    eta_min = eta_ranges[eta_ite]
                    eta_max = eta_ranges[eta_ite+1]
	    	    
	        prof_effic_pt_dir_ass_npu+= Create1DProf(prof_effic_merged[dir_ite][ass_ite], eta_min,eta_max, npu_min,npu_max,"y"),
	    	    
	    prof_effic_pt_dir_ass+= prof_effic_pt_dir_ass_npu,

        ##npu
        for eta_ite in range(len(eta_ranges)):
	    if eta_ite==(len(eta_ranges)-1):
		eta_min = eta_ranges[0]
		eta_max = eta_ranges[len(eta_ranges)-1]
            else:
		eta_min = eta_ranges[eta_ite]
		eta_max = eta_ranges[eta_ite+1]
			                
            prof_effic_npu_dir_ass_eta = ()	 
  
            for pt_ite in range(len(pt_ranges)):
                if pt_ite==(len(pt_ranges)-1):
                    pt_min = pt_ranges[0]
                    pt_max = pt_ranges[len(pt_ranges)-1]
                else:
                    pt_min = pt_ranges[pt_ite]
                    pt_max = pt_ranges[pt_ite+1]
	    	    
	        prof_effic_npu_dir_ass_eta+= Create1DProf(prof_effic_merged[dir_ite][ass_ite], eta_min,eta_max, pt_min,pt_max,"z"),
	    	    
	    prof_effic_npu_dir_ass+= prof_effic_npu_dir_ass_eta,
	    	    
	prof_effic_eta_dir+= prof_effic_eta_dir_ass,
	prof_effic_pt_dir+= prof_effic_pt_dir_ass,
	prof_effic_npu_dir+= prof_effic_npu_dir_ass,
	    	    
    prof_effic_eta+= prof_effic_eta_dir,
    prof_effic_pt+= prof_effic_pt_dir,
    prof_effic_npu+= prof_effic_npu_dir,
    
    
print " done"
print ""
    

#_____________________
#
# CHANGE ORDER	
#_____________________

# old order is directory - association - axis1 - axis2
# new order is directory - axis1 - axis2 - association 

print "Change order again...",

##eta
prof_effic_eta_final = ()

for dir_ite in range(dir_number):
    
    prof_effic_eta_final_dir = () 
       
    for axis1_ite in range(len(prof_effic_eta[0][0])):

        prof_effic_eta_final_dir_axis1 = ()
       
        for axis2_ite in range(len(prof_effic_eta[0][0][axis1_ite])):

            prof_effic_eta_final_dir_axis1_axis2 = ()
       
            for ass_ite in range(ass_number):

                prof_effic_eta_final_dir_axis1_axis2+= prof_effic_eta[dir_ite][ass_ite][axis1_ite][axis2_ite],

            prof_effic_eta_final_dir_axis1+= prof_effic_eta_final_dir_axis1_axis2,

        prof_effic_eta_final_dir+= prof_effic_eta_final_dir_axis1,

    prof_effic_eta_final+= prof_effic_eta_final_dir,

##pt
prof_effic_pt_final = ()

for dir_ite in range(dir_number):
    
    prof_effic_pt_final_dir = () 
       
    for axis1_ite in range(len(prof_effic_pt[0][0])):

        prof_effic_pt_final_dir_axis1 = ()
       
        for axis2_ite in range(len(prof_effic_pt[0][0][axis1_ite])):

            prof_effic_pt_final_dir_axis1_axis2 = ()
       
            for ass_ite in range(ass_number):

                prof_effic_pt_final_dir_axis1_axis2+= prof_effic_pt[dir_ite][ass_ite][axis1_ite][axis2_ite],

            prof_effic_pt_final_dir_axis1+= prof_effic_pt_final_dir_axis1_axis2,

        prof_effic_pt_final_dir+= prof_effic_pt_final_dir_axis1,

    prof_effic_pt_final+= prof_effic_pt_final_dir,

##npu
prof_effic_npu_final = ()

for dir_ite in range(dir_number):
    
    prof_effic_npu_final_dir = () 
       
    for axis1_ite in range(len(prof_effic_npu[0][0])):

        prof_effic_npu_final_dir_axis1 = ()
       
        for axis2_ite in range(len(prof_effic_npu[0][0][axis1_ite])):

            prof_effic_npu_final_dir_axis1_axis2 = ()
       
            for ass_ite in range(ass_number):

                prof_effic_npu_final_dir_axis1_axis2+= prof_effic_npu[dir_ite][ass_ite][axis1_ite][axis2_ite],

            prof_effic_npu_final_dir_axis1+= prof_effic_npu_final_dir_axis1_axis2,

        prof_effic_npu_final_dir+= prof_effic_npu_final_dir_axis1,

    prof_effic_npu_final+= prof_effic_npu_final_dir,
    
print " done"
print ""
    

#________________________
#
# PLOT THEM ALL	
#________________________

# the order is directory - axis1 - axis2 - association

print "Plot them all ..."

##eta	    
  
x_bin_number = prof_effic_eta_final[0][0][0][0].GetNbinsX()
x_min = prof_effic_eta_final[0][0][0][0].GetXaxis().GetBinLowEdge(1)
x_max = prof_effic_eta_final[0][0][0][0].GetXaxis().GetBinUpEdge(x_bin_number)
  
y_bin_number = 1
y_min = 0.
y_max = 1.05 

histo_help_eta = ROOT.TH2F("histo_help_eta", "", x_bin_number, x_min, x_max, y_bin_number, y_min, y_max)
       
for axis1_ite in range(len(prof_effic_eta_final[0])):
       
    for axis2_ite in range(len(prof_effic_eta_final[0][0])):
				
	histo_title = "Npu " + npu_namesP[axis1_ite] + " - Pt " + pt_namesP[axis2_ite] + "; eta; efficiency"
	histo_help_eta.SetTitle(histo_title)
	
        savename = "eta/JetValidator_efficiency_VsEta_Npu" + npu_namesP[axis1_ite] + "_Pt" +pt_namesP[axis2_ite]
        DrawProfs2D(prof_effic_eta_final[0][axis1_ite][axis2_ite],histo_help_eta,savename)

##pt    
  
x_bin_number = prof_effic_pt_final[0][0][0][0].GetNbinsX()
x_min = prof_effic_pt_final[0][0][0][0].GetXaxis().GetBinLowEdge(1)
x_max = prof_effic_pt_final[0][0][0][0].GetXaxis().GetBinUpEdge(x_bin_number)
  
y_bin_number = 1
y_min = 0.
y_max = 1.05 

histo_help_pt = ROOT.TH2F("histo_help_pt", "", x_bin_number, x_min, x_max, y_bin_number, y_min, y_max)
       
for axis1_ite in range(len(prof_effic_pt_final[0])):

    histo_help_pt_axis1 = ()
       
    for axis2_ite in range(len(prof_effic_pt_final[0][0])):
				
	histo_title = "Npu " + npu_namesP[axis1_ite] + " - Eta " + eta_namesP[axis2_ite] + "; p_{t}^{Gen} / GeV; efficiency"
	histo_help_pt.SetTitle(histo_title)
	
        savename = "pt/JetValidator_efficiency_VsPt_Npu" + npu_namesP[axis1_ite] + "_Eta" +eta_namesP[axis2_ite]
        DrawProfs2D(prof_effic_pt_final[0][axis1_ite][axis2_ite],histo_help_pt,savename)

##npu	    
  
x_bin_number = prof_effic_npu_final[0][0][0][0].GetNbinsX()
x_min = prof_effic_npu_final[0][0][0][0].GetXaxis().GetBinLowEdge(1)
x_max = prof_effic_npu_final[0][0][0][0].GetXaxis().GetBinUpEdge(x_bin_number)
  
y_bin_number = 1
y_min = 0.
y_max = 1.05 

histo_help_npu = ROOT.TH2F("histo_help_npu", "", x_bin_number, x_min, x_max, y_bin_number, y_min, y_max)
       
for axis1_ite in range(len(prof_effic_npu_final[0])):

    histo_help_npu_axis1 = ()
       
    for axis2_ite in range(len(prof_effic_npu_final[0][0])):
				
	histo_title = "Eta " + eta_namesP[axis1_ite] + " - Pt " + pt_namesP[axis2_ite] + "; number of pileup interactions; efficiency"
	histo_help_npu.SetTitle(histo_title)
	
        savename = "npu/JetValidator_efficiency_VsNpu_Eta" + eta_namesP[axis1_ite] + "_Pt" + pt_namesP[axis2_ite]
        DrawProfs2D(prof_effic_npu_final[0][axis1_ite][axis2_ite],histo_help_npu,savename)