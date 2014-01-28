#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetListOfSubdirectories, GetMergedProf3DRek, npu_ranges, pt_rangesTRACK, eta_rangesTRACK, Create1DProf, npu_namesP, pt_namesPTRACK, eta_namesPTRACK, DrawProfs2DWithHisto, colorsMap, my_colorsP, my_colorsF, GetProjection, UpperEdges, LowerEdges, CreateFractionProfsRek


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:k:m:d:s:c:'
keywords = ['path','keyword','mother','directory','subdirectory','collections']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)


path=""
mother=""
directory=""
histo_name=""
subdirectory=""
collections= ()

for o,p in opts: 
    if o in ['-p','--path']:
        path = p 
    if o in ['-k','--keyword']:
        histo_name = p
    if o in ['-m','--mother']:
        mother = p
    if o in ['-d','--directory']:
        directory = p
    if o in ['-s','--subdirectory']:
        subdirectory = p
    if o in ['-c','--collections']:
        collections = p.split(",")
	
label = histo_name.split("_")[1]	
output_name = 'DaughtersFinder_' + label

useDefaultOnly = False
if "Q0" in collections[0]:
    useDefaultOnly = True
    

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
dir_names = GetListOfSubdirectories(file_names[0],"",directory)
dir_number = len(dir_names)

ass_names = ()
ass_number = len(collections)

for dir_ite in range(dir_number):
	
    ass_names_help = ()
	
    for coll_ite in range(ass_number):
		
        ass_names_help+= collections[coll_ite],
	
    ass_names+= ass_names_help,

if ass_number<1:
    print "Too less associations"
    sys.exit()
	
print ""

File_ref = ()
prof_queue_in = ()

for file_ite in range(file_number):
	
    File_ref+= ROOT.TFile.Open(file_names[file_ite]),    
    print "File " + str(file_ite+1)+ ": " + str(file_names[file_ite]) + " ...",
    
    prof_file = () 
       
    for dir_ite in range(dir_number):

        prof_file_dir = ()
       
        for ass_ite in range(ass_number):
	
	    help = int(random.random() *1000)    

            prof_help = File_ref[file_ite].Get(ass_names[dir_ite][ass_ite]+histo_name)
            prof_help.SetName(histo_name+"_"+str(help))
            prof_file_dir+= prof_help,

        prof_file+= prof_file_dir,

    prof_queue_in+= prof_file,
    
    print " done"
    
print ""
    

#_____________________
#
# CHANGE ORDER	
#_____________________

# old order is file - directory - association
# new order is directory - association - file

print "Change order ...",

prof_queue = ()

for dir_ite in range(dir_number):
    
    prof_dir = () 
       
    for ass_ite in range(ass_number):

        prof_dir_ass = ()
       
        for file_ite in range(file_number):

            prof_dir_ass+= prof_queue_in[file_ite][dir_ite][ass_ite],

        prof_dir+= prof_dir_ass,

    prof_queue+= prof_dir,
    
print " done"
print ""
    

#_____________________
#
# Merge files	
#_____________________

# old order is directory - association - file
# new order is directory - association

print "Merge files ...",

prof_merged_pre = GetMergedProf3DRek(prof_queue)
    
print " done"
print ""


#________________________
#
# DO ADDITIONAL STUFF	
#________________________

# the order is directory - association

print "Create additional fraction profs ... ",

if label == "trackEfficiency":
    prof_merged = CreateFractionProfsRek(prof_merged_pre)
else:
    prof_merged = prof_merged_pre
    
print " done"
print ""
    

#________________________
#
# CREATE SEPARATED PLOTS	
#________________________

# old order is directory - association
# new order is directory - association - axis1 - axis2

print "Create separated plots ...",

prof_eta = ()
prof_pt  = ()
prof_npu = ()

for dir_ite in range(dir_number):
    
    prof_eta_dir = () 
    prof_pt_dir = ()
    prof_npu_dir = ()  
       
    for ass_ite in range(len(prof_merged[dir_ite])):
    
        prof_eta_dir_ass = ()
        prof_pt_dir_ass = ()
        prof_npu_dir_ass = ()

        ##eta && pt
        for npu_ite in range(len(npu_ranges)):
            if npu_ite==(len(npu_ranges)-1):
                npu_min = npu_ranges[0]
                npu_max = npu_ranges[len(npu_ranges)-1]
            else:
                npu_min = npu_ranges[npu_ite]
                npu_max = npu_ranges[npu_ite+1]
			
            ##eta    
            prof_eta_dir_ass_npu = ()	 
  
            for pt_ite in range(len(pt_rangesTRACK)):
                if pt_ite==(len(pt_rangesTRACK)-1):
                    pt_min = pt_rangesTRACK[0]
                    pt_max = pt_rangesTRACK[len(pt_rangesTRACK)-1]
                else:
                    pt_min = pt_rangesTRACK[pt_ite]
                    pt_max = pt_rangesTRACK[pt_ite+1]
	    	    
	        prof_eta_dir_ass_npu+= Create1DProf(prof_merged[dir_ite][ass_ite], pt_min,pt_max, npu_min,npu_max,"x"),
	    	    
	    prof_eta_dir_ass+= prof_eta_dir_ass_npu,
    
            ##pt
            prof_pt_dir_ass_npu = ()	 
  
            for eta_ite in range(len(eta_rangesTRACK)):
                if eta_ite==(len(eta_rangesTRACK)-1):
                    eta_min = eta_rangesTRACK[0]
                    eta_max = eta_rangesTRACK[len(eta_rangesTRACK)-1]
                else:
                    eta_min = eta_rangesTRACK[eta_ite]
                    eta_max = eta_rangesTRACK[eta_ite+1]
	    	    
	        prof_pt_dir_ass_npu+= Create1DProf(prof_merged[dir_ite][ass_ite], eta_min,eta_max, npu_min,npu_max,"y"),
	    	    
	    prof_pt_dir_ass+= prof_pt_dir_ass_npu,

        ##npu
        for eta_ite in range(len(eta_rangesTRACK)):
	    if eta_ite==(len(eta_rangesTRACK)-1):
		eta_min = eta_rangesTRACK[0]
		eta_max = eta_rangesTRACK[len(eta_rangesTRACK)-1]
            else:
		eta_min = eta_rangesTRACK[eta_ite]
		eta_max = eta_rangesTRACK[eta_ite+1]
			                
            prof_npu_dir_ass_eta = ()	 
  
            for pt_ite in range(len(pt_rangesTRACK)):
                if pt_ite==(len(pt_rangesTRACK)-1):
                    pt_min = pt_rangesTRACK[0]
                    pt_max = pt_rangesTRACK[len(pt_rangesTRACK)-1]
                else:
                    pt_min = pt_rangesTRACK[pt_ite]
                    pt_max = pt_rangesTRACK[pt_ite+1]
	    	    
	        prof_npu_dir_ass_eta+= Create1DProf(prof_merged[dir_ite][ass_ite], eta_min,eta_max, pt_min,pt_max,"z"),
	    	    
	    prof_npu_dir_ass+= prof_npu_dir_ass_eta,
	    	    
	prof_eta_dir+= prof_eta_dir_ass,
	prof_pt_dir+= prof_pt_dir_ass,
	prof_npu_dir+= prof_npu_dir_ass,
	    	    
    prof_eta+= prof_eta_dir,
    prof_pt+= prof_pt_dir,
    prof_npu+= prof_npu_dir,
    
    
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
prof_eta_final = ()

for dir_ite in range(dir_number):
    
    prof_eta_final_dir = () 
       
    for axis1_ite in range(len(prof_eta[0][0])):

        prof_eta_final_dir_axis1 = ()
       
        for axis2_ite in range(len(prof_eta[0][0][axis1_ite])):

            prof_eta_final_dir_axis1_axis2 = ()
       
            for ass_ite in range(len(prof_eta[0])):

                prof_eta_final_dir_axis1_axis2+= prof_eta[dir_ite][ass_ite][axis1_ite][axis2_ite],

            prof_eta_final_dir_axis1+= prof_eta_final_dir_axis1_axis2,

        prof_eta_final_dir+= prof_eta_final_dir_axis1,

    prof_eta_final+= prof_eta_final_dir,

##pt
prof_pt_final = ()

for dir_ite in range(dir_number):
    
    prof_pt_final_dir = () 
       
    for axis1_ite in range(len(prof_pt[0][0])):

        prof_pt_final_dir_axis1 = ()
       
        for axis2_ite in range(len(prof_pt[0][0][axis1_ite])):

            prof_pt_final_dir_axis1_axis2 = ()
       
            for ass_ite in range(len(prof_pt[0])):

                prof_pt_final_dir_axis1_axis2+= prof_pt[dir_ite][ass_ite][axis1_ite][axis2_ite],

            prof_pt_final_dir_axis1+= prof_pt_final_dir_axis1_axis2,

        prof_pt_final_dir+= prof_pt_final_dir_axis1,

    prof_pt_final+= prof_pt_final_dir,

##npu
prof_npu_final = ()

for dir_ite in range(dir_number):
    
    prof_npu_final_dir = () 
       
    for axis1_ite in range(len(prof_npu[0][0])):

        prof_npu_final_dir_axis1 = ()
       
        for axis2_ite in range(len(prof_npu[0][0][axis1_ite])):

            prof_npu_final_dir_axis1_axis2 = ()
       
            for ass_ite in range(len(prof_npu[0])):

                prof_npu_final_dir_axis1_axis2+= prof_npu[dir_ite][ass_ite][axis1_ite][axis2_ite],

            prof_npu_final_dir_axis1+= prof_npu_final_dir_axis1_axis2,

        prof_npu_final_dir+= prof_npu_final_dir_axis1,

    prof_npu_final+= prof_npu_final_dir,
    
print " done"
print ""
    

#________________________
#
# PLOT THEM ALL	
#________________________

# the order is directory - axis1 - axis2 - association

print "Plot them all ... "
print ""

for ass_ite in range(ass_number-1): 
	
    color = my_colorsP[ass_ite]
    if useDefaultOnly:
        color = my_colorsF[ass_ite]   
    print ass_names[0][ass_ite].split("/")[1], "is drawn in", colorsMap[color]
    
print "" 
 
y_bin_number = 1
y_min = LowerEdges[label]
y_max = UpperEdges[label]

##eta	    
  
x_bin_number = prof_eta_final[0][0][0][0].GetNbinsX()
x_min = prof_eta_final[0][0][0][0].GetXaxis().GetBinLowEdge(1)
x_max = prof_eta_final[0][0][0][0].GetXaxis().GetBinUpEdge(x_bin_number)

histo_help_eta = ROOT.TH2F("histo_help_eta", "", x_bin_number, x_min, x_max, y_bin_number, y_min, y_max)
       
for axis1_ite in range(len(prof_eta_final[0])):
	
    if axis1_ite==(len(npu_ranges)-1):
        npu_min = npu_ranges[0]
        npu_max = npu_ranges[len(npu_ranges)-1]
    else:
        npu_min = npu_ranges[axis1_ite]
        npu_max = npu_ranges[axis1_ite+1]
       
    for axis2_ite in range(len(prof_eta_final[0][0])):
	
        if axis2_ite==(len(pt_rangesTRACK)-1):
            pt_min = pt_rangesTRACK[0]
            pt_max = pt_rangesTRACK[len(pt_rangesTRACK)-1]
        else:
            pt_min = pt_rangesTRACK[axis2_ite]
            pt_max = pt_rangesTRACK[axis2_ite+1]
	    
	histo_proj = GetProjection(prof_merged[0], histo_help_eta.GetYaxis(), pt_min,pt_max, npu_min, npu_max,  "x")
				
	histo_title = 'Npu ' + npu_namesP[axis1_ite] + ' - Pt ' + pt_namesPTRACK[axis2_ite] +'; eta;' + label
	histo_help_eta.SetTitle(histo_title)
	
        savename = mother + '/eta/' + subdirectory + '/' + label + '/' + output_name + 'VsEta_Npu' + npu_namesP[axis1_ite] + '_Pt' +pt_namesPTRACK[axis2_ite]
        DrawProfs2DWithHisto(prof_eta_final[0][axis1_ite][axis2_ite],histo_help_eta,histo_proj,savename,useDefaultOnly)
	

##pt    
  
x_bin_number = prof_pt_final[0][0][0][0].GetNbinsX()
x_min = prof_pt_final[0][0][0][0].GetXaxis().GetBinLowEdge(1)
x_max = prof_pt_final[0][0][0][0].GetXaxis().GetBinUpEdge(x_bin_number)

histo_help_pt = ROOT.TH2F("histo_help_pt", "", x_bin_number, x_min, x_max, y_bin_number, y_min, y_max)
       
for axis1_ite in range(len(prof_pt_final[0])):
	
    if axis1_ite==(len(npu_ranges)-1):
        npu_min = npu_ranges[0]
        npu_max = npu_ranges[len(npu_ranges)-1]
    else:
        npu_min = npu_ranges[axis1_ite]
        npu_max = npu_ranges[axis1_ite+1]
       
    for axis2_ite in range(len(prof_pt_final[0][0])):
	
        if axis2_ite==(len(eta_rangesTRACK)-1):
            eta_min = eta_rangesTRACK[0]
            eta_max = eta_rangesTRACK[len(eta_rangesTRACK)-1]
        else:
            eta_min = eta_rangesTRACK[axis2_ite]
            eta_max = eta_rangesTRACK[axis2_ite+1]
	    
	histo_proj = GetProjection(prof_merged[0], histo_help_pt.GetYaxis(), eta_min,eta_max, npu_min, npu_max,  "y")
							
	histo_title = 'Npu ' + npu_namesP[axis1_ite] + ' - Eta ' + eta_namesPTRACK[axis2_ite] + '; p_{t} / GeV;' + label
	histo_help_pt.SetTitle(histo_title)
	
        savename = mother + '/pt/' + subdirectory + '/' + label + '/' + output_name + 'VsPt_Npu' + npu_namesP[axis1_ite] + '_Eta' +eta_namesPTRACK[axis2_ite]
        DrawProfs2DWithHisto(prof_pt_final[0][axis1_ite][axis2_ite],histo_help_pt,histo_proj,savename,useDefaultOnly)

##npu	    
  
x_bin_number = prof_npu_final[0][0][0][0].GetNbinsX()
x_min = prof_npu_final[0][0][0][0].GetXaxis().GetBinLowEdge(1)
x_max = prof_npu_final[0][0][0][0].GetXaxis().GetBinUpEdge(x_bin_number)  

histo_help_npu = ROOT.TH2F("histo_help_npu", "", x_bin_number, x_min, x_max, y_bin_number, y_min, y_max)
       
for axis1_ite in range(len(prof_npu_final[0])):
	
    if axis1_ite==(len(eta_rangesTRACK)-1):
        eta_min = eta_rangesTRACK[0]
        eta_max = eta_rangesTRACK[len(eta_rangesTRACK)-1]
    else:
        eta_min = eta_rangesTRACK[axis1_ite]
        eta_max = eta_rangesTRACK[axis1_ite+1]
       
    for axis2_ite in range(len(prof_npu_final[0][0])):
	
        if axis2_ite==(len(pt_rangesTRACK)-1):
            pt_min = pt_rangesTRACK[0]
            pt_max = pt_rangesTRACK[len(pt_rangesTRACK)-1]
        else:
            pt_min = pt_rangesTRACK[axis2_ite]
            pt_max = pt_rangesTRACK[axis2_ite+1]
	    
	histo_proj = GetProjection(prof_merged[0], histo_help_npu.GetYaxis(), eta_min,eta_max, pt_min, pt_max,  "z")
	    			
	histo_title = 'Eta ' + eta_namesPTRACK[axis1_ite] + ' - Pt ' + pt_namesPTRACK[axis2_ite] + '; #PU;' + label
	histo_help_npu.SetTitle(histo_title)
	
        savename = mother + '/npu/' + subdirectory + '/' + label + '/' + output_name + 'VsNpu_Eta' + eta_namesPTRACK[axis1_ite] + '_Pt' + pt_namesPTRACK[axis2_ite]
        DrawProfs2DWithHisto(prof_npu_final[0][axis1_ite][axis2_ite],histo_help_npu,histo_proj,savename,useDefaultOnly)
