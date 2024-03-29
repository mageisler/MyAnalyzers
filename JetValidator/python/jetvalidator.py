#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetListOfSubdirectories, GetMergedProf3DRek, npu_ranges, pt_rangesJET, eta_rangesJET,  Create1DProf, npu_namesP, pt_namesPJET, eta_namesPJET, DrawProfs2DWithHisto, colorsMap, my_colorsP, GetProjection, UpperEdges, LowerEdges


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:k:'
keywords = ['path','keyword']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)


path=""
histo_name=""

for o,p in opts: 
    if o in ['-p','--path']:
        path = p 
    if o in ['-k','--keyword']:
        histo_name = p
	
label = histo_name.split("_")[1]
output_name = 'JetValidator_' + label
    

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
	
    ass_names+=GetListOfSubdirectories(file_names[0],dir_names[i],"patJets"),
    
ass_number = len(ass_names[0])

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

prof_merged = GetMergedProf3DRek(prof_queue)
    
print " done"
print ""
    

#________________________
#
# CREATE SEPERATED PLOTS	
#________________________

# old order is directory - association
# new order is directory - association - axis1 - axis2

print "Create seperated plots ...",

prof_eta = ()
prof_pt  = ()
prof_npu = ()

for dir_ite in range(dir_number):
    
    prof_eta_dir = () 
    prof_pt_dir = ()
    prof_npu_dir = ()  
       
    for ass_ite in range(ass_number):
    
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
  
            for pt_ite in range(len(pt_rangesJET)):
                if pt_ite==(len(pt_rangesJET)-1):
                    pt_min = pt_rangesJET[0]
                    pt_max = pt_rangesJET[len(pt_rangesJET)-1]
                else:
                    pt_min = pt_rangesJET[pt_ite]
                    pt_max = pt_rangesJET[pt_ite+1]
	    	    
	        prof_eta_dir_ass_npu+= Create1DProf(prof_merged[dir_ite][ass_ite], pt_min,pt_max, npu_min,npu_max,"x"),
	    	    
	    prof_eta_dir_ass+= prof_eta_dir_ass_npu,
    
            ##pt
            prof_pt_dir_ass_npu = ()	 
  
            for eta_ite in range(len(eta_rangesJET)):
                if eta_ite==(len(eta_rangesJET)-1):
                    eta_min = eta_rangesJET[0]
                    eta_max = eta_rangesJET[len(eta_rangesJET)-1]
                else:
                    eta_min = eta_rangesJET[eta_ite]
                    eta_max = eta_rangesJET[eta_ite+1]
	    	    
	        prof_pt_dir_ass_npu+= Create1DProf(prof_merged[dir_ite][ass_ite], eta_min,eta_max, npu_min,npu_max,"y"),
	    	    
	    prof_pt_dir_ass+= prof_pt_dir_ass_npu,

        ##npu
        for eta_ite in range(len(eta_rangesJET)):
	    if eta_ite==(len(eta_rangesJET)-1):
		eta_min = eta_rangesJET[0]
		eta_max = eta_rangesJET[len(eta_rangesJET)-1]
            else:
		eta_min = eta_rangesJET[eta_ite]
		eta_max = eta_rangesJET[eta_ite+1]
			                
            prof_npu_dir_ass_eta = ()	 
  
            for pt_ite in range(len(pt_rangesJET)):
                if pt_ite==(len(pt_rangesJET)-1):
                    pt_min = pt_rangesJET[0]
                    pt_max = pt_rangesJET[len(pt_rangesJET)-1]
                else:
                    pt_min = pt_rangesJET[pt_ite]
                    pt_max = pt_rangesJET[pt_ite+1]
	    	    
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
       
            for ass_ite in range(ass_number):

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
       
            for ass_ite in range(ass_number):

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
       
            for ass_ite in range(ass_number):

                prof_npu_final_dir_axis1_axis2+= prof_npu[dir_ite][ass_ite][axis1_ite][axis2_ite],

            prof_npu_final_dir_axis1+= prof_npu_final_dir_axis1_axis2,

        prof_npu_final_dir+= prof_npu_final_dir_axis1,

    prof_npu_final+= prof_npu_final_dir,
    
print " done"
print ""
												
sn = output_name.split("_")[0] + '_' + label

out = ROOT.TFile(sn+".root","RECREATE")


for ite_4 in range(len(prof_eta_final[0][0][0])):
	prof_eta_final[0][0][0][ite_4].Write("eta_" + ass_names[0][ite_4].split("/")[1])
	prof_pt_final[0][0][0][ite_4].Write("pt_" + ass_names[0][ite_4].split("/")[1])
	prof_npu_final[0][0][0][ite_4].Write("npu_" + ass_names[0][ite_4].split("/")[1])

out.Close()

sys.exit()
    

#________________________
#
# PLOT THEM ALL	
#________________________

# the order is directory - axis1 - axis2 - association

print "Plot them all ... "
print ""

for ass_ite in range(ass_number): 
    print ass_names[0][ass_ite].split("/")[1], "is drawn in", colorsMap[my_colorsP[ass_ite]]
    
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
			
			if axis2_ite==(len(pt_rangesJET)-1):
				pt_min = pt_rangesJET[0]
				pt_max = pt_rangesJET[len(pt_rangesJET)-1]
			else:
				pt_min = pt_rangesJET[axis2_ite]
				pt_max = pt_rangesJET[axis2_ite+1]
				
			histo_proj = GetProjection(prof_merged[0], histo_help_eta.GetYaxis(), pt_min,pt_max, npu_min, npu_max,  "x")
			
			histo_title = 'Npu ' + npu_namesP[axis1_ite] + ' - Pt ' + pt_namesPJET[axis2_ite] +'; eta;' + label
			histo_help_eta.SetTitle(histo_title)
			
			savename = 'eta/' + label + '/' + output_name + 'VsEta_Npu' + npu_namesP[axis1_ite] + '_Pt' +pt_namesPJET[axis2_ite]
			DrawProfs2DWithHisto(prof_eta_final[0][axis1_ite][axis2_ite],histo_help_eta,histo_proj,savename)
	

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
	
        if axis2_ite==(len(eta_rangesJET)-1):
            eta_min = eta_rangesJET[0]
            eta_max = eta_rangesJET[len(eta_rangesJET)-1]
        else:
            eta_min = eta_rangesJET[axis2_ite]
            eta_max = eta_rangesJET[axis2_ite+1]
	    
	histo_proj = GetProjection(prof_merged[0], histo_help_pt.GetYaxis(), eta_min,eta_max, npu_min, npu_max,  "y")
							
	histo_title = 'Npu ' + npu_namesP[axis1_ite] + ' - Eta ' + eta_namesPJET[axis2_ite] + '; p_{t} / GeV;' + label
	histo_help_pt.SetTitle(histo_title)
	
        savename = 'pt/' + label + '/' + output_name + 'VsPt_Npu' + npu_namesP[axis1_ite] + '_Eta' +eta_namesPJET[axis2_ite]
        DrawProfs2DWithHisto(prof_pt_final[0][axis1_ite][axis2_ite],histo_help_pt,histo_proj,savename)

##npu	    
  
x_bin_number = prof_npu_final[0][0][0][0].GetNbinsX()
x_min = prof_npu_final[0][0][0][0].GetXaxis().GetBinLowEdge(1)
x_max = prof_npu_final[0][0][0][0].GetXaxis().GetBinUpEdge(x_bin_number)  

histo_help_npu = ROOT.TH2F("histo_help_npu", "", x_bin_number, x_min, x_max, y_bin_number, y_min, y_max)
       
for axis1_ite in range(len(prof_npu_final[0])):
	
    if axis1_ite==(len(eta_rangesJET)-1):
        eta_min = eta_rangesJET[0]
        eta_max = eta_rangesJET[len(eta_rangesJET)-1]
    else:
        eta_min = eta_rangesJET[axis1_ite]
        eta_max = eta_rangesJET[axis1_ite+1]
       
    for axis2_ite in range(len(prof_npu_final[0][0])):
	
        if axis2_ite==(len(pt_rangesJET)-1):
            pt_min = pt_rangesJET[0]
            pt_max = pt_rangesJET[len(pt_rangesJET)-1]
        else:
            pt_min = pt_rangesJET[axis2_ite]
            pt_max = pt_rangesJET[axis2_ite+1]
	    
	histo_proj = GetProjection(prof_merged[0], histo_help_npu.GetYaxis(), eta_min,eta_max, pt_min, pt_max,  "z")
	    			
	histo_title = 'Eta ' + eta_namesPJET[axis1_ite] + ' - Pt ' + pt_namesPJET[axis2_ite] + '; #PU;' + label
	histo_help_npu.SetTitle(histo_title)
	
        savename = 'npu/' + label + '/' + output_name + 'VsNpu_Eta' + eta_namesPJET[axis1_ite] + '_Pt' + pt_namesPJET[axis2_ite]
        DrawProfs2DWithHisto(prof_npu_final[0][axis1_ite][axis2_ite],histo_help_npu,histo_proj,savename)