#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetListOfSubdirectories, npu_ranges, npu_namesP, GetCombinedTH2FsWeightRek, GetBinNumber, DrawTH1Fs1D, colorsMap, my_colorsP, DrawProf

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
	
output_name = 'METValidator_'
    

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
    
dir_names = GetListOfSubdirectories(file_names[0],"","metvalidator")
dir_number = len(dir_names)

ass_names = ()

for i in range(dir_number):
	
    ass_names+=GetListOfSubdirectories(file_names[0],dir_names[i],""),
    
#ass_number = len(ass_names[0]) - 1
ass_number = len(ass_names[0])
	
print ""
    
print ""

File_ref = ()
histo_queue_in = ()

for file_ite in range(file_number):
	
    File_ref+= ROOT.TFile.Open(file_names[file_ite]),    
    print "File " + str(file_ite+1)+ ": " + str(file_names[file_ite]) + " ...",
    
    histo_file = () 
       
    for dir_ite in range(dir_number):

        histo_file_dir = ()
       
        for ass_ite in range(ass_number):
	
	    help = int(random.random() *1000)    

            histo_help = File_ref[file_ite].Get(ass_names[dir_ite][ass_ite]+histo_name)
            histo_help.SetName(histo_name+"_"+str(help))
            histo_file_dir+= histo_help,

        histo_file+= histo_file_dir,

    histo_queue_in+= histo_file,
    
    print " done"
    
print ""
    

#_____________________
#
# CHANGE ORDER	
#_____________________

# old order is file - directory - association
# new order is directory - association - file

print "Change order ...",

histo_queue = ()

for dir_ite in range(dir_number):
    
    histo_dir = () 
       
    for ass_ite in range(ass_number):

        histo_dir_ass = ()
       
        for file_ite in range(file_number):

            histo_dir_ass+= histo_queue_in[file_ite][dir_ite][ass_ite],

        histo_dir+= histo_dir_ass,

    histo_queue+= histo_dir,
    
print " done"
print ""
    

#_____________________
#
# Merge files	
#_____________________

# old order is directory - association - file
# new order is directory - association

print "Merge files ...",

histo_merged = GetCombinedTH2FsWeightRek(histo_queue)
    
print " done"
print ""
    

#________________________
#
# CREATE SEPERATED PLOTS	
#________________________

# old order is directory - association
# new order is directory - association - axis1

print "Create seperated plots ...",

histo_npu = ()

for dir_ite in range(dir_number):
	histo_npu_dir = ()  
	
	for ass_ite in range(len(histo_merged[dir_ite])):
		histo_npu_dir_ass = ()
		
		##npu
		for npu_ite in range(len(npu_ranges)):
			
			if npu_ite==(len(npu_ranges)-1):
				npu_min = npu_ranges[0]
				npu_max = npu_ranges[len(npu_ranges)-1]
			else:
				npu_min = npu_ranges[npu_ite]
				npu_max = npu_ranges[npu_ite+1]	 
				
			min_bin = GetBinNumber(histo_merged[dir_ite][ass_ite].GetYaxis(),npu_min)
			max_bin = GetBinNumber(histo_merged[dir_ite][ass_ite].GetYaxis(),npu_max)
			
			help_name = histo_merged[dir_ite][ass_ite].GetName() + "_Npu" + npu_namesP[npu_ite]
			
			histo_npu_dir_ass+= histo_merged[dir_ite][ass_ite].ProjectionX(help_name,min_bin,max_bin),
			
		histo_npu_dir+= histo_npu_dir_ass,
		
	histo_npu+= histo_npu_dir,
    
    
print " done"
print ""
    

#________________________
#
# ADDITIONAL PLOTS	
#________________________

# old order is directory - association
# new order is directory - association

print "Create additional plots ...",

profs = ()

for dir_ite in range(dir_number):
    profs_dir = ()  
       
    for ass_ite in range(len(histo_merged[dir_ite])):
        
		profs_dir+= histo_merged[dir_ite][ass_ite].ProfileY(),
	
    profs+= profs_dir,
    
print " done"
print ""
	
    

#_____________________
#
# CHANGE ORDER	
#_____________________

# old order is directory - association - axis1
# new order is directory - axis1 - association 

print "Change order again...",

##npu
histo_npu_final = ()

for dir_ite in range(dir_number):
    
    histo_npu_final_dir = () 
       
    for axis1_ite in range(len(histo_npu[0][0])):

        histo_npu_final_dir_axis1 = ()
       
        for ass_ite in range(len(histo_npu[0])):

            histo_npu_final_dir_axis1+= histo_npu[dir_ite][ass_ite][axis1_ite],

        histo_npu_final_dir+= histo_npu_final_dir_axis1,

    histo_npu_final+= histo_npu_final_dir,
    
print " done"
print ""

out = ROOT.TFile("METValidator_" + histo_name + ".root","RECREATE")
		
for ite_2 in range(len(profs[0])):
	histo_npu_final[0][0][ite_2].Write(ass_names[0][ite_2].split("/")[1])
	profs[0][ite_2].Write("prof_" + ass_names[0][ite_2].split("/")[1])

out.Close()

sys.exit()

#________________________
#
# PLOT THEM ALL	
#________________________

# the order is directory - axis1 - association

print "Plot them all ... "
print ""

for ass_ite in range(ass_number): 
    print ass_names[0][ass_ite].split("/")[1], "is drawn in", colorsMap[my_colorsP[ass_ite]]
    
print "" 

##npu
       
for npu_ite in range(len(histo_npu_final[0])):
	
    savename = 'METValidator_Npu' + npu_namesP[npu_ite] 
    
    histo_npu_final[0][npu_ite][0].SetTitle(savename.split("_")[0] + " - " + savename.split("_")[1])
    histo_npu_final[0][npu_ite][0].GetXaxis().SetTitle("MET / GeV")
    histo_npu_final[0][npu_ite][0].GetYaxis().SetTitle("Entries")
    
    DrawTH1Fs1D(histo_npu_final[0][npu_ite],savename)
	    
	    
##additional

x_bin_number = histo_merged[0][0].GetNbinsX()
x_min = histo_merged[0][0].GetXaxis().GetBinLowEdge(1)
x_max = histo_merged[0][0].GetXaxis().GetBinUpEdge(x_bin_number)
	    
y_bin_number = histo_merged[0][0].GetNbinsY()
y_min = histo_merged[0][0].GetYaxis().GetBinLowEdge(1)
y_max = histo_merged[0][0].GetYaxis().GetBinUpEdge(y_bin_number)

histo_help = ROOT.TH2F("histo_help"," ; number of pileup interactions; <MET> / GeV", y_bin_number,y_min,y_max,x_bin_number,x_min,x_max)

DrawProf(profs[0],histo_help,"METValidator_Prof")    