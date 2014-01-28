#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os
sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfFilesWithEnding, GetListOfSubdirectories, my_colorsP, colorsMap, my_markers, markersMap, my_lines, linesMap, DrawProf, DrawTH1Fs1D


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:v:k:d:'
keywords = ['path','versions','keywords','directories']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)


path_tag=""
versions_tag= ()
keywords_tag= ()
diretory_tag= ()

for o,p in opts: 
    if o in ['-p','--path']:
        path_tag = p 
    if o in ['-v','--versions']:
        versions_tag = p.split(",")
    if o in ['-k','--keywords']:
        keywords_tag = p.split(",")
    if o in ['-d','--directories']:
        diretory_tag = p.split(",")

print "" 

num_ver = len(versions_tag)
if num_ver<1:
    print "Too less versions"
    print ""
    sys.exit()
else:
    print "You have chosen", num_ver, "versions:"
    for ver_ite in range(num_ver):
        print " ", versions_tag[ver_ite]
    print "" 

num_keys= len(keywords_tag)
if num_keys<1:
    print "Too less keywords"
    print ""
    sys.exit()
else:
    print "You have chosen", num_keys, "keywords:"
    for key_ite in range(num_keys):
        print " ", keywords_tag[key_ite] 
    print "" 

num_dirs= len(diretory_tag)
if num_dirs<1:
    print "Too less directories"
    print ""
    sys.exit()
else:
    print "You have chosen", num_dirs, "directories:"
    for dir_ite in range(num_dirs):
        print " ", diretory_tag[dir_ite] 
    print "" 
			
output_name = 'METValidatorData_'
    

#_____________________
#
# READ INPUT FILES
#_____________________

# order is version(MC/Data) - directory - histograms

print ""
print "Read input file"
print ""

prof_input = ()
files_refs = ()

for v_ite in range(num_ver):
    
	prof_input_ver= () 
	
	version_path = path_tag + versions_tag[v_ite] + "/" 

	file_names = GetListOfFilesWithEnding(version_path,".root")
	num_file = len(file_names)

	if num_file<1:
		print "Too less files found in", version_path 
		sys.exit()

	files_refs+= ROOT.TFile.Open(file_names[0]),
	print "File " + str(v_ite+1)+ ": " + str(file_names[0]) + " ...",
	
	for dir_ite in range(num_dirs):
    
		prof_input_ver_dir = () 
				
		for key_ite in range(num_keys):

			prof_help = files_refs[v_ite].Get("metvalidator/" + diretory_tag[dir_ite] + "/" + keywords_tag[key_ite])
			prof_help.SetName(diretory_tag[dir_ite] + "_" + keywords_tag[key_ite])
					    
			prof_input_ver_dir+= prof_help,
	            
		prof_input_ver+= prof_input_ver_dir,
	
	print " done"
	
	prof_input+= prof_input_ver,
    
print ""
    

#_____________________
#
# CHANGE ORDER	
#_____________________

# old order is version(MC/Data) - directory - histograms
# new order is histograms - directory - version(MC/Data)

print "Change order ...",

prof_queue = ()

for key_ite in range(num_keys):
    
    prof_key = () 
    
    for dir_ite in range(num_dirs):

		prof_key_dir = ()
		
		for ver_ite in range(num_ver):

			prof_key_dir+= prof_input[ver_ite][dir_ite][key_ite],

		prof_key+= prof_key_dir,

    prof_queue+= prof_key,
    
print " done"
print ""
    

#________________________
#
# CREATE SEPERATED PLOTS	
#________________________

# old order is histograms - directory - version(MC/Data)
# new order is histograms - directory - version(MC/Data) - axis1

print "Create seperated plots ...",

histo_npu = ()

for key_ite in range(num_keys):
	
	histo_npu_key = () 
	
	for dir_ite in range(num_dirs):

		histo_npu_key_dir = ()
		
		for ver_ite in range(num_ver):

			histo_npu_key_dir_ver = ()
		
			##npu
			for npu_ite in range(len(npu_ranges)):
			
				if npu_ite==(len(npu_ranges)-1):
					npu_min = npu_ranges[0]
					npu_max = npu_ranges[len(npu_ranges)-1]
				else:
					npu_min = npu_ranges[npu_ite]
					npu_max = npu_ranges[npu_ite+1]	 
				
				min_bin = GetBinNumber(prof_queue[key_ite][dir_ite][ver_ite].GetYaxis(),npu_min)
				max_bin = GetBinNumber(prof_queue[key_ite][dir_ite][ver_ite].GetYaxis(),npu_max)
				
				help_name = prof_queue[key_ite][dir_ite][ver_ite].GetName() + "_Npu" + npu_namesP[npu_ite]
				
				histo_npu_key_dir_ver+= prof_queue[key_ite][dir_ite][ver_ite].ProjectionX(help_name,min_bin,max_bin),
			
			histo_npu_key_dir+= histo_npu_key_dir_ver,
			
		histo_npu_key+= histo_npu_key_dir,
		
	histo_npu+= histo_npu_key,
    
    
print " done"
print ""
    

#________________________
#
# ADDITIONAL PLOTS	
#________________________

# old order is histograms - directory - version(MC/Data)
# new order is histograms - directory - version(MC/Data)

print "Create additional plots ...",

profs = ()

for key_ite in range(num_keys):
    
    profs_key = () 
    
    for dir_ite in range(num_dirs):
		
		profs_key_dir = ()
		
		for ver_ite in range(num_ver):

			profs_key_dir+= prof_queue[key_ite][dir_ite][ver_ite].ProfileY(),
			
		profs_key+= profs_key_dir,
	
    profs+= profs_dir,
    
print " done"
print ""
	
    

#_____________________
#
# CHANGE ORDER	
#_____________________

# old order is histograms - directory - version(MC/Data) - axis1
# new order is histograms - directory - axis1 - version(MC/Data) 

print "Change order again...",

##npu
histo_npu_final = ()

for key_ite in range(num_keys):
    
    histo_npu_final_key = () 
    
    for dir_ite in range(num_dirs):
		
		histo_npu_final_key_dir = ()
		
		for npu_ite in range(len(npu_ranges)):
			
			histo_npu_final_key_dir_npu = ()
		
			for ver_ite in range(num_ver):
				
				histo_npu_final_key_dir_npu+= histo_npu[key_ite][dir_ite][ver_ite][npu_ite],

        	histo_npu_final_key_dir+= histo_npu_final_key_dir_npu,
			
		histo_npu_final_key+= histo_npu_final_key_dir,

    histo_npu_final+= histo_npu_final_key,
    
print " done"
print ""

#________________________
#
# PLOT THEM ALL	
#________________________

# the order is histograms - directory - axis1 - version(MC/Data) 

print "Plot them all ... "
print ""

##npu
       
for npu_ite in range(len(histo_npu_final[0][0])):
			
	savenameData = 'METValidatorData_Npu' + npu_namesP[npu_ite]
	
	histo_npu_final[0][0][npu_ite][0].SetTitle(savenameData.split("_")[0] + " - " + savenameData.split("_")[1])
	histo_npu_final[0][0][npu_ite][0].GetXaxis().SetTitle("MET / GeV")
	histo_npu_final[0][0][npu_ite][0].GetYaxis().SetTitle("Entries")
	
	DrawTH1Fs1D(histo_npu_final[0][0][npu_ite][0],savenameData)
	

##additional	
	
x_bin_number = histo_merged[0][0].GetNbinsX()
x_min = histo_merged[0][0].GetXaxis().GetBinLowEdge(1)
x_max = histo_merged[0][0].GetXaxis().GetBinUpEdge(x_bin_number)
	    
y_bin_number = histo_merged[0][0].GetNbinsY()
y_min = histo_merged[0][0].GetYaxis().GetBinLowEdge(1)
y_max = histo_merged[0][0].GetYaxis().GetBinUpEdge(y_bin_number)

histo_help = ROOT.TH2F("histo_help"," ; number of pileup interactions; <MET> / GeV", y_bin_number,y_min,y_max,x_bin_number,x_min,x_max)

DrawProf(profs[0][0],histo_help,"METValidatorData_Prof")    
