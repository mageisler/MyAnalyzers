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
			
output_name = 'TrackValidatorData_'
    

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

			prof_help = files_refs[v_ite].Get("trackvalidatordata/" + diretory_tag[dir_ite] + "/" + keywords_tag[key_ite])
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
# PLOT THEM ALL	
#________________________

# the order is histograms - directory - version(MC/Data)

print "Plot them all ... "
print ""

for ver_ite in range(num_ver):
    
    print versions_tag[ver_ite], "is drawn in", colorsMap[my_colorsP[ver_ite]]
	
print ""

for key_ite in range(num_keys):
    
    for dir_ite in range(num_dirs):
			
		dir_tag = diretory_tag[dir_ite]
				
		if 'p_' in keywords_tag[key_ite]:

			l_help = ()
			  
			x_bin_number = prof_queue[key_ite][dir_ite][0].GetNbinsX()
			x_min = prof_queue[key_ite][dir_ite][0].GetXaxis().GetBinLowEdge(1)
			x_max = prof_queue[key_ite][dir_ite][0].GetXaxis().GetBinUpEdge(x_bin_number)
			
			y_bin_number = 1
			mini = prof_queue[key_ite][dir_ite][0].GetMinimum()
			maxi = prof_queue[key_ite][dir_ite][0].GetMaximum()
				
			for ver_ite in range(num_ver):    
				if ( prof_queue[key_ite][dir_ite][ver_ite].GetMaximum()>maxi ): 
					maxi = prof_queue[key_ite][dir_ite][ver_ite].GetMaximum()
				if ( prof_queue[key_ite][dir_ite][ver_ite].GetMinimum()<mini ): 
					mini = prof_queue[key_ite][dir_ite][ver_ite].GetMinimum()
							
							
			if ( mini>2*math.sqrt(mini) ):
				y_min = mini - 2*math.sqrt(mini) 
			else:
				y_min = 0
					
			y_max = maxi + 1.1*math.sqrt(maxi)	
				
				
			for version in versions_tag:    
				if version == "MC":
					l_help+= "Simulation",
				else:
					l_help+= version,  
				
			titles = ';'
			titles+= prof_queue[key_ite][dir_ite][0].GetXaxis().GetTitle() + ';'
			titles+= prof_queue[key_ite][dir_ite][0].GetYaxis().GetTitle()
			
			help = str(random.random()*1000)
			
			histo_help = ROOT.TH2F("histo_help"+help,titles,x_bin_number,x_min,x_max,y_bin_number,y_min,y_max)
			
			savename = output_name + keywords_tag[key_ite].split("_")[1] + "_Vs_"  + keywords_tag[key_ite].split("_")[2] + "_" + dir_tag
			DrawProf(prof_queue[key_ite][dir_ite],histo_help,savename,changeMarkerStyle=True,legend=l_help)
			
		else:
			
			for ver_ite in range(num_ver):
				if versions_tag[ver_ite] == "MC":
					prof_queue[key_ite][dir_ite][ver_ite].SetName("Simulation")
				else:
					prof_queue[key_ite][dir_ite][ver_ite].SetName(versions_tag[ver_ite])
			
			prof_queue[key_ite][dir_ite][0].SetMinimum(0.)
			
			prof_queue[key_ite][dir_ite][0].SetTitle("")
			
			savename = output_name + keywords_tag[key_ite].split("_")[1] 
			if len(keywords_tag[key_ite].split("_"))>2:
				savename+= "_Vs_"  + keywords_tag[key_ite].split("_")[2]
				
			savename+= "_" + dir_tag
			
			DrawTH1Fs1D(prof_queue[key_ite][dir_ite],savename,scale=True,legend=True,changeLineStyle=True)