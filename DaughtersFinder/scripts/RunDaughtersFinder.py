import os, fcntl, fcntl, select, sys, getopt, time

sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfKeyWords, GetListOfFilesWithEnding, GetListOfSubdirectories

motherDic = {"Z":"zdf","W":"wdf","B":"bdf","T":"tdf"}

print "\n This script should create the histogram based on the DaughtersFinder \n" 


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:m:s:'
keywords = ['path','mother','subdirectories']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

path=''
mother=''
subdirectories=''

for o,p in opts: 
    if o in ['-p','--path']:
	path = p
    if o in ['-m','--mother']:
	mother = p
    if o in ['-s','--subdirectories']:
	subdirectories = p.split(",")
	
directory = motherDic[mother]


#_____________________
#
# Check Input Parameters
#_____________________

PathGiven = not (path=='')

if not PathGiven:
    print " No Input path given! =( \n"
    sys.exit()


##____________________
##
## PYTHON
##____________________
        

if PathGiven:
	
    print "\n Will create the histograms using the files from", path

    file_names = GetListOfFilesWithEnding(path,".root")
    file_number = len(file_names)

    if file_number<1:
        print "Too less files"
        sys.exit()
	
    dir_names = GetListOfSubdirectories(file_names[0],"",directory)
    dir_number = len(dir_names)

    if dir_number<1:
        print "Too less directories"
        sys.exit()
    
    for subdir in subdirectories:
		
		ass_names = ()
		
		for i in range(dir_number):
			ass_names+= GetListOfSubdirectories(file_names[0],dir_names[i],subdir),
		
		ass_number = len(ass_names[0])
		collections = ''
		for ass_ite in range(ass_number):
			collections+= ass_names[0][ass_ite] + ","
			
		collections+= dir_names[0] + "cutsRecoTracks/"
		
		if ass_number<1:
			print "Too less associations"
			continue
		
		keywords = GetListOfKeyWords(file_names[0],ass_names[0][0],"p_")
		
		for keyword in keywords:
	     
			label = keyword.split("_")[1]
		
			makePicEtaDir = 'mkdir -p ../pictures/' + mother + '/eta/' + subdir + '/' + label + '/'    
			os.system(makePicEtaDir)
			makePicPtDir = 'mkdir -p ../pictures/' + mother + '/pt/' + subdir + '/' + label + '/'   
			os.system(makePicPtDir)
			makePicNpuDir = 'mkdir -p ../pictures/' + mother + '/npu/' + subdir + '/' + label + '/'    
			os.system(makePicNpuDir)
	
			makeRootEtaDir = 'mkdir -p ../root/' + mother + '/eta/' + subdir + '/' + label + '/'    
			os.system(makeRootEtaDir)
			makeRootPtDir = 'mkdir -p ../root/' + mother + '/pt/' + subdir + '/' + label + '/'    
			os.system(makeRootPtDir)
			makeRootNpuDir = 'mkdir -p ../root/' + mother + '/npu/' + subdir + '/' + label + '/'    
			os.system(makeRootNpuDir)
	
			makePdfEtaDir = 'mkdir -p ../pdfs/' + mother + '/eta/' + subdir + '/' + label + '/'    
			os.system(makePdfEtaDir)
			makePdfPtDir = 'mkdir -p ../pdfs/' + mother + '/pt/' + subdir + '/' + label + '/'    
			os.system(makePdfPtDir)
			makePdfNpuDir = 'mkdir -p ../pdfs/' + mother + '/npu/' + subdir + '/' + label + '/'    
			os.system(makePdfNpuDir)
		
			pyRun = 'python -u ../python/daughtersfinder.py -p ' + path + ' -k ' + keyword + ' -m ' + mother + ' -d ' + directory + ' -c ' + collections + ' -s ' + subdir + ' &>../logs/df_' + label + '_log.txt \n'
			print "", pyRun
	
			os.system(pyRun)
