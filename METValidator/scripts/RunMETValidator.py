import os, fcntl, fcntl, select, sys, getopt, time

sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfKeyWords, GetListOfFilesWithEnding, GetListOfSubdirectories

print "\n This script should create the histogram based on the METValidator \n" 


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'n:r:p:o:'
keywords = ['number','run','path','outputPath']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

NumberOfEvents=-1
RunMETValidator=False
path=''
outputPath=''

for o,p in opts: 
	if o in ['-r','--run']:
		if p == "True":
			RunMETValidator = True
	if o in ['-n','--number']:
		NumberOfEvents = int(p)
	if o in ['-p','--path']:
		path = p
	if o in ['-o','--outputPath']:
		outputPath = p


#_____________________
#
# Check Input Parameters
#_____________________

PathGiven = not (path=='')

if (RunMETValidator) and (PathGiven):
	print " Please give either a path or run the METValidator. \n"
	sys.exit()

if RunMETValidator:
	print " Will run the METValidator on", NumberOfEvents, "events first. This may take some time!"


##____________________
##
## CMSRUN
#_____________________

if RunMETValidator:		
    
	OutputPathGiven = not (outputPath=='')
	if not OutputPathGiven: 
		date = time.localtime()[0:3]
		folder = "%02i"% date[2] + "%02i"% date[1] + str(date[0]) + "/"
		outputPath = '/user/geisler/METValidator/' + folder
	
	path = outputPath
	
	if NumberOfEvents == -1:
		cmsRun = 'cmsRun ../test/metvalidator_cfg.py &> ../logs/JV_log.txt \n'
	else:
		cmsRun = 'cmsRun ../test/metvalidator_cfg.py ' + str(NumberOfEvents) + ' &>../logs/MV_log.txt \n'
	print "", cmsRun
	os.system(cmsRun)
    
	print " Will move the output file to", outputPath
	makeDir = 'mkdir -p ' + outputPath
	moveFile = 'mv metvalidator.root ' + outputPath + 'metvalidator.root'
	print "", moveFile
	os.system(makeDir)
	os.system(moveFile)


##____________________
##
## PYTHON
##____________________
        
PathGiven = not (path=='')

if PathGiven:
	
	print "\n Will create the histograms using the files from", path


	file_names = GetListOfFilesWithEnding(path,".root")
	file_number = len(file_names)

	if file_number<1:
		print "Too less files"
		sys.exit()
	
	dir_names = GetListOfSubdirectories(file_names[0],"","metvalidator")
	dir_number = len(dir_names)

	if dir_number<1:
		print "Too less directories"
		sys.exit()

	ass_names = ()

	for i in range(dir_number):
	
		ass_names+= GetListOfSubdirectories(file_names[0],dir_names[i],"CorrPFMET"),
    
	ass_number = len(ass_names[0])

	if ass_number<1:
		print "Too less associations"
		sys.exit()
	
	keywords = GetListOfKeyWords(file_names[0],ass_names[0][0],"h_")
              
	for keyword in keywords:
	    
		label = keyword.split("_")[1]
	
		#makePicDir = 'mkdir -p ../pictures/'  
        #os.system(makePicDir)
	
        #makeRootDir = 'mkdir -p ../root/'    
        #os.system(makeRootDir)
	
        #makePdfDir = 'mkdir -p ../pdfs/'
        #os.system(makePdfDir)
	    
		pyRun = 'python -u ../python/METValidator.py -p ' + path + ' -k ' + keyword + ' &>../logs/mv_' + label + '_log.txt \n'
		print "", pyRun
    
        #os.system(pyRun)