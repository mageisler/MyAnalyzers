import os, fcntl, fcntl, select, sys, getopt, time

sys.path.append('/user/geisler/CMSSW/Helpers/')
from AnalysisTools import GetListOfKeyWords, GetListOfFilesWithEnding, GetListOfSubdirectories

print "\n This script should create the histogram based on the JetValidator \n" 


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'n:r:p:o:'
keywords = ['number','run','path','outputPath']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

NumberOfEvents=-1
RunJetValidator=False
path=''
outputPath=''

for o,p in opts: 
	if o in ['-r','--run']:
		if p == "True":
			RunJetValidator = True
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

if (RunJetValidator) and (PathGiven):
	print " Please give either a path or run the JetValidator. \n"
	sys.exit()

if RunJetValidator:
	print " Will run the JetValidator on", NumberOfEvents, "events first. This may take some time!"


##____________________
##
## CMSRUN
#_____________________

if RunJetValidator:		
    
	OutputPathGiven = not (outputPath=='')
	if not OutputPathGiven: 
		date = time.localtime()[0:3]
		folder = "%02i"% date[2] + "%02i"% date[1] + str(date[0]) + "/"
		outputPath = '/user/geisler/JetValidator/' + folder
	
	path = outputPath
	
	if NumberOfEvents == -1:
		cmsRun = 'cmsRun ../test/jetvalidator_cfg.py &> ../logs/JV_log.txt \n'
	else:
		cmsRun = 'cmsRun ../test/jetvalidator_cfg.py ' + str(NumberOfEvents) + ' &>../logs/JV_log.txt \n'
	print "", cmsRun
	os.system(cmsRun)
    
	print " Will move the output file to", outputPath
	makeDir = 'mkdir -p ' + outputPath
	moveFile = 'mv jetvalidator.root ' + outputPath + 'jetvalidator.root'
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
	
	dir_names = GetListOfSubdirectories(file_names[0],"","jetvalidator")
	dir_number = len(dir_names)

	if dir_number<1:
		print "Too less directories"
		sys.exit()

	ass_names = ()

	for i in range(dir_number):
	
		ass_names+= GetListOfSubdirectories(file_names[0],dir_names[i],"patJets"),
    
	ass_number = len(ass_names[0])

	if ass_number<1:
		print "Too less associations"
		sys.exit()
	
	keywords = GetListOfKeyWords(file_names[0],ass_names[0][0],"p_ptResponse")
              
	for keyword in keywords:
	    
		label = keyword.split("_")[1]
	
        #makePicEtaDir = 'mkdir -p ../pictures/eta/' + label + '/'    
        #os.system(makePicEtaDir)
        #makePicPtDir = 'mkdir -p ../pictures/pt/' + label + '/'    
        #os.system(makePicPtDir)
        #makePicNpuDir = 'mkdir -p ../pictures/npu/' + label + '/'    
        #os.system(makePicNpuDir)
	
        #makeRootEtaDir = 'mkdir -p ../root/eta/' + label + '/'    
        #os.system(makeRootEtaDir)
        #makeRootPtDir = 'mkdir -p ../root/pt/' + label + '/'    
        #os.system(makeRootPtDir)
        #makeRootNpuDir = 'mkdir -p ../root/npu/' + label + '/'    
        #os.system(makeRootNpuDir)
	
        #makePdfEtaDir = 'mkdir -p ../pdfs/eta/' + label + '/'    
        #os.system(makePdfEtaDir)
        #makePdfPtDir = 'mkdir -p ../pdfs/pt/' + label + '/'    
        #os.system(makePdfPtDir)
        #makePdfNpuDir = 'mkdir -p ../pdfs/npu/' + label + '/'    
        #os.system(makePdfNpuDir)
	    
		pyRun = 'python -u ../python/jetvalidator.py -p ' + path + ' -k ' + keyword + ' &>../logs/jv_' + label + '_log.txt \n'
		print "", pyRun
    
        #os.system(pyRun)