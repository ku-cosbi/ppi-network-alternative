import os,sys,time
import multiprocessing
from multiprocessing import Queue, Process
import subprocess
import socket

def TMAlignResultReader(outputMessage, structureName_A, structureName_B):
	splittedOutputMessage = outputMessage.split('\n')
	structureLength_A = -1
	structureLength_B = -1
	alignedLength = -1
	RMSD = -1
	TMScoreNormalized_A = -1
	TMScoreNormalized_B = -1
	for outputLine in splittedOutputMessage:
		if 'Length of structure A:' in outputLine:
			try:
				structureLength_A = int(outputLine.replace('Length of structure A:', '').replace('residues','').strip())
			except:
				structureLength_A = -1
		if 'Length of structure B:' in outputLine:
			try:
				structureLength_B = int(outputLine.replace('Length of structure B:', '').replace('residues','').strip())
			except:
				structureLength_B = -1
		if 'Aligned length=' in outputLine:
			splittedOutputLine = outputLine.split(',')
			try:
				alignedLength = int(splittedOutputLine[0].replace('Aligned length=', '').strip())
			except:
				alignedLength = -1
			try:
				RMSD = float(splittedOutputLine[1].replace('RMSD=', '').strip())
			except:
				RMSD = -1
			
		if 'TM-score=' in outputLine:
			if 'normalized by length of structure A' in outputLine:
				try:
					TMScoreNormalized_A = float(outputLine.replace('TM-score=','').split('(if')[0].strip())
				except:
					TMScoreNormalized_A = -1
		if 'TM-score=' in outputLine:
			if 'normalized by length of structure B' in outputLine:
				try:
					TMScoreNormalized_B = float(outputLine.replace('TM-score=','').split('(if')[0].strip())
				except:
					TMScoreNormalized_B = -1
	return structureLength_A, structureLength_B, alignedLength, RMSD, TMScoreNormalized_A, TMScoreNormalized_B

def TMAlignRunWork(taskQueue_TMAlignDict, TMAlignGeneIDResultsFileDirectory, TMAlignResultsFileDirectory, TMAlignRunFileDirectory, generatedPDBFilesDirectory, TMAlignErrorQueue):
	while True:
		geneID, pdbIDList = taskQueue_TMAlignDict.get()
		taskQueue_TMAlignDict.task_done()
		if geneID is None:
			break
		TMAlignGeneIDResultFileDirectory = '%s/%s_TMAlignResults.txt' %(TMAlignGeneIDResultsFileDirectory, geneID)
		tempTMAlignGeneIDResultFileDirectory = '%s/%s_TMAlignResults.txt_temp' %(TMAlignGeneIDResultsFileDirectory, geneID)
		if os.path.exists(tempTMAlignGeneIDResultFileDirectory):
			os.system('rm %s' %(tempTMAlignGeneIDResultFileDirectory))
		alignmentDict = {}
		if os.path.exists(TMAlignGeneIDResultFileDirectory):
			TMAlignGeneIDResultFile = open(TMAlignGeneIDResultFileDirectory, 'r')
			tempTMAlignGeneIDResultFile = open(tempTMAlignGeneIDResultFileDirectory, 'w')
			for alignmentLine in TMAlignGeneIDResultFile:
				tempTMAlignGeneIDResultFile.write(alignmentLine)
				splittedAlignmentLine = alignmentLine.strip().split('\t')
				alignmentKey = '%s_vs_%s' %(splittedAlignmentLine[0], splittedAlignmentLine[1])
				alignmentDict[alignmentKey] = alignmentLine
			TMAlignGeneIDResultFile.close()
			tempTMAlignGeneIDResultFile.close()
		pdbIDList = sorted(pdbIDList)
		numberOfPDBIDs = len(pdbIDList)
		tempTMAlignGeneIDResultFile = open(tempTMAlignGeneIDResultFileDirectory, 'a')
		for i in range(numberOfPDBIDs):
			pdbID_1 = pdbIDList[i]
			pdbFileDirectory_1 = '%s/%s.pdb' %(generatedPDBFilesDirectory, pdbID_1)
			if os.path.exists(pdbFileDirectory_1):
				for j in range(i+1, numberOfPDBIDs):
					pdbID_2 = pdbIDList[j]
					pdbFileDirectory_2 = '%s/%s.pdb' %(generatedPDBFilesDirectory, pdbID_2)
					if os.path.exists(pdbFileDirectory_2):
						alignmentKey = '%s_vs_%s' %(pdbID_1, pdbID_2)
						if not alignmentKey in alignmentDict:
							outputFileDirectory = '%s/%s_TMAlignOutput.txt' %(TMAlignResultsFileDirectory, alignmentKey)
							outputFileDirectory_all = '%s/%s_TMAlignOutput.txt_all' %(TMAlignResultsFileDirectory, alignmentKey)
							stderr = socket.socketpair()
							stderr[0].settimeout(0.01)
							stdout = socket.socketpair()
							stdout[0].settimeout(0.01)
							proc = subprocess.Popen([TMAlignRunFileDirectory, '-A', pdbFileDirectory_1, '-B', pdbFileDirectory_2, '-o', outputFileDirectory], stdout=stdout[1], stderr=stderr[1], close_fds=True)
							errorMessage = u''
							outputMessage = u''
							while True:
								proc.poll()
								try:
									outtmp = stdout[0].recv(4096)
								except socket.timeout as exc:
									outtmp = ''
								outputMessage = outputMessage + outtmp
								try:
									errtmp = stderr[0].recv(4096)
								except socket.timeout as exc:
									errtmp = ''
								errorMessage = errorMessage + errtmp
								if len(errorMessage) > 4096:
									proc.kill()
									proc.wait()
								if proc.returncode != None:
									returnCode = proc.returncode
									break
							os.system('rm %s' %(outputFileDirectory))
							os.system('rm %s' %(outputFileDirectory_all))
							if errorMessage:
								TMAlignErrorText = 'Error: TMAlign: %s %s.\n' %(alignmentKey, errorMessage.strip())
								TMAlignErrorQueue.put(TMAlignErrorText)
							if outputMessage:
								#resultFileDirectory = '%s/%s_TMAlignResult.txt'  %(TMAlignResultsFileDirectory, alignmentKey)
								#resultFile = open(resultFileDirectory, 'w')
								#resultFile.write(outputMessage)
								#resultFile.close()
								structureLength_A, structureLength_B, alignedLength, RMSD, TMScoreNormalized_A, TMScoreNormalized_B = TMAlignResultReader(outputMessage, pdbID_1, pdbID_2)
								minStructureLength = structureLength_A
								if minStructureLength > structureLength_B:
									minStructureLength = structureLength_B
								alignmentRatio = float(alignedLength)/minStructureLength
								tempTMAlignGeneIDResultFile.write('%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.5f\t%.5f\n' %(pdbID_1, pdbID_2, structureLength_A, structureLength_B, alignedLength, alignmentRatio, RMSD, TMScoreNormalized_A, TMScoreNormalized_B))
								
			else:
				TMAlignErrorText = 'Error: %s does not exist.\n' %(pdbFileDirectory_1)
				TMAlignErrorQueue.put(TMAlignErrorText)
		tempTMAlignGeneIDResultFile.close()
		os.system('mv %s %s' %(tempTMAlignGeneIDResultFileDirectory, TMAlignGeneIDResultFileDirectory))

def mainTMAlignDriverForStructures(geneIDToPDBMappingFileDirectory, TMAlignGeneIDResultsFileDirectory, TMAlignResultsFileDirectory, TMAlignRunFileDirectory, generatedPDBFilesDirectory, TMAlignDriverForStructuresLogFileDirectory, numberOfProcesses):
	print('\n* TMAlign DRIVER FOR STRUCTURES STARTED *\n')
	print('Time stamp: %s' %(time.asctime()))
	t1 = time.time()
	
	if not os.path.exists(geneIDToPDBMappingFileDirectory):
		sys.exit('\nThe %s does not exist.\n' %(geneIDToPDBMappingFileDirectory))
	if not os.path.exists(TMAlignGeneIDResultsFileDirectory):
		os.system('mkdir %s' %(TMAlignGeneIDResultsFileDirectory))
	if not os.path.exists(TMAlignResultsFileDirectory):
		os.system('mkdir %s' %(TMAlignResultsFileDirectory))
	if not os.path.exists(TMAlignRunFileDirectory):
		sys.exit('\nThe TMAlign run file does not exist in %s' %(TMAlignRunFileDirectory))
	
	geneIDDict = {}
	geneIDToPDBMappingFile = open(geneIDToPDBMappingFileDirectory, 'r')
	for geneIDLine in geneIDToPDBMappingFile:
		splittedGeneIDLine = geneIDLine.strip().split('\t')
		if len(splittedGeneIDLine) > 1:
			geneID = splittedGeneIDLine[0]
			splittedPDBEntries = splittedGeneIDLine[1].split(';')
			
			if not geneID in geneIDDict: 
				geneIDDict[geneID] = []
			for pdbEntry in splittedPDBEntries:
				pdbEntryList = pdbEntry.strip().split(':')
				if len(pdbEntryList) > 1:
					pdbName = '%s_%s' %(pdbEntryList[0], pdbEntryList[1])
				else:
					pdbName = '%s' %(pdbEntryList[0])
				geneIDDict[geneID].append(pdbName)
	geneIDToPDBMappingFile.close()

	taskQueue_TMAlignDict = multiprocessing.JoinableQueue()
	TMAlignErrorQueue = multiprocessing.Manager().Queue()
	TMAlignRunWorker = [Process(target=TMAlignRunWork, args=(taskQueue_TMAlignDict, TMAlignGeneIDResultsFileDirectory, TMAlignResultsFileDirectory, TMAlignRunFileDirectory, generatedPDBFilesDirectory, TMAlignErrorQueue)) for i in range(numberOfProcesses)]
	for tempWorkers in TMAlignRunWorker:
		tempWorkers.start()
	numberOfTMAlignRuns = 0
	for geneID in geneIDDict:
		numberOfPDBIDsForGeneID = len(geneIDDict[geneID])
		if numberOfPDBIDsForGeneID > 1:
			taskQueue_TMAlignDict.put([geneID, geneIDDict[geneID]])
			numberOfTMAlignRuns = numberOfTMAlignRuns + (numberOfPDBIDsForGeneID * (numberOfPDBIDsForGeneID + 1))/2 
	for i in range(numberOfProcesses):
		taskQueue_TMAlignDict.put([None, None])
	print('Number of TMAlign runs = %d' %(numberOfTMAlignRuns))
	taskQueue_TMAlignDict.join()

	TMAlignErrorQueue.put('Done')
	TMAlignDriverForStructuresLogFile = open(TMAlignDriverForStructuresLogFileDirectory, 'w')
	while True:
		errorEntry = TMAlignErrorQueue.get()
		if errorEntry == 'Done':
			break
		TMAlignDriverForStructuresLogFile.write('%s' %(errorEntry))
	TMAlignDriverForStructuresLogFile.close()

	t2 = time.time()
	print('Elapsed time = %f' %(t2-t1))
	print('Time stamp: %s' %(time.asctime()))
	print('\n* TMAlign DRIVER FOR STRUCTURES COMPLETED *\n')
