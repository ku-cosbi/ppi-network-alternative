#################
## Created by Engin Cukuroglu
#################

import os,sys,time
from multiprocessing import Queue, Process
import multiprocessing

def generateMonomerPDBFilesWork(taskQueue_pdbDict, allPDBFilesDirectory, generateMonomerPDBFilesErrorQueue):
	while True:
		pdbName, pdbDict = taskQueue_pdbDict.get()
		taskQueue_pdbDict.task_done()
		if pdbName is None:
			break
		pdbFileDirectory = '%s/%s.pdb' %(allPDBFilesDirectory, pdbName)
		if os.path.exists(pdbFileDirectory):
			createPDBFileList = []
			for chainInfo in pdbDict:
				createPDBFileList.append([chainInfo, '', pdbDict[chainInfo][1]])
			pdbFile = open(pdbFileDirectory, 'r')
			pdbResidueDict = {}
			pdbAtomDict = {}
			pdbChainIDDict = {}
			while 1:
				pdbLine = pdbFile.readline()
				if pdbLine == '' or pdbLine[0:3] == 'END':
					break
				if pdbLine[0:4] == 'ATOM':
					chainID = pdbLine[21]
					pdbChainIDDict[chainID] = 1
					resNo = pdbLine[22:26].strip()
					resICode = pdbLine[26]
					atomAlternativeLocationIndicator = pdbLine[16]
					resType = pdbLine[17:20].strip()
					atomType = pdbLine[12:16].strip()
					resDictKey = '%s_%s' %(resNo, chainID)
					resDictValue = '%s_%s_%s' %(resType, resNo, chainID)
					atomDictKey = '%s_%s_%s_%s' %(resType, resNo, chainID, atomType)
					if not resDictKey in pdbResidueDict:
						pdbResidueDict[resDictKey] = [resDictValue, resICode, atomAlternativeLocationIndicator]
						pdbAtomDict[atomDictKey] = 1
						for tempPDBFileProperty in createPDBFileList:
							if chainID == 'all':
								tempPDBFileProperty[1] = '%s%s' %(tempPDBFileProperty[1], pdbLine)
							elif chainID in tempPDBFileProperty[0]:
								tempPDBFileProperty[1] = '%s%s' %(tempPDBFileProperty[1], pdbLine)
					else:
						if pdbResidueDict[resDictKey][0] == resDictValue:
							if pdbResidueDict[resDictKey][1] == resICode:
								if not atomDictKey in pdbAtomDict:
									pdbAtomDict[atomDictKey] = 1
									for tempPDBFileProperty in createPDBFileList:
										if chainID in tempPDBFileProperty[0]:
											tempPDBFileProperty[1] = '%s%s' %(tempPDBFileProperty[1], pdbLine)
			pdbFile.close()
			for tempPDBFileProperty in createPDBFileList:
				chainInfo = tempPDBFileProperty[0]
				chainExistenceStatus = 1
				for chainID in chainInfo:
					if chainID == 'all':
						if not len(pdbChainIDDict) > 0:
							chainExistenceStatus = 0
							errorText = '%s\t2\t%s\tChain info does not present in PDB file.\n' %(pdbDict[chainID][0], pdbName)
							generateMonomerPDBFilesErrorQueue.put(errorText)
							break
					elif not chainID in pdbChainIDDict:
						chainExistenceStatus = 0
						#errorText = '%s\t2\t%s\tChain %s does not present in PDB file.\n' %(pdbDict[chainID][0], pdbName, chainID)
						errorText = chainID + ' does not present in ' + pdbName + ' PDB file.\n'
						generateMonomerPDBFilesErrorQueue.put(errorText)
						break
				if chainExistenceStatus == 1:
					tempPDBFile = open(tempPDBFileProperty[2], 'w')
					tempPDBFile.write(tempPDBFileProperty[1])
					tempPDBFile.close()
		else:
			for chainID in pdbDict:
				errorText = '%s\t2\t%s\tPDB file does not present in %s.\n' %(pdbDict[chainID][0], pdbName, allPDBFilesDirectory)
				generateMonomerPDBFilesErrorQueue.put(errorText)
		
def generateMonomerPDBFiles(allPDBFilesDirectory, generatedPDBFilesDirectory, pdbDict, generateMonomerPDBFilesLogFileDirectory, numberOfProcesses):
	if not os.path.exists(allPDBFilesDirectory):
		sys.exit('\nThe %s does not exist.\n' %(allPDBFilesDirectory))
	if not os.path.exists(generatedPDBFilesDirectory):
		os.system('mkdir %s' %(generatedPDBFilesDirectory))

	taskDict = {}
	for pdbName in pdbDict:
		if len(pdbDict[pdbName]) > 0:
			if not pdbName in taskDict:
				taskDict[pdbName] = {}
			for pdbChain in pdbDict[pdbName]:
				if pdbChain == 'all':
					pdbMonomer = '%s' %(pdbName)
				else:
					pdbMonomer = '%s_%s' %(pdbName, pdbChain)
				pdbMonomerFileDirectory = '%s/%s.pdb' %(generatedPDBFilesDirectory, pdbMonomer)
				if not os.path.exists(pdbMonomerFileDirectory):
					taskDict[pdbName][pdbChain] = [pdbMonomer, pdbMonomerFileDirectory]
	
	taskQueue_pdbDict = multiprocessing.JoinableQueue()
	generateMonomerPDBFilesErrorQueue = multiprocessing.Manager().Queue()
	#generateMonomerPDBFilesErrorQueue = multiprocessing.freeze_support()
	generateMonomerPDBFilesWorker = [Process(target=generateMonomerPDBFilesWork, args=(taskQueue_pdbDict, allPDBFilesDirectory, generateMonomerPDBFilesErrorQueue)) for i in range(numberOfProcesses)]
	for tempWorkers in generateMonomerPDBFilesWorker:
		tempWorkers.start()

	for pdbName in taskDict:
		if len(taskDict[pdbName]) > 0:
			taskQueue_pdbDict.put([pdbName, taskDict[pdbName]])
	
	for i in range(numberOfProcesses):
		taskQueue_pdbDict.put([None, None])

	taskQueue_pdbDict.join()
	generateMonomerPDBFilesErrorQueue.put('Done')
	generateMonomerPDBFilesLogFile = open(generateMonomerPDBFilesLogFileDirectory, 'w')
	while True:
		errorEntry = generateMonomerPDBFilesErrorQueue.get()
		if errorEntry == 'Done':
			break
		generateMonomerPDBFilesLogFile.write('%s' %(errorEntry))
	generateMonomerPDBFilesLogFile.close()

def mainPDBMonomerExtractor(allPDBFilesDirectory, generatedPDBFilesDirectory, geneIDToPDBMappingFileDirectory, generateMonomerPDBFilesLogFileDirectory, numberOfProcesses):
	print('\n* PDB MONOMER EXTRACTOR STARTED*\n')
	print('Time stamp: %s' %(time.asctime()))
	t1 = time.time()
	
	if not os.path.exists(geneIDToPDBMappingFileDirectory):
		sys.exit('\nThe %s does not exist.\n' %(geneIDToPDBMappingFileDirectory))
	
	pdbDict = {}
	geneIDToPDBMappingFile = open(geneIDToPDBMappingFileDirectory, 'r')
	for geneIDLine in geneIDToPDBMappingFile:
		splittedGeneIDLine = geneIDLine.strip().split('\t')
		if len(splittedGeneIDLine) > 1:
			splittedPDBEntries = splittedGeneIDLine[1].split(';')
			for pdbEntry in splittedPDBEntries:
				pdbEntryList = pdbEntry.strip().split(':')
				if len(pdbEntryList) > 1:
					chainID = pdbEntryList[1]
				else:
					chainID = 'all'
				if not pdbEntryList[0] in pdbDict:
					pdbDict[pdbEntryList[0]] = {}
					pdbDict[pdbEntryList[0]][chainID] = 1
				else:
					pdbDict[pdbEntryList[0]][chainID] = 1
					
	geneIDToPDBMappingFile.close()
			
	generateMonomerPDBFiles(allPDBFilesDirectory, generatedPDBFilesDirectory, pdbDict, generateMonomerPDBFilesLogFileDirectory, numberOfProcesses)
			
	t2 = time.time()
	print('\nElapsed time = %f seconds\n' %(t2-t1))
	print('Time stamp: %s' %(time.asctime()))
	print('\n* PDB MONOMER EXTRACTOR COMPLETED*\n')
	
