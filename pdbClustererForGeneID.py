#################
## Created by Engin Cukuroglu
#################

import os,sys,time
import multiprocessing
from multiprocessing import Queue, Process
import networkx as nx

def pdbSimilarityNetworkBuilder(pdbMonomerList, similarityDict):
	G = nx.Graph()
	for pdbMonomer in pdbMonomerList:	
		G.add_node(pdbMonomer)
	for i in range(len(pdbMonomerList)):
		pdb_1 = pdbMonomerList[i]
		for j in range(i+1,len(pdbMonomerList)):
			pdb_2 = pdbMonomerList[j]
			if pdb_1 in similarityDict:
				if pdb_2 in similarityDict[pdb_1]:
					G.add_edge(pdb_1, pdb_2, weight=similarityDict[pdb_1][pdb_2])
	return G
	

def geneIDTOPDBMappingFileReader(geneIDToPDBMappingFileDirectory):
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
	return geneIDDict

def TMAlignGeneIDResultFileReader(TMAlignGeneIDResultFileDirectory, maxRMSD, minSimilarity):
	similarityDict = {}
	TMAlignGeneIDResultFile = open(TMAlignGeneIDResultFileDirectory, 'r')
	for pdbInfoLine in TMAlignGeneIDResultFile:
		splittedPDBInfoLine = pdbInfoLine.strip().split('\t')
		RMSD = float(splittedPDBInfoLine[6])
		if RMSD <= maxRMSD:
			similarity = float(splittedPDBInfoLine[5])
			if similarity >= minSimilarity:
				pdb_1 = splittedPDBInfoLine[0]
				pdb_2 = splittedPDBInfoLine[1]
				if not pdb_1 in similarityDict:
					similarityDict[pdb_1] = {}
					similarityDict[pdb_1][pdb_2] = similarity
				else:
					similarityDict[pdb_1][pdb_2] = similarity
				if not pdb_2 in similarityDict:
					similarityDict[pdb_2] = {}
					similarityDict[pdb_2][pdb_1] = similarity
				else:
					similarityDict[pdb_2][pdb_1] = similarity
				
	TMAlignGeneIDResultFile.close()
	return similarityDict
	


def pdbClustererForGeneIDRunWork(taskQueue_geneID, geneIDResultQueue, generatedPDBFilesDirectory, minResidueSize, maxRMSD, minSimilarity):
	while True:
		geneID, pdbIDList, TMAlignGeneIDResultFileDirectory = taskQueue_geneID.get()
		taskQueue_geneID.task_done()
		if geneID is None:
			break
		pdbMonomerSizeDict = {}
		pdbMonomerList = []
		for pdbID in pdbIDList:
			tempPDBMonomerFileDirectory = '%s/%s.pdb' %(generatedPDBFilesDirectory, pdbID)
			if os.path.exists(tempPDBMonomerFileDirectory):
				tempPDBSize = 0
				tempPDBMonomerFile = open(tempPDBMonomerFileDirectory, 'r')
				for pdbLine in tempPDBMonomerFile:
					if pdbLine[0:4] == 'ATOM':
						if pdbLine[13:15] == 'CA':
							tempPDBSize = tempPDBSize + 1
				tempPDBMonomerFile.close()
				if tempPDBSize >= minResidueSize:
					pdbMonomerSizeDict[pdbID] = tempPDBSize
					pdbMonomerList.append(pdbID)
		candidateMonomerNumber = len(pdbMonomerSizeDict)
		if candidateMonomerNumber != 0:
			if candidateMonomerNumber == 1:
				geneIDResultText = '%s\t%d' %(geneID, candidateMonomerNumber)
				for pdbMonomer in pdbMonomerSizeDict:
					geneIDResultText = '%s\t%s' %(geneIDResultText, pdbMonomer)
				geneIDResultText = '%s\n' %(geneIDResultText)
				geneIDResultQueue.put(geneIDResultText)
			elif not os.path.exists(TMAlignGeneIDResultFileDirectory):
				geneIDResultText = '%s\t%d' %(geneID, candidateMonomerNumber)
				pdbMonomerCounter = 0
				for pdbMonomer in pdbMonomerSizeDict:
					if pdbMonomerCounter == 0:
						geneIDResultText = '%s\t%s' %(geneIDResultText, pdbMonomer)
					else:
						geneIDResultText = '%s,%s' %(geneIDResultText, pdbMonomer)
					pdbMonomerCounter = pdbMonomerCounter + 1
				geneIDResultText = '%s\n' %(geneIDResultText)
				geneIDResultQueue.put(geneIDResultText)
			else:
				similarityDict = TMAlignGeneIDResultFileReader(TMAlignGeneIDResultFileDirectory, maxRMSD, minSimilarity)
				pdbSimilarityNetwork = pdbSimilarityNetworkBuilder(pdbMonomerList, similarityDict)
				nodesToCombine = 1
				while nodesToCombine > 0:
					maxEdgeWeight = -1
					maxEdgeNodes = []
					nodesToCombine = 0
					for edge in pdbSimilarityNetwork.edges_iter():
						edgeWeight = pdbSimilarityNetwork[edge[0]][edge[1]]['weight']
						if edgeWeight == 1:
							maxEdgeWeight = edgeWeight
							maxEdgeNodes = [edge[0], edge[1]]
							nodesToCombine = 1
							break
						if edgeWeight >= minSimilarity:
							if edgeWeight > maxEdgeWeight:
								maxEdgeWeight = edgeWeight
								maxEdgeNodes = [edge[0], edge[1]]
								nodesToCombine = 1
					
					if nodesToCombine > 0:
						pdb_1 = maxEdgeNodes[0]	
						pdb_2 = maxEdgeNodes[1]
						selectedMonomer = pdb_1
						discardedMonomer = pdb_2
						if pdbMonomerSizeDict[pdb_1] < pdbMonomerSizeDict[pdb_2]:	
							selectedMonomer = pdb_2
							discardedMonomer = pdb_1
						selectedMonomer_neighborDict = {}
						discardedMonomer_neighborDict = {}
						for neighbor in pdbSimilarityNetwork.neighbors(selectedMonomer):
							if not neighbor == discardedMonomer:
								selectedMonomer_neighborDict[neighbor] = pdbSimilarityNetwork[selectedMonomer][neighbor]['weight']
						for neighbor in pdbSimilarityNetwork.neighbors(discardedMonomer):
							if not neighbor == selectedMonomer:
								discardedMonomer_neighborDict[neighbor] = pdbSimilarityNetwork[discardedMonomer][neighbor]['weight']
						for pdbMonomer in discardedMonomer_neighborDict:
							tempWeight = discardedMonomer_neighborDict[pdbMonomer]
							if pdbMonomer in selectedMonomer_neighborDict:
								if selectedMonomer_neighborDict[pdbMonomer] < discardedMonomer_neighborDict[pdbMonomer]:
									selectedMonomer_neighborDict[pdbMonomer] = discardedMonomer_neighborDict[pdbMonomer]
									pdbSimilarityNetwork.add_edge(selectedMonomer, pdbMonomer, weight=discardedMonomer_neighborDict[pdbMonomer])
								else:
									pdbSimilarityNetwork.add_edge(selectedMonomer, pdbMonomer, weight=discardedMonomer_neighborDict[pdbMonomer])
						pdbSimilarityNetwork.remove_node(discardedMonomer)	
				geneIDResultText = '%s\t%d' %(geneID, pdbSimilarityNetwork.number_of_nodes())
				pdbMonomerCounter = 0
				for pdbNode in pdbSimilarityNetwork.nodes_iter():
					if pdbMonomerCounter == 0:
						geneIDResultText = '%s\t%s' %(geneIDResultText, pdbNode)
					else:
						geneIDResultText = '%s,%s' %(geneIDResultText, pdbNode)
					pdbMonomerCounter = pdbMonomerCounter + 1
				geneIDResultText = '%s\n' %(geneIDResultText)
				geneIDResultQueue.put(geneIDResultText)
				pdbSimilarityNetwork.clear()

def mainPDBClustererForGeneID(geneIDToPDBMappingFileDirectory, TMAlignGeneIDResultsFileDirectory, clusteredGeneIDToPDBMappingFileDirectory, generatedPDBFilesDirectory, minResidueSize, maxRMSD, minSimilarity, numberOfProcesses):
	print('\n* PDB CLUSTERER FOR GENE ID STARTED *\n')
	print('Time stamp: %s' %(time.asctime()))
	t1 = time.time()

	if not os.path.exists(geneIDToPDBMappingFileDirectory):
		sys.exit('\nThe %s does not exist.' %(geneIDToPDBMappingFileDirectory))

	geneIDDict = geneIDTOPDBMappingFileReader(geneIDToPDBMappingFileDirectory)
	
	taskQueue_geneID = multiprocessing.JoinableQueue()
	geneIDResultQueue = multiprocessing.Queue()
	pdbClustererRunWorker = [Process(target=pdbClustererForGeneIDRunWork, args=(taskQueue_geneID, geneIDResultQueue, generatedPDBFilesDirectory, minResidueSize, maxRMSD, minSimilarity)) for i in range(numberOfProcesses)]

	for tempWorkers in pdbClustererRunWorker:
		tempWorkers.start()
	
	for geneID in geneIDDict:
		numberOfPDBIDsForGeneID = len(geneIDDict[geneID])
		if numberOfPDBIDsForGeneID > 0:
			TMAlignGeneIDResultFileDirectory = '%s/%s_TMAlignResults.txt' %(TMAlignGeneIDResultsFileDirectory, geneID)
			taskQueue_geneID.put([geneID, geneIDDict[geneID], TMAlignGeneIDResultFileDirectory])

	for i in range(numberOfProcesses):
		taskQueue_geneID.put([None, None, None])

	taskQueue_geneID.join()	
	
	geneIDResultQueue.put('Done')
	clusteredGeneIDToPDBMappingFile = open(clusteredGeneIDToPDBMappingFileDirectory, 'w')
	while True:
		geneIDEntry = geneIDResultQueue.get()
		if geneIDEntry == 'Done':
			break
		clusteredGeneIDToPDBMappingFile.write('%s' %(geneIDEntry))
	clusteredGeneIDToPDBMappingFile.close()
	
	t2 = time.time()
	print('Elapsed time = %f' %(t2-t1))
	print('Time stamp: %s' %(time.asctime()))
	print('\n* PDB CLUSTERER FOR GENE ID COMPLETED *\n')
