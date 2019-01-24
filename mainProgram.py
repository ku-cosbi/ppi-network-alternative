#################
## Created by Engin Cukuroglu
#################

from pdbMonomerExtractor import mainPDBMonomerExtractor
from TMAlignDriverForStructures import mainTMAlignDriverForStructures
from pdbClustererForGeneID import mainPDBClustererForGeneID
import os,sys,time

print('\n MAIN PROGRAM STARTED')
print('Time stamp : %s' %(time.asctime()))
t1 = time.time()

timeStamp = '%s' %(time.strftime('%Y_%B_%d'))

archivedFilesDirectory = 'archivedFiles'
fullListFileDirectory = 'fullLists'
allPDBFilesDirectory = 'pdbFiles'
mainLogFileDirectory = 'logFiles'
logFileDirectory = '%s/log_%s' %(mainLogFileDirectory, timeStamp)
archivedPDBFileDirectory = 'archivedPDBFiles'

if not os.path.exists(allPDBFilesDirectory):
	os.system('mkdir %s' %(allPDBFilesDirectory))

if not os.path.exists(archivedPDBFileDirectory):
	os.system('mkdir %s' %(archivedPDBFileDirectory))

if not os.path.exists(archivedFilesDirectory):
	os.system('mkdir %s' %(archivedFilesDirectory))

if not os.path.exists(fullListFileDirectory):
	os.system('mkdir %s' %(fullListFileDirectory))

if not os.path.exists(mainLogFileDirectory):
	os.system('mkdir %s' %(mainLogFileDirectory))

if not os.path.exists(logFileDirectory):
	os.system('mkdir %s' %(logFileDirectory))


generatedPDBFilesDirectory = 'generatedPDBFiles'
geneIDToPDBMappingFileDirectory = 'usedGeneIDToPDBMapping.txt'
generateMonomerPDBFilesLogFileDirectory = '%s/generateMonomerPDBFilesLogs_%s.txt' %(logFileDirectory, timeStamp)
numberOfProcesses = 20
mainPDBMonomerExtractor(allPDBFilesDirectory, generatedPDBFilesDirectory, geneIDToPDBMappingFileDirectory, generateMonomerPDBFilesLogFileDirectory, numberOfProcesses)
print('\n Monomer extraction completed.')

TMAlignGeneIDResultsFileDirectory = 'TMAlignGeneIDResults'
TMAlignResultsFileDirectory = 'TMAlignResults'
TMAlignRunFileDirectory = 'externalTools/TMalign/TMalign'
TMAlignDriverForStructuresLogFileDirectory = '%s/TMAlignDriverForStructuresLogs_%s.txt' %(logFileDirectory, timeStamp)
mainTMAlignDriverForStructures(geneIDToPDBMappingFileDirectory, TMAlignGeneIDResultsFileDirectory, TMAlignResultsFileDirectory, TMAlignRunFileDirectory, generatedPDBFilesDirectory, TMAlignDriverForStructuresLogFileDirectory, numberOfProcesses)
print('\n TMalign step completed.')

#geneIDToPDBMappingFileDirectory = 'usedGeneIDToPDBMappingDeneme.txt'
clusteredGeneIDToPDBMappingFileDirectory = 'clusteredGeneIDToPDBMapping.txt'
minResidueSize = 30
maxRMSD = 2
minSimilarity = 0.95
mainPDBClustererForGeneID(geneIDToPDBMappingFileDirectory, TMAlignGeneIDResultsFileDirectory, clusteredGeneIDToPDBMappingFileDirectory, generatedPDBFilesDirectory, minResidueSize, maxRMSD, minSimilarity, numberOfProcesses)
print('\n Clustering completed.')

t2 = time.time()
print('\nElapsed time = %f seconds\n' %(t2-t1))
print('Time stamp : %s' %(time.asctime()))
print('\n MAIN PROGRAM COMPLETED')
