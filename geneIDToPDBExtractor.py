import os,sys, time

t1 = time.time()

uniprotMappingFileDirectory = 'idmapping_selected.tab'
uniprotMappingFile = open(uniprotMappingFileDirectory, 'r')

geneIDToPDBMapFileDirectory = 'geneIDToPDBMapping.txt'
geneIDToPDBMapFile = open(geneIDToPDBMapFileDirectory, 'w')

for uniprotEntry in uniprotMappingFile:
	splittedUniprotEntry = uniprotEntry.strip().split('\t')
	if len(splittedUniprotEntry[2]) > 0:
		if len(splittedUniprotEntry[5]) > 0:
			splittedGeneID = splittedUniprotEntry[2].split(';')
			if len(splittedGeneID) > 1:
				for geneID in splittedGeneID:
					geneIDToPDBMapFile.write('%s\t%s\n' %(geneID.strip(), splittedUniprotEntry[5]))
			else:
				geneIDToPDBMapFile.write('%s\t%s\n' %(splittedUniprotEntry[2], splittedUniprotEntry[5]))
		
		

uniprotMappingFile.close()
geneIDToPDBMapFile.close()
t2 = time.time()
print('Elapsed time = %f' %(t2-t1))
