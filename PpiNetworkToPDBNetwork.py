import time


t1 = time.time()

geneIdNetworkFile = open('PPI_Network.txt', 'r')
PdbNetworkFile = open('PDB_Network.txt', 'w')
alternativeConformationsFile = open('clusteredGeneIDToPDBMapping.txt', 'r')
alternativeConformations = alternativeConformationsFile.readlines();
alternativeConformationsFile.close()


for PPI in geneIdNetworkFile:
        
        splittedPPI = PPI.split(' ')
        protein1 = splittedPPI[0]
        protein2 = splittedPPI[2][:-1]
        #print(protein1)
        #print(protein2)
        alternativeConfs1 = []
        alternativeConfs2 = []
        found1 = False
        found2 = False

        for geneId in alternativeConformations:
                
                splittedGeneId = geneId.split('\t')
                #print(splittedGeneId[1])
                alternativeConfs = splittedGeneId[2].split(',')
                
                if splittedGeneId[0] == protein1:
                        for i in range(0, int(splittedGeneId[1])):
                                alternativeConfs1.append(alternativeConfs[i][:4] + alternativeConfs[i][5:6])
                        found1 = True
                        
                elif splittedGeneId[0] == protein2:
                        for i in range(0, int(splittedGeneId[1])):
                                alternativeConfs2.append(alternativeConfs[i][:4] + alternativeConfs[i][5:6])
                        found2 = True
                        
                if found1 and found2:
                        for alternativeConf1 in alternativeConfs1:
                                for alternativeConf2 in alternativeConfs2:
                                        PdbNetworkFile.write(alternativeConf1 + ' ' + alternativeConf2 + '\n')
                                        #print(alternativeConf1 + ' ' + alternativeConf2 + '\n')
                        break

PdbNetworkFile.close()
geneIdNetworkFile.close()
t2 = time.time()
print('Elapsed time = %f' %(t2-t1))
