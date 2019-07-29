import time


t1 = time.time()

geneIdNetworkFile = open('PPI_Network.txt', 'r')
PpiNetworkPredictionsFile = open('PPI_Network_Prism_Predictions.txt', 'w')
prismPredictionsFile = open('PrismPredictions.txt', 'r')
prismPredictions = prismPredictionsFile.readlines();
prismPredictionsFile.close()
alternativeConformationsFile = open('clusteredGeneIDToPDBMapping.txt', 'r')
alternativeConformations = alternativeConformationsFile.readlines();
alternativeConformationsFile.close()



for PPI in geneIdNetworkFile:
        
        splittedPPI = PPI.split(' ')
        protein1 = splittedPPI[0]
        protein2 = splittedPPI[2][:-2]
        #print(protein1)
        #print(protein2)
        alternativeConfs1 = []
        alternativeConfs2 = []
        found1 = False
        found2 = False

        for geneId in alternativeConformations:
                
                splittedGeneId = geneId.split('\t')
                #print(splittedGeneId[1])
                
                if len(splittedGeneId) > 2:
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
                        
                        PpiNetworkPredictionsFile.write(protein1 + '-' + protein2 + '\n')
                        for alternativeConf1 in alternativeConfs1:
                                PpiNetworkPredictionsFile.write(alternativeConf1 + '\t')
                        PpiNetworkPredictionsFile.write('\n')
                        for alternativeConf2 in alternativeConfs2:
                                PpiNetworkPredictionsFile.write(alternativeConf2 + '\t')
                        PpiNetworkPredictionsFile.write('\n')
                        
                                        
                        for alternativeConf1 in alternativeConfs1:
                                for alternativeConf2 in alternativeConfs2:
                                        for index in range(4, len(prismPredictions)):
                                                splittedPrediction = prismPredictions[index].split()
                                                #print(splittedPrediction[0])
                                                #print(splittedPrediction[1])
                                                #print(alternativeConf1)
                                                #print(alternativeConf2)
                                                if splittedPrediction[0].upper() == alternativeConf1.upper() and splittedPrediction[1].upper() == alternativeConf2.upper():
                                                        PpiNetworkPredictionsFile.write(prismPredictions[index])
                        break

geneIdNetworkFile.close()
PpiNetworkPredictionsFile.close()
t2 = time.time()
print('Elapsed time = %f' %(t2-t1))
