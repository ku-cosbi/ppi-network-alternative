import time


t1 = time.time()


BianaNetworkFile = open('subnetwork.sif.5', 'r')
geneIdNetworkFile = open('PPI_Network.txt', 'w')
nodeGeneIdsFile = open('Network_Gene_Ids.txt', 'w')
scoresFile = open('guild_scores.txt', 'r')
geneScores = scoresFile.readlines();
scoresFile.close()
geneIds = set()


for edge in BianaNetworkFile:
        node1 = 0
        node2 = 0
        found1 = False
        found2 = False
        splittedEdge = edge.split(' ')

        for geneInfo in geneScores:
                
                splittedGeneInfo = geneInfo.split('\t')
                
                if splittedGeneInfo[0] == splittedEdge[0]:
                        node1 = splittedGeneInfo[5]
                        found1 = True
                        
                elif splittedGeneInfo[0] == splittedEdge[2][:-1]:
                        node2 = splittedGeneInfo[5]
                        found2 = True
                        
                if found1 and found2:
                        break

        if node1.isdigit() and node2.isdigit():
                geneIdNetworkFile.write(str(node1) + ' interaction ' + str(node2) + '\n')
                geneIds.add(node1)
                geneIds.add(node2)
		
for geneId in geneIds:
        nodeGeneIdsFile.write(geneId + '\n')
		

BianaNetworkFile.close()
geneIdNetworkFile.close()
nodeGeneIdsFile.close()
t2 = time.time()
print('Elapsed time = %f' %(t2-t1))
