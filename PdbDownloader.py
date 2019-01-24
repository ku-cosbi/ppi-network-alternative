import time
import os.path
import requests
from multiprocessing.pool import ThreadPool


#------------------------------------------------------------
def fetch_url(entry):
        pdb, path, uri = entry
	result = ""
        if not os.path.exists(path):
                r = requests.get(uri, stream=True)
                if r.status_code == 200:
			result = "Downloading " + pdb + ".pdb"
                        with open(path, 'wb') as f:
                                for chunk in r:
                                        f.write(chunk)
		elif r.status_code == 404:
			result = "The requested pdb file (" + pdb + ") was not found on the server."
		
	return result



#--------------------- Main Procedure ------------------------

t1 = time.time()

nodeGeneIdsFile = open('Network_Gene_Ids.txt', 'r')
idMappingFile = open('geneIDToPDBMapping.txt', 'r')
idMapping = idMappingFile.readlines();
idMappingFile.close()
usedIdsFile = open('usedGeneIDToPDBMapping.txt', 'w')
urls = []

for geneId in nodeGeneIdsFile:
        
        geneId = geneId[:-2]
                        
        for geneProducts in idMapping:
                
                splittedgeneProducts = geneProducts.split()
                
                if splittedgeneProducts[0] == geneId:
                    
                        usedIdsFile.write(geneProducts)
 
                        for i in range(1, len(splittedgeneProducts)):
                               
                                pdbId = splittedgeneProducts[i][:4]  
                                pdbPath = 'pdbFiles/' + pdbId + '.pdb'
                                file_url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=" + pdbId
                                          
                                urls.append((pdbId, pdbPath, file_url))

                        break

nodeGeneIdsFile.close()
usedIdsFile.close()

'''for entry in urls:
        fetch_url(entry)'''

results = ThreadPool(8).imap_unordered(fetch_url, urls)

for result in results:
        print(result)
    

t2 = time.time()
print('Elapsed time = %f' %(t2-t1))
