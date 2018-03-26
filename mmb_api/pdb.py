import requests
import os
import shutil
import json

class MmbPdb(object):
    """Wrapper class for the MMB group PDB REST API.
    This class is a wrapper for the PDB (http://www.rcsb.org/pdb/home/home.do)
    mirror of the MMB group REST API (http://mmb.irbbarcelona.org/api/)
    """
    def __init__(self, url=None):
        self.url = url if url else "http://mmb.irbbarcelona.org/api/pdb/"

    def get_pdb(self, pdb_code, output_pdb_path, filt="filter=/1&group=ATOM"):
        """
        Writes the PDB file content of pdb_code
        to output_pdb_path
        """
        if os.path.isfile(pdb_code):
            shutil.copy(pdb_code, output_pdb_path)
        else:
            url = self.url+pdb_code.lower()+"/coords/?"+filt
            pdb_string = requests.get(url).content
            with open(output_pdb_path, 'w') as pdb_file:
                pdb_file.write(pdb_string)

    def get_cluster_pdb_codes(self, pdb_code, cluster="cl-90"):
        """
        Returns the list of pdb_codes of the selected cluster. 'pdb_code' format: [pdbCode] 
        """
        pdb_codes = set()
        url = self.url+pdb_code.lower()+'/clusters/'+cluster+".json"
        cluster = json.loads(requests.get(url).content)['clusterMembers']
        for elem in cluster:
            pdb_codes.add(elem['_id'].lower())
        return pdb_codes

    def get_cluster_pdb_code_chains(self, pdb_code, cluster="cl-90"):
        """
        Returns the list of PDB codes and chains of the selected cluster. 'pdb_code' format: [pdbCode](_[chain])
        """
        pdb_code_chains = set()
        url = self.url+pdb_code+'/clusters/'+cluster+".json"
        content = json.loads(requests.get(url).content)
	try:
	    cluster = json.loads(requests.get(url).content)['clusterMembers']
	except:
	    return pdb_code_chains
	
        for elem in cluster:
            for chain in elem['chainIds']:
                pdb_code_chains.add(chain)
        return pdb_code_chains
