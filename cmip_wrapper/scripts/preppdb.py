#!/usr/bin/python
#  CMIP utility to prepare PDB files for COORFMT=1 or 2
#

import argparse
import re
import sys
from cmip_wrapper.ResLib import ResiduesDataLib

def main():
    parser = argparse.ArgumentParser(
        prog='preppdb',
        description='Utility to prepare PDB files for CMIP COORFMT=2'
    )
    
    parser.add_argument (
        '--noheaders',
        action='store_true',
        dest='nohead',
        help='Remove headers'
    )
    
    parser.add_argument (
        '--nohetats',
        action='store_true',
        dest='nohetats',
        help='Remove HETATM'
    )
    
    parser.add_argument('res_lib',help='CMIP Residue library')
    parser.add_argument('pdb_in',help='Input PDB File')
    parser.add_argument('pdb_out',help='Output PDB File')
    
    args = parser.parse_args()
    
    print ("Settings")
    print ("--------")
    for k,v in vars(args).items():
        print ('{:10}:'.format(k),v)
        
    aaLib = ResiduesDataLib(args.res_lib)
    print ("{} residue/atom pairs loaded from {}".format(aaLib.nres,args.res_lib))

    try:
        inPDB = open(args.pdb_in,'r')
    except OSError:
        print ("#ERROR: loading PDB from {}".format(args.pdb_in))
        sys.exit(1)
        
    try:
        outPDB = open(args.pdb_out,'w')
    except OSError:
        print ("#ERROR: Could not open file {} for writing ",format(args.pdb_out))
        sys.exit(1)

    for line in inPDB:
        line = line.rstrip()
        if not re.match('^ATOM',line) and not re.match('^HETATM',line):
            if not args.nohead:
                outPDB.write(line+"\n")
            continue
        
        if args.nohetats and re.match('^HETATM',line):
            continue

        pre = line[:54]
        nomat = line[12:16]
        if re.match('^[1-9]',nomat):
            nomat = nomat[1:4]+nomat[:1] 
        nomat = nomat.replace(' ','')
        nomr = line[17:20].replace(' ','')
        parms = aaLib.getParams(nomr,nomat)
        if parms.atType == 'X':
            print ("Unknown:",pre)
        outPDB.write("{}{:8.4f}  {}\n".format(pre,parms.charg,parms.atType))
    inPDB.close()
    outPDB.close()

if __name__== '__main__':
    main()