#!/usr/bin/env python 
import argparse
import pandas as pd
from ete3 import NCBITaxa

ncbi = NCBITaxa()

def make_new_row(taxid,columns):
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    lineagenames = ncbi.get_taxid_translator(lineage)
    new_dict = {"taxID":taxid}
    for c in columns:
        new_dict[c] = "missing"
    for tid in lineage:
        new_dict[lineage2ranks[tid]] = lineagenames[tid]

    new_row = {k: new_dict[k] for k in ["taxID"] + columns}
    return(new_row)


##------------------------------------------
# COMMAND LINE ARGUMENTS\

parser = argparse.ArgumentParser(description='Extend bracken file')

parser.add_argument('--bracken_file', dest='bracken_file', action='store', required=True,help='bracken output file to extend (usually a .backen file)')

parser.add_argument('--output', dest='bracken_out', action='store', required=True,help='output file of this script - extended bracken file')    

args = parser.parse_args()

bracken_file = args.bracken_file
barcken_out = args.bracken_out

##------------------------------------------
df = pd.read_csv(bracken_file,sep="\t")

columns = ["superkingdom","phylum","class","order","family","genus","species"]

tax_rows = []

for taxid in df["taxonomy_id"]:
    try:
        tax_rows.append(make_new_row(taxid,columns))
    except:
        print("Failed for :"+str(taxid))

    
tax_df = pd.DataFrame.from_dict(tax_rows,orient="columns")

mdf = df.merge(tax_df,left_on = "taxonomy_id",right_on="taxID")

mdf = mdf.drop(['taxID',"name","taxonomy_lvl"],axis=1)

mdf.to_csv(barcken_out,sep="\t",index=False)