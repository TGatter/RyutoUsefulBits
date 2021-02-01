#!/usr/bin/env python

#call GTF_File ID_List

import sys, ast, os
import gffutils
import random
	
def transform_func(x):
    # adds some text to the end of transcript IDs
    if 'transcript_id' in x.attributes:
        x.attributes['transcript_id'][0] += '_t'
        return x

#transform=transform_func	

fn = sys.argv[1]
db = gffutils.create_db(fn, dbfn= fn+'.db', force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True, disable_infer_transcripts=True, id_spec={'gene': 'gene_id', 'transcript': "transcript_id"})

#db = gffutils.FeatureDB(fn+'.db', keep_order=True)

keys = [t.id for t in db.features_of_type('gene')]


items = []
for line in open(sys.argv[2]):
    items.append(line.strip().split()[0])
 
for i in items:
    for e in db.children(db[i], level=1):
        print e
            

