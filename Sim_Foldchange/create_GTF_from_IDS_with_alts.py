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

items = []
for line in open(sys.argv[2]):
    items.append(line.strip().split()[1])

for i in items:
    for e in db.children(db[i], level=1):
        print e


modifyable = []
for i in items:
    if len(list(db.children(db[i], level=1))) > 3:
        modifyable.append(i)

modify = random.sample(modifyable, int(len(items)*float(sys.argv[3])))
            
for i, t in enumerate(modify):
    collect = list(db.children(db[t], level=1))
    collect[0][8]['transcript_id'] = "synthetic_"+str(i)
    
    print collect[0]
    collect.pop(0)

    rm_items = random.randint(1, min(max(2, int(len(collect) * 0.1)) , len(collect) - 2))
    for m in range(rm_items):
        collect.remove(random.choice(collect[1:-1]))

    collect.sort(key=lambda x: int(x[3]) ) 

    mod_items = random.randint(0, min(1, len(collect) - 1))
 
    for m in range(mod_items):
        dist = 0
        attempts = 0
        while dist < 20 and attempts < 20:
            attempts = attempts + 1
            pos = random.randint(0, len(collect) - 2)
            dist = collect[pos + 1][3] - collect[pos][4]
        
        if dist < 20:
            continue

        flip = random.randint(0, 1)
        #print "---1", collect[pos]
        #print "---2", collect[pos+1]
        #print >> sys.stderr, "Here", dist
        if flip:
            collect[pos][4] = collect[pos][4] + random.randint(20, dist)
        else:
            collect[pos+1][3] = collect[pos+1][3] - random.randint(20, dist)

    for e in collect:
        e[8]['transcript_id'] = "synthetic_"+str(i)
        print e





