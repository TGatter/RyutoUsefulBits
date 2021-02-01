#!/usr/bin/env python

import sys, ast, os


tmp = []
all_trans = []

for line in open(sys.argv[1]):
	
	split=line.strip().split()

	if split[2] == "transcript":
		s = int(split[3])
		e = int(split[4])
		all_trans.append(( split[0], min(s,e), max(s,e), []) )
	
	all_trans[-1][3].append(line)

all_trans.sort(key=lambda e: (e[0], e[1], e[2]))


for at in all_trans:
	for l in at[3]:
		print l, 
