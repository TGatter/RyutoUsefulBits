#!/usr/bin/env python

import sys, ast, os

left_t = 0
right_t = 0
chrom_t = ""
temp = []


chrom = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

with open(sys.argv[1]) as f:
	for line in f:

		split=line.strip().split()
		if split[2] == "transcript":
			if chrom_t == chrom and start_t <= end and end_t >= start:
				for t in temp:
					print t

			temp = []
			chrom_t = split[0]
			start_t = int(split[3])
			end_t = int(split[4])
	
		temp.append(line.strip())

		
