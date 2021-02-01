#!/usr/bin/env python

import sys, ast, os


bam_name="split" 

# 2 Bam1
# 3 Bam2
# 4 Bam3
# 5 Bam4

# 6 target_folder
# 7 truth

with open(sys.argv[1]) as f:
	for line in f:

		split=line.strip().split()

		chrom = split[0]
		start = int(split[1])
		end = int(split[2])

		if end - start < 250:
			continue 

		os.system('~/Downloads/samtools-1.5/samtools view -b '+ sys.argv[2] +' "'+chrom+':'+str(start)+'-'+str(end)+'" > '+ sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align1.bam")
		os.system('~/Downloads/samtools-1.5/samtools index  ' + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align1.bam " + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align1.bai")			

		os.system('~/Downloads/samtools-1.5/samtools view -b '+ sys.argv[3] +' "'+chrom+':'+str(start)+'-'+str(end)+'" > '+ sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align2.bam")
		os.system('~/Downloads/samtools-1.5/samtools index  ' + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align2.bam " + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align2.bai")			

		os.system('~/Downloads/samtools-1.5/samtools view -b '+ sys.argv[4] +' "'+chrom+':'+str(start)+'-'+str(end)+'" > '+ sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align3.bam")
		os.system('~/Downloads/samtools-1.5/samtools index  ' + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align3.bam " + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align3.bai")	

		os.system('~/Downloads/samtools-1.5/samtools view -b '+ sys.argv[5] +' "'+chrom+':'+str(start)+'-'+str(end)+'" > '+ sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align4.bam")
		os.system('~/Downloads/samtools-1.5/samtools index  ' + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align4.bam " + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".align4.bai")	

		os.system("./search_gtf.py "+ sys.argv[7] + " " + chrom + " " + str(start) + " " + str(end) + " > " + sys.argv[6]+"/"+bam_name+"."+chrom+"."+str(start)+"."+str(end)+".truth.gtf")
		
