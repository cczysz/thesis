#!/bin/python3
from getSNAPResults import snapresults

# snapresults(list of rsids,population of interest, r2 cutoff)

# snp_list = []

out_file = open('/home/t.cczysz/broad_results.txt','w+')

with open('/home/t.cczysz/rsids.txt') as f:
	snp_list = f.readlines()

sub_lists = [snp_list[i:i+100] for i in range(0,len(snp_list),100)]

for sub_list in sub_lists:
	#sub_list = [snp.strip() for snp in sub_list]
		# snapresults(snp,'CEU','0.8')
	# snp_sub_list = snp_list[sub_list]

	snapresults(sub_list,'CEU','0.8')

	with open('SNAPResults_0.8.txt') as f:
		results = f.readlines()

	for line in results:

		if "SNP" in line:
			continue
		if "Error" in line:
			print("error")
			break
		if "WARNING" in line:
			continue
		# line = line.strip().split('\t')
		# rsid_1 , rsid_2 = line[0] , line[1]
		out_file.write(line)
out_file.close()

