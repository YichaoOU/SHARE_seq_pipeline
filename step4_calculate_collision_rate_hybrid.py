#!/usr/bin/env python
import subprocess
import sys
# from Levenshtein import distance
import gzip
# from multiprocessing import Pool
from functools import partial
import os
import argparse
import pandas as pd
import subprocess
import glob
from copy import deepcopy as dp
'''

In the working dir, run python step3.py {{jid}}

python ../../../../src/step4_calculate_collision_rate.py --table KOK676_S3.total_number_reads.tsv --human human/KOK676_S3.R1.bed --mouse mouse/KOK676_S3.R1.bed

'''





def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	mainParser.add_argument("--table",  help="barcode count table (input)", required=True)
	mainParser.add_argument("--reads",  help="human_mouse bamtobed r1 (input)", required=True)
	mainParser.add_argument("--threshold",  help="threshold for collision", default=0.8,type=float)

	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def get_file(x):
	return glob.glob(x)[0]

def unique_dict(d):
	for k in d:
		d[k] = len(set(d[k]))
	return d
def parse_file(x,total_human,total_mouse,total_reads):
	h=0
	m=0
	t=0
	
	with open(x) as f:
		for line in f:
			t+=1
			chr = line.split()[0]
			temp = line.split()[3].split("_")
			label = temp[1].replace("/1","")
			read_name = temp[0]
			# print (read_name)
			if "hg" in chr:
				h+=1
				total_human[label].append(read_name)
			if "mm" in chr:
				m+=1
				total_mouse[label].append(read_name)
			total_reads[label].append(read_name)
	# print (h,m,t)
	return unique_dict(total_human),unique_dict(total_mouse),unique_dict(total_reads)



def main():

	args = my_args()
	
	barcode = pd.read_csv(args.table,sep="\t")
	barcode = barcode[barcode['total_reads']>0]
	sample = args.table.split(".")[0]
	barcode['name'] = barcode['barcode_1']+barcode['barcode_2']+barcode['barcode_3']
	barcode.index = barcode['name']
	barcode['init'] = [[] for _ in range(len(barcode))]
	
	total_human = dp(barcode['init'].to_dict())
	total_mouse = dp(barcode['init'].to_dict())
	total_reads = dp(barcode['init'].to_dict())

	total_human,total_mouse,total_reads = parse_file(args.reads,total_human,total_mouse,total_reads)

	df1 = pd.DataFrame.from_dict(total_human,orient='index')
	# print (df1.sort_values(0).tail())
	df2 = pd.DataFrame.from_dict(total_mouse,orient='index')
	df3 = pd.DataFrame.from_dict(total_reads,orient='index')
	barcode['human_count']  = df1[0]
	barcode['mouse_count']  = df2[0]
	barcode['total_reads'] = df3[0]
	barcode['human_percent'] =  barcode['human_count'] / barcode['total_reads']
	barcode['mouse_percent'] =  barcode['mouse_count'] / barcode['total_reads']
	barcode['max_mapping_rate'] =  barcode[['human_percent','mouse_percent']].max(axis=1)
	barcode['collision'] = barcode['max_mapping_rate']<args.threshold
	# print (barcode.head())
	print ("collision rate for sample %s is %s "%(sample,float(barcode['collision'].sum())/barcode.shape[0]))
	barcode.to_csv(f'{sample}.mapping_rate.tsv',sep="\t")
	barcode[['human_count','mouse_count','total_reads','human_percent','mouse_percent']].to_csv(f"{sample}.for_collision_plot.tsv",sep=" ",index=True,header=False)
	

	
if __name__ == "__main__":
	main()
	
		