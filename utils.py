
import sys
import os
import argparse
import gzip
from Bio.Seq import Seq
from Bio import SeqIO
import logging
import Colorer
import itertools
import string
import pandas as pd
from Levenshtein import distance
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns


def collision_boxplot(input,label,cutoff=10,threshold=0.8):

	df = pd.read_csv(input,sep=" ",header=None)
	fig, ax =plt.subplots(1,2,figsize=(10,5))
	sns.boxplot(y=df[3], ax=ax[0])
	ax[0].set_yscale('log')
	ax[0].set_title(f'total cells: {df.shape[0]}')
	ax[0].set_ylabel('#reads per cell')
	df = df[df[3]>=cutoff]
	sns.boxplot(y=df[3], ax=ax[1])
	ax[1].set_yscale('log')
	ax[1].set_title(f'total cells (after filter): {df.shape[0]}')
	ax[1].set_ylabel('#reads per cell')
	df['max_mapping_rate'] =  df[[4,5]].max(axis=1)
	df['collision'] = df['max_mapping_rate']<threshold
	fig.suptitle(f"collision rate for sample {label} is {df['collision'].sum()/df.shape[0]:.2f}", fontsize=16)
	plt.savefig(f"{label}.collision.boxplot.png",bbox_inches='tight')




def revcomp(seq):
	try: ## python2
		tab = string.maketrans(b"ACTG", b"TGAC")
	except:  ## python3
		tab = bytes.maketrans(b"ACTG", b"TGAC")
	return seq.translate(tab)[::-1]

def read_fasta(f,rc=False):
	my_dict = {}
	for r in SeqIO.parse(f, "fasta"):
		seq = str(r.seq).upper()
		if rc:
			seq = revcomp(seq)
		my_dict[r.id] = seq
	return my_dict	
	
	
def read_barcode_table(f,rc=False):
	df = pd.read_csv(f,sep="\t",header=None)
	# print (df)
	df.index = df[0].to_list()
	if rc:
		df[2] = [revcomp(x) for x in df[1]]
		return df[2].to_dict()
	return df[1].to_dict()
	

def mismatch_sequence(seq,n):
	bases=['A','T','G','C']
	k_mer = [''.join(p) for p in itertools.product(bases, repeat=len(seq))]
	result = []
	for k in k_mer:
		if hamming_distance(seq,k)<=n:
			result.append(k)
	return result

def hamming_distance(s1, s2):
	return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))

def setup_custom_logger(jid):
	logFormatter = logging.Formatter('%(asctime)s - %(levelname)s - %(funcName)s - %(message)s',"%Y-%m-%d %H:%M:%S")
	rootLogger = logging.getLogger("root")
	fileHandler = logging.FileHandler(jid+".log")
	fileHandler.setFormatter(logFormatter)
	rootLogger.addHandler(fileHandler)
	consoleHandler = logging.StreamHandler()
	consoleHandler.setFormatter(logFormatter)
	rootLogger.addHandler(consoleHandler)
	rootLogger.setLevel(logging.INFO)
	return rootLogger


def dict3d_to_df(x):
	bc3 = []
	bc2 = []
	bc1 = []
	counts = []
	for i in x:
		for j in x[i]:
			for k in x[i][j]:
				bc1.append(i)
				bc2.append(j)
				bc3.append(k)
				counts.append(x[i][j][k])
	df = pd.DataFrame()
	df['barcode_1'] = bc1
	df['barcode_2'] = bc2
	df['barcode_3'] = bc3
	df['total_reads'] = counts
	df = df.sort_values("total_reads",ascending=False)
	return df
	

def k_mer_distance(k,error,barcode_list=None):
	out_dict={}
	bases=['A','T','G','C']
	if error != 0:
		k_mer = [''.join(p) for p in itertools.product(bases, repeat=k)]
	else:
		k_mer = barcode_list
	for ii in range(len(k_mer)):
		i = k_mer[ii]
		out_dict[i]={}
		for jj in range(ii,len(k_mer)):
			j = k_mer[jj]
			dist = distance(i, j)
			# dist = 0
			# if i != j:
				# dist = 1
			if dist <= error:
				out_dict[i][j] = dist
				if not j in out_dict:
					out_dict[j]={}
				out_dict[j][i] = dist
	return out_dict

def get_barcode_dict(bc1):
	out = {}
	for i in bc1:
		out[i]={}
		for j in bc1:
			out[i][j]={}
			for k in bc1:
				out[i][j][k]=0
	return out

def find_dist(kmer_dict,target_barcode,barcode_list):
	for x in barcode_list:
		try:
			kmer_dict[x][target_barcode] 
			return True,x
		except:
			pass
	return False,None
	
