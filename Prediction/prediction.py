# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import argparse
import datetime
from math import log,exp
import _pickle as cpickle

def timestamp(name):
	print("{0} : Timestamp: {1:%Y-%m-%d %H:%M:%S}".format(name, datetime.datetime.now()))
def argsparse():
	parser=argparse.ArgumentParser(description='Group ethnic prediction with vcf file from individual.')
	parser.add_argument('-d','--pickle', action='store', dest='pickle', metavar='pickle', help='dictionary of GnomAD frequency',required=True)
	parser.add_argument('-i','--input', action='store', dest='input', metavar='file', help='vcf file',required=True)
	parser.add_argument('-o','--output', action='store', dest='output', metavar='file', help='output file, default=output.out',default='output.out')
	parser.add_argument('-m','--multiple', action='store_true',dest='multiple', help='multiple input file, default=False')
	parser.add_argument('-p','--prior', action='store_true',dest='prior', help='use priors in the determination of ethnicity, default=False)')
	parser.add_argument('-oth','--other', action='store_true',dest='oth', help='add the oth population to analysis, default=False')
	args=parser.parse_args()
	return args

def unpickle (name):
	timestamp('unpickling...')
	with open(name,'rb') as file:
		pickler=cpickle.Unpickler(file)
		objects=pickler.load()
	return objects

def pickle(name,objects):
	timestamp('pickling...')
	with open(name,'wb') as file:
		pickler = cpickle.Pickler(file)
		pickler.dump(objects)

def reading(file): # vcf file reading
	Snps=[]
	with open(file,'r') as f:
		for line in f:
			if not line.startswith('##'):
				if line.startswith('#'):
					name=line.replace("\n", "").split('\t')[9]
				else:
					line=line.replace("\n", "").split('\t')
					for a,alt in enumerate(line[4].split(',')):
						snp=line[0].replace('chr','')+'_'+line[1]+'_'+line[3]+'_'+alt
						Snps.append(snp)
	return Snps,name

def analysis(population,dictionary,Snps,prior,use_prior,use_oth) :
	prob=[[]for i in population]
	for snp in Snps:
		if snp in dictionary : 
			x=[prob[p].append(log(dictionary[snp][p])) for p in range(len(population))] 

	if use_prior:
		total=[(sum(prob[p])+log(prior[p])) for p in range(len(population))]
	else:
		total=[(sum(prob[p])) for p in range(len(population))]
	Total=[(exp(total[p]-max(total))/sum([exp(total[i]-max(total)) for i in range(len(population))])) for p in range(len(population))]

	if use_oth:
		ethnicity=population[total.index(max(total))]
	else:
		ethnicity=population[total.index(max(total))] if population[total.index(max(total))] != 'OTH' else population[total.index((sorted(total))[-2])]

	score=round(Total[population.index(ethnicity)], 4)
	return ethnicity,score

def main():
	args=argsparse()
	population=['AFR','AMR','ASJ','EAS','FIN','NFE','OTH','SAS']
	prior=[0.0001,0.0001,0.0001,0.0001,0.0001,1.0000,0.0001,0.0001]
	dictionary=unpickle(args.pickle)

	o=open(args.output,'w')
	o.write('#ID'+'\t'+'PREDICTION'+'\t'+'SCORE')

	if args.multiple :
		files=[line.replace('\n','') for line in open(args.input,'r') ]
		for file in files:
			Snps,ID=reading(file)
			ethnicity,score=analysis(population,dictionary,Snps,prior,args.prior,args.oth)
			o.write('\n'+ID+'\t'+ethnicity+'\t'+str(score))

	else:
		Snps,ID=reading(args.input)
		ethnicity,score=analysis(population,dictionary,Snps,prior,args.prior,args.oth)
		o.write('\n'+ID+'\t'+ethnicity+'\t'+str(score))

	o.close()

if "__main__" == "__main__":
	timestamp('STARTING')
	main()
	timestamp('DONE')