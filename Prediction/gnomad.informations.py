#python3
#python64/3.5.2
"""
Usage : python3 gnomad.informations.py gnomad_data_base_files
"""
import os
import sys
import csv
import datetime
import pandas as pd
import _pickle as cpickle
from multiprocessing import Pool

path_output='/home/jacques_group/Frederique/LevesqueS/databases/'

def timestamp(name):
	print("{0} : Timestamp: {1:%Y-%m-%d %H:%M:%S}".format(name, datetime.datetime.now()))

def pickle(name,objects) :
	with open(name,'wb') as file:
		pickler = cpickle.Pickler(file)
		pickler.dump(objects)

def reading(file,population):	
	timestamp(file)										
	Snps = []
	with open(file,'r') as f :
			for line in f :	
				if not line.startswith('#'):
					line = line.replace("\n","").replace(";","\t").replace("=","\t").split('\t')
					info = line[6:]
					for idx,ac in enumerate(info[info.index('AC')+1].split(',')) :# split into line snp with multiple ALT allele
						if int(ac) != 0: 
							snp = line[0] +"_"+ line[1] +"_"+ line[3]+'_'+line[4].split(',')[idx]
							AC=[0 if (pop=='SAS') and ('AN_SAS' not in info) else int((info[info.index('AC_'+str(pop))+1]).split(',')[idx]) for pop in population ]
							AN=[0 if (pop=='SAS') and ('AN_SAS' not in info) else int(info[info.index('AN_'+str(pop))+1]) for pop in population ]
							Snps.append(([snp]+AC+AN))
	return Snps

def writing(Snps):
	timestamp("Writing")
	with open('gnomad.informations.raw.csv', 'w', newline='') as csvfile:
			table = csv.writer(csvfile, delimiter='\t',quotechar='|', quoting=csv.QUOTE_MINIMAL)
			for snp in Snps :
				table.writerow(snp)

def concatenate(Snps):
	timestamp("concatenate")
	labels = ['#SNP','AC_AFR','AC_AMR','AC_ASJ','AC_EAS','AC_FIN','AC_NFE','AC_OTH','AC_SAS','AN_AFR','AN_AMR','AN_ASJ','AN_EAS','AN_FIN','AN_NFE','AN_OTH','AN_SAS']
	df = pd.DataFrame.from_records(Snps, columns=labels)
	cat = df.groupby(['#SNP']).sum()
	cat.to_csv(path_output+'gnomad.coding.cat.csv', sep ="\t")
	cat=cat.reset_index()
	cat=cat.values.tolist()
	return cat

def notAjusted(cat,population):
	timestamp("ajusting")
	Informations={}
	for line in cat:
		snp=line[0]
		AC=line[1:9]
		AN=line[9:]
		AF=[(float(AC[p]))/(float(AN[p])) if float(AN[p])!=0 else 0.0 for p in range(len(population))]
		Informations[snp]=[AF,AC,AN]
	pickle(path_output+'gnomad.coding.notAjusted.pickle',Informations)

def ajusted(cat,population):
	timestamp("ajusting")
	Informations={}

	with open(path_output+'gnomad.coding.cat.csv','r') as cat:
		for line in cat:
			if not line.startswith('#'):
				line=line.replace('\n','').split('\t')
				snp=line[0]
				AC=line[1:9]
				AN=line[9:]
				if int(AN[-1]) != 0:	# only common variants btw exomes and genomes
					AF=[(float(AC[p])+1)/(float(AN[p])+4) if float(AN[p])!=0 else 0 for p in range(len(population))]
					Informations[snp]=AF					
	pickle(path_output+'gnomad.coding.ajusted.AF.pickle',Informations)

	"""
	for line in cat:
		snp=line[0]
		AC=line[1:9]
		AN=line[9:]
		if int(AN[-1]) != 0:	# only common variants btw exomes and genomes
			AF=[(float(AC[p])+1)/(float(AN[p])+4) if float(AN[p])!=0 else 0 for p in range(len(population))]
			Informations[snp]=[AF,AC,AN]
	pickle(path_output+'gnomad.coding.ajusted.pickle',Informations)
	"""

def main():
	
	population = ['AFR','AMR','ASJ','EAS','FIN','NFE','OTH','SAS']
	
	files = sys.argv[1]
	All_Snps = []
	with open(files,'r') as f :				
		for file in f :
			file = file.replace("\n", "")
			Snps = reading(file,population)	# reading GNOMAD files
			All_Snps = All_Snps+Snps

	cat = concatenate(All_Snps)		
#	notAjusted(population)
	ajusted(cat,population)

#	writing(All_Snps)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")














"""
def ajusted_tmp(cat,population):
	timestamp("ajusting")
	with open(path_output+'gnomad.informations.cat.csv','r') as f :
		Frequency={}
		SNP=[{} for i in range(len(population))]
		for line in f :
			if not line.startswith('#'):
				line=line.replace("\n","").split('\t')
				snp=line[0]
				AF=[]
				print (line[9:17])
				for idx,pop in enumerate(population) :
					AC = float(line[1+idx])
					AN = float(line[9+idx])
					if not '0' in line[9:17]:
						SNP[idx][snp] = ((float(AC)+1)/(float(AN)+4))
						AF.append(AC/AN)
					else :
						AF.append(0.0)
				Frequency[snp]=(AF)

	with open('frequency.pickle','wb') as file:
		pickler = cpickle.Pickler(file)
		pickler.dump(Frequency)
	
	timestamp("pickling")
	with Pool(8) as P :
		P.starmap(pickle, zip(population,SNP)) 
"""