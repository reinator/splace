'''
Author: Renato Oliveira
Date: 05-04-2022
Version: 1.2


1.Separa mesmos genes de diferentes organismos
2.Alinha separadamente cada gene
3.Concatena alinhamentos desde que os ids das sequencias tenham o mesmo nome


Usage: python /bio_temp/share_bio/softwares/SPLACE/SPLACE_v1.2.py <Fasta files list>.txt <num_threads> <output_name>
Ex: python /bio_temp/share_bio/softwares/SPLACE/SPLACE_v1.2.py filesList.txt 24 Glomeridesmus

Ex: Fasta_files_list.txt

--------------------------------
								|
bases_genes_ITV1046I2.fasta 	|
bases_genes_ITV8918.fasta 		|
								|
--------------------------------


Exigencias: O nome dos genes devem ser os mesmos nos diferentes arquivos fasta dos organismos, sendo necessario o uso do separador "_" 
entre o nome do gene e qualquer outra especificacao do mesmo.

Ex: Arquivo bases_genes_ITV1046I2.fasta:

>atp8_ITV1046I2 atp8 ATP synthase F0 subunit 8 7816:7956 forward
attccacaaatatctccaataaattgagaaataatattcttaacttctattttaattctt
ttaataatttcaattattattcatcaaaactcaaattttaaattatctaaaagaaaaaaa
attccttcaaaaatttattaa
>atp6_ITV1046I2 atp6 ATP synthase F0 subunit 6 7964:8638 forward
atgataacaaatttattctcaatttttgatccttcttctattaccccaatttcactaaat
tgattaagtataattttaattataatttttataaatttaactttctggttattcaagtca
aaaaatcaaattattattaataatctaaattctatcattattaaggaattaaaaacaaca
ttaaaaacaagaaactatcctaactttattttaattttactaactttatttattttaatt
ttaattaataatttaataggtttattcccctatattttcacaagtacaagccacataata
attactttatcattagctttacctctatgatttatatccattataatattaataacaata
aacacaataaattttttagctcatctagttcctcaaagaactccttcttacttaatatcc
tttatagttttaattgaatctattagaaatattattcgaccaataacattagcaattcga
ctaacggctaatataattgctggacatcttttaatcactcttttaagatcaataaatgaa
aaaataaatattttttcatcaattattattattttatcttcaacaactctcataatctta
gaattagctgtagcaattattcaagcctacgtatttataactctattatcactatattta
agagaaattaattaa


Ex: Arquivo bases_genes_ITV8918.fasta

>atp8_ITV8918 atp8 atp8 7851:7991 forward
attccacaaatatccccaataaattgagaaataatatttttaacttctattttaattctt
ttaataacttcaattattatccatcaaaattctaattttcgattatctaaaaaacaaaaa
attccttcaaaaatttattaa
>atp6_ITV8918 atp6 atp6 7999:8673 forward
atgataacaaatttattctcaatttttgatccttcttctattacctcaatctcattaaat
tgattaagtatacttttaattataatctttataaatttaactttttgattatttaaatca
aaaaatcaaattcttattaataacctaaattctattattaataaagaattaaaaacaaca
ttaaaaacaagaaactaccctaattttattttaattttattaactttattcactctaatt
ttaattaataacttgataggcttatttccctacatctttacaagtacaagtcacataata
attactttatcattagctttacccctatgacttatatctattataatattaataataata
aatacaataaattttctagctcatctagtccctcaaagaactccttcatacttaatatcc
tttatagttttaattgaatctattagaaatattattcgaccaataacattagcaattcga
ttaactgctaatataattgccggacatcttttaattacccttttaagatcaataaatgaa
aaaataaatattttctcatcaattattattattttatcctcaacaatccttataattcta
gaattagctgtagcaatcattcaagcctacgtatttataactctactatcattatattta
agagaaatcaattaa


Os resultados serao salvos na pasta <output_name>_results criada.
Na pasta <output_name>_genes serao criados arquivos fasta para cada gene encontrado, assim como o resultado dos alinhamentos

'''

import sys
import os


def readFilesList(arqList):
	filesList = []
	print("Organisms gene files:\n")
	for line in arqList:
		line=line.split("\n")
		line=line[0]
		print(line)
		filesList.append(line)

	return (filesList)

def readGenesList(genes_list_file):
	genes_list = []
	print("Genes list:\n")
	for line in genes_list_file:
		line=line.split("\n")
		line=line[0]
		print(line)
		genes_list.append(line)

	return (genes_list)


def run_shared_genes(output_name, genes_keys, gene_index_list, lines):
	gene_arq_name=[]
	gene_order = dict()

	sharedFile = open(output_name+"_sharedGenes.txt","w")
	sharedTxt=""


	for gk in genes_keys:
		gene_order[gene_index_list[0][gk]]=gk[1:]

	gene_name_index = dict()

	for g in genes_keys:
		#print(g)
		sharedTxt+=g+"\n"
		gene=""
		i = gene_index_list[0][g]
		line = lines[0][i]
		gene+=line
		i+=1
		line = lines[0][i]
		while(line[0]!=">"):
			gene+=line
			i+=1
			if(i<len(lines[0])):
				line = lines[0][i]
			else:
				break

		for l in range(1,len(lines)):
			#print(filesList[l])
			i = gene_index_list[l][g]
			line = lines[l][i]
			gene+=line
			i+=1
			line = lines[l][i]
			while(line[0]!=">"):
				gene+=line
				i+=1
				if(i<len(lines[l])):
					line = lines[l][i]
				else:
					break

		gene_arq_name.append(output_name+"_genes/"+g[1:]+".fasta")

		gene_name_index[g[1:]] = gene_arq_name[-1]
		gene_arq = open(gene_arq_name[-1],"w")
		gene_arq.write(gene)
		gene_arq.close()

	sharedFile.write(sharedTxt)
	sharedFile.close()

	return (gene_arq_name, gene_name_index, gene_order)

def run_list_genes(output_name, genes_keys, gene_index_list, genes_list, lines):
	gene_arq_name=[]
	gene_order = dict()

	listedFile = open(output_name+"_listedGenes.txt","w")
	listedTxt=""

	i=0
	for gk in genes_list:
		gene_order[i]=gk.strip("\r").strip("\n")
		i+=1

	gene_name_index = dict()

	for g in genes_list:
		#print(g)
		listedTxt+=g+"\n"
		gene=""
		if(">"+g in gene_index_list[0].keys()):
			i = gene_index_list[0][">"+g]
			line = lines[0][i]
			gene+=line
			i+=1
			line = lines[0][i]
			while(line[0]!=">"):
				gene+=line
				i+=1
				if(i<len(lines[0])):
					line = lines[0][i]
				else:
					break
		else:
			gene+=">"+g+"_MISSING\n"
			gene+="??????????\n"

		for l in range(1,len(lines)):
			#print(filesList[l])
			if(">"+g in gene_index_list[l].keys()):
				i = gene_index_list[l][">"+g]
				line = lines[l][i]
				gene+=line
				i+=1
				line = lines[l][i]
				while(line[0]!=">"):
					gene+=line
					i+=1
					if(i<len(lines[l])):
						line = lines[l][i]
					else:
						break
			else:
				gene+=">"+g+"_MISSING\n"
				gene+="??????????\n"

		gene_arq_name.append(output_name+"_genes/"+g+".fasta")

		gene_name_index[g] = gene_arq_name[-1]
		gene_arq = open(gene_arq_name[-1],"w")
		gene_arq.write(gene)
		gene_arq.close()

	listedFile.write(listedTxt)
	listedFile.close()

	return (gene_arq_name, gene_name_index, gene_order)



def splitGenes(filesList, output_name, genes_list):
	os.system("mkdir -p "+output_name+"_genes")
	os.system("cd "+output_name+"_genes")

	genes_and_org = open(output_name+"_genes_and_orgs.tsv", "w")
	genes_and_org_text = ""

	num_org = len(filesList)

	gene_index_list=[]
	
	gene=""
	gene_name=""
	lines=[] #Linhas dos arquivos de genes de cada organismo
	for n in range(num_org):
		arq_genes = open(filesList[n],"r")
		#print(filesList[n])
		genes_and_org_text+=filesList[n]+"\t"
		lines.append(arq_genes.readlines()) #linhas de cada arquivo de genes guardado em uma lista
		num_lines=len(lines[n])
		i=0
		gene_index = dict()
		while(i < num_lines):
			line=lines[n][i]
			if(line[0]==">"):
				if(line.find("_") != -1):
					gene_name=line[:line.find("_")].strip("\r").strip("\n")
				elif(line.find(" ") != -1):
					gene_name=line[:line.find(" ")].strip("\r").strip("\n")
				elif(line.find("\t") != -1):
					gene_name=line[:line.find("\t")].strip("\r").strip("\n")
				elif(line.find("-") != -1):
					gene_name=line[:line.find("-")].strip("\r").strip("\n")
				else:
					gene_name = line.strip("\r").strip("\n")
					#print(gene_name)
				gene_index[gene_name]=i
				#print(gene_name)
				
			i+=1
		gene_index_list.append(gene_index)

		#write the genes of an organism sorted alphabetically
		genes_of_organism = gene_index.keys()
		for gene_name in sorted(genes_of_organism):
			genes_and_org_text+=gene_name[1:]+"\t"

		genes_and_org_text=genes_and_org_text.strip("\t")
		genes_and_org_text+="\n"

	genes_and_org.write(genes_and_org_text)
	genes_and_org.close()

	print(str(len(gene_index_list))+" organisms")
	genes_keys = sorted(gene_index_list[0].keys())
	genes_keys_union = sorted(gene_index_list[0].keys())
	#shared_genes = gene_index_list[0].keys()
	for n in range(1, num_org):
		genes_names = gene_index_list[n].keys()
		genes_keys = list(set(genes_names) & set(genes_keys))
		genes_keys_union = list(set(genes_names) | set(genes_keys_union))
		#print(sorted(genes_keys_union))
		#print(n)

	num_genes = len(genes_keys)
	print(str(num_genes)+" shared genes")
	print(str(len(genes_keys_union))+" total genes")


	if(genes_list == None):
		gene_arq_name, gene_name_index, gene_order = run_shared_genes(output_name, genes_keys, gene_index_list, lines)
	else:
		gene_arq_name, gene_name_index, gene_order = run_list_genes(output_name, genes_keys_union, gene_index_list, genes_list, lines)
	
	return (gene_arq_name, gene_name_index, gene_order, num_org)


	
def concatenateGenes(output_name,gene_name_index,gene_order, num_org, filesList):
	os.system("mkdir -p "+output_name+"_result")
	conc_gene_arq = open(output_name+"_result/concat_aligned_"+output_name+".fasta", "w")

	conc_gene = [] #lista onde cada item refere a concatenacao dos genes de um organismo

	for n in range(num_org):
		conc_gene.append("")


	gene_order_keys = gene_order.keys()
	sorted(gene_order_keys)
	print("Shared genes among organisms")
	for o in gene_order_keys:
		gene_name = gene_order[o]
		print(gene_name)
		gene_name_dir = gene_name_index[gene_name]
		ponto = gene_name_dir.rfind(".")
		gene_name_dir = gene_name_dir[:ponto]+"_mafft.fasta"

		gene_name_arq = open(gene_name_dir, "r")
		#print(gene_name_dir)
		gene_name_text = gene_name_arq.readlines()
		#print(gene_name_text)
		
		i=-1
		j=0
		
		while(j<len(gene_name_text)):
			is_missing = False
			line = gene_name_text[j]
			if("MISSING" in line):
				is_missing = True
			if(line[0]==">"):
				i+=1
				#conc_gene[i]+=line
				j+=1
				line = gene_name_text[j]
				if(is_missing):
					line = line.replace("-","?")
				
			while(line[0]!=">"):
				conc_gene[i]+=line
				j+=1
				if(j<len(gene_name_text)):
					line = gene_name_text[j]
					if(is_missing):
						line = line.replace("-","?")
				else:
					break
		#print(conc_gene)
		#break

	#print(conc_gene)
	norg=0

	for string in conc_gene:
		string = string.split("\n")
		final_string=""
		for s in string:
			final_string+=s

		#print(final_string)
		name_org=filesList[norg]
		ponto=name_org.rfind(".")
		name_org=name_org[:ponto]

		conc_gene_arq.write(">"+name_org+"\n")
		i=0
		for c in final_string:
			conc_gene_arq.write(c)
			#i+=1
			#if(i==59):
				#conc_gene_arq.write("\n")
				#i=0

		conc_gene_arq.write("\n")
		norg+=1

	conc_gene_arq.close()

	conc_gene_arq = open(output_name+"_result/concat_aligned_"+output_name+".fasta", "r")

	readFasta(conc_gene_arq, output_name)

	conc_gene_arq.close()

def readFasta(arq, output_name):
	novoFasta=""
	novoFasta+=arq.readline()

	for line in arq:
		
		if(line[0]!=">"):
			begin=0
			while(line[begin]=="-"):
				novoFasta+="?"
				begin+=1

			end=len(line)-2
			while(line[end]=="-"):
				end-=1
			end+=1

			for i in range(begin, end):
				novoFasta+=line[i]

			i = end

			while(i<len(line)-2):
				novoFasta+="?"
				i+=1

		else:
			novoFasta+="\n"+line
				

	novoArq = open(output_name+"_result/final_concat_aligned_"+output_name+".fasta", "w")
	novoArq.write(novoFasta)

	novoArq.close()

	
def runMafft(file, num_threads):
	ponto = file.rfind(".")
	command="mafft --thread "+str(num_threads)+" --auto "+file+" > "+file[:ponto]+"_mafft.fasta"
	os.system(command)


arqList = open(sys.argv[1], "r")
num_threads = sys.argv[2]
output_name = sys.argv[3]

try:
    sys.argv[4]
except IndexError:
    genes_list_file = None
else:
    genes_list_file = open(sys.argv[4], "r")

genes_list = readGenesList(genes_list_file)
filesList = readFilesList(arqList)

gene_arq_name,gene_name_index,gene_order,num_org = splitGenes(filesList, output_name, genes_list)

print("Running MAFFT\n")
for i in gene_arq_name:
	runMafft(i, num_threads)

concatenateGenes(output_name, gene_name_index, gene_order, num_org, filesList)

arqList.close()
