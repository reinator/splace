# SPLACE
A tool to SPLit, Align and ConcatenatE genes sequences for phylogenetic inference

## How to run SPLACE

SPLACE will:

1.Split the same genes from different organisms;\
2.Align each gene group separately;\
3.Concatenate the aligned genes originated from the same organism.


Usage:

~~~
./splace.sh -l <Fasta files list>.txt -t <num_threads> -o <output_name> -g <genes_list>
~~~

#-l <Fasta files list> = Text file with a list of all fasta files to be included in the result.\
#-t <num_threads> = Number or threads to use in the blast step. Default is 1.\
#-o <output_name> = Name of the result file to be generated.\
#-g <genes_list> = If specified, any organism that do not have such gene, will have "?" in the supermatrix. If no specified, analysis is carried out only with the shared genes.


~~~
./splace.sh -l filesList.txt -t 24 -o Glomeridesmus - g genes_list.txt
~~~~
  
Example of a Fasta_files_list.txt:

~~~
bases_genes_ITV1046I2.fasta 	
bases_genes_ITV8918.fasta 	
~~~ 

Example of a genes_list.txt:
~~~
nad1
cox1
nad2
cox2
atp6
~~~

### Requisites:

 The gene names must be the same in the fasta file of the different organisms, being necessary the use of the "_" separator 
between the gene name and any other gene specification or none separator and specification at all.

Ex: File bases_genes_ITV1046I2.fasta:

~~~
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
~~~

Ex: File bases_genes_ITV8918.fasta

~~~
>atp8
attccacaaatatccccaataaattgagaaataatatttttaacttctattttaattctt
ttaataacttcaattattatccatcaaaattctaattttcgattatctaaaaaacaaaaa
attccttcaaaaatttattaa
>atp6
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
~~~

Results will be saved at the created folder <output_name>_results.
In the folder <output_name>_genes will be created fasta file for each gene and its alignment.
  
## Workflow:
![](https://github.com/reinator/splace/blob/main/workflow.tif?raw=true)
  
A - To generate the supermatrix of n organisms, SPLACE will need n fasta files, each one containing all the g genes from a particular organism;\
B - First, SPLACE splits the genes from an organism, gathering the genes that have the same name from the n organisms into a single fasta file, therefore generating g new fasta files, each one containing the same gene from different organisms (Figure 1B);\
C - Then, SPLACE aligns each one of the g fasta files using the MAFFT aligner (with default parameter –auto), generating g’ new aligned fasta files;\
D - Finally, the genes in the g’ aligned fasta files that came from the same organism are concatenated into a single sequence, generating a single fasta file with the supermatrix containing n sequences, each representing one of the n organisms;\
E - Phylogeny can then be reconstructed using the supermatrix fasta file and the method of choice by the user.
