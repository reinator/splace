'''
Author: Renato Oliveira
Date: 07-03-2018
Version: 1.2


1.Split the same genes from different organisms;
2.Align each gene group separately;
3.Concatenate the aligned genes originated from the same organism.


Usage: ./splace.sh -l <Fasta files list>.txt -t <num_threads> -o <output_name> -g <genes_list>

#-l <Fasta files list> = Text file with a list of all fasta files to be included in the result.
#-t <num_threads> = Number or threads to use in the blast step. Default is 1.
#-o <output_name> = Name of the result file to be generated.
#-g <genes_list> = If specified, any organism that do not have such gene, will have "?" in the supermatrix. If no specified, analysis is carried out only with the shared genes.

Ex: ./splace.sh -l filesList.txt -t 24 -o Glomeridesmus -g genes_list.txt

Ex: Fasta_files_list.txt

--------------------------------
								|
bases_genes_ITV1046I2.fasta 	|
bases_genes_ITV8918.fasta 		|
								|
--------------------------------

Ex: genes_list.txt

--------------------------------
								|
nad1							|
cox1							|
nad2							|
cox3							|
atp6 							|
								|
--------------------------------

Requisites: The gene names must be the same in the fasta file of the different organisms, being necessary the use of the "_" separator 
between the gene name and any other gene specification or none specification at all.

Ex: File bases_genes_ITV1046I2.fasta:

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


Ex: File bases_genes_ITV8918.fasta

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


Results will be saved at the created folder <output_name>_results.
In the folder <output_name>_genes will be created fasta file for each gene and its alignment.
'''

THREADS=1
OUTPUT="output"
GENES_LIST="no"
#DB_FILE=/bio/pimba_metabarcoding/databases.txt

while getopts "l:o:t:g:" opt; do
	case $opt in
		l) FILELIST="$OPTARG"
		;;
		o) OUTPUT="$OPTARG"
		;;
		t) THREADS="$OPTARG"
        ;;
        g) GENES_LIST="$OPTARG"
        ;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

CURRENT_PATH=$(pwd)

DIR_NAME_FILE=$(dirname $FILELIST)
cd $DIR_NAME_FILE
FULL_PATH_FILE=$(pwd)
cd $CURRENT_PATH

COMMON_PATH=$({ echo $FULL_PATH_FILE; echo $CURRENT_PATH;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')
FILELIST=$(echo ${FULL_PATH_FILE#"$COMMON_PATH"})/$(basename $FILELIST)

echo "Creating a SPLACE Container: "
docker run -id -v $COMMON_PATH:/common/ -v $CURRENT_PATH:/output/ --name splace itvdsbioinfo/splace:latest

if [ $GENES_LIST = "no" ];
then

	echo "Running the SPLACE Container: "
	docker exec -i splace /bin/bash -c 'cd /output/; \
		python3 /splace/SPLACE_v2.py /common/'$FILELIST' '$THREADS' '$OUTPUT'; \
		chmod -R 777 /output/'$OUTPUT'_*;'

else
	echo "Running the SPLACE Container: "
	docker exec -i splace /bin/bash -c 'cd /output/; \
		python3 /splace/SPLACE_v2.py /common/'$FILELIST' '$THREADS' '$OUTPUT' '$GENES_LIST'; \
		chmod -R 777 /output/'$OUTPUT'_*;'
fi

#conda deactivate

echo "Stopping Containeres: "
docker stop splace

echo "Removing Containeres: "
docker rm splace