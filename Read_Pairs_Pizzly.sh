#!/bin/bash
###MAINTAINER: SARA LATOUR
###CONTACT: saralatour@outlook.com
###UPDATED:21/09/17


#Customization:

#Path to Pizzly 
pizzlypath=/home/saral/Documents/Fusioncallers/pizzly/pizzly

#Path to Kallisto
kallistopath=/home/saral/Documents/Fusioncallers/pizzly/kallisto_linux-v0.43.1

#Path to Reference cDNA Fasta File
cdna=/home/saral/Documents/Fusioncallers/reference_genome/Homo_sapiens.GRCh38.cdna.all.fa.gz  

#Path to Reference GTF File
GTF=/home/saral/Documents/Fusioncallers/reference_genome/Homo_sapiens.GRCh38.90.gtf.gz 

####NOTHING BELOW THIS LINE REQUIRES CHANGING####

cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF
#Maintainer:
echo "\nPizzly Fusion Caller for RNA-seq Read Pairs:"
echo "\nMAINTAINER:\nSara Latour"
echo "\nCONTACT:\nsaralatour@outlook.com"
#Date:
now="$(date +'%d_%m_%Y')"
echo "\nDATE(dd_mm_yyyy):\n$now\n"
cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF
#CMD Line Arguments:
inputdir=$1

if [ $# -ne 1 ]
  then
    echo "\nERROR MESSAGE:\nPlease ensure that the full path to the Fastq Directory is provided:\n"
	exit 1
fi

echo "\nChecking if Input directory is valid..."
if   [ -d "${inputdir}" ]
then echo "✓"
else echo "Input is not a valid directory... exiting now\n";
     exit 1
fi

#Input:
echo "The Directory Selected is:$inputdir\n"

#Create Output File
echo "Creating Output Directory (Output Files Location:$inputdir/Output)\n"
rm -Rf $inputdir/Output
mkdir $inputdir/Output
outputdir=$inputdir/Output
mkdir -p $inputdir/Output/Fastas
mkdir -p $inputdir/Output/Json

#Create File Lists
ls $inputdir > $outputdir/filenames.txt 
cd $outputdir 

cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF

#RUN Index
echo "\nBeginning Index Step for Kallisto...\n"
$kallistopath/kallisto index -i $outputdir/index.idx -k 31 $cdna 
echo "Completed Index Step!\n"

#RUN Kallisto Quant & Pizzly
cat filenames.txt | while read line; do sed 's/R[1,2]_001.fastq.gz//g' ; done > $outputdir/kallisto_list.txt
sort -u $outputdir/kallisto_list.txt > $outputdir/kallisto_unique_list.txt

echo "Running Pizzly Fusion Caller on the following RNA-seq Pairs IDs...\n"
cat $outputdir/kallisto_unique_list.txt
echo "\nStarting Kallisto and Pizzly runs for files from $inputdir directory now...\n"

cat $outputdir/kallisto_unique_list.txt | while read line; do echo "$line"; $kallistopath/kallisto quant -i $outputdir/index.idx --fusion -o $outputdir/"$line"_pizzly_out $inputdir/"$line"R1_001.fastq.gz $inputdir/"$line"R2_001.fastq.gz ; $pizzlypath/build/pizzly -k 31 --gtf $GTF --cache index.cache.txt --align-score 2 --insert-size 400 --fasta $cdna --output "$line"_pizzly_out $outputdir/"$line"_pizzly_out/fusion.txt ; $pizzlypath/scripts/flatten_json.py "$line"_pizzly_out.json > "$line"_"$now".txt;tr ' ' '\t'<"$line"_"$now".txt > "$line"_"$now"_final.txt ; awk '{print $1,$2,$3,$4,$5,$6}' "$line"_"$now"_final.txt > "$line"_"$now"no_filter.txt ; awk '{if ($5&&$6 != 0 ) print $1,$2,$3,$4,$5,$6;}' "$line"_"$now"_final.txt > "$line"_"$now"_sc_pc_filter.txt ; rm "$line"_"$now".txt; done

cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF

#Organize Files
echo "Removing Intermediate Files and Refining Output...\n"
cd $outputdir
cp *.fasta Fastas/
cp *.json Json/
rm *.json
rm *.fasta
rm kallisto_*
rm filenames.txt
rm *_final.txt


cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF
echo "Pizzly Fusion Caller Run Complete for $inputdir!\n"
echo "Exiting now..."
cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF
