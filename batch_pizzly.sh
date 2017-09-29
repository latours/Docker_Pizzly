#!/bin/bash
###CREATOR: SARA LATOUR
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
echo -e "\nPizzly Fusion Caller for RNA-seq Read Pairs:"
echo -e "\nMAINTAINER:\nSara Latour"
echo -e "\nCONTACT:\nsaralatour@outlook.com"
#Date:
now="$(date +'%d_%m_%Y')"
echo -e "\nDATE(dd_mm_yyyy):\n$now\n"
cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF

#CMD Line Arguments:
inputdir=$1
outputdir=$2
if [ $# -ne 2 ]
  then
    echo -e "\nERROR MESSAGE:\nPlease ensure that the full path to the Input and Output Directories are provided:\n"
	exit 1
fi

echo -e "\nChecking if Input directory is valid..."
if   [ -d "${inputdir}" ]
then echo -e "✓"
else echo -e "Input is not a valid directory... exiting now\n";
     exit 1
fi


echo -e "Checking if Output directory is valid..."
if   [ -d "${outputdir}" ]
then echo -e "✓"
else echo -e "Output is not a valid directory... exiting now\n";
     exit 1
fi


#Input:
echo -e "The Directory Selected is:$inputdir\n"

#Create Output File
echo -e "Creating Sub-Directories for Output..."
mkdir -p $outputdir/Fastas
mkdir -p $outputdir/Json
mkdir -p $outputdir/Pizzly_Run_Information
mkdir -p $outputdir/FusionCaller_Output_"$now"

#Create File Lists - Note to self: Change it to find later.
cd $inputdir
find *."fastq.gz" > $outputdir/filenames.txt 
cd $outputdir

cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF

#RUN Index
echo -e "\nBeginning Index Step for Kallisto...\n"
$kallistopath/kallisto index -i $outputdir/index.idx -k 31 $cdna 
echo -e "Completed Index Step!\n"

#RUN Kallisto Quant & Pizzly
cat filenames.txt | while read line; do sed 's/R[1,2]_001.fastq.gz//g' ; done > $outputdir/kallisto_list.txt
sort -u $outputdir/kallisto_list.txt > $outputdir/kallisto_unique_list.txt
echo -e "Running Pizzly Fusion Caller on the following RNA-seq Pairs IDs...\n"
cat $outputdir/kallisto_unique_list.txt
echo -e "\nStarting Kallisto and Pizzly runs for files from $inputdir directory now...\n"
cat $outputdir/kallisto_unique_list.txt | while read line; do echo "$line"; $kallistopath/kallisto quant -i $outputdir/index.idx --fusion -o $outputdir/"$line"pizzly_out $inputdir/"$line"R1_001.fastq.gz $inputdir/"$line"R2_001.fastq.gz ; $pizzlypath/build/pizzly -k 31 --gtf $GTF --cache index.cache.txt --align-score 2 --insert-size 400 --fasta $cdna --output "$line"pizzly_out $outputdir/"$line"pizzly_out/fusion.txt ; $pizzlypath/scripts/flatten_json.py "$line"pizzly_out.json > "$line""$now".txt;tr ' ' '\t'<"$line""$now".txt > "$line""$now"_final.txt ; awk '{print $1,$2,$3,$4,$5,$6}' "$line""$now"_final.txt > "$line""$now"_no_filter.txt ; awk '{if ($5&&$6 != 0 ) print $1,$2,$3,$4,$5,$6;}' "$line""$now"_final.txt > "$line""$now"_sc_pc_filter.txt ; rm "$line""$now".txt; done
cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF

#Organize Files
echo -e "Removing Intermediate Files and Organizing Output...\n"
cd $outputdir
cp *.fasta Fastas/
cp *.json Json/
cp -R *pizzly_out Pizzly_Run_Information/
rm -R *pizzly_out 
cp *filter.txt FusionCaller_Output_"$now"/
rm *filter.txt
rm *.json
rm *.fasta
rm kallisto_*
rm filenames.txt
rm *_final.txt
echo -e "Final Output Files are located in the FusionCaller_Output_"$now"/ Directory" > FinalOutputLocation.txt

cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF
echo -e "Pizzly Fusion Caller Run Complete for $inputdir!\n"
echo -e "Exiting now..."
cat << "EOF"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
EOF
