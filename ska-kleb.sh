
ls *.fasta| cut -f 1 -d "." > id.txt
ls ${PWD}/*.fasta > path.txt
paste id.txt path.txt > skalist.txt




paste list.txt ls.txt > skalist.txt # creates a list tab separated for -f option name of sample and fasta name, ska also takes fastqs
ska build -o kleb -f ska_filelist.txt 
ska distance -o kleb_only_distance -m 1 --allow-ambiguous kleb.skf  --threads 10 # kleb.skf is a file from build
ska align --min-freq 1 --filter no-filter -o klebska_alignment  kleb.skf
raxmlHPC -s bcc_align -n bcc -m GTRCAT -T 10 -f a -x 1234 -N 100 -p 456


ska build -o ecoli -f ska_filelist.txt 
ska align --min-freq 1 --filter no-filter -o ecoli_alignment ecoli.skf
raxmlHPC -s ecoli_alignment -n ecoli -m GTRCAT -T 10 -f a -x 1234 -N 100 -p 456



for sample in `cat ../ecoli_list`; do amrfinder --name ${sample} --organism Escherichia --nucleotide ${sample}.fasta -o ./amr_finder/${sample}.txt; done


for sample in `cat ../kleb_list.txt`; do amrfinder --name ${sample} --organism Klebsiella_pneumoniae --nucleotide ${sample}.fasta -o ./amr_finder/${sample}.txt; done

