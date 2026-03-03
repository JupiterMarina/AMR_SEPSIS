mash sketch -o assemblies_sketch *.fasta
mash dist /home/kanye/Downloads/refseq.genomes.k21s1000.msh assemblies_sketch.msh > mash.txt
##refseq in /home/kanye/marina/julius-pt/run2/missing
awk '{print $2, $1, $3, $4, $5, $6}' mash.txt | sort -k1,1 -k3,3n | awk '!seen[$1]++' > best_hits.txt
cat best_hits.txt | cut -f 2 -d " " | cut -d "A" -f 1 | sed 's/_$//g' >>ass2.txt
conda activate mlst 
for ass in `cat ass2.txt`; do esearch -db assembly -query ${ass} |  esummary | xtract -pattern DocumentSummary -element AssemblyAccession,Organism >> hits.txt; done
