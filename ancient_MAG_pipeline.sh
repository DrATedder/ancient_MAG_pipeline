#!/bin/sh
### SLURM submission parameters should be updated if you use this system, or removed/substituted if not.

#SBATCH --job-name=MAG_pipeline
#SBATCH --account= ### Details needed
#SBATCH --partition= ### Details needed
#SBATCH --cpus-per-task=20  ### Adjust as required
#SBATCH --time=24:00:00 ### Adjust as required



### 1. Define all of the separate location paths required.

DECONTAM_READS= ### This location can be changed to where your samples are stored
MEGAHIT= ### Needs creating
METABAT= ### Needs creating
MAXBIN= ### Needs creating
CONCOCT= ### Needs creating
DAStool= ### Needs creating
CHECKM= ### Needs creating
DASOUT= ### Needs creating
GTDBTK= ### Needs creating


### 2. Perform key actions - listed below

for i in ${DECONTAM_READS}/*_R1.fastq.gz ### Needs correcting for actual sequence names
  do
  name=$(basename ${i} _R1.fastq.gz)
  echo -e '\n***** Beginning sample ${name} analysis. *****\n'
  ### MEGAHIT
  echo -e '\n***Begin MEGAHIT***\n'
  anvio megahit -1 ${DECONTAM_READS}/${name}_R1.fastq.gz -2 ${DECONTAM_READS}/${name}_R2.fastq.gz -o ${MEGAHIT}/${name}_coassembly -t ${SLURM_CPUS_PER_TASK} --kmin-1pass --no-mercy
  anvio anvi-script-reformat-fasta ${MEGAHIT}/${name}_coassembly/final.contigs.fa -o ${MEGAHIT}/${name}_coassembly/final.contigs_simplified_2500.fa -l 2500 --simplify-names
  ### bowtie2
  echo -e '\n***Bowtie2 Index***\n'
  anvio bowtie2-build ${MEGAHIT}/${name}_coassembly/final.contigs_simplified_2500.fa ${MEGAHIT}/${name}_coassembly/${name}
  echo -e '\n***Bowtie2 assembly***\n'
  anvio bowtie2 -x  ${MEGAHIT}/${name}_coassembly/${name} -q -1 ${DECONTAM_READS}/${name}_R1.fastq.gz -2 ${DECONTAM_READS}/${name}_R2.fastq.gz --no-unal -p ${SLURM_CPUS_PER_TASK} -S ${MEGAHIT}/${name}_coassembly/${name}_coassembly.sam
  echo -e '\n***Samtools view***\n'
  samtools view -b -o ${MEGAHIT}/${name}_coassembly/${name}_coassembly-raw.bam ${MEGAHIT}/${name}_coassembly/${name}_coassembly.sam
  echo -e '\n***Samtools sort***\n'
  samtools sort -o ${MEGAHIT}/${name}_coassembly/${name}_coassembly.bam ${MEGAHIT}/${name}_coassembly/${name}_coassembly-raw.bam
  echo -e '\n***Samtools Index***\n'
  samtools index ${MEGAHIT}/${name}_coassembly/${name}_coassembly.bam
  ### clean up intermediate files
  rm ${MEGAHIT}/${name}_coassembly/${name}_coassembly-raw.bam
  rm ${MEGAHIT}/${name}_coassembly/${name}_coassembly.sam
  ### CONCOCT
  echo -e '\n***CONCOCT***\n'
  anvio cut_up_fasta.py ${MEGAHIT}/${name}_coassembly/final.contigs_simplified_2500.fa -c 10000 -o 0 --merge_last -b ${CONCOCT}/${name}_final.contigs_simplified_10k.bed > ${CONCOCT}/${name}_final.contigs_simplified_10k.fa
  anvio concoct_coverage_table.py ${CONCOCT}/${name}_final.contigs_simplified_10k.bed ${MEGAHIT}/${name}_coassembly/${name}_coassembly.bam > ${CONCOCT}/${name}_final.contigs_simplified_concoct_coverage_table.tsv
  anvio concoct --composition_file ${CONCOCT}/${name}_final.contigs_simplified_10k.fa --coverage_file ${CONCOCT}/${name}_final.contigs_simplified_concoct_coverage_table.tsv -b ${CONCOCT}/${name}_final.contigs -t ${SLURM_CPUS_PER_TASK}
  anvio merge_cutup_clustering.py ${CONCOCT}/${name}_final.contigs_clustering_gt1000.csv > ${CONCOCT}/${name}_final.contigs_clustering_merged.csv
  anvio extract_fasta_bins.py ${MEGAHIT}/${name}_coassembly/final.contigs_simplified_2500.fa ${CONCOCT}/${name}_final.contigs_clustering_merged.csv --output_path ${CONCOCT}/${name}/${name}_concoct_final.contigs/
  echo -e '\n***METABAT***\n'
  anvio jgi_summarize_bam_contig_depths --outputDepth ${MEGAHIT}/${name}_coassembly/${name}_mapped_depth.txt --pairedContigs ${MEGAHIT}/${name}_coassembly/${name}_paired.txt ${MEGAHIT}/${name}_coassembly/${name}_coassembly.bam
  anvio metabat -i ${MEGAHIT}/${name}_coassembly/final.contigs_simplified_2500.fa -a ${MEGAHIT}/${name}_coassembly/${name}_mapped_depth.txt -m 2000 --saveCls -t ${SLURM_CPUS_PER_TASK} -o ${METABAT}/${name}/${name}_output-file
  echo -e '\n***MAXBIN***\n'
  pileup.sh in=${MEGAHIT}/${name}_coassembly/${name}_coassembly.bam out=${MEGAHIT}/${name}_coassembly/${name}_cov.txt
  awk '{print$1"\t"$5}' ${MEGAHIT}/${name}_coassembly/${name}_cov.txt | grep -v '^#' > ${MEGAHIT}/${name}_coassembly/${name}_abundance.txt
  perl /storage02/or-microbio/Maxbin/MaxBin-2.2.7/run_MaxBin.pl -thread ${SLURM_CPUS_PER_TASK} -contig ${MEGAHIT}/${name}_coassembly/final.contigs_simplified_2500.fa -out ${MAXBIN}/${name}/${name}_Maxbin2 -abund ${MEGAHIT}/${name}_coassembly/${name}_abundance.txt
  echo -e '\n***DAS_Tool prep.***\n'
  sh /storage02/or-microbio/tool_scripts/Fasta_to_Contigs2Bin.sh -i ${METABAT}/${name}/ -e fa >${DAStool}/${name}_MetaBAT2_bins.tsv
  sh /storage02/or-microbio/tool_scripts/Fasta_to_Contigs2Bin.sh -i ${CONCOCT}/${name}/${name}_concoct_final.contigs/ -e fa >${DAStool}/${name}_concoct_bins.tsv
  sh /storage02/or-microbio/tool_scripts/Fasta_to_Contigs2Bin.sh -i ${MAXBIN}/${name}/ -e fasta >${DAStool}/${name}_Maxbin_bins.tsv
  echo -e '\n***CHECKM***\n'
  checkm lineage_wf ${CONCOCT}/${name}/${name}_concoct_final.contigs/ ${CONCOCT}/${name}/${name}_concoct_final.contigs/checkm_concoct -x .fa -t ${SLURM_CPUS_PER_TASK} --reduced_tree > ${CHECKM}/${name}_concoct_output.txt
  checkm lineage_wf ${METABAT}/${name}/ ${METABAT}/${name}/checkm_metabat -x .fa -t ${SLURM_CPUS_PER_TASK} --reduced_tree > ${CHECKM}/${name}_metabat_output.txt
  checkm lineage_wf ${MAXBIN}/${name}/ ${MAXBIN}/${name}/checkm_maxbin -x .fasta -t ${SLURM_CPUS_PER_TASK} --reduced_tree > ${CHECKM}/${name}_maxbin_output.txt
  ### DAStool
  echo -e '\n***DAS_Tool***\n'
  if [ ! -d ${DASOUT}/${name} ]; then
    mkdir ${DASOUT}/${name}
  fi
  DAS_Tool -i ${DAStool}/${name}_concoct_bins.tsv,${DAStool}/${name}_Maxbin_bins.tsv,${DAStool}/${name}_MetaBAT2_bins.tsv -l concoct,maxbin,metabat -c ${MEGAHIT}/${name}_coassembly/final.contigs_simplified_2500.fa -o ${DASOUT}/${name}/${name} -t ${SLURM_CPUS_PER_TASK} --write_bin_evals
  ### Change maxbin fasta to fa in DASOUT
  echo -e '\n***Prep. maxbin outputs for GTDBTK***\n'
  ### Move all fasta files from the binning tools into the correct place
  if [ ! -d ${DASOUT}/${name}/fasta_files ]; then
    mkdir ${DASOUT}/${name}/fasta_files
  fi
  for file in ${MAXBIN}/${name}/*.fasta
  do
    tmp_name=$(basename ${file} .fasta)
    mv ${file} ${DASOUT}/${name}/fasta_files/${tmp_name}.fa
  done
  for file in ${CONCOCT}/${name}/${name}_concoct_final.contigs/*.fa
  do
    tmp_name=$(basename ${file} .fa)
    mv ${file} ${DASOUT}/${name}/fasta_files/concoct_${tmp_name}.fa
  done
  for file in ${CONCOCT}/${name}/${name}_concoct_final.contigs/*.fa
  do
    tmp_name=$(basename ${file} .fa)
    mv ${file} ${DASOUT}/${name}/fasta_files/concoct_${tmp_name}.fa
  done
    for file in ${METABAT}/${name}/*.fa
  do
    tmp_name=$(basename ${file} .fa)
    mv ${file} ${DASOUT}/${name}/fasta_files/metabat_${tmp_name}.fa
  done
  echo -e '\n***GTDBTK***\n'
  if [ ! -d ${GTDBTK}/${name}/${name}_gtdbtk_output ]; then
    mkdir ${GTDBTK}/${name}/${name}_gtdbtk_output
  fi
  gtdbtk classify_wf --genome_dir ${DASOUT}/${name}/fasta_files --out_dir ${GTDBTK}/${name}/${name}_gtdbtk_output  --cpus ${SLURM_CPUS_PER_TASK} --extension fa
done
