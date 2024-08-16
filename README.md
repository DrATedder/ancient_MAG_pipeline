# ancient_MAG_pipeline

A shell script pipeline intended to aid in the creation of **M**etagenomic **A**ssembled **G**enomes (**MAGs**) from ancient microbiome samples. The paper associated with this data is available here: [Standeven *et al*.]().

## Pipeline software pre-requisites

- **Megahit v1.2.9**: A tool used to assemble large and complex *de novo* datasets [(Li *et al*. 2015)](https://pubmed.ncbi.nlm.nih.gov/25609793/). Can be accessed/downloaded [here](https://github.com/voutcn/megahit).

- **MetaBAT v2.12.1**: A tool used for contig binning [(Kang *et al*. 2015)](https://peerj.com/articles/1165/). Can be accessed/downloaded [here](https://bitbucket.org/berkeleylab/metabat/src/master/).

- **Maxbin2 v2.2.7**: A tool used for contig binning [(Wu *et al*. 2016)](https://academic.oup.com/bioinformatics/article/32/4/605/1744462). Can be accessed/downloaded [here](https://sourceforge.net/projects/maxbin2/).

- **CONCOCT v1.1.0**: A contig binning tool [(Alneberg *et al*. 2014)](https://www.nature.com/articles/nmeth.3103). Can be accessed/downloaded [here](https://github.com/BinPro/CONCOCT).

- **DAStool v1.1.6**: A tool to select high qualtiy bins [(Sieber *et al*. 2018)](https://www.nature.com/articles/s41564-018-0171-1). Can be accessed/downloaded [here](https://github.com/cmks/DAS_Tool).

- **CheckM v1.2.2**: A tool to assess completeness and contamination of MAGs [(Parks *et al*. 2015)](http://genome.cshlp.org/content/25/7/1043.short). Can be accessed/downloaded [here](https://github.com/Ecogenomics/CheckM).

- **GTDB-Tk v1.0.2**: A Genome Taxonomy Database Toolkit to classify MAGs [(Chaumeil *et al*. 2019)](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz848/5626182). Can be accessed/downloaded [here](https://github.com/Ecogenomics/GTDBTk).

- **Bowtie2 v2.5.0**: A tool for read mapping sequences to a reference genome [(Langmead and Salzberg, 2012)](https://academic.oup.com/bioinformatics/article/35/3/421/5055585?login=false). Can be accessed/downloaded [here](https://github.com/BenLangmead/bowtie2).

- **SAMtools v1.12**: A tool for interacting with high-throughput sequencing data [(Danecek *et al*. 2021)](https://pubmed.ncbi.nlm.nih.gov/33590861). Can be accessed/downloaded [here](http://www.htslib.org/).

## Considerations before beginning

- **Quality filtering**: This pipeline won't run any quality control steps (data cleaning, filtering, adaptor trimming etc...). We assume that you will have performed suitable steps for your own data type.

- **Sequence type(s)**: This pipeline assumes your data is in the form of paired-end (PE) illumina reads. While the steps *may* work on other data types, it is likely that you will need to tweak the script based on reviewing the documentation whihc packages with each of the software employed here.

- **Ancient vs. Modern**: While this pipeline is intended for use with aDNA, it will work similarly with modern sequence reads.

- **Job submission system**: Given the likely size of data files being run through this pipeline, we assume you will be using a batch job submission system. This pipeline assumes you will be using [SLURM](https://slurm.schedmd.com/documentation.html). If your system is different, please do take a look through the script to adjust relevant sections. This is mainly relevant for threading of various pipeline elements. For example:

  ```megahit -1 ${DECONTAM_READS}/${name}_R1.fastq.gz -2 ${DECONTAM_READS}/${name}_R2.fastq.gz -o ${MEGAHIT}/${name}_coassembly -t ${SLURM_CPUS_PER_TASK} --kmin-1pass --no-mercy```

