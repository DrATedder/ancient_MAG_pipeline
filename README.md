# ancient_MAG_pipeline

A shell script pipeline intended to aid in the creation of **M**etagenomic **A**ssembled **G**enomes (**MAGs**) from ancient microbiome samples. The paper associated with this data is available here: [Standeven *et al*.]().

## Pipeline software pre-requisites

- `Megahit v1.2.9`: A tool used to assemble large and complex *de novo* datasets [(Li *et al*. 2015)](https://pubmed.ncbi.nlm.nih.gov/25609793/). Can be accessed/downloaded [here](https://github.com/voutcn/megahit).

- `MetaBAT v2.12.1`: A tool used for contig binning [(Kang *et al*. 2015)](https://peerj.com/articles/1165/). Can be accessed/downloaded [here](https://bitbucket.org/berkeleylab/metabat/src/master/).

- `Maxbin2 v2.2.7`: A tool used for contig binning [(Wu *et al*. 2016)](https://academic.oup.com/bioinformatics/article/32/4/605/1744462). Can be accessed/downloaded [here](https://sourceforge.net/projects/maxbin2/).

- `CONCOCT v1.1.0`: A contig binning tool [(Alneberg *et al*. 2014)](https://www.nature.com/articles/nmeth.3103). Can be accessed/downloaded [here](https://github.com/BinPro/CONCOCT).

- `DAStool v1.1.6`: A tool to select high qualtiy bins [(Sieber *et al*. 2018)](https://www.nature.com/articles/s41564-018-0171-1). Can be accessed/downloaded [here](https://github.com/cmks/DAS_Tool).

- `CheckM v1.2.2`: A tool to assess completeness and contamination of MAGs [(Parks *et al*. 2015)](http://genome.cshlp.org/content/25/7/1043.short). Can be accessed/downloaded [here](https://github.com/Ecogenomics/CheckM).

- `GTDB-Tk v1.0.2`: A Genome Taxonomy Database Toolkit to classify MAGs [(Chaumeil *et al*. 2019)](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz848/5626182). Can be accessed/downloaded [here](https://github.com/Ecogenomics/GTDBTk).

- `Bowtie2 v2.5.0`: A tool for read mapping sequences to a reference genome [(Langmead and Salzberg, 2012)](https://academic.oup.com/bioinformatics/article/35/3/421/5055585?login=false). Can be accessed/downloaded [here](https://github.com/BenLangmead/bowtie2).

- `SAMtools v1.12`: A tool for interacting with high-throughput sequencing data [(Danecek *et al*. 2021)](https://pubmed.ncbi.nlm.nih.gov/33590861). Can be accessed/downloaded [here](http://www.htslib.org/).

- `Anvi'o`: A multi-omics platform for microbial genomics. Can be accessed/downloaded [here](https://anvio.org/).

## Considerations before beginning

- **Anvi'o tool use**: For ease of use, some of the tools used in this pipeline were employed through `Anvi'o`, including: `Bowtie2`, `CONCOCT` and `MetaBAT`. It is likely possible to perform this analysis without the use of `Anvi'o`, but the pipeline will need some tweeks.

- **Quality filtering**: This pipeline won't run any quality control steps (data cleaning, filtering, adaptor trimming etc...). We assume that you will have performed suitable steps for your own data type.

- **Sequence type(s)**: This pipeline assumes your data is in the form of paired-end (PE) illumina reads. While the steps *may* work on other data types, it is likely that you will need to tweak the script based on reviewing the documentation whihc packages with each of the software employed here.

- **Ancient vs. Modern**: While this pipeline is intended for use with aDNA, it will work similarly with modern sequence reads.

- **Job submission system**: Given the likely size of data files being run through this pipeline, we assume you will be using a batch job submission system. This pipeline assumes you will be using [SLURM](https://slurm.schedmd.com/documentation.html). If your system is different, please do take a look through the script to adjust relevant sections. This is mainly relevant for threading of various pipeline elements. For example: `-t ${SLURM_CPUS_PER_TASK}` If you are not using a job submission system, it *may* be possible to simply replace these incidences with an appropriate number of threads linked to your own system.

- **PATH**: Throughout the pipeline, we assume all software packages are installed in your `PATH`. If you do not have administrator access for install, it would be best to give explicit paths to each install.

- **Sample grouping**: Typically, one would create MAGs from groups of samples (grouped based on pre-defined criteria). `ancient_MAG_pipeline.sh` performs the pipeline analysis on single (paired-end) samples. This *may* reduce the number of MAGs you retrieve. It would likely be worth thinking about meaningful sample groupings (i.e. time-points, geographic locations etc...). The script can be tweaked to accomodate this.

## Directory creation and file path locations

Prior to running the script for the first time, you will need to input the file path for each of the required input and output directories (shown below). All output directories only require creating once (unless you specifically require different locations for different projects). Each sample you run through the pipleine will get it's own directory within each of these, so overwriting of directories/files between projects *should* be avoided.

```
### 1. Define all of the separate location paths required.

MEGAHIT=/path/to/desired/megahit/output
METABAT=/path/to/desired/metabat/output
MAXBIN=/path/to/desired/maxbin/output
CONCOCT=/path/to/desired/concoct/output
DAStool=/path/to/desired/dastool/output
CHECKM=/path/to/desired/checkm/output
DASOUT=/path/to/desired/dasout/output
GTDBTK=/path/to/desired/gtdbtk/output
```

## Cite ancient_MAG_pipeline

If you use `ancient_MAG_pipeline` in your study, please reference it using the details below. Please note, this pipeline leverages a lot of existing software (listed and linked above), please ensure you also give appropriate credit to the authors of those software.

Francesca J. Standeven, Gwyn Dahlquist-Axe, Camilla F. Speller, Conor J. Meehan and Andrew Tedder (XXXX) **An efficient pipeline for creating metagenomic-assembled genomes from ancient DNA in oral microbiome samples.**
