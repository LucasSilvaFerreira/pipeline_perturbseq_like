# pipeline_perturbseq_like
---
## Running the Pipeline:  


### config file 
Create a perturb.config files with the following requirements    
---
```
params.GTF_GZ_LINK = 'http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz'
params.TRANSCRIPTOME_REFERENCE = "human"
params.KALLISTO_BIN = '/home/lf114/miniconda3/envs/perturbseq_pipeline/bin/kallisto'
params.GENOME = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
params.GUIDE_FEATURES = 'sgRNA_guide_metainfo_modified.xlsx'
params.CHEMISTRY = '0,0,16:0,16,26:1,0,0'
params.THREADS = 15
params.DISTANCE_NEIGHBORS = 1000000
params.IN_TRANS = "FALSE"
params.FASTQ_FILES_TRANSCRIPTS = ['path/FOR/you/fasta/scRNAseq_reads_R1.fastq.gz path/FOR/you/fasta/scRNAseq_reads_R2.fastq.gz']
params.FASTQ_NAMES_TRANSCRIPTS = ['S1_L1']
params.FASTQ_FILES_GUIDES = ['path/FOR/you/fasta/scGUIDE_reads_R1.fastq.gz path/FOR/you/fasta/scGUIDE_reads_R2.fastq.gz'] 
params.FASTQ_NAMES_GUIDES = ['S1_L1']
params.CREATE_REF = false
```


### GUIDE_FEATURES format



### Requirements  
Creating the conda env

```

conda create -n perturbseq_like_pipeline python=3.8
conda activate perturbseq_like_pipeline
conda install -c conda-forge mamba -y
```

installs:

Following this order
```
pip uninstall GTFProcessing -y
pip install git+https://github.com/LucasSilvaFerreira/GTFProcessing.git
pip install gtfparse==1.3.0
pip install git+https://github.com/LucasSilvaFerreira/Perturb_Loader.git
mamba install -c bioconda nextflow -y
mamba install -c bioconda kallisto -y
pip install --quiet kb-python
mamba install -c anaconda openpyxl -y 
mamba install -c conda-forge r-base -y 
Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
Rscript -e 'BiocManager::install("Rhdf5lib")'
Rscript -e 'BiocManager::install("rhdf5")'
Rscript -e 'install.packages("doParallel", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'
mamba install -c conda-forge r-gert -y 
mamba install -c conda-forge r-ragg -y
Rscript -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'
mamba install -c conda-forge r-ggplot2
Rscript -e 'devtools::install_github("katsevich-lab/sceptre")'
```

---

## TODO PIPELINE



