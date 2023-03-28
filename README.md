# Pipeline single Cell Perturb-seq like
---
## Running the Pipeline:  


```
ex: 
nextflow run main.nf -c perturb.config -with-timeline timeline_OUT_exp_1 -resume -w  perturbseq_OUT_exp_1
```
```
add the full path for your:
  - main.nf
  - perturb.config
  - timeline_OUT_exp_1  
  - perturbseq_OUT_exp_1
```

```
PARAMETERS:  
  -with-timeline timeline_OUT_exp_1  # This will generate a report with the resourcers for each step in a file timeline_OUT_exp_1  
  -resume  case the pipeline is interrupted return from the no concluded steps 
  - perturbseq_OUT_exp_1   this will be the outdirectory 
  
```



## config file 
### Create a perturb.config files with the following requirements    
---
```
params.GTF_GZ_LINK = 'http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz'
params.TRANSCRIPTOME_REFERENCE = "human"
params.KALLISTO_BIN = '/home/lf114/miniconda3/envs/perturbseq_pipeline/bin/kallisto'
params.GENOME = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
params.GUIDE_FEATURES = 'sgRNA_guide_metainfo_modified.xlsx'
params.CHEMISTRY = '10XV3'
params.THREADS = 15
params.DISTANCE_NEIGHBORS = 1000000
params.IN_TRANS = "FALSE"
params.FASTQ_FILES_TRANSCRIPTS = ['path/FOR/you/fasta/scRNAseq_reads_R1.fastq.gz path/FOR/you/fasta/scRNAseq_reads_R2.fastq.gz']
params.FASTQ_NAMES_TRANSCRIPTS = ['S1_L1']
params.FASTQ_FILES_GUIDES = ['path/FOR/you/fasta/scGUIDE_reads_R1.fastq.gz path/FOR/you/fasta/scGUIDE_reads_R2.fastq.gz'] 
params.FASTQ_NAMES_GUIDES = ['S1_L1']
params.CREATE_REF = false
```

__params.GTF_GZ_LINK:__  
The link to download the gzipped GTF file containing the genome annotation for ex: Homo sapiens (GRCh38.106 release).  
ex: 'http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz'  
__params.TRANSCRIPTOME_REFERENCE:__ 
The name of the transcriptome reference used for the alignment (ex: "human").  
This will be used to compute the mitochondrial genes aligments and use it as a covariate in the perturbation analysis.   
__params.KALLISTO_BIN:__   
The path to the kallisto binary used for the quantification of transcript abundance.(Find where kallisto is installed in your system)  
__params.GENOME:__   
The URL to download the gzipped reference genome file in FASTA format for hg38 (UCSC goldenPath).  
ex: 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'  
__params.GUIDE_FEATURES:__  
The filename of the Excel file containing the metadata information of the sgRNA guides used in the experiment. (see next section for more details)  
__params.CHEMISTRY:__  
The type of chemistry used for the single cell experiment (10XV3 in this case).  
- 10XV3
- 10VX2
- '0,0,16:0,16,26:1,0,0' for 5' PE  
To custom extract your barcode, UMI, guide follow the rule: The first part of your string (0,0,16) indicates that the Cell Barcode is in the first fast file (R1) and starts at position 1 (0 often in computing) and goes until position 16. The UMI is the (read1, position 16 until 26), The sequence is in the R2 file (1) and starts from the first position ultil the end...you can see the information here and learn about  different chemistry specification: <link>https://pachterlab.github.io/kallisto/manual</link>

__params.THREADS:__  
The number of threads used for parallel execution of the pipeline.  
__params.DISTANCE_NEIGHBORS:__   
The maximum distance (in bp) used to define the neighborhood of a gene for in-cis analysis.    
__params.IN_TRANS:__  
Case __"TRUE"__ it will test the differential perturbation between each guide and all the genes in the scRNAseq dataset.  
Case __"FALSE"__ it will test the differential perturbation between each guide and the genes in the neighborhood of the guide defined in the DISTANCE_NEIGHBORS parameter.  
__params.FASTQ_FILES_TRANSCRIPTS:__      
A list of input FASTQ files containing the scRNA-seq reads (R1 and R2) for the transcript-level analysis. 
The format should be ['file1_path_R1.fastq.gz file1_path_R2.fastq.gz', 'file2_path_R1.fastq.gz file2_path_R2.fastq.gz', ...].  
__params.FASTQ_NAMES_TRANSCRIPTS:__     
For each FASTQ file in FASTQ_FILES_TRANSCRIPTS, a sample name should be provided in FASTQ_NAMES_TRANSCRIPTS. following the structure
['S1_L1', 'S2_L1'']. Case the same sample came from same sample, but in a different lane, the sample name should be the same,
 but the lane should be different. ex: ['S1_L1', 'S2_L1'].  
__params.FASTQ_FILES_GUIDES:__     
a list of input FASTQ files containing the scGUIDE-seq reads (R1 and R2) for the perturbation analysis.
The format should be ['guideFile1_path_R1.fastq.gz guideFile1_path_R2.fastq.gz', 'guideFile2_path_R1.fastq.gz guideFile2_path_R2.fastq.gz', ...].
__params.FASTQ_NAMES_GUIDES:__     
For each FASTQ file in FASTQ_FILES_GUIDES, a sample name should be provided in FASTQ_NAMES_GUIDES. following the structure
['S1_L1', 'S2_L1'']. Case the same sample came from same sample, but in a different lane, the sample name should be the same,
 but the lane should be different. ex: ['S1_L1', 'S2_L1'].





### GUIDE_FEATURES format

<image src="https://raw.githubusercontent.com/LucasSilvaFerreira/pipeline_perturbseq_like/master/image/feature_example.png">  
Example Guide features file:<link>https://github.com/LucasSilvaFerreira/pipeline_perturbseq_like/blob/master/feature_example.xlsx</link>  

 
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
pip install nextflow
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
pip install muon
```

---

## TODO  JAMBOREE 

- [x] Create a docker image  
- [ ] Create a conda yml  
- [x] Decide on the experimental demo dataset (Maybe subsampling the gasperini 2019)  

## TODO pipeline

- [X]  checkedRUN AGAIN WITH A MODIFIED EXCEL (CHANGING THE NAME OF THE NEGATIVE CONTROL GENES)
- [X] Eliminate guides present with less than 30 cells
- [X] Fix the single batch problem 
- [X]  SOLVE THE PROBLEM TO RUN AGAINST ALL GENES (TRANS ANALYSIS) (Should be configured by guide in the future?)
- [x]  custom whitelist
- [ ]  ADD PARAMETERS TO SOME OF THE FUNCTIONS (Some ARE HARDCODED NOW ..EX: merge..number of cells and etc)
- [ ]  Create a nice description to the input files
- [ ]  fix the Rbase version to 4.2.2
- [ ]  change threads to cpus  (taks.cpu)
- [ ]  probably I can provide a list of elements and ask nextflow to handle the paralelizing (receiving the elements names in a channel?)
- [ ]  Create easing options to add 3'and 5'chemistry
- [ ]  Cell hashing and Multi-seq ? https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=719672

