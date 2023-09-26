// Declare syntax version
nextflow.enable.dsl=2
// Script parameters
//params-ALIGN_MEMORY = 20G
//params.GTF_GZ_LINK = 'http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz'
//params.TRANSCRIPTOME_REFERENCE = "human"
//params.KALLISTO_BIN = '/home/lf114/miniconda3/envs/perturbseq_pipeline/bin/kallisto'
//params.GENOME = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
//params.GUIDE_FEATURES = '5p27sgRNA_guide_metainfo_modified.xlsx'
//params.CHEMISTRY = '0,0,16:0,16,26:0,26,0,1,0,0'
//params.THREADS = 15
//params.DISTANCE_NEIGHBORS = 1000000
//params.IN_TRANS = "FALSE"
//params.FASTQ_FILES_TRANSCRIPTS = ['5p27sgRNAGex_02KRWD_11408_S3_L001_R1_001.fastq.gz 5p27sgRNAGex_02KRWD_11408_S3_L001_R2_001.fastq.gz']
//params.FASTQ_NAMES_TRANSCRIPTS = ['S1_L1']
//params.FASTQ_FILES_GUIDES = ['5p27sgRNAsgRNA_02KRWK_11408_S15_L001_R1_001.fastq.gz 5p27sgRNAsgRNA_02KRWK_11408_S15_L001_R2_001.fastq.gz' ]
//params.FASTQ_NAMES_GUIDES = ['S1_L1']
//params.CREATE_REF = false
//params.DIRECTION  = "left"  "right"  "both"
//params.RUN_MULTISEQ = true
//params.R1_MULTI =                     '/n/data1/bch/hemonc/bauer/archana/NovaSeq5_second_release/01.RawData/MSeqhighMOI_RPI_3/MSeqhighMOI_RPI_3_CKDL230003506-1A_HT2MJDSX5_L4_1.fq.gz' 
//params.R2_MULTI =                     '/n/data1/bch/hemonc/bauer/archana/NovaSeq5_second_release/01.RawData/MSeqhighMOI_RPI_3/MSeqhighMOI_RPI_3_CKDL230003506-1A_HT2MJDSX5_L4_2.fq.gz'
//params.BARCODES_MULTIBAR_LIST_MULTI = '/n/data1/bch/hemonc/bauer/lucassilva/yanhua_perturb/high_moi_multi_barcode.csv'
//params.BAR_MULTI= [1,16]         
//params.UMI_MULTI= [17,28]         
//params.R2_MULTI_TAG = [1,8]          



workflow {

   
   dir_images_composition_scrna =  compositionREADSscRNA(Channel.from(params.FASTQ_NAMES_TRANSCRIPTS), Channel.from(params.FASTQ_FILES_TRANSCRIPTS))
   dir_images_composition_guides = compositionREADSGuides(Channel.from(params.FASTQ_NAMES_GUIDES), Channel.from(params.FASTQ_FILES_GUIDES))

   gtf_out = downloadGTF(params.GTF_GZ_LINK)
   downloadReference(Channel.of(params.TRANSCRIPTOME_REFERENCE), Channel.of(params.KALLISTO_BIN) )
   downloadGenome (Channel.of(params.GENOME))
   guide_feature_preprocessed = guidePreprocessing(params.GUIDE_FEATURES)
   creatingGuideRef( downloadGenome.out.genome, Channel.of(params.KALLISTO_BIN), guide_feature_preprocessed.guide_features, params.CREATE_REF )
    
    
   map_rna = mappingscRNA (
                 Channel.from(params.FASTQ_NAMES_TRANSCRIPTS),
                 Channel.from(params.FASTQ_FILES_TRANSCRIPTS),
                 downloadReference.out.transcriptome_idx.collect(),
                 downloadReference.out.t2t_transcriptome_index.collect(),
                 Channel.of(params.KALLISTO_BIN).collect(),
                 Channel.of(params.CHEMISTRY).collect(),
                 Channel.of(params.THREADS).collect(),
		 Channel.of(params.WHITELIST).collect(),
                )
    
   map_guide = mappingGuide (
                 Channel.from(params.FASTQ_NAMES_GUIDES),
                 Channel.from(params.FASTQ_FILES_GUIDES), 
                 creatingGuideRef.out.guide_index.collect(),
                 creatingGuideRef.out.t2tguide_index.collect(),
                 Channel.of(params.KALLISTO_BIN).collect(),
                 Channel.of(params.CHEMISTRY).collect(),
                 Channel.of(params.THREADS).collect(),
		 Channel.of(params.WHITELIST).collect())
    
    
    dir_count_files = map_guide.ks_guide_out.join(map_rna.ks_transcripts_out, remainder: true).flatten().toList().view()
    //dir_count_files.view()
    
    df_initialized = preprocessing(dir_count_files)
    

    
    //out_processed = filtering(df_initialized.df_initial_files, dir_count_files )
    
    concact_prefiltering_out = concact_prefiltering(df_initialized.df_initial_files,
                                                    dir_count_files,
                                                params.EXPECTED_CELL_NUMBER,
                                                params.MITO_SPECIE,
                                                params.MITO_EXPECTED_PERCENTAGE,
                                                params.PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS,
                                                params.TRANSCRIPTS_UMI_TRHESHOLD,
                                                )
    
    

    //concact_prefiltering_out.guide_ann
    //concact_prefiltering_out.transcripts_ann
    
    
        
    moun_raw_creation = moun_raw_creation(  concact_prefiltering_out.guide_ann,
                                            concact_prefiltering_out.transcripts_ann,
                                            downloadReference.out.t2t_transcriptome_index,     
                                            gtf_out.gtf )
    
    
    
    
    
    merge_bin_and_muon_out =  merge_bin_and_muon( moun_raw_creation.raw_moun_data, params.GUIDE_UMI_LIMIT, params.MERGE )



    out_pre_bar_multiseq = preprocess_bar_multiseq(merge_bin_and_muon_out.muon_processed)
    out_pre_bar_multiseq.view()
    
    
    multi_out = MultiSeq(params.R1_MULTI,
                params.R2_MULTI,
                out_pre_bar_multiseq.cell_barcode_to_multi,
                params.BARCODES_MULTIBAR_LIST_MULTI ,
                params.BAR_MULTI,
                params.UMI_MULTI,
                params.R2_MULTI_TAG,
                merge_bin_and_muon_out.muon_processed)



    muon_final =  select_final_muon(merge_bin_and_muon_out.muon_processed, multi_out.muon_with_multiseq)

    pert_loader_out = PerturbLoaderGeneration(muon_final.final_muon , gtf_out.gtf , params.DISTANCE_NEIGHBORS, params.IN_TRANS, params.ADDGENENAMES )
    runSceptre_out = runSceptre(pert_loader_out.perturb_piclke, params.DIRECTION)
    create_anndata_from_sceptre_out = create_anndata_from_sceptre(runSceptre_out.sceptre_out_dir, muon_final.final_muon )
    
    
}


process downloadGTF {
    input:
        val gtf_gz_path
    output:
        path "transcripts.gtf", emit: gtf
    script:
        if (params.CUSTOM_REFERENCE == false)
        """    
        wget -O - $gtf_gz_path | gunzip -c > transcripts.gtf
        """ 
        if (params.CUSTOM_REFERENCE == true)
        """ 
        cp $params.CUSTOM_GTF_PATH transcripts.gtf
        """ 
}

process downloadReference {
    input:
        val ref_name 
        path k_bin 
    
    output:
        path "transcriptome_index.idx" , emit: transcriptome_idx
        path "transcriptome_t2g.txt"   , emit: t2t_transcriptome_index
    script:
        if (params.CUSTOM_REFERENCE == false)
            """
                kb ref -d $ref_name -i transcriptome_index.idx -g transcriptome_t2g.txt --kallisto ${k_bin}
            """
        if (params.CUSTOM_REFERENCE == true)
            """
                cp $params.CUSTOM_REFERENCE_IDX  transcriptome_index.idx
                cp $params.CUSTOM_REFERENCE_T2T  transcriptome_t2g.txt
            """
}

process downloadGenome {
    input:
    val genome_path
    output:
    path "genome.fa.gz" , emit: genome
    script:
    """
        wget $genome_path -O genome.fa.gz 
    """


}


process guidePreprocessing {
    cache 'lenient'
    debug true
    input:
    path (guide_input_table)
    output:
    path "guide_features.txt" , emit: guide_features
    script:
    """    
    guide_table_processing.py  $guide_input_table
    """
}





process creatingGuideRef {
    cache 'lenient'
    input:
    val genome_path
    path k_bin
    path guide_features
    val create_ref
    output:
    path "guide_index.idx" ,  emit: guide_index
    path "t2guide.txt" , emit: t2tguide_index
    script:
    
    """
        kb ref -i guide_index.idx -f1 $genome_path -g t2guide.txt --kallisto $k_bin  --workflow kite $guide_features 

    """


    
}

process compositionREADSscRNA {
    debug true
    input:
    tuple val(out_name_dir)
    tuple val(string_fastqz)
    output:
    path ("transcript_${out_name_dir}_composition"),  emit: composition_plot_dir

    script:
    
        """
        mkdir transcript_"${out_name_dir}_composition"
        fq_composition.py $string_fastqz transcript_${out_name_dir}                                                                
        """
} 


process compositionREADSGuides {
    debug true
    input:
    tuple val(out_name_dir)
    tuple val(string_fastqz)
    output:
    path ("guide_${out_name_dir}_composition"),  emit: composition_plot_dir

    script:
    
        """
        mkdir "guide_${out_name_dir}_composition"
        fq_composition.py $string_fastqz guide_${out_name_dir}                                                                
        """
} 






process mappingscRNA {
    cache 'lenient'
    cpus 10
    debug true
    input:
    tuple val(out_name_dir)
    tuple val(string_fastqz)
    tuple path(transcriptome_idx)
    tuple path(t2t_transcriptome_index)
    tuple path(k_bin)
    tuple val(chemistry)
    tuple val(threads)
    tuple val(whitelist)
    
    output:
    path ("${out_name_dir}_ks_transcripts_out"),  emit: ks_transcripts_out

    script:
    
        """
	echo "kb count -i $transcriptome_idx -g  $t2t_transcriptome_index --verbose --workflow kite -w $whitelist --h5ad --kallisto $k_bin -x $chemistry -o ${out_name_dir}_ks_transcripts_out -t $threads $string_fastqz  --overwrite" -m 48G
        kb count -i $transcriptome_idx -g  $t2t_transcriptome_index --verbose --workflow kite -w $whitelist --h5ad --kallisto $k_bin -x $chemistry -o ${out_name_dir}_ks_transcripts_out -t ${task.cpus} $string_fastqz  --overwrite    -m 48G                                                                
        """
} 

process mappingGuide {
    cache 'lenient'
    cpus 10
    debug true
    input:
    tuple val(out_name_dir)
    tuple val(string_fastqz)
    tuple (path guide_index)
    tuple (path t2tguide_index)
    tuple path(k_bin)
    tuple val(chemistry)
    tuple val(threads)
    tuple val(whitelist)

    output:
    path ("${out_name_dir}_ks_guide_out"),  emit: ks_guide_out
    script:
        """
        kb count -i $guide_index       -g  $t2tguide_index --verbose   --report  --workflow kite -w $whitelist  --h5ad --kallisto $k_bin -x $chemistry  -o ${out_name_dir}_ks_guide_out -t ${task.cpus} $string_fastqz --overwrite -m 48G
        """
} 


process capture_variables_and_save_list{
    debug true
    input:
    val received
    val out_name
    output:
    path "${out_name}.txt",  emit: out_file
    script:
    """
    echo  '${received}' > ${out_name}.txt 
    """
}


process preprocessing {
    debug true
    input:
    path (count_list)
    output:
    path 'initial_preprocessing_file_names.txt', emit: df_initial_files
    script:
    """    
    preprocessing.py  ${count_list} 

    """   
    
}


process filtering{
    debug true
    input:
    path (path_df)
    path (all_files_context)
    output:
    path 'results_per_lane/processed_anndata_guides_data.h5ad', emit: guide_ann
    path 'results_per_lane/processed_anndata_transcripts_data.h5ad',  emit: transcripts_ann
    script:
    """
    #use -merge to merge the guides
    # I need to add these parameters to the pipeline config  mito...cellnumber...merge...guide_limit
    filtering_and_lane_merging.py --path ${path_df} --expected_cell_number 8000 --mito_specie hsapiens --mito_expected_percentage 0.2 --percentage_of_cells_to_include_transcript 0.2  --guide_umi_limit 5

    """
}




process concact_prefiltering{
    debug true
    input:
    path (path_df)
    path (dir_count_context)
    val  (expected_cell_number)
    val  (mito_specie)
    val  (mito_expected_percentage)
    val (percentage_of_cells_to_include_transcript)
    val  (transcripts_umi_treshold)
    output:
    path 'results_per_lane/full_raw_guide_ann_data.h5ad', emit: guide_ann
    path 'results_per_lane/full_raw_scrna_ann_data.h5ad',  emit: transcripts_ann
    
    
    
    
    script:
    """
    #chmod 700 /n/scratch3/users/l/lf114/pipeline_perturbseq_like/bin/concact_and_pre_filtering.py; 
    concact_and_pre_filtering.py --path ${path_df} --expected_cell_number $expected_cell_number --mito_specie $mito_specie --mito_expected_percentage $mito_expected_percentage --percentage_of_cells_to_include_transcript  $percentage_of_cells_to_include_transcript --transcripts_umi_treshold $transcripts_umi_treshold 
    """
}




process moun_raw_creation{
    debug true
    input:
    path (ann_guide)
    path (ann_exp)
    path (transcript_file)
    path (gtf_in)
    output:
    path 'raw_mudata_guide_and_transcripts.h5mu', emit: raw_moun_data

    
    script:
    """
    #chmod 700 /n/scratch3/users/l/lf114/pipeline_perturbseq_like/bin/muon_creation.py; 
    muon_creation.py --ann_guide $ann_guide --ann_exp $ann_exp --gtf_in $gtf_in --transcript_file $transcript_file

    
    """
}

process merge_bin_and_muon {
    debug true
    input:
    path (muon_data)
    val  (guide_umi_limit)
    val (merge)
    output:
    path 'processed_mudata_guide_and_transcripts.h5mu', emit: muon_processed
    script:
    """
    echo "merge option will  work in  future versions";
    merge_bin_and_muon.py --muon_data  $muon_data  --guide_umi_limit $guide_umi_limit --merge $merge
    """
}


process PerturbLoaderGeneration {
    debug true
    input:
    path (muon_data)
    path (gtf_in)
    val (distance_from_guide)
    val (in_trans)
    val (addgenes)
    output:
    path 'perturbdata.pkl', emit: perturb_piclke
    
   """ 
   PerturbLoader_generation.py --muon_data $muon_data --gtf_in $gtf_in  --distance_from_guide $distance_from_guide --in_trans $in_trans --add_gene_names $addgenes 

    """   
}


process runSceptre {
    debug true
    input:
    path (perturbloader_pickle)
    val (direction)
    output:
    path 'sceptre_out' , emit: sceptre_out_dir
    
   """ 
   mkdir sceptre_out  2>/dev/null
   mv $perturbloader_pickle sceptre_out
   cd sceptre_out
   runSceptre.py $perturbloader_pickle  $direction
   """   
}

process create_anndata_from_sceptre {
    debug true
    input:
    path (sceptre_results_dir)
    path (mudata_processed)
    
    output:
    path 'mudata_results.h5mu' , emit:  mudata_perturbation_results
    
   """ 
   sceptre_anndata_creation.py $sceptre_results_dir $mudata_processed
   """   


}




process preprocess_bar_multiseq {
    debug true
    input:
        path (MUON_DATA)
    output:
        path 'cell_barcode_capturing.csv',     emit: cell_barcode_to_multi
    script:
        if (params.RUN_MULTISEQ)
            """ 
                preprocessing_muon_multi.py $MUON_DATA
            """   
        if (params.RUN_MULTISEQ == false)

            """
                echo 'skipping preprocessing muon multi'
                touch cell_barcode_capturing.csv
            """
    
}

process  MultiSeq {
    debug true
    input:
        path (R1_MULTI)
        path (R2_MULTI)
        path (BARCODES_CELL_LIST_MULTI)
        path (BARCODES_MULTIBAR_LIST_MULTI)
        val (BAR_MULTI)
        val (UMI_MULTI)
        val (R2_MULTI_TAG)
        path (MUON_DATA)
        
    output:
        path 'final_class.csv',               emit: multi_class
        path 'final_class_cell_barcode.csv', emit: multi_class_barcode
        path  'bar_table.csv' ,         emit: bar_count
        path  'processed_mudata_guide_and_transcripts_multiseq_filtered.h5mu', emit: muon_with_multiseq

    script:
        BAR_MULTI_0 = BAR_MULTI[0]
        UMI_MULTI_0 = UMI_MULTI[0]
        R2_MULTI_TAG_0 = R2_MULTI_TAG[0]
        BAR_MULTI_1 =  BAR_MULTI[1]
        UMI_MULTI_1 =  UMI_MULTI[1]
        R2_MULTI_TAG_1 =  R2_MULTI_TAG[1]     
        if (params.RUN_MULTISEQ)
            
            """
                echo $R1_MULTI
                multiseq.py $R1_MULTI $R2_MULTI $BARCODES_CELL_LIST_MULTI $BARCODES_MULTIBAR_LIST_MULTI $BAR_MULTI_0 $BAR_MULTI_1 $UMI_MULTI_0 $UMI_MULTI_1 $R2_MULTI_TAG_0 $R2_MULTI_TAG_1 $MUON_DATA
            """ 
        
        if (params.RUN_MULTISEQ == false)
            """
            echo 'skiping muon'
            touch final_class.csv
            touch final_class_cell_barcode.csv
            touch bar_table.csv
            touch processed_mudata_guide_and_transcripts_multiseq_filtered.h5mu
            """
       
}






process select_final_muon{
    debug true
    input:
        path(SCRNA_GUIDE)
        path(SCRA_GUIDE_MULTI)
    output:
        path 'muon_selected_final.h5mu', emit:  final_muon
    script:
        if (params.RUN_MULTISEQ)   
            """
            cp $SCRA_GUIDE_MULTI muon_selected_final.h5mu
            """
        
        if (params.RUN_MULTISEQ == false)

            """
            cp $SCRA_GUIDE muon_selected_final.h5mu
            """

}
