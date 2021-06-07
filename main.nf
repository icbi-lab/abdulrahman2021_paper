#!/usr/bin/env nextflow

RES_DIR = params.resultsDir


process p01_process_data {
    def id = "01_process_counts"
    cpus = 16
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_scanpy.sif"
    cache 'lenient'
    publishDir "$RES_DIR/01_process_data", mode: params.publishDirMode

    input:
        // it would be better to include all input files explicitly.
        // However not that easy due to the fact that file identifier
        // is its parent folder.
        file 'data' from Channel.fromPath("data")
        file 'sample_sheet.csv' from Channel.fromPath("tables/vanderburg_01_samples.csv")
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")

    output:
        file "adata.h5ad" into process_data_adata
        file "${id}.html" into process_data_html

    """
    execute_notebook.sh ${id} ${task.cpus} notebook.Rmd \\
       "-r sample_sheet sample_sheet.csv -r output_file adata.h5ad -r data_dir data -r n_cpus ${task.cpus}"
    """
}


process p02_filter_data {
    def id = "02_filter_data"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_scanpy.sif"
    cpus = 16
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/jupytertools.py")
        file 'tables/*' from Channel.fromPath(
            "tables/{mitochondrial_genes,biomart,ribosomal_genes}.tsv"
        ).collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from process_data_adata

    output:
        file "adata.h5ad" into filter_data_adata_1, filter_data_adata_2
        file "${id}.html" into filter_data_html

    """
    execute_notebook.sh ${id} ${task.cpus} notebook.Rmd \\
       "-r input_file input_adata.h5ad -r output_file adata.h5ad -r table_dir tables"
    """
}

/**
 * Use the pre-computed doublets from the `tables` directory
 * that was generated using this process instead.
 * This process takes some time to run and is not numerically stable,
 * i.e. the result is slightly different every time.
 */

// process p02b_doublet_detection {
//     def id = "02b_doublet_detection"
//     conda "/home/sturm/.conda/envs/vanderburg_scanpy"
//     cpus = 16
//     clusterOptions '-V -S /bin/bash -l gpu=1 -q all.q'
//     publishDir "$RES_DIR/$id", mode: params.publishDirMode

//     input:
//         file "input_adata.h5ad" from filter_data_adata_1
//         file "model.json" from Channel.fromPath("tables/solo_model.json")

//     output:
//         file "out/is_doublet.npy" into doublet_detection_is_doublet
//         file "out/*.pdf"

//     """
//     export OPENBLAS_NUM_THREADS=${task.cpus}
//     export OMP_NUM_THREADS=${task.cpus}
//     export MKL_NUM_THREADS=${task.cpus}
//     export OMP_NUM_cpus=${task.cpus}
//     export MKL_NUM_cpus=${task.cpus}
//     export OPENBLAS_NUM_cpus=${task.cpus}
//     solo -o out -p model.json input_adata.h5ad
//     """
// }


process p03_normalize {
    def id = "03_normalize"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_scanpy.sif"
    cpus = 16
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file "is_doublet.npy" from Channel.fromPath("tables/is_doublet.npy")
        file 'lib/*' from Channel.fromPath("lib/{jupytertools,scio,scpp}.py").collect()
        file 'tables/*' from Channel.fromPath(
            "tables/{biomart.tsv,cell_cycle_regev.tsv,adata_pca.pkl.gz}"
        ).collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from filter_data_adata_2

    output:
        file "adata.h5ad" into correct_data_adata
        file "${id}.html" into correct_data_html

    """
    execute_notebook.sh ${id} ${task.cpus} notebook.Rmd \\
       "-r input_file input_adata.h5ad -r output_file adata.h5ad -r tables_dir tables -r doublet_file is_doublet.npy"
    """
}


process p04_annotate_cell_types {
    def id = "04_annotate_cell_types"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_scanpy.sif"
    cpus = 16
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/jupytertools.py")
        file 'tables/*' from Channel.fromPath(
            "tables/cell_type_markers.csv"
        ).collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from correct_data_adata

    output:
        file "adata.h5ad" into annotate_cell_types_adata
        file "${id}.html" into annotate_cell_types_html
        file "cell_types_per_sample.csv"

    """
    execute_notebook.sh ${id} ${task.cpus} notebook.Rmd \\
       "-r input_file input_adata.h5ad -r output_dir . -r table_dir tables"
    """

}


process p05_prepare_adata_t_nk {
    def id = "05_prepare_adata_nk_t"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_scanpy.sif"
    cpus 1
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/jupytertools.py")
        file 'tables/*' from Channel.fromPath(
            "tables/{cell_type_markers.csv,adata_pca*.pkl.gz}"
        ).collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from annotate_cell_types_adata

    output:
        file "adata.h5ad" into prepare_adata_t_nk, prepare_adata_t_nk_3, prepare_adata_t_nk_6
        file "${id}.html" into prepare_adata_t_nk_html
        file "adata_obs.tsv" into prepare_adata_t_nk_obs,
           prepare_adata_t_nk_obs_2
        file "norm_counts.tsv" into prepare_adata_t_nk_norm_counts
    """
    execute_notebook.sh ${id} ${task.cpus} notebook.Rmd \\
       "-r input_file input_adata.h5ad -r output_file adata.h5ad -r output_file_obs adata_obs.tsv -r output_file_norm_counts norm_counts.tsv -r table_dir tables -r cpus ${task.cpus} -r results_dir ."
    """
}

process p20_prepare_cluster_de_analysis {
    def id = "20_prepare_cluster_de_analysis"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_de_results.sif"
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'obs.tsv' from prepare_adata_t_nk_obs_2
        file 'counts.tsv' from prepare_adata_t_nk_norm_counts

    output:
        file "*.rda" into prepare_cluster_de_analysis_rda
        file "${id}.html" into prepare_cluster_de_analysis_html

    """
    reportsrender notebook.Rmd \
        ${id}.html \
        --cpus=${task.cpus} \
        --params="input_obs=obs.tsv \
                  input_counts=counts.tsv \
                  output_dir='.'"
    """
}

process p21_run_de_analysis_clusters {
    def id = "21_run_de_analysis_clusters"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_edger.sif"
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    cpus 6

    input:
        file input_data from prepare_cluster_de_analysis_rda.flatten()

    output:
        file "${input_data}.res.tsv" into run_de_analysis_clusters_results
        file "${input_data}.res.xlsx" into run_de_analysis_clusters_results_xlsx

    """
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
            MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus} \
            MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus} \
            MKL_THREADING_LAYER=GNU
    run_de.R ${input_data} ${input_data}.res.tsv \
        --cpus=${task.cpus} \
        --excel=${input_data}.res.xlsx
    """
}

process p22_cluster_de_analysis {
    def id = "22_cluster_de_analysis"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_de_results.sif"
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file "*" from run_de_analysis_clusters_results_xlsx.collect()
        file "*" from run_de_analysis_clusters_results.collect()

    output:
        file "${id}.html" into cluster_de_analysis_html
        file "*.zip" into cluster_de_analysis_zip

    """
    # use python, zip not available in container
    python -m zipfile -c ${id}.zip *.xlsx
    reportsrender notebook.Rmd \
        ${id}.html \
        --cpus=${task.cpus} \
        --params="de_dir='.'"
    """
}

process p60_tcr_analysis {
    def id = "60_tcr_analysis"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_scanpy.sif"
    cpus 42
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from prepare_adata_t_nk_3

    output:
        file "${id}.html" into tcr_analysis_html
        file "${id}.zip" into tcr_analysis_tsv

    """
    execute_notebook.sh ${id} ${task.cpus} notebook.Rmd \\
        "-r input_file input_adata.h5ad -r n_cpus ${task.cpus} -r output_dir . -r tcr_dir cellranger"
    python -m zipfile -c ${id}.zip *.tsv
    """
}

process p61_cluster_analysis {
    def id = "61_cluster_analysis"
    container "https://github.com/icbi-lab/abdulrahman2021_paper/releases/download/containers-0.1.0/vanderburg_scanpy.sif"
    cpus 1
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from prepare_adata_t_nk_6

    output:
        file "${id}.zip" into cluster_analysis_figures
        file "${id}.html" into cluster_analysis_html
    """
    execute_notebook.sh ${id} ${task.cpus} notebook.Rmd \\
        "-r input_file input_adata.h5ad"
    python -m zipfile -c ${id}.zip figures/*.pdf
    """
}




process deploy {
    publishDir "${params.deployDir}", mode: "copy"
    executor "local"

    input:
        file "input/*" from Channel.from().mix(
            process_data_html,
            filter_data_html,
            correct_data_html,
            annotate_cell_types_html,
            prepare_adata_t_nk_html,
            cluster_de_analysis_html,
            cluster_de_analysis_zip,
            tcr_analysis_html,
            tcr_analysis_tsv,
            cluster_analysis_html,
            cluster_analysis_figures,
        ).collect()

    output:
        file "*.html"
        file "*.zip"

    """
    cp input/*.{html,zip} .
    """
}

