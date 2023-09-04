#Calculate Points of Departure from BMDExpress2 - paths generated for bmd_gene.txt and bmd_path.txt from calculate_pods.R script. BMDExpress input directory added
# along with output directories for pathways and genes. PoD's calculated and csv files created with all PoD's of interest added.

calculate_PoDs_from_BMDExpress2(gene_bmd_file_path = "bmd_gene.txt",
                                path_bmd_file_path = "bmd_pathway.txt",
                                bmdExpress_input_dir = "BMD_input",
                                gene_csv_output_dir = "BMD_output_filtered/gene_bmds_filtered",
                                path_csv_output_dir = "BMD_output_filtered/pathway_bmds_filtered",
                                pod_csv_output_prefix = "pods",
                                pod_csv_output_dir = ".")
