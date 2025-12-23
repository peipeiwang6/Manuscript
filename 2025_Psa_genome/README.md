##########################################################################################
#                        Scripts for allelic group classification,                       #
#     and Ka and Ks calculation for syntenic and reciprocal best matched gene pairs      #
#                                Peipei Wang, 20250231                                   #
##########################################################################################

# parse the blast results and get the tandem duplicates
00_parse_Psa_allbyall_blast_results.py

# parse the original McScanX and tandem results
01_parse_McScanX_and_tandem_information.py

# get allelic genes in syntenic, inversion, translocation and duplication regions, based on Syri information
02_get_genes_in_syntenic_inversion_translocation_regions_based_on_syri.py

# infer allelic genes based on Syri and GFF information first, than fill out missing genes using McScanX information
03_get_allelic_genes_based_on_GFF_Syri_McScanX_Tandem.py

# get reciprocal best blast hits for genes between set3 and set4 chromosomes
04_get_reciprocal_blast_hits_for_genes_on_setA_chrs.py

# get reciprocal blast hits for type 5 (conserved) and type 2 (with heterozyous loss) allelic groups
05_get_reciprocal_blast_hits_for_conserved_and_biased_loss.py

# calculate Ka and Ks between reciprocal best matched gene pairs between set3 and set4, and within set3 and set4
06_KaKs_for_RBM.py

# convert protein alignment to CDS alignment, this script is used when calculating the Ka and Ks values
07_convert_pep_alignment_to_CDS_alignment.py

# calculate Ka and Ks for syntenic and non-syntenic gene pairs
08_calculate_syntenic_and_non-syntenic_gene_pairs_Ks.py

# draw the Ks distribution for syntenic and RBM gene pairs
09_compare_Ks_between_allelic_syntenic_blocks.r

# draw the locations of allelic groups on chromosomes
10_draw_the_allele_location_for_targe_genes.r