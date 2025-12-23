###########################################################
#         get the allelic gene groups for Psa             #
#             Peipei Wang, 20250231                       #
###########################################################

# parse the blast results and get the tandem duplicates
00_parse_Psa_allbyall_blast_results.py

# parse the original McScanX and tandem results
01_parse_McScanX_and_tandem_information.py

# get allelic genes in syntenic, inversion, translocation and duplication regions, based on Syri information
02_get_genes_in_syntenic_inversion_translocation_regions_based_on_syri.py

# infer allelic genes based on Syri and GFF information first, than fill out missing genes using McScanX information
03_get_allelic_genes_based_on_GFF_Syri_McScanX_Tandem.py

# get reciprocal blast hits for genes between set3 and set4 chromosomes
04_get_reciprocal_blast_hits_for_genes_on_setA_chrs.py

# get reciprocal blast hits for type 5 (conserved) and type 2 (with heterozyous loss) allelic groups
05_get_reciprocal_blast_hits_for_conserved_and_biased_loss.py


# calculate KaKs between genes in allelic collinear block and other blocks
04_compare_KaKs_between_allelic_syntenic_blocks.r
04_KaKs_for_other_syntenic_blocks_jobs.py
13_convert_pep_alignment_to_CDS_alignment.py
14_compare_KaKs_and_Exp_dist_between_allelic_syntenic_blocks.r

# draw the locations of allelic groups on chromosomes
15_draw_the_allele_location_for_targe_genes.r