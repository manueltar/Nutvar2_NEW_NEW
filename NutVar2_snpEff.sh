# Define time and stamp it in every result directory

echo `ulimit -a`

Softwaredir=$1
# ~/Downloads/nutvar2-master
vcfinput=$2
# USER INPUT p.e. vep_example.vcf
output=$3

# Get date and time to add to the results folder

now="$(date +'%d_%m_%Y')"
now2="$(date +"%T")"
now4="_"
now3=$now$now4$now2
#~ foldername="$(date +"%T")"
foldername="$now3"
echo "$now3"


snpEFFdir=SOFTWARE/snpEff
VEPdir=SOFTWARE/ensembl-tools-release-75/scripts/variant_effect_predictor

bindir1=bin/shared
bindir2=bin/snpEff
bindir3=bin/VEP
datadir1=data/intermediate
datadir2=data/build_tables
datadir3=data/external

datadir4=$output

#~ echo "$datadir4"

#~ exit

mypwd=$(pwd)

#~ exit;

	# Open file and eliminate header lines. Here we should check in the input vcf has the appropriate fields. ISSUE.

cat $2 | perl -ne 'chomp;unless($_=~/^##/){$_=~s/^[Cc]hr//;print "$_\n";}' > ${datadir1}/vcfinput.vcf

	# Run the minimal representation script
	
	echo "1-Transforming_vcf_to_minimal_representation"

	perl $1/${bindir1}/2_Script_minimal_representation_vcf_7.0.pl ${datadir1}/vcfinput.vcf $1/${datadir1}/vcfinput_mr.vcf

	# Ask the user whether she wants to run SnpEff, VEP or both #ISSUE

	# NOTICE GRCh37.75 database for snpEff and the ENSEMBL version of genome 37.75 are installed through the install script


	# Run SnpEff

echo "2-Runing_snpEff"

java -Xmx4g -jar $1/${snpEFFdir}/snpEff.jar eff -c $1/${snpEFFdir}/snpEff.config -v GRCh37.75 -lof -csvStats -nextProt -sequenceOntology $1/${datadir1}/vcfinput_mr.vcf > $1/${datadir1}/vcfinput_mr_eff.vcf

mv snpEff_genes.txt $1/${datadir1}/snpEff_genes.txt
mv snpEff_summary.csv $1/${datadir1}/snpEff_summary.csv

#Parsing the results of snpEff

# ISSUE There are two more scripts in this folder /${bindir2}/ 24 and 25 ---> Erase them?

echo "3-Parsing_snpEff_results_and_calculating_NMD"

perl $1/${bindir2}/24_snpEff_parser_def_minus_heather_3.0.pl ${datadir1}/vcfinput_mr_eff.vcf ${datadir2}/gtf_INTRONS_Plus_NMD_threshold.txt ${datadir1}/out_snpeff_parsed.txt ${datadir1}/snpeff_NMD.txt

echo "snpEff_done.Runing_NUTVAR2"

# Here there is a chance to MPI as script 25 is the longest the others can take place untill 27 while 25 is executing

echo "4-Calculating_downstream_pTCs_in _frameshift_variants"

perl $1/${bindir1}/25_Downstream_frameshift_API_independent_5.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/CDS_genomic_coordinates_full_compresed.txt ${datadir1}/snpeff_derived_PTCS_API_independent.txt

echo "5-Calculating_percentage_of_sequence_affected"

perl $1/${bindir1}/26_NEW_EXTRA_key_%_sequence_2.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir1}/snpeff_percentage.txt

	
# Here
	
	#~ DOCUMENT!
	
	# Solved in MPZ we lost track of 
	
	#~ retained_intron	exon	ENST00000488271
	#~ processed_transcript	ENST00000526189
	#~ nonsense_mediated_decay	ENST00000463290
	#~ 
	#~ because they are non protein coding transcripts an NutVar is for protein coding transcripts!
	
	# Document; 
	
	#this is the script of the warnings in the splicin in last component and splice in monoexonic
	
	
echo "6-Calculating_percentage_of_domains_and_sites_affected I"
 
perl $1/${bindir1}/32_key_PROTEINS_8.0_GLOBAL_3.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir2}/ALL_ISOFORMS_PROTEIN_table_full.txt ${datadir1}/snpeff_detailed_ProtAndSite_Pre_step.txt

#this is the script of the warnings in the splicin in last component and splice in monoexonic

sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 ${datadir1}/snpeff_detailed_ProtAndSite_Pre_step.txt > ${datadir1}/snpeff_detailed_ProtAndSite_Pre_step_ordered.txt

perl $1/${bindir1}/32_key_PROTEINS_8.0_GLOBAL_ParteII.pl ${datadir1}/snpeff_detailed_ProtAndSite_Pre_step_ordered.txt ${datadir1}/snpeff_detailed_ProtAndSite_Post_step.txt

echo "7-Calculating_percentage_of_domains_and_sites_affected_II"

perl $1/${bindir1}/32_key_PROTEINS_8.0_GLOBAL_3.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir2}/ALL_ISOFORMS_DOMAIN_table_full.txt ${datadir1}/snpeff_DOMAINS_Pre_step.txt

sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 ${datadir1}/snpeff_DOMAINS_Pre_step.txt > ${datadir1}/snpeff_DOMAINS_Pre_step_ordered.txt

perl $1/${bindir1}/32_key_PROTEINS_8.0_GLOBAL_ParteII.pl ${datadir1}/snpeff_DOMAINS_Pre_step_ordered.txt ${datadir1}/snpeff_DOMAINS_Post_step.txt

echo "8-Calculating_NMD_for_frameshifts_with_derived_pTCs"

perl $1/${bindir1}/27_key_NMD_5.0_DERIVED_STOPS_3.0.pl ${datadir2}/gtf_INTRONS_Plus_NMD_threshold.txt ${datadir1}/snpeff_derived_PTCS_API_independent.txt ${datadir1}/out_snpeff_parsed.txt ${datadir1}/snpeff_derived_NMD.txt
			
echo "9-Creating_first_matrix"

perl $1/${bindir1}/38_global_feature_table_1_4_paralell.pl ${datadir1}/snpeff_NMD.txt ${datadir1}/snpeff_derived_NMD.txt ${datadir1}/snpeff_detailed_ProtAndSite_Post_step.txt ${datadir1}/snpeff_DOMAINS_Post_step.txt ${datadir1}/snpeff_percentage.txt ${datadir1}/snpeff_first_table.txt

echo "10-Creating_matrix _def.Adding_Pincipal_Isoform_information"

perl $1/${bindir1}/40_tabla_PEJMAN_16_def_5.0.pl ${datadir1}/snpeff_first_table.txt ${datadir1}/snpeff_NMD.txt ${datadir1}/snpeff_derived_NMD.txt ${datadir2}/gtf_output_ENST.txt ${datadir2}/gtf_output_ENSG.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir3}/appris_principal_isoform_gencode_19_15_10_2014.txt ${datadir3}/Pervasive.txt ${datadir1}/Matrix_snpeff.txt

## Document NMD is frameshift derived always in frameshifts

echo "11-Creating_first_matrix_I_CCDS"

perl $1/${bindir1}/41_CCDS_collapser_3.0.pl ${datadir2}/gtf_tabladef_sorted_by_SYMBOL.txt ${datadir1}/snpeff_NMD.txt ${datadir1}/snpeff_NMD_CCDS.txt ${datadir1}/snpeff_derived_NMD.txt ${datadir1}/snpeff_derived_NMD_CCDS.txt ${datadir2}/gtf_output_ENSG.txt ${datadir2}/gtf_output_ENSG_CCDS.txt ${datadir2}/gtf_output_ENST.txt ${datadir2}/gtf_output_ENST_CCDS.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir2}/ENST_table_full_condensed_CCDS.txt ${datadir3}/appris_principal_isoform_gencode_19_15_10_2014.txt ${datadir3}/appris_principal_isoform_gencode_19_15_10_2014_CCDS.txt ${datadir3}/Pervasive.txt ${datadir3}/Pervasive_CCDS.txt ${datadir1}/snpeff_first_table.txt ${datadir1}/snpeff_first_table_CCDS.txt

echo "12-Creating_matrix_def_CCDS"

perl $1/${bindir1}/42_tabla_PEJMAN_15.0_version_paralel_5.0.pl ${datadir1}/snpeff_first_table_CCDS.txt ${datadir1}/snpeff_NMD_CCDS.txt ${datadir1}/snpeff_derived_NMD_CCDS.txt ${datadir2}/gtf_output_ENST_CCDS.txt ${datadir2}/gtf_output_ENSG_CCDS.txt ${datadir2}/ENST_table_full_condensed_CCDS.txt ${datadir3}/appris_principal_isoform_gencode_19_15_10_2014_CCDS.txt ${datadir3}/Pervasive_CCDS.txt ${datadir1}/Matrix_snpeff_CCDS.txt

echo "Printing_results_in_$datadir4"

mkdir -p ${datadir4}/"$foldername"

echo "13-Adding_Gene_based_features"

perl $1/${bindir1}/53BIS_Fuse_Matrix\&Gene_based.pl ${datadir1}/Matrix_snpeff.txt ${datadir3}/pRDG2.txt ${datadir3}/Genes_AllInnateImmunity.txt ${datadir3}/Genes_Antiviral.txt ${datadir3}/Genes_ISGs.txt ${datadir3}/Genes_OMIMrecessive.txt ${datadir3}/RVIS2.txt ${datadir4}/"$foldername"/Matrix_snpeff_added_gene_based_scores.txt

perl $1/${bindir1}/53BIS_Fuse_Matrix\&Gene_based.pl ${datadir1}/Matrix_snpeff_CCDS.txt ${datadir3}/pRDG2.txt ${datadir3}/Genes_AllInnateImmunity.txt ${datadir3}/Genes_Antiviral.txt ${datadir3}/Genes_ISGs.txt ${datadir3}/Genes_OMIMrecessive.txt ${datadir3}/RVIS2.txt ${datadir4}/"$foldername"/Matrix_snpeff_CCDS_added_gene_based_scores.txt

echo "snpEff_data_matrix_generated"


