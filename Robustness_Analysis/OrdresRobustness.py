# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 12:04:54 2022

@author: benja
"""

from os import listdir
import os
SamplesList = listdir("/home/benja/robustness-master_REV/TSV")

for sample in SamplesList:
    sampleName = sample.split('_')[1]
    individual = sampleName[0:3]
    timePoint = sampleName[3:5]
    ordre1 = 'src/generate_perturbations.R'+' '+ '/home/benja/robustness-master_REV/TSV/16S_'+str(individual)+str(timePoint)+'.tsv'+' '+ str(individual)+str(timePoint)+'_perturbed_taxonomic_compositions'+' '+ str(individual)+str(timePoint)+'_perturbed_function_compositions'+' '+str(individual)+str(timePoint)+'_Unifrac'+ ' --genome_content_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/KO_predicted.tsv.gz --copy_number_norm_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/marker_predicted_and_nsti.tsv.gz --tree /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/out.tre'
    print ordre1
    os.system(ordre1)
    ordre2 = 'src/summarize_perturbations.R'+' '+ str(individual)+ str(timePoint)+'_perturbed_taxonomic_compositions'+' '+ str(individual)+str(timePoint)+'_perturbed_function_compositions'+' '+str(individual)+str(timePoint)+'_Unifrac'+' '+ str(individual)+str(timePoint)+'_function_distribution_features' +' '+ str(individual)+str(timePoint)+'_perturbation_distances'+' '+str(individual)+str(timePoint)+'_function_specific_differences'+' '+'--genome_content_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/KO_predicted.tsv.gz'
    print ordre2    
    os.system(ordre2)
    ordre3 = 'src/calculate_robustness_factors.R'+' '+ str(individual)+ str(timePoint)+'_perturbation_distances'+' '+ str(individual)+str(timePoint)+'_robustness_factors'
    print ordre3    
    os.system(ordre3)
    ordre4 = 'src/calculate_function_specific_robustness_factors.R'+' '+ str(individual)+ str(timePoint)+'_function_specific_differences'+' '+ str(individual)+ str(timePoint)+ '_function_specific_robustness_factors'
    print ordre4    
    os.system(ordre4)


#src/generate_perturbations.R /home/benja/robustness-master_REV/TSV/16S_A01T2.tsv A1T2_perturbed_taxonomic_compositions A1T2_perturbed_function_compositions A1T2_Unifrac --genome_content_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/KO_predicted.tsv.gz --copy_number_norm_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/marker_predicted_and_nsti.tsv.gz --tree /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/out.tre
#src/generate_perturbations.R /home/benja/robustness-master_REV/TSV/16S_A01T1.tsv A01T1_perturbed_taxonomic_compositions A01T1_perturbed_function_compositions A01T1_Unifrac--genome_content_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/KO_predicted.tsv.gz --copy_number_norm_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/marker_predicted_and_nsti.tsv.gz --tree /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/out.tre


#src/summarize_perturbations.R A1T2_perturbed_taxonomic_compositions A1T2_perturbed_function_compositions A1T2_Unifrac    A1T2_function_distribution_features A1T2_perturbation_distances A1T2_function_specific_differences --genome_content_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/KO_predicted.tsv.gz
#src/summarize_perturbations.R A01T1_perturbed_taxonomic_compositions A01T1_perturbed_function_compositions A01T1_Unifrac A01T1_function_distribution_features A01T1_perturbation_distances A01T1_function_specific_differences --genome_content_table /home/benja/robustness-master_REV/PICRUSt2_4robustness/picrust2_results/KO_predicted.tsv.gz


#src/calculate_robustness_factors.R A1T2_perturbation_distances A01T2_robustness_factors
#src/calculate_robustness_factors.R A01T1_perturbation_distances A01T1_robustness_factors'

#src/calculate_function_specific_robustness_factors.R A1T2_function_specific_differences A1T2_function_specific_robustness_factors
#src/calculate_function_specific_robustness_factors.R A01T1_function_specific_differences A01T1_function_specific_robustness_factors