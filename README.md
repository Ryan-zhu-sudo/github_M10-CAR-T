# github_M10 CAR-T
Picbio sequel Third-generation amplicon sequencing of HIV before and after M10 CAR-T treatment. Scripts were further subjected to genetic evolutionary analysis.

##Motivation
Chimeric antigen receptor T (CAR-T) cells, long proposed for HIV-1 treatment, has yet to achieve desirable therapeutic efficacy. In our study, the PCR products were sequenced on Third-generation PacBio sequel platform by Genewiz to evaluate the genetic change (Genewiz, Suzhou, China) in 8 HIV patients treated by M10 CAR-T immunotherapy and 2 controls only applied with HAART therapy. High-throughput sequencing of partial pol gene and env gene amplified from whole blood samples at baseline or final visit.The results demonstrated that several participants showed a genetic shift in the composition of viral reservoir populations after treatment. Additionally, in some participants, the remaining proviruses showed reduced overall diversity in the virus pol gene. The results also reveal the selective pressure mainly targeted the env gene. This comprehensive survey supported the potential of M10 CAR-T cells as a new, safe, and effective therapeutic option toward HIV-1 eradication in AIDS patients. 

##Description

![](https://github.com/Ryan-zhu-sudo/github_M10-CAR-T/blob/main/Workflow.png)

 + Scripts in `01_extract_homologous_seqs` was applied to extract the homologous sequences from AD/ID (protease/reverse) and Env region.
    The corresponding results of the script run are stored in the `05_Homologous_nucleic_sequence.zip` file on the figshare website about this project.

 + Scripts in `02_LowFreq_Mutation_filter` was applied to keep only the major viral genotypes for subsequent analysis without changing the trend.
   The corresponding results of the script run are stored in the `06_Homologous_protein_sequence.zip` file on the figshare website about this project.

 + Scripts in `03_Data_stats` was applied to count the change of HIV viral reservoir in Pol and Env region before and after M10 CART treatment.
	The corresponding results of the script run are stored in the `07_CDhit_cluster.zip` and `08_Data_stats.zip` file on the figshare website about this project.
	
 + `HIV(HXB2)_Ref-seq.fasta` is the complete HIV genome file.

 + protease/reverse_intergrate_Reference.fasta provides the sequence of AD/ID regions (protease/reverse and integrase) used in this project, as well as Env region.

 + The `M10 CART treatment sample grouping label.xlsx` records whether all data are grouped before or after M10 CART treatment. The homologous sequences were renamed according to the grouping after translation into protein, and the pre-treatment label is SHAI, and the post-treatment label is AFTER.The renaming results are stored in fishare website under 06_Homologous_protein_sequence of this project in `03_AD_protein-rename`,`03_ID_protein-rename`,`03_Env_protein-rename`.

 + `The Third-generation amplicon sequencing analysis report.pdf` file illustrates the library construction method used in this study and the processing of the raw sequencing data after off-boarding. 


##Notes：
+ The AD group represents the protease/reverse region of HIV.
+ The ID group represents the integrase region of HIV.
+ The Env group represents the Env region of HIV.
+ The original down-sequencing files ending in .fastq in this study are stored in NCBI Bioproject (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA903552). The intermediate processing files are located in the fishare platform (https://figshare.com/articles/online_resource/figshare_M10_CAR-T/21671042).
+ All data under the folders are for testing purposes, for reproduction of all results, the complete data is stored in can be obtained through the figshare platform.

## Contact
Author:
Name：Yang Liuke<br>
Email：yangliuke816@163.com<br>
Phone：+86-18761863703

Name：Youwei Zhu<br>
Email：ywzhu22@m.fudan.edu.com<br>
Phone：+86-19145632696
