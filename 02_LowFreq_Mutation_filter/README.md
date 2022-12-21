# LowFreq Mutation filter

## Usage
For counting the frequency of amino acid types at each position of the amino acid sequence and comparing the
differences in amino acid types and frequencies at each position between two sets of data.Finally, the sequences with low frequency mutations are filtered by the s2_filter_low_freq.py script

## Requirements
Python3.6+ is required

## Running
+ run`python b4_filter_translations.py test_AD`, output `test_AD-filter`<br>to Remove sequences with * in the protein translation fragment .The ID area and Env area are handled in the same way as the AD area.

+ run`python /p1_stat_seq.py [input_dir][start_pos][end_pos][output_dir]`<br>to count the amino acid type and frequency at each position of the sequence.
  This python script can be skipped.
  
+ run`python p2_compare.py [input_file1][input_file2]`to compare two sets of data.
  This python script can be skipped.
  
+ run`python s2_filter_low_freq.py [input_dir][start][end][min_freq][max_lows]`<br>to filtering sequences containinglow-frequency amino acids.

  for example:<br>
  run`python b4_filter_translations.py test_AD-filter 0 356 0.01 5`,output
  `test_AD-filter_filter_low_freq`.
  
   run`python b4_filter_translations.py test_ID-filter 0 356 0.01 7`,output
  `test_ID-filter_filter_low_freq`.

  The following parameters were set to ensure that the subsequent CD-HIT clustering and tree building clusters were in the range of 5-10 for AD and ID group. These parameters are used to keep only the major viral genotypes for subsequent analysis without changing the trend.
  
 1. For the AD group `TD004, TD005, TD008 (AD3_AD24 samples)`, the filtering parameter was set to `0.01 5`. These samples were done in Picbio Sequel I sequencing platform, which has a larger data volume. The other AD groupings were set with a filtering parameter of `0.01 10` and these samples were done on the Picbio Sequel II sequencing platform, which has a smaller data volume.
 
 2. For the ID group `TD004, TD005, TD006 (ID4_ID25 samples)`, the filtering parameter was set to `0.01 7`. These samples were done in Picbio Sequel I sequencing platform, which has a larger data volume. The other ID groupings were set with a filtering parameter of `0.01 20` and these samples were done on the Picbio Sequel II sequencing platform, which has a smaller data volume.

## Notes：
  + The AD group represents the protease/reverse region of HIV.
  + The ID group represents the integrase region of HIV.
  + The Env group represents the Env region of HIV.The Env region is not filtered for low-frequency sequences because the number of translated homologous protein sequences after extraction is lower compared to AD and ID, which does not affect the subsequent cluster clustering.
  
## Contact
Author:<br>
Name：Yang Liuke<br>
Email：yangliuke816@163.com<br>
Phone：+86-18761863703

Name：Youwei Zhu<br>
Email：ywzhu22@m.fudan.edu.com<br>
Phone：+86-19145632696

