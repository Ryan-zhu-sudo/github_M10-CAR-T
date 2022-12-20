###########################################################################################
For counting the frequency of amino acid types at each position of the amino acid sequence 
and comparing the differences in amino acid types and frequencies at each position between
 two sets of data. Plot amino acid sequence logos. Finally, the sequences with low frequency 
mutations are filtered by the s2_filter_low_freq.py script
###########################################################################################

# Use
1. Run p1_stat_seq.py
    to count the amino acid type and frequency at each position of the sequence
    python . /p1_stat_seq.py [input_dir] [start_pos] [end_pos] [output_dir]
    input_dir: home directory containing amino acid sequences, with txt suffix, cannot contain extraneous files
    start_pos: the starting position to be compared, e.g. 0
    end_pos: the end position to be counted, cannot be larger than the length of the sequence
    Example: python . /p1_stat_seq.py . /test_before 11 100 . /output
    Output: . /output/test_before-seqstats.csv, and similarly for the directory test_after
    get . /output/test_after-seqstats.csv

2.Run p2_compare.py
    to compare two sets of data
    python . /p2_compare.py [input_file1] [input_file2]
    input_file1: the file output by the first step, such as the before group output file
    input_file2: the file output by the first step, such as the after group output file
    Example: python . /p2_compare.py . /output/test_before-seqstats.csv . /output/test_after-seqstats.csv
    Output: . /output/amino_compare.csv
    At this point, you have to analyze the file amino_compare.csv to determine where you need the final statistics and plotting, and also combine the structural domain area selection, for example
    eventually found that the base difference from position 20 to position 30 is relatively large, then 20 and 30 will be used as sequence selection parameters

3.Run p3_plot.py
    Plot according to the position determined in step 2
    python . /p3_plot.py [input_file] [start_pos] [end_pos]
    Example: python . /p3_plot.py . /output/test_before-seqstats.csv 20 30. Similarly for after group operations
    output pdf image file
    Rscript needs to exist in the environment variable, you need to install the motifStack package

4.s2_filter_low_freq.py
    Filtering sequences containing low-frequency amino acids
    Using.
    python s2_filter_low_freq.py [input_dir] [start] [end] [min_freq] [max_lows]
    input_dir: directory containing fasta files (.fas or .txt endings)
    start: Please set to 0
    end: integer, AD: 346, ID: 309, this series of scripts starts at 0, AD length 347, ID length 310
    min_freq: float, the minimum amino acid frequency for filtering, e.g. 0.01
    max_lows: integer, the number of amino acids containing low frequencies (<min_freq), anything greater than that will be filtered
    Note: This script will reference p1_stat_seq.py, you need to put them in the same directory

Author:
Name：Yang Liuke
Email：yangliuke816@163.com
Address：Institutes of Biomedical Sciences, Fudan University, Shanghai, China
Phone：+86-18761863703

Name：Youwei Zhu
Email：ywzhu22@m.fudan.edu.cn
Address：Shanghai Public Health Clinical Center & Institutes of Biomedical Sciences, Fudan University, Shanghai, China
Phone：+86-19145632696

