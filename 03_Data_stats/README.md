# Data_stats

## Usage
Statistics for before and after clustering viral  reservoir change . When using this method of analysis, statistics can only be a sequence type, such as AD, ID and Env. you cannot put the two sequences into a directory.

## Requirements
Python3.6+ is required.

## Running
cd 	`test_AD`

+ run`for file in $(ls *.fas)` to change the same source protein sequence files ending with .fas in the `pass` folder of the `test_AD-filter_filter_low_freq` result file in the previous step to files ending with .txt. Then put these files into the `test_AD`, `test_ID`, `test_Env` folders under `03_Data_stats`.The ID area and Env area are handled in the same way as the AD area.


+ run`for i in $(ls *.txt);do cat ${i%%_*}_SHAI_AD.txt  ${i%%_*}_AFTER_AD.txt > ../sub_combine_filter/${i%%_*}_combine_AD.txt;done`to Get the sequences with high similarity before and after M10 CART treatment, and the post-treatment sequences belonging to the same cluster may be variants of the pre-treatment sequences.output `test_AD_combine`.The ID area are handled in the same way as the AD area.

+ run`for file in $(ls *.txt); do cd-hit -i $file  -o ../test_AD_combine_cluster/$file.0.95.txt -c 0.95 -n 5 -g 1 -d 0 -M 16000 -T 8; done`command to cluster homologous protein sequences into different Clusters,output `test_AD_combine`, the folder to include the files ending with .clstr. The ID area and Env area are handled in the same way as the AD area.

+ run `python stat_Pol.py ./test_AD_combine_cluster ./test_AD_combine` to get Statistical data for before and after clusters in AD/ID(protease/reverse) region.The ID area are handled in the same way as the AD area. This script is specifically designed to handle the AD/ID (protease/reverse) region.

+ run `python stat_Env.py ./test_AD_combine_cluster ./test_AD_combine` to get Statistical data for before and after clusters in env region. This script is specifically designed to handle the Env region.

## Notes：
+ The AD group represents the protease/reverse region of HIV.
+ The ID group represents the integrase region of HIV.
+ The Env group represents the Env region of HIV.
+ The Env group does not perform protein sequence combine before and after M10 CART treatment. Since Env sequence exactly clusters into two complete clusters before and after treatment without data combine shown in `test_Env_cluster`.Therefore, `stat_Env.py` was designed to separately count the change of virus reservoir in Env region before and after treatment.

## Contact
Author:<br>
Name：Yang Liuke<br>
Email：yangliuke816@163.com<br>
Phone：+86-18761863703

Name：Youwei Zhu<br>
Email：ywzhu22@m.fudan.edu.com<br>
Phone：+86-19145632696
