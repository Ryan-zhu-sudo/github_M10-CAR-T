# Extract homologous seqs

## Usage

Used for extracting target homologous sequences from the raw sequenceing data.
The target homologous sequences were marked in `HIV-genome.docx`, the first region was tagged as `AD` and the second region was tagged as `ID`. Raw sequenceing data was in the package `03_SeqCount.zip`.

## Requirements

Python3.6+ is required, other packages need to be installed were listed in `requirements.txt`.<br>
NCBI-BLAST need to be installed.<br>
EMBOSS:transeq need to be installed.

## Running

`cd ./scripts`

+ run `b0_make_blastdb.sh` 
Create a blastn database for every fasta in `03_SeqCount`, such as `03_SeqCount/AD24/blastdb`.

+ run `b1_create_query.py`  
Generate query fasta for `AD` and `ID`, `output `query/AD-query.fa`` and `query/ID-query.fa`.

+ run `b2_run_blast.sh`  
`query/AD-query.fa` and `query/ID-query.fa` will blast to each db in `03_SeqCount`, 
result files saved in `./output`.

+ run `b3_extract_seqs.py`  
It will parses the blast results and filters sequences, outputed fasta files saved in `result`.

`cd ./result`

+ run`for i in $(ls *.fa); do transeq ./${i%.*}.fa ../protein/${i%.*}.txt -frame=6;done`  
It will yield six homologous protein sequences translated from both the forward and reverse strands of each homologous sequence. Subsequently, the translated proteins are grouped according to the AD, ID, and Env regions, and the test files are placed in the `test_AD`,`test_ID`,`test_Env` folders.

## Notes：
+ The AD group represents the protease/reverse region of HIV.
+ The ID group represents the integrase region of HIV.
+ The Env group represents the integrase region of HIV.

## Contact
Author:
Name：Yang Liuke<br>
Email：yangliuke816@163.com<br>
Phone：+86-18761863703

Name：Youwei Zhu<br>
Email：ywzhu22@m.fudan.edu.com<br>
Phone：+86-19145632696
