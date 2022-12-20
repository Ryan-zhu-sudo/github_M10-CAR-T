#!/usr/bin/env bash

fa_ad="../query/AD-query.fa"
fa_id="../query/ID-query.fa"

fdir="../03_SeqCount"

outdir="../output"
if [[ ! -d ${outdir} ]];then
  mkdir ${outdir}
fi

dirs=`ls ${fdir}`

for dr in ${dirs};do
  echo ${dr}

  blastdb=${fdir}/${dr}/blastdb/${dr}

  output_ad=${outdir}/${dr}.AD.blast
  output_id=${outdir}/${dr}.ID.blast

  blastn -query ${fa_ad} -db ${blastdb} -out ${output_ad} -task blastn-short -word_size 10 -max_target_seqs 100000 -outfmt 6 -num_threads 4

  blastn -query ${fa_id} -db ${blastdb} -out ${output_id} -task blastn-short -word_size 10 -max_target_seqs 100000 -outfmt 6 -num_threads 4

done

