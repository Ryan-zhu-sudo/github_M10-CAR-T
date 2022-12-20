#!/usr/bin/env bash

fdir="../03_SeqCount"

dirs=`ls ${fdir}`

for dr in ${dirs};do
  echo ${dr}
  db_dir=${fdir}/${dr}/blastdb
  fa_in=${fdir}/${dr}/${dr}.target.fa

  if [[ ! -d ${db_dir} ]];then
    nkdir ${db_dir}
  fi

  makeblastdb -in ${fa_in} -dbtype nucl -parse_seqids -out ${db_dir}/${dr}
done

