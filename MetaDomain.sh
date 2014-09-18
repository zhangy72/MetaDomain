#!/bin/bash
# Copyright (c) 2011 Yuan Zhang, Yanni Sun.
# You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
# pipeline of MetaDomain. 

# input: 
# -m: hmm file;
# -f: fasta file;
# -g: gamma (alignment score threshold rate), default: 0.6;
# -n: aligned read number threshold, default: 20;
# -c: domain coverage threshold, default: 0.3.

# output:
# -o: domain expression result in the following format:
# [domain name] [domain length] [aligned read number] [covering rate]
# [transcribed/nontranscribed].

usage() {
  echo "MetaDomain.sh -m <HMMER3 HMM file> -f <fasta file> -o <output folder> [options] 
  Options:
    -h:  show this message
    -g:  gamma (alignment score threshold rate, in the range of [0,1], default: 0.6)
    -n:  aligned read number threshold (default: 20)
    -c:  domain coverage threshold (default: 0.3)" >&2
}

hmm=
fasta=
gamma=0.6
num=20
coverage=0.3
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while getopts "hm:f:g:n:c:o:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    m)
      hmm=$OPTARG
      ;;
    f)
      fasta=$OPTARG
      ;;
    g)
      gamma=$OPTARG
      ;;
    n)
      num=$OPTARG
      ;;
    c)
      coverage=$OPTARG
      ;;
    o)
      out=$OPTARG     
      ;;
    esac
done

if [ "$hmm" == "" ];then
  echo "Please specify the hmm file."
  usage
  exit
fi

if [ "$fasta" == "" ];then
  echo "Please specify the input fasta file."
  usage
  exit
fi

if [ "$out" == "" ];then
  echo "Please specify the output folder."
  usage
  exit
fi


# generate a temporary folder.
if [ ! -d $out/ ];then
  mkdir $out/
else
  echo 'Output folder exists. Please specify another output folder.'
  exit
fi
tmp="$(cd $out && pwd)"
base_fasta=`echo $fasta | awk '{split($1,a,"/"); print a[length(a)]}'`

python $DIR/parse_hmm_files.py $hmm $tmp/HMMs
cat /dev/null >$out/${base_fasta}.metadomain

ls $tmp/HMMs | while read line
do
  hmm=$tmp/HMMs/$line
  hmm_acc=`echo $line | awk '{print substr($1,1,7)}'`
  hmm_length=`grep "^LENG" $tmp/HMMs/$line | awk '{print $2}'`
  echo "${hmm_acc} ${hmm_length}" >$tmp/${hmm_acc}_length.list
  $DIR/hmmscore $hmm $fasta $gamma 0 >$tmp/${fasta_base}_${hmm_acc}.hmmscore
  $DIR/generate_domain_expression $tmp/${hmm_acc}_length.list $tmp/${fasta_base}_${hmm_acc}.hmmscore $num $coverage $tmp/temp.out.metadomain
  cat $tmp/temp.out.metadomain >>$out/${base_fasta}.metadomain
  rm $tmp/temp.out.metadomain
  rm $tmp/${fasta_base}_${hmm_acc}.hmmscore
  rm $tmp/${hmm_acc}_length.list
done
rm -r $tmp/HMMs
