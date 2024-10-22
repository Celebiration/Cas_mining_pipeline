#!/usr/bin/bash

Usage()
{
    echo -e "本程序用于提取单个per_seq.txt文件的结果\nUsage: $0 -f contig_file -i per_seq_file -c cut_off=50"
}
cut_off=50
lines1='【------------------------------'
lines2='------------------------------】'
nametype="contig_prodigal"
while getopts ':f:i:c:r:' OPT; do
    case $OPT in
        f) contig="$OPTARG";;
        i) per_seq="$OPTARG";;
        c) cut_off="$OPTARG";;
        r) nametype="$OPTARG";;
        *) Usage; exit 1;;
    esac
done
if [ -z $contig ];then Usage; exit 1; fi
if [ -z $per_seq ];then Usage; exit 1; fi

outdir=${per_seq%/*}
bb=${per_seq##*/}
b=${bb::-12}

grep -P "^[^#]" $per_seq > ${outdir}/tmp
echo -n > ${outdir}/${b}_cas12i_tmp.txt
echo -n > ${outdir}/hitted_contigs.fa
echo "${lines1}开始写入文件：${outdir}/${b}_cas12i_tmp.txt，${outdir}/hitted_contigs.fa${lines2}"
while read line
do
    l=${line%%\#*}
    r=${line#*\#}
    hmm_start=`echo ${r}|awk -v FS="[ #]*" '{print $1}'`
    hmm_end=`echo ${r}|awk -v FS="[ #]*" '{print $2}'`
    hmm_strand=`echo ${r}|awk -v FS="[ #]*" '{print $3}'`
    line1="${l}\"#${r}\""
    i=${line1%% *}
    j=${line1#* }
    ii=${i%_*}
    ##JGI
    if [ $nametype != "contig_prodigal" ];then
    ii=${ii%_*}
    hmm_start=`echo ${r}|awk -v FS="[ #]*" '{print $3}'`
    hmm_end=`echo ${r}|awk -v FS="[ #]*" '{print $4}'`
    hmm_strand=`echo ${r}|awk -v FS="[ #]*" '{print $5}'`
    fi
    #没有contig_lengths_dir时：
    contig_seq=`seqkit grep -I -w 0 -p $ii ${contig}|tail -1`
    protein_seq=`seqkit grep -I -w 0 -p $i ${outdir}/${b}.fa|tail -1`
    echo "$i $hmm_start $hmm_end $hmm_strand ${#contig_seq} $j $contig_seq $protein_seq" >> ${outdir}/${b}_cas12i_tmp.txt
    echo ">${i}" >> ${outdir}/hitted_contigs.fa
    echo $contig_seq >> ${outdir}/hitted_contigs.fa
done < ${outdir}/tmp
seqkit rmdup -n ${outdir}/hitted_contigs.fa -o ${outdir}/hitted_contigs.uniq.fa
echo "完成。"
rm ${outdir}/tmp