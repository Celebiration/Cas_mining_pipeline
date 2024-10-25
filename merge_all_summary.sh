#!/usr/bin/bash
Usage()
{
    echo -e "本程序用于合并snakemake的结果。\nUsage: $0 -o outputfile"
}
outputfile=merged_all_summary.csv
while getopts ':o:' OPT; do
    case $OPT in
        o) outputfile="$OPTARG";;
        *) Usage; exit 1;;
    esac
done
echo '"file","id","hmm_start","hmm_end","hmm_strand","aa_length","contig_len","hmm_query_name","full_E-value","full_score","full_bias","best_E-value","best_score","best_bias","exp","reg","clu","ov","env","dom","rep","inc","target_description","hmm_seq","protein_seq","integraity","Name","Start","End","DR_Consensus","Repeat_ID","DR_Length","Spacers","Potential_Orientation","CRISPRDirection","Evidence_Level","Conservation_DRs","Conservation_Spacers"' > ${outputfile}
#for i in `ls -d *_out`;do file=${i::-4};if [ -s ${i}/all_summary.csv  ];then sed '1d' ${i}/all_summary.csv | awk -v ff="$file" -v OFS="," '{print "\""ff"\"",$0}'>> ${outputfile};fi;done
for i in `find . -type f -name "all_summary.csv"`;do ii=${i%/*};iii=${ii##*/};iii=${ii::-4};cat $i|sed '1d'|awk -v ff="$iii" -v OFS="," '{print "\""ff"\"",$0}' >> ${outputfile};done
