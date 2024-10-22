#!/usr/bin/bash
#参数为输出文件名
output=$1
echo '"sample","id","hmm_start","hmm_end","hmm_strand","contig_len","hmm_query_name","full_E-value","full_score","full_bias","best_E-value","best_score","best_bias","exp","reg","clu","ov","env","dom","rep","inc","target_description","hmm_seq","protein_seq","integraity","Name","Start","End","DR_Consensus","Repeat_ID","DR_Length","Spacers","Potential_Orientation","CRISPRDirection","Evidence_Level","Conservation_DRs","Conservation_Spacers"' > ${output};for i in `find . -type f -name "all_summary.csv"`;do j=${i%_out*};j=${j##*/};j=${j%.*};sed '1d' $i | awk -v jj="$j" -v q='"' -v OFS="," '{print q jj q,$0}' >> ${output};done
