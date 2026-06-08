#!/bin/bash
# 06_download_atac.sh
# Downloads wild-type ATAC-seq peak files from DANIO-CODE

set -e

echo "=== Downloading ATAC-seq peaks ==="

cd peaks

# Zhang Lab - Early stages (wild-type)
echo "Downloading Zhang Lab samples..."
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007335BS/ATAC-seq_Zhang_Lab_0001AS.DCD003127SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007338BS/ATAC-seq_Zhang_Lab_0001AS.DCD003152SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007338BS/ATAC-seq_Zhang_Lab_0001AS.DCD003157SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007337BS/ATAC-seq_Zhang_Lab_0001AS.DCD003146SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007336BS/ATAC-seq_Zhang_Lab_0001AS.DCD003131SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007336BS/ATAC-seq_Zhang_Lab_0001AS.DCD003141SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007334BS/ATAC-seq_Zhang_Lab_0001AS.DCD003126SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007334BS/ATAC-seq_Zhang_Lab_0001AS.DCD003151SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz

# Skarmeta Lab - Multiple stages (wild-type)
echo "Downloading Skarmeta Lab samples..."
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007265BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003082SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007265BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003084SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007266BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003087SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007266BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003090SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007267BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003086SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007267BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003089SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007268BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003092SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007268BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003096SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007269BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003079SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007269BS/ATAC-seq_Skarmeta_Lab_0002AS.DCD003088SQ.danRer11.USERmickael.dong.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD000152BS/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.danRer11.USERpanosfirbas.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz
curl -O https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD000152BS/ATAC-seq_Skarmeta_Lab_0001AS.DCD000395SQ.danRer11.USERpanosfirbas.R1.merged.filt.tn5.pval0.01.300K.narrowPeak.gz

echo ""
echo "=== Download complete ==="
ls -lh *.gz | wc -l
echo "files downloaded"
ls -lh *.gz
