#!/usr/bin/env python3
"""
06b_normalize_scores.py

Quantile-normalize ATAC-seq signal scores across samples to make them comparable
before aggregating into consensus peaks. Raw narrowPeak scores span ~85–27,000
and are not comparable across samples/labs.

Outputs normalized score BED files and provides metadata for replicate collapsing.
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def quantile_normalize(scores_array):
    """
    Quantile normalize a 1D array of scores.
    Returns normalized scores (0–1 range).
    """
    if len(scores_array) == 0:
        return scores_array
    
    # Sort scores and get ranks
    sorted_indices = np.argsort(scores_array)
    ranks = np.empty_like(sorted_indices)
    ranks[sorted_indices] = np.arange(len(scores_array))
    
    # Convert ranks to quantiles (0–1)
    normalized = ranks / (len(scores_array) - 1) if len(scores_array) > 1 else np.ones_like(ranks)
    return normalized

def rank_scores(scores_array):
    """
    Rank scores (1 = lowest, N = highest).
    Returns ranks normalized to 0–1 range.
    """
    if len(scores_array) == 0:
        return scores_array
    
    sorted_indices = np.argsort(scores_array)
    ranks = np.empty_like(sorted_indices, dtype=float)
    ranks[sorted_indices] = np.arange(1, len(scores_array) + 1)
    
    # Normalize to 0–1
    normalized = (ranks - 1) / (len(scores_array) - 1) if len(scores_array) > 1 else np.ones_like(ranks)
    return normalized

def normalize_sample_scores(narrowpeak_file, method="quantile"):
    """
    Load narrowPeak file, normalize col-5 (score) using specified method,
    return normalized BED.
    """
    try:
        df = pd.read_csv(
            narrowpeak_file,
            sep='\t',
            header=None,
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            dtype={0: str, 1: int, 2: int, 3: str, 4: int, 5: str, 
                   6: float, 7: float, 8: float, 9: int}
        )
        df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                      'signal_val', 'pvalue', 'qvalue', 'summit']
    except Exception as e:
        print(f"ERROR reading {narrowpeak_file}: {e}")
        return None
    
    if len(df) == 0:
        return df
    
    # Extract and normalize scores
    scores = df['score'].values
    
    if method == "quantile":
        normalized_scores = quantile_normalize(scores)
    elif method == "rank":
        normalized_scores = rank_scores(scores)
    else:  # fallback: max normalization
        normalized_scores = scores / scores.max() if scores.max() > 0 else scores
    
    # Replace score column with normalized values
    df['score'] = normalized_scores
    
    return df

def main():
    os.chdir('peaks')
    
    # Load configuration
    sys.path.insert(0, '..')
    # We'll read config via subprocess since it's bash
    import subprocess
    config_cmd = "source ../config.sh && echo $SCORE_NORMALIZATION_METHOD"
    result = subprocess.run(config_cmd, shell=True, capture_output=True, text=True, executable='/bin/bash')
    normalize_method = result.stdout.strip() or "quantile"
    
    print(f"=== Normalizing ATAC-seq Scores (method: {normalize_method}) ===")
    
    # Find all narrowPeak files
    narrowpeak_files = sorted(Path('.').glob('*.narrowPeak'))
    
    if not narrowpeak_files:
        print("ERROR: No .narrowPeak files found in current directory")
        sys.exit(1)
    
    print(f"Found {len(narrowpeak_files)} narrowPeak files")
    
    # Create normalized output directory
    os.makedirs('normalized_scores', exist_ok=True)
    
    # Normalize each file
    print("\nNormalizing individual sample scores...")
    normalization_log = []
    
    for narrowpeak_file in narrowpeak_files:
        sample_name = narrowpeak_file.stem
        print(f"  Processing {sample_name}...")
        
        df = normalize_sample_scores(str(narrowpeak_file), method=normalize_method)
        
        if df is not None and len(df) > 0:
            # Save normalized BED (chrom, start, end, name, normalized_score, strand, summit)
            output_file = f'normalized_scores/{sample_name}.normalized.bed'
            df[['chrom', 'start', 'end', 'name', 'score', 'strand', 'summit']].to_csv(
                output_file,
                sep='\t',
                header=False,
                index=False
            )
            
            # Log statistics
            normalization_log.append({
                'sample': sample_name,
                'num_peaks': len(df),
                'original_score_min': narrowpeak_file,  # placeholder
                'normalized_score_mean': f"{df['score'].mean():.4f}",
                'normalized_score_median': f"{df['score'].median():.4f}"
            })
    
    # Save normalization log
    log_df = pd.DataFrame(normalization_log)
    log_df.to_csv('normalized_scores/normalization_log.txt', sep='\t', index=False)
    
    print(f"\n=== Normalization complete ===")
    print(f"Normalized scores saved to: normalized_scores/")
    print(f"Normalization log: normalized_scores/normalization_log.txt")
    print("\nNext step: Collapse replicates by biosample, then by stage")
    print("See DATA_NEEDED.md for DANIO_CODE_METADATA (SQ → biosample → stage mapping)")

if __name__ == '__main__':
    main()
