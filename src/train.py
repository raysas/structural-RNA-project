#!/usr/bin/env python3
"""
RNA Scoring Function Training Script

Computes observed and reference frequencies from distance histograms,
then calculates pseudo-energy scores using the log-ratio formula.

Usage:
  python train.py --input-dir dist_data --output-dir training_output
"""

import argparse
import os
import numpy as np
import json
from pathlib import Path


def load_histogram(filepath):
    """Load a histogram file and return counts as numpy array."""
    if not os.path.exists(filepath):
        return None
    return np.loadtxt(filepath, dtype=float)


def compute_frequency(counts, pseudocount=1e-6):
    """
    Compute frequency from counts with pseudocount smoothing.
    
    Args:
        counts: Array of counts per bin
        pseudocount: Small value to avoid division by zero
    
    Returns:
        Frequency array (normalized counts)
    """
    total = np.sum(counts)
    if total == 0:
        return np.zeros_like(counts, dtype=float)
    
    # Add pseudocount to avoid zero frequencies
    smoothed_counts = counts + pseudocount
    smoothed_total = total + pseudocount * len(counts)
    
    return smoothed_counts / smoothed_total


def compute_score(obs_freq, ref_freq, max_score=10.0):
    """
    Compute pseudo-energy score: -log(f_obs / f_ref)
    
    Args:
        obs_freq: Observed frequency
        ref_freq: Reference frequency
        max_score: Maximum allowed score value
    
    Returns:
        Score array, capped at max_score
    """
    # Avoid division by zero
    safe_obs = np.maximum(obs_freq, 1e-10)
    safe_ref = np.maximum(ref_freq, 1e-10)
    
    score = -np.log(safe_obs / safe_ref)
    
    # Cap at maximum score
    score = np.minimum(score, max_score)
    
    return score


def main():
    parser = argparse.ArgumentParser(
        description="Train RNA scoring function from distance histograms"
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        default='dist_data',
        help='Directory containing histogram files from extract_distances.py'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='training_output',
        help='Directory for output score and frequency tables'
    )
    parser.add_argument(
        '--max-score',
        type=float,
        default=10.0,
        help='Maximum score value (default: 10.0)'
    )
    parser.add_argument(
        '--pseudocount',
        type=float,
        default=1e-6,
        help='Pseudocount for smoothing (default: 1e-6)'
    )
    parser.add_argument(
        '--cutoff',
        type=float,
        default=20.0,
        help='Maximum distance in Angstroms (default: 20.0)'
    )
    parser.add_argument(
        '--bin-width',
        type=float,
        default=1.0,
        help='Bin width in Angstroms (default: 1.0)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Base pairs to process
    pairs = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'GG', 'CU', 'GU', 'UU']
    
    # Load reference histogram (XX - all pairs)
    ref_file = os.path.join(args.input_dir, 'XX_histogram.txt')
    if not os.path.exists(ref_file):
        print(f"ERROR: Reference file not found: {ref_file}")
        print("Run extract_distances.py first to generate histograms.")
        return 1
    
    ref_counts = load_histogram(ref_file)
    ref_freq = compute_frequency(ref_counts, args.pseudocount)
    
    # Generate distance bin information
    n_bins = len(ref_counts)
    bin_edges = np.linspace(0, args.cutoff, n_bins + 1)
    bin_min = bin_edges[:-1]
    bin_max = bin_edges[1:]
    bin_mid = (bin_min + bin_max) / 2
    
    print(f"Training scoring function...")
    print(f"  Input directory: {args.input_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Number of bins: {n_bins}")
    print(f"  Distance range: 0-{args.cutoff} Å")
    print(f"  Max score: {args.max_score}")
    
    # Save reference frequency table
    ref_table = np.column_stack([bin_min, bin_max, bin_mid, ref_counts, ref_freq])
    ref_path = os.path.join(args.output_dir, 'freq_XX.csv')
    np.savetxt(
        ref_path,
        ref_table,
        delimiter=',',
        header='Distance_Min,Distance_Max,Distance_Mid,Count,Frequency',
        comments='',
        fmt='%.6f'
    )
    
    # Process each base pair
    processed_pairs = []
    for pair in pairs:
        hist_file = os.path.join(args.input_dir, f'{pair}_histogram.txt')
        
        if not os.path.exists(hist_file):
            print(f"  WARNING: Histogram not found for {pair}, skipping...")
            continue
        
        # Load observed counts
        obs_counts = load_histogram(hist_file)
        
        if len(obs_counts) != len(ref_counts):
            print(f"  WARNING: {pair} has different bin count, skipping...")
            continue
        
        # Compute observed frequency
        obs_freq = compute_frequency(obs_counts, args.pseudocount)
        
        # Compute score
        score = compute_score(obs_freq, ref_freq, args.max_score)
        
        # Save frequency table
        freq_table = np.column_stack([bin_min, bin_max, bin_mid, obs_freq])
        freq_path = os.path.join(args.output_dir, f'freq_{pair}.csv')
        np.savetxt(
            freq_path,
            freq_table,
            delimiter=',',
            header='Distance_Min,Distance_Max,Distance_Mid,Frequency',
            comments='',
            fmt='%.6f'
        )
        
        # Save score table
        score_table = np.column_stack([bin_min, bin_max, bin_mid, score])
        score_path = os.path.join(args.output_dir, f'score_{pair}.csv')
        np.savetxt(
            score_path,
            score_table,
            delimiter=',',
            header='Distance_Min,Distance_Max,Distance_Mid,Score',
            comments='',
            fmt='%.6f'
        )
        
        processed_pairs.append(pair)
        print(f"  ✓ Processed {pair}")
    
    # Save metadata
    metadata = {
        'input_dir': args.input_dir,
        'output_dir': args.output_dir,
        'max_score': args.max_score,
        'pseudocount': args.pseudocount,
        'cutoff': args.cutoff,
        'bin_width': args.bin_width,
        'n_bins': int(n_bins),
        'pairs': processed_pairs
    }
    
    metadata_path = os.path.join(args.output_dir, 'metadata.json')
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\n✓ Training complete!")
    print(f"  Processed {len(processed_pairs)} base pairs")
    print(f"  Frequency tables: freq_*.csv")
    print(f"  Score tables: score_*.csv")
    print(f"  Metadata: metadata.json")
    
    return 0


if __name__ == '__main__':
    exit(main())
