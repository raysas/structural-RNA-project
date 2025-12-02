"""
RNA Distance Extraction Pipeline

This script processes RNA structures (PDB or mmCIF files), extracts atom coordinates,
computes pairwise distances between selected atoms, and outputs either counts per distance
bin (histogram) or raw distances (kde) for nucleotide pair types. Optional detailed logging 
of all interactions is supported. Parallel processing is implemented for efficiency.

Usage
-----
1. Process a single PDB file:
    python script.py --pdb 1ABC --chains A B --format pdb --atom-mode C3' --dist-mode intra \
                     --cutoff 20.0 --seq-sep 4 --bin-width 1.0 --method histogram \
                     --out-dir dist_data --save-detailed

2. Process a list of PDB IDs from a text file:
    python script.py --list pdb_list.txt --cores 4 --format mmcif --atom-mode C3' C4' \
                     --dist-mode inter --cutoff 25.0 --seq-sep 3 --method kde \
                     --out-dir output_dir --save-detailed

Parameters
----------
--pdb : str
    Single PDB ID or filename to process.
--list : str
    Text file with a list of PDB IDs or filenames to process.
--chains : list of str
    Specific RNA chains to process (optional).
--format : str, default='pdb'
    Input file format: 'pdb' or 'mmcif'.
--all-models : bool
    Process all NMR models if True (default False).
--cores : int
    Number of CPU cores for parallel processing.
--atom-mode : str or list of str, default="C3'"
    Atom selection mode for RNA distance calculations. Options:
      - Single atom name (e.g., "C3'", "C4'")
      - "centroid" to use the centroid of all atoms in a nucleotide
      - "all" to include all atoms
      - List of atom names (e.g., ["C3'", "C4'"]) to include multiple specific atoms
--dist-mode : str, default='intra'
    Distance computation mode: 'intra' (within chain) or 'inter' (between chains).
--cutoff : float, default=20.0
    Maximum distance to consider between atoms (Angstroms).
--seq-sep : int, default=4
    Minimum sequence separation for intra-chain distances.
--bin-width : float, default=1.0
    Width of distance bins when counting interactions.
--method : str, default='histogram'
    Output format: 'histogram' (counts per bin) or 'kde' (raw distances).
--out-dir : str, default='dist_data'
    Directory to save outputs.
--save-detailed : bool
    Save a detailed CSV log of all interactions.

Outputs
-------
- Histogram (counts per bin) of distances for nucleotide pair types: AA, AC, AG, AU, CC, CG, CU, GG, GU, UU, XX
- Raw distances for nucleotide pair types (if method='kde')
- Optional CSV file with detailed interaction data (if --save-detailed)
"""

import argparse
import os
import numpy as np
import pandas as pd
import time
from concurrent.futures import ProcessPoolExecutor
from core import FastParser, OnlineFetcher, DistanceComputer

_WORKER_CONFIG = {}

def parse_input_file(filename):
    """
    Parse a text file listing RNA structures to process.

    Each line should contain:
        PDB_ID [CHAIN1 CHAIN2 ...]

    Parameters
    ----------
    filename : str
        Path to the text file containing PDB IDs or local filenames.

    Returns
    -------
    list of dict
        Each dict contains:
        - 'id': PDB ID or filename (str)
        - 'chains': list of chain identifiers to process (or None)
        - 'is_local': bool indicating if the file exists locally
    """
    targets = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts: 
                continue
            
            pdb_id = parts[0]
            chains = parts[1:] if len(parts) > 1 else None
            is_file = os.path.exists(pdb_id)
            
            targets.append({
                'id': pdb_id,
                'chains': chains,
                'is_local': is_file
            })
    return targets

def worker_process(target):
    """
    Process a single RNA structure: load, parse, and compute distances.

    Parameters
    ----------
    target : dict
        Dictionary with keys:
        - 'id': PDB ID or local filename
        - 'chains': list of chains to process (or None)
        - 'is_local': bool indicating if file exists locally

    Returns
    -------
    tuple
        - dict: distance lists grouped by nucleotide pair type (AA, AC, ..., XX)
        - pandas.DataFrame: detailed interaction log (if requested)
    """
    try:
        cfg = _WORKER_CONFIG
        ident = target['id']
        
        content = None
        if target['is_local']:
            try:
                with open(ident, 'r') as f:
                    content = f.read().splitlines()
            except:
                return None, None
        else:
            content = OnlineFetcher.get_content(ident, file_format=cfg['fmt'])
            
        if not content:
            return None, None

        parser = FastParser(atom_mode=cfg['atom_mode'])
        data = parser.parse(content, file_format=cfg['fmt'], pdb_name=ident, 
                            chains_filter=target['chains'], use_all_models=cfg['all_models'])
        
        if data is None:
            return None, None

        computer = DistanceComputer(cutoff=cfg['cutoff'], seq_separation=cfg['seq_sep'])
        dists, detailed_df = computer.compute(data, dist_mode=cfg['dist_mode'], detailed_output=cfg['detailed'])
        
        return dists, detailed_df
        
    except Exception:
        return None, None

def process_dataset_parallel(targets, config, max_workers=None):
    """
    Process multiple RNA structures in parallel.

    Parameters
    ----------
    targets : list of dict
        Each dict specifies a structure (from parse_input_file or manually constructed).
    config : dict
        Configuration options:
        - 'atom_mode': str or list, which atoms to include
        - 'dist_mode': 'intra' or 'inter'
        - 'cutoff': float, maximum distance to consider
        - 'seq_sep': int, minimum sequence separation for intra-chain distances
        - 'fmt': 'pdb' or 'mmcif'
        - 'all_models': bool, process all NMR models if True
        - 'detailed': bool, include detailed DataFrame output
    max_workers : int, optional
        Maximum number of parallel worker processes. Defaults to None (auto).

    Returns
    -------
    tuple
        - dict: aggregated distance lists per nucleotide pair type
        - list of pandas.DataFrame: detailed logs from all structures (if requested)
    """
    global _WORKER_CONFIG
    _WORKER_CONFIG = config
    
    keys = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'CU', 'GG', 'GU', 'UU', 'XX']
    final_results = {k: [] for k in keys}
    all_detailed_dfs = []
    
    atom_mode_print = config['atom_mode']
    if isinstance(atom_mode_print, list):
        atom_mode_print = "+".join(atom_mode_print)

    print(f"Starting Parallel Processing (Mode: {atom_mode_print})...")
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(worker_process, targets))
        
        for i, (dists, df) in enumerate(results):
            if dists is None:
                continue 
                
            for k, v in dists.items():
                if k in final_results:
                    final_results[k].extend(v)
            
            if df is not None:
                all_detailed_dfs.append(df)
            
            if (i+1) % 10 == 0:
                print(f"  ...processed {i+1}/{len(targets)} structures")

    elapsed = time.time() - start_time
    print(f"Batch processing complete in {elapsed:.2f} seconds.")
    return final_results, all_detailed_dfs

def main():
    """
    Command-line interface for RNA distance extraction.

    Parameters
    ----------
    Uses command-line arguments; see argparse definitions.

    Outputs
    -------
    - Counts per distance bin for nucleotide pair types (if method='histogram')
    - Raw distances for nucleotide pair types (if method='kde')
    - Optional CSV log of all interactions (if --save-detailed)
    """
    parser = argparse.ArgumentParser(description="RNA Distance Extractor")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pdb', type=str, help="Single PDB ID or filename.")
    group.add_argument('--list', type=str, help="Text file containing a list of PDB IDs.")
    
    parser.add_argument('--chains', nargs='+', help="Specific chains to process.")
    parser.add_argument('--format', choices=['pdb', 'mmcif'], default='pdb', help="Input file format.")
    parser.add_argument('--all-models', action='store_true', help="Process all NMR models.")
    
    parser.add_argument('--cores', type=int, default=None, help="Number of CPU cores to use.")
    
    parser.add_argument('--atom-mode', nargs='+', default=["C3'"], 
                        help="Atom selection mode (e.g., C3', centroid, all, or list).")
    
    parser.add_argument('--dist-mode', choices=['intra', 'inter'], default='intra',
                        help="Interaction mode.")
    parser.add_argument('--cutoff', type=float, default=20.0, help="Interaction cutoff distance.")
    parser.add_argument('--seq-sep', type=int, default=4, 
                        help="Minimum sequence separation.")
    parser.add_argument('--bin-width', type=float, default=1.0, help="Width of histogram bins.")
    
    parser.add_argument('--method', choices=['histogram', 'kde'], default='histogram', 
                        help="Output format: counts (histogram) or raw distances (kde).")
    parser.add_argument('--out-dir', type=str, default='dist_data', help="Directory for output files.")
    parser.add_argument('--save-detailed', action='store_true', help="Export detailed interaction log.")
    
    args = parser.parse_args()

    mode_arg = args.atom_mode
    if len(mode_arg) == 1:
        mode_arg = mode_arg[0]

    config = {
        'atom_mode': mode_arg,
        'dist_mode': args.dist_mode,
        'cutoff': args.cutoff,
        'seq_sep': args.seq_sep,
        'fmt': args.format,
        'all_models': args.all_models,
        'detailed': args.save_detailed
    }

    targets = []
    if args.pdb:
        is_local = os.path.exists(args.pdb)
        targets.append({'id': args.pdb, 'chains': args.chains, 'is_local': is_local})
        final_results, all_detailed_dfs = process_dataset_parallel(targets, config, max_workers=1)
    elif args.list:
        targets = parse_input_file(args.list)
        final_results, all_detailed_dfs = process_dataset_parallel(targets, config, max_workers=args.cores)

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"\nWriting results to '{args.out_dir}/' (Format: {args.method.upper()})...")
    
    keys = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'CU', 'GG', 'GU', 'UU', 'XX']

    if args.method == 'histogram':
        bin_edges = np.arange(0, args.cutoff + args.bin_width, args.bin_width)
        for k in keys:
            vals = final_results[k]
            if not vals:
                hist = np.zeros(len(bin_edges) - 1, dtype=int)
            else:
                hist, _ = np.histogram(vals, bins=bin_edges)
            np.savetxt(f"{args.out_dir}/{k}_histogram.txt", hist, fmt='%d')
            
    elif args.method == 'kde':
        for k in keys:
            vals = final_results[k]
            if not vals:
                open(f"{args.out_dir}/{k}_kde_raw.txt", 'w').close()
            else:
                np.savetxt(f"{args.out_dir}/{k}_kde_raw.txt", vals, fmt='%.3f')

    if args.save_detailed and all_detailed_dfs:
        print("Writing detailed interaction log...")
        final_df = pd.concat(all_detailed_dfs, ignore_index=True)
        csv_path = os.path.join(args.out_dir, "detailed_interactions.csv")
        final_df.to_csv(csv_path, index=False)
        print(f"Log saved to {csv_path}")

    print("Pipeline completed successfully.")

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    main()
