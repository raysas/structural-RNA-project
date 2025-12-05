#!/usr/bin/env python3
"""
RNA Scoring Profile Plotter

Visualizes pseudo-energy scores as a function of distance for each base pair.

Usage:
  python plot_scores.py --input-dir training_output --output-dir plots
    python plot_scores.py --input-dir training_output --output-dir plots --combined --plotly
"""

import argparse
import os
from pathlib import Path
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize


def load_score_table(filepath):
    """Load a score table CSV file"""
    if not os.path.exists(filepath):
        return None
    
    if 'kde' in filepath:
        # Load KDE data
        data = np.loadtxt(filepath, delimiter=',', skiprows=1)
        print()
        return {
            'distance': data[:, 0],
            'score': data[:, 1]
        }
    
    data = np.loadtxt(filepath, delimiter=',', skiprows=1)
    return {
        'distance_min': data[:, 0],
        'distance_max': data[:, 1],
        'distance': data[:, 2],
        'score': data[:, 3]
    }

def plot_single_profile(pair, score_data, output_path):
    sns.set_theme(style="whitegrid")
    x = score_data['distance']
    y = score_data['score']
    cmap = plt.get_cmap('plasma')
    norm = Normalize(vmin=np.min(y), vmax=np.max(y))

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    segment_values = (y[:-1] + y[1:]) / 2

    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=3)
    lc.set_array(segment_values)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.add_collection(lc)
    sc = ax.scatter(x, y, c=y, cmap=cmap, norm=norm, s=40, zorder=3)

    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y) - 0.1*np.ptp(y), np.max(y) + 0.1*np.ptp(y))
    ax.set_xlabel('Distance (Å)', fontsize=12)
    ax.set_ylabel('Pseudo-energy Score', fontsize=12)
    ax.set_title(f'Scoring Profile: {pair}', fontsize=14, fontweight='bold')
    ax.axhline(0, linestyle='--', color='gray', alpha=0.5, linewidth=1)
    ax.grid(alpha=0.3)

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label('Score', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    

def plot_single_profile_plotly(pair, score_data, output_path):
    x = score_data['distance']
    y = score_data['score']
    y_norm = (y - np.min(y)) / (np.max(y) - np.min(y))
    segment_colors = sample_colorscale('Plasma', y_norm)
    traces = []
    for i in range(len(x)-1):
        trace = go.Scatter(
            x=x[i:i+2],
            y=y[i:i+2],
            mode='lines+markers',
            line=dict(color=segment_colors[i], width=3),
            marker=dict(color=segment_colors[i], size=6),
            showlegend=False,
            hoverinfo='x+y'
        )
        traces.append(trace)
    fig = go.Figure(data=traces)
    fig.update_layout(
        title=f'Scoring Profile: {pair}',
        xaxis_title='Distance (Å)',
        yaxis_title='Pseudo-energy Score',
        template='plotly_white',
        hovermode='closest'
    )
    fig.write_html(output_path)
    return fig



def plot_combined_profiles(pairs_data, output_path):
    '''Put each pair in its own axis within a grid'''
    sns.set_theme(style="whitegrid")
    
    n_pairs = len(pairs_data)
    n_cols = 5  # Adjust columns as desired
    n_rows =2
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows), squeeze=False)
    
    # Flatten axes array for easy indexing
    axes_flat = axes.flatten()
    
    for i, (pair, score_data) in enumerate(pairs_data.items()):
        ax = axes_flat[i]
        x = score_data['distance']
        y = score_data['score']
        
        cmap = plt.get_cmap('plasma')
        norm = Normalize(vmin=np.min(y), vmax=np.max(y))
        
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        segment_values = (y[:-1] + y[1:]) / 2
        
        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=2)
        lc.set_array(segment_values)
        
        ax.add_collection(lc)
        ax.scatter(x, y, c=y, cmap=cmap, norm=norm, s=20, zorder=3)
        
        ax.set_xlim(0, np.max(x))
        ax.set_ylim(np.min(y)-0.5, np.max(y)+0.5)
        ax.set_xlabel('Distance (Å)', fontsize=10)
        ax.set_ylabel('Pseudo-energy Score', fontsize=10)
        ax.set_title(f'Pair: {pair}', fontsize=12, fontweight='bold')
        ax.axhline(0, linestyle='--', color='gray', alpha=0.5, linewidth=1)
        ax.grid(alpha=0.3)
    
    # Hide unused axes if number of pairs < n_rows * n_cols
    for j in range(i+1, n_rows * n_cols):
        axes_flat[j].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_combined_profiles_plotly(pairs_data, output_path):
    """
    Plot multiple base pair profiles in a single Plotly figure with colored segments.
    
    Args:
        pairs_data (dict): Dictionary with pair names as keys and score_data (DataFrame) as values.
        output_path (str): Path to save the HTML figure.
    
    Returns:
        fig (go.Figure): Plotly figure object.
    """
    fig = go.Figure()
    trace_indices = {}
    
    # Add traces for each pair
    for pair_idx, (pair, score_data) in enumerate(pairs_data.items()):
        x = score_data['distance']
        y = score_data['score']
        y_norm = (y - np.min(y)) / (np.max(y) - np.min(y))
        segment_colors = sample_colorscale('Plasma', y_norm)
        
        trace_indices[pair] = []

        is_first_pair = (pair_idx == 0)
        for i in range(len(x) - 1):
            trace = go.Scatter(
                x=x[i:i+2],
                y=y[i:i+2],
                mode='lines+markers',
                line=dict(color=segment_colors[i], width=3),
                marker=dict(color=segment_colors[i], size=6),
                name=pair if i == 0 else None,  
                hoverinfo='x+y+name',
                visible=is_first_pair
            )
            fig.add_trace(trace)
            trace_indices[pair].append(len(fig.data)-1)

    buttons = []
    first_pair = list(trace_indices.keys())[0]
    for pair, indices in trace_indices.items():
        visibility = [False] * len(fig.data)
        for idx in indices:
            visibility[idx] = True
        buttons.append(dict(
            label=pair,
            method='update',
            args=[{'visible': visibility},
                  {'title': f'Scoring Profile: {pair}'}]
        ))

    buttons.append(dict(
        label='All',
        method='update',
        args=[{'visible': [True]*len(fig.data)},
              {'title': 'Combined Scoring Profiles'}]
    ))

    fig.update_layout(
        title=f'Scoring Profile: {first_pair}',
        xaxis_title='Distance (Å)',
        yaxis_title='Pseudo-energy Score',
        template='plotly_white',
        hovermode='closest',
        updatemenus=[dict(active=0, buttons=buttons, x=1.1, y=0.9)],
        showlegend=False
    )

    fig.write_html(output_path)
    return fig



def main():
    parser = argparse.ArgumentParser(
        description="Plot RNA scoring profiles from trained tables"
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        default='training_output',
        help='Directory containing score tables from train.py'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='plots',
        help='Directory for output plot files'
    )
    parser.add_argument(
        '--combined',
        action='store_true',
        help='Generate a combined plot with all profiles'
    )
    parser.add_argument(
        '--plotly',
        action='store_true',
        help='Use Plotly for interactive plots'
    )
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, 'png'), exist_ok=True)
    
    pairs = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'GG', 'UC', 'UG', 'UU']
    
    print(f"Plotting scoring profiles...")
    print(f"  Input directory: {args.input_dir}")
    print(f"  Output directory: {args.output_dir}")
    
    pairs_data = {}
    for pair in pairs:
        score_file = os.path.join(args.input_dir, f'score_{pair}.csv')
        score_data = load_score_table(score_file)
        
        if score_data is None:
            print(f"  WARNING: Score table not found for {pair}, skipping...")
            continue
        
        pairs_data[pair] = score_data
        
        output_path = os.path.join(args.output_dir, 'png', f'{pair}_profile.png')
        plot_single_profile(pair, score_data, output_path)
        print(f"  ✓ Plotted {pair}")
    
    
    if args.combined and pairs_data:
        combined_path = os.path.join(args.output_dir,'png', 'all_profiles.png')
        plot_combined_profiles(pairs_data, combined_path)
        print(f"  ✓ Generated combined plot")

    if args.plotly:
        os.makedirs(os.path.join(args.output_dir, 'html'), exist_ok=True)
        for pair, score_data in pairs_data.items():
            output_path = os.path.join(args.output_dir, 'html', f'{pair}_profile.html')
            plot_single_profile_plotly(pair, score_data, output_path)
            print(f"  ✓ Plotted {pair} (Plotly)")
        if args.combined:
            combined_path = os.path.join(args.output_dir, 'html', 'all_profiles.html')
            plot_combined_profiles_plotly(pairs_data, combined_path)
            print(f"  ✓ Generated combined plot (Plotly)")
    
    print(f"\n✓ Plotting complete!")
    print(f"  Generated {len(pairs_data)} individual plots")
    if args.combined:
        print(f"  Combined plot: all_profiles.png")
    
    return 0


if __name__ == '__main__':
    exit(main())
