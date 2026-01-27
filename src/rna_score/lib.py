import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Union
import os
import pandas as pd
import glob
import numpy as np

from pathlib import Path
from typing import List, Dict, Optional, Union

from access_rna_structures import RNAStructureDownloader

try:
    import plotly.graph_objects as go
    from plotly.colors import sample_colorscale
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


def download_rna_structures(
    output_dir: Union[str, Path] = "rna_structures",
    max_structures: int = 10,
    download_all: bool = False,
    rna_only: bool = False,
    formats: List[str] = ["pdb", "cif"],
    workers: int = 5,
    show_info: bool = False,
    list_only: bool = False,
) -> Dict[str, Union[List[str], Path, List[Dict]]]:
    """
    Download RNA structures from the RCSB PDB database.

    This function is a programmatic wrapper around the RNA structure
    downloader pipeline, equivalent to the CLI command `rna-score access`.

    Parameters
    ----------
    output_dir : str or Path, default "rna_structures"
        Directory where structures and metadata will be saved.
    max_structures : int, default 10
        Maximum number of structures to retrieve.
        Ignored if `download_all=True`.
    download_all : bool, default False
        If True, download all available RNA structures.
    rna_only : bool, default False
        If True, restrict results to RNA-only structures
        (no protein or DNA).
    formats : list of {"pdb", "cif"}, default ["pdb", "cif"]
        File formats to download.
    workers : int, default 5
        Number of parallel download workers.
    show_info : bool, default False
        If True, fetch and return metadata for a subset of structures
        before downloading.
    list_only : bool, default False
        If True, only retrieve and save PDB IDs without downloading files.

    Returns
    -------
    dict
        Dictionary with the following keys:

        - ``pdb_ids`` : list of str  
          Retrieved PDB identifiers.
        - ``output_dir`` : Path  
          Base output directory.
        - ``pdb_dir`` : Path  
          Directory containing downloaded PDB files.
        - ``mmcif_dir`` : Path  
          Directory containing downloaded mmCIF files.
        - ``info`` : list of dict or None  
          Structure metadata if `show_info=True`.

    Notes
    -----
    This function mirrors the CLI behavior exactly but is suitable
    for notebooks, pipelines, and Python APIs.
    """
    output_dir = Path(output_dir)

    downloader = RNAStructureDownloader(
        output_dir=str(output_dir),
        max_workers=workers,
    )

    max_results = None if download_all else max_structures

    # --- Search phase ---
    if rna_only:
        pdb_ids = downloader.get_rna_only_structures(max_results)
    else:
        pdb_ids = downloader.search_rna_structures(max_results)

    if not pdb_ids:
        raise RuntimeError("No RNA structures found.")

    # Save ID list
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "pdb_ids.txt").write_text("\n".join(pdb_ids))

    info = None
    if show_info:
        info = downloader.get_structure_info(pdb_ids)

    if list_only:
        return {
            "pdb_ids": pdb_ids,
            "output_dir": output_dir,
            "pdb_dir": output_dir / "pdb",
            "mmcif_dir": output_dir / "mmcif",
            "info": info,
        }

    # --- Download phase ---
    downloader.download_batch(pdb_ids, formats=formats)

    (output_dir / "downloaded_ids.txt").write_text("\n".join(pdb_ids))

    return {
        "pdb_ids": pdb_ids,
        "output_dir": output_dir,
        "pdb_dir": output_dir / "pdb",
        "mmcif_dir": output_dir / "mmcif",
        "info": info,
    }



class RNAScorer:
    """
    Orchestrates RNA distance extraction, scoring, KDE training,
    and full workflow execution, with memoization of both paths and data tables.
    
    Parameters
    ----------
    base_dir : str
        Base directory containing the `src/` folder with scripts.
    """

    def __init__(self, base_dir: str = "."):
        self.base_dir = Path(base_dir)
        
        # Try to find scripts directory - check both development and installed locations
        if (self.base_dir / "src").exists():
            self.scripts_dir = self.base_dir / "src"
        elif (Path(__file__).parent.parent).exists():
            # When installed, scripts are siblings to rna_score package
            self.scripts_dir = Path(__file__).parent.parent
        else:
            self.scripts_dir = self.base_dir / "src"

        # Paths
        self.distances_dir: Optional[Path] = None
        self.training_dir: Optional[Path] = None
        self.kde_dir: Optional[Path] = None
        self.scoring_tables_dir: Optional[Path] = None
        self.scores_file: Optional[Path] = None

        # Config tracking
        self.extraction_method: Optional[str] = None  # Track extraction method

        # Loaded tables
        self.histograms: dict[str, pd.DataFrame] = {}
        self.scoring_tables: dict[str, pd.DataFrame] = {}
        self.kde_tables: dict[str, pd.DataFrame] = {}
        self.scores: Optional[pd.DataFrame] = None


    def extract_distances(self, **kwargs) -> Path:
        """
        Extract interatomic distances from RNA structures.

        Computes pairwise atom-atom or residue-residue distances from PDB/mmCIF files
        and exports binned histograms or raw distances for KDE training. Equivalent to
        the CLI command `rna-score extract`.

        Parameters
        ----------
        pdb : str, optional
            Path to a single PDB/mmCIF file to process.
        pdb_list : str, optional
            Path to a text file containing PDB IDs and optional chain IDs.
            Format: `<PDB_ID> [CHAIN_ID1 CHAIN_ID2 ...]`
            Example lines: `1EHZ A`, `1Y26 B C`, `2OEU`
        folder : str, optional
            Directory containing multiple PDB/mmCIF files to process.
        fmt : str, default="pdb"
            Structure format. Options: "pdb" or "mmcif". Auto-detected for local files.
        atom_mode : list of str, default=["C3'"]
            Atom selection mode(s). Options:
            - "C3'" (default): C3' atoms only
            - "centroid": geometric center of each residue
            - "all": all atoms
            - Multiple atom names: e.g., ["P", "C4"]
        chains : list of str, optional
            Specific chain IDs to process. If not provided, all chains are used.
        dist_mode : str, default="intra"
            Distance calculation mode:
            - "intra": within a single chain (default)
            - "inter": between distinct chains
        cutoff : float, default=20.0
            Maximum atom-atom distance (Å) counted as a contact.
        seq_sep : int, default=4
            Minimum sequence separation for intra-chain contacts (residues).
            Ignored in "inter" mode. Distances considered from residue i to i+seq_sep.
        bin_width : float, default=1.0
            Bin width (Å) for histogram discretization (histogram method only).
        cores : int, optional
            Number of CPU cores to use for parallel processing. Default: all available.
        out_dir : str or Path, default="dist_data"
            Output directory where results will be written.
        save_detailed : bool, default=False
            If True, saves a detailed CSV log with full information for every measured
            distance including PDB ID, model, chain IDs, residue IDs, atom names,
            B-factors, altlocs, distance, and pair type.
        all_models : bool, default=False
            If True, processes all models in NMR structures. Default: first model only.
        method : str, default="histogram"
            Extraction method:
            - "histogram": binned counts for direct scoring
            - "kde": raw distances for kernel density estimation

        Returns
        -------
        Path
            Path to the output directory containing extracted distance distributions.

        Raises
        ------
        ValueError
            If none of `pdb`, `pdb_list`, or `folder` is provided.
        FileNotFoundError
            If the extract_distances.py script is not found.

        Examples
        --------
        Extract distances from a single PDB file:
        
        >>> scorer = RNAScorer()
        >>> scorer.extract_distances(pdb="structure.pdb", out_dir="distances")
        
        Extract from a list of PDB IDs with specific chains:
        
        >>> scorer.extract_distances(
        ...     pdb_list="structures.txt",
        ...     chains=["A", "B"],
        ...     dist_mode="inter",
        ...     save_detailed=True
        ... )
        
        Extract using multiple atom types for KDE training:
        
        >>> scorer.extract_distances(
        ...     folder="rna_structures/mmcif",
        ...     fmt="mmcif",
        ...     atom_mode=["P", "C4"],
        ...     method="kde"
        ... )

        Notes
        -----
        - The method caches extracted data in `self.distances_dir` and loads histograms
          into `self.histograms` for downstream use.
        - For large datasets, use the `cores` parameter to parallelize processing.
        - The `save_detailed` option generates very large files for extensive datasets.
        """
        out_dir_path = Path(kwargs.get("out_dir", "dist_data"))
        script_path = self.scripts_dir / "extract_distances.py"
        if not script_path.exists():
            raise FileNotFoundError(f"extract_distances.py not found at {script_path}")

        cmd = [sys.executable, str(script_path)]
        for key in ["pdb", "pdb_list", "folder"]:
            if kwargs.get(key):
                cmd += [f"--{key.replace('_','-')}", kwargs[key]]
                break
        else:
            raise ValueError("Must provide one of pdb, pdb_list, or folder")

        if kwargs.get("chains"):
            cmd += ["--chains"] + kwargs["chains"]

        cmd += ["--format", kwargs.get("fmt", "pdb")]
        if kwargs.get("all_models"):
            cmd.append("--all-models")
        if kwargs.get("cores"):
            cmd += ["--cores", str(kwargs["cores"])]
        cmd += ["--atom-mode"] + kwargs.get("atom_mode", ["C3'"])
        cmd += ["--dist-mode", kwargs.get("dist_mode", "intra")]
        cmd += ["--cutoff", str(kwargs.get("cutoff", 20.0))]
        cmd += ["--seq-sep", str(kwargs.get("seq_sep", 4))]
        cmd += ["--bin-width", str(kwargs.get("bin_width", 1.0))]
        cmd += ["--method", kwargs.get("method", "histogram")]
        cmd += ["--out-dir", str(out_dir_path)]
        if kwargs.get("save_detailed"):
            cmd.append("--save-detailed")

        subprocess.run(cmd, check=True)
        self.distances_dir = out_dir_path

        # -- load extracted histograms based on method
        self.histograms = {}
        method = kwargs.get("method", "histogram")        
        self.extraction_method = method  # remember the method used        
        if method == "kde":
            # -- load KDE raw distance files
            kde_files = glob.glob(str(out_dir_path / "*_kde_raw.txt"))
            for f in kde_files:
                try:
                    df = pd.read_csv(f, header=None, names=["distance"])
                    self.histograms[Path(f).stem] = df
                except Exception as e:
                    print(f"Warning: Could not load KDE file {Path(f).name}: {e}")
        else:
            # -- load histogram files
            hist_files = glob.glob(str(out_dir_path / "*_histogram.txt"))
            for f in hist_files:
                try:
                    df = pd.read_csv(f, header=None, names=["count"])
                    self.histograms[Path(f).stem] = df
                except Exception as e:
                    print(f"Warning: Could not load histogram {Path(f).name}: {e}")

        if not self.histograms:
            expected_pattern = "*_kde_raw.txt" if method == "kde" else "*_histogram.txt"
            raise FileNotFoundError(
                f"No histogram files found in {out_dir_path}. Expected pattern: {expected_pattern}"
            )

        return out_dir_path

    def train_scoring(self, **kwargs) -> Path:
        """
        Train scoring tables from extracted distance distributions.

        Converts raw distance histograms or KDE data into scoring tables by computing
        statistical potentials or probability densities. These tables are used to score
        new RNA structures. Equivalent to the CLI command `rna-score train`.

        Parameters
        ----------
        input_dir : str or Path, optional
            Directory containing histogram or KDE distance files (typically output from
            `extract_distances`). If not provided, uses `self.distances_dir` from a
            previous extraction step.
        output_dir : str or Path, default="training_output"
            Output directory where scoring tables and metadata will be written.
        max_score : float, default=10.0
            Maximum score cap applied to individual distance scores to prevent
            extreme outliers from dominating the total score.
        pseudocount : float, default=0.0
            Pseudocount added to histogram bins to avoid division by zero and
            stabilize scoring for rare distance bins. Higher values increase
            smoothing but may reduce discrimination.
        cutoff : float, default=20.0
            Maximum distance (Å) to consider in scoring table generation.
            Must match the cutoff used during distance extraction.
        bin_width : float, default=1.0
            Bin width (Å) for histogram discretization. Must match the bin width
            used during distance extraction.
        method : str, default="histogram"
            Training method:
            - "histogram": converts binned counts to statistical potentials
            - "kde": trains kernel density estimation models from raw distances

        Returns
        -------
        Path
            Path to the output directory containing scoring tables.

        Raises
        ------
        ValueError
            If no input directory is provided and distances have not been extracted.
        FileNotFoundError
            If the train.py script is not found.

        Examples
        --------
        Train scoring tables after extracting distances:
        
        >>> scorer = RNAScorer()
        >>> scorer.extract_distances(folder="structures/", out_dir="distances")
        >>> scorer.train_scoring(output_dir="scoring_tables")
        
        Train with custom parameters and pseudocount for smoothing:
        
        >>> scorer.train_scoring(
        ...     input_dir="distances",
        ...     output_dir="training_output",
        ...     max_score=15.0,
        ...     pseudocount=1.0,
        ...     method="histogram"
        ... )
        
        Train KDE-based scoring from raw distance data:
        
        >>> scorer.train_scoring(
        ...     input_dir="kde_distances",
        ...     output_dir="kde_scoring",
        ...     method="kde"
        ... )

        Notes
        -----
        - The method caches the output directory in `self.training_dir` and
          `self.scoring_tables_dir`, and loads scoring tables into
          `self.scoring_tables` for downstream use.
        - Scoring tables are saved as CSV files named `score_<PAIR_TYPE>.csv`
          (e.g., `score_AA.csv`, `score_CG.csv`).
        - A metadata JSON file is also generated containing training parameters.
        - Ensure `cutoff` and `bin_width` match those used in distance extraction.
        """
        input_dir = kwargs.get("input_dir") or self.distances_dir
        if input_dir is None:
            raise ValueError("No input_dir provided and distances not extracted yet.")
        
        # Resolve input_dir to absolute path
        input_dir = Path(input_dir)
        if not input_dir.is_absolute():
            # Try relative to base_dir first
            candidate = self.base_dir / input_dir
            if candidate.exists():
                input_dir = candidate.resolve()
            # Try relative to package root
            elif not input_dir.exists():
                pkg_root = Path(__file__).resolve().parents[2]
                candidate = pkg_root / input_dir
                if candidate.exists():
                    input_dir = candidate.resolve()
            else:
                input_dir = input_dir.resolve()

        output_dir_path = Path(kwargs.get("output_dir", "training_output"))
        script_path = self.scripts_dir / "train.py"
        if not script_path.exists():
            raise FileNotFoundError(f"train.py not found at {script_path}")

        # Use extraction method if available and not overridden
        method = kwargs.get("method") or self.extraction_method or "histogram"
        
        cmd = [
            sys.executable, str(script_path),
            "--input-dir", str(input_dir),
            "--output-dir", str(output_dir_path),
            "--max-score", str(kwargs.get("max_score", 10.0)),
            "--pseudocount", str(kwargs.get("pseudocount", 0.0)),
            "--cutoff", str(kwargs.get("cutoff", 20.0)),
            "--bin-width", str(kwargs.get("bin_width", 1.0)),
            "--method", method
        ]
        
        result = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if result.returncode != 0:
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
            
        self.training_dir = output_dir_path
        self.scoring_tables_dir = output_dir_path

        # Load scoring tables
        self.scoring_tables = {}
        for f in glob.glob(str(output_dir_path / "score_*.csv")):
            self.scoring_tables[Path(f).stem] = pd.read_csv(f)
        return output_dir_path


    def save_plot_scores(self, **kwargs) -> Path:
        """
        Generate scoring profile plots from trained score tables.

        Wraps the CLI command ``rna-score plot`` to render individual pair-type
        profiles and (optionally) a combined multi-panel plot.

        Parameters
        ----------
        input_dir : str or Path, optional
            Directory containing score tables (e.g., outputs from
            :meth:`train_scoring`). If not provided, falls back to
            ``self.scoring_tables_dir``.
        output_dir : str or Path, default="plots"
            Directory where plot images will be saved (PNG/HTML as produced by
            the underlying plotting script).
        combined : bool, default=False
            If True, also generates a combined plot that includes all pair-type
            profiles in one figure.

        Returns
        -------
        Path
            Path to the directory containing the generated plots.

        Raises
        ------
        ValueError
            If no input directory is provided and ``self.scoring_tables_dir`` is
            not set.
        FileNotFoundError
            If the ``plot_scores.py`` script cannot be located.

        Examples
        --------
        Plot profiles after training (default individual plots):

        >>> scorer = RNAScorer()
        >>> scorer.train_scoring(input_dir="distances", output_dir="training_output")
        >>> scorer.plot_scores(output_dir="plots")

        Plot with a combined figure:

        >>> scorer.plot_scores(
        ...     input_dir="training_output",
        ...     output_dir="plots",
        ...     combined=True,
        ... )

        Notes
        -----
        - The plotting behavior mirrors the CLI ``rna-score plot`` command.
        - Generated files depend on the plotting script (e.g., separate PNGs per
          pair type and an optional combined plot).
        """
        input_dir = kwargs.get("input_dir") or self.scoring_tables_dir
        if input_dir is None:
            raise ValueError("No input_dir provided and scoring tables not loaded yet.")

        output_dir_path = Path(kwargs.get("output_dir", "plots"))
        script_path = self.scripts_dir / "plot_scores.py"
        if not script_path.exists():
            raise FileNotFoundError(f"plot_scores.py not found at {script_path}")

        cmd = [
            sys.executable, str(script_path),
            "--input-dir", str(input_dir),
            "--output-dir", str(output_dir_path),
        ]
        if kwargs.get("combined"):
            cmd.append("--combined")

        subprocess.run(cmd, check=True)
        return output_dir_path
    
    def plot_scores(self, **kwargs):
        """
        Generate an interactive Plotly figure of scoring profiles.

        Creates a combined interactive plot with dropdown buttons to switch between
        different base-pair types (AA, AC, AG, etc.) and view all profiles together.

        Parameters
        ----------
        input_dir : str or Path, optional
            Directory containing score tables (e.g., outputs from
            :meth:`train_scoring`). If not provided, falls back to
            ``self.scoring_tables_dir``.

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive Plotly figure with dropdown menu to select pair types.

        Raises
        ------
        ImportError
            If Plotly is not installed.
        ValueError
            If no input directory is provided and ``self.scoring_tables_dir`` is not set.
        FileNotFoundError
            If no score table files are found in the input directory.

        Examples
        --------
        Create and display an interactive plot:

        >>> scorer = RNAScorer()
        >>> scorer.train_scoring(input_dir="distances", output_dir="training_output")
        >>> fig = scorer.plot_scores()
        >>> fig.show()

        Save the plot to HTML:

        >>> fig = scorer.plot_scores(input_dir="training_output")
        >>> fig.write_html("interactive_plot.html")

        Use in Jupyter notebook:

        >>> fig = scorer.plot_scores()
        >>> fig  # displays inline

        Notes
        -----
        - The figure includes an interactive dropdown menu to switch between pair types.
        - Each profile is color-coded using the Plasma colorscale based on score values.
        - Use the "All" button to view all profiles simultaneously.
        """
        if not PLOTLY_AVAILABLE:
            raise ImportError(
                "Plotly is required for interactive plotting. "
                "Install it with: pip install plotly"
            )

        # -- resolve input directory robustly (handles notebook cwd vs project root)
        input_dir = kwargs.get("input_dir") or self.scoring_tables_dir
        if input_dir is None:
            raise ValueError("No input_dir provided and scoring tables not loaded yet.")

        input_dir = Path(input_dir)
        if not input_dir.exists():
            # Try resolving relative to provided base_dir
            if hasattr(self, "base_dir"):
                candidate = Path(self.base_dir) / input_dir
                if candidate.exists():
                    input_dir = candidate
            # Try resolving relative to package root (src/..)
            if not input_dir.exists():
                pkg_root = Path(__file__).resolve().parents[2]
                candidate = pkg_root / input_dir
                if candidate.exists():
                    input_dir = candidate

        pairs = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'CU', 'GG', 'GU', 'UU']
        pairs_data = {}

        # Load score tables
        for pair in pairs:
            path = input_dir / f"score_{pair}.csv"
            if not path.exists():
                continue

            try:
                df = pd.read_csv(path)

                # Normalize column names for flexible loading
                colmap = {c.lower(): c for c in df.columns}

                # Identify distance column or derive it
                distance = None
                if 'distance' in colmap:
                    distance = df[colmap['distance']].values
                elif 'distance_mid' in colmap:
                    distance = df[colmap['distance_mid']].values
                elif 'distance_min' in colmap and 'distance_max' in colmap:
                    distance = (df[colmap['distance_min']].values + df[colmap['distance_max']].values) / 2.0

                # Identify score column
                score_col = None
                for key in ['score', 'log_score']:
                    if key in colmap:
                        score_col = colmap[key]
                        break
                if score_col is None and 'Score' in df.columns:
                    score_col = 'Score'

                if distance is not None and score_col is not None:
                    pairs_data[pair] = pd.DataFrame({
                        'distance': distance,
                        'score': df[score_col].values,
                    })
                else:
                    # Skip files that don't have usable columns
                    continue
            except Exception as e:
                print(f"Warning: Could not load {pair}: {e}")
                continue

        if not pairs_data:
            # Fallback: show what files are present to aid debugging
            available = list(input_dir.glob("score_*.csv"))
            raise FileNotFoundError(
                f"No valid score tables found in {input_dir}. "
                "Expected files like score_AA.csv, score_CG.csv, etc. "
                f"Found: {[p.name for p in available]}"
            )

        # Create the figure
        fig = go.Figure()
        trace_map = {}

        # Add traces for each pair
        for idx, (pair, score_df) in enumerate(pairs_data.items()):
            x = score_df['distance'].values
            y = score_df['score'].values
            
            # Normalize y values for color mapping
            y_norm = (y - np.min(y)) / (np.max(y) - np.min(y)) if np.max(y) != np.min(y) else np.zeros_like(y)
            colors = sample_colorscale("Plasma", y_norm)

            trace_map[pair] = []
            visible = (idx == 0)  # Only first pair visible initially

            # Create line segments with individual colors
            for i in range(len(x) - 1):
                fig.add_trace(go.Scatter(
                    x=x[i:i+2],
                    y=y[i:i+2],
                    mode="lines+markers",
                    name=pair if i == 0 else None,
                    line=dict(color=colors[i], width=3),
                    marker=dict(color=colors[i], size=6),
                    visible=visible,
                    showlegend=False
                ))
                trace_map[pair].append(len(fig.data) - 1)

        # Create dropdown buttons
        buttons = []
        all_visible = [True] * len(fig.data)

        for pair, indices in trace_map.items():
            visibility = [False] * len(fig.data)
            for i in indices:
                visibility[i] = True

            buttons.append({
                "label": pair,
                "method": "update",
                "args": [{"visible": visibility}, {"title": f"Scoring Profile: {pair}"}]
            })

        buttons.append({
            "label": "All",
            "method": "update",
            "args": [{"visible": all_visible}, {"title": "Combined Scoring Profiles"}]
        })

        # Update layout
        fig.update_layout(
            title=f"Scoring Profile: {list(pairs_data.keys())[0]}",
            xaxis_title="Distance (Å)",
            yaxis_title="Pseudo-energy Score",
            template="plotly_white",
            updatemenus=[{
                "buttons": buttons,
                "direction": "down",
                "showactive": True,
                "x": 1.15,
                "xanchor": "left",
                "y": 0.95,
                "yanchor": "top"
            }],
            showlegend=False
        )

        return fig



    def score_structure(self, **kwargs) -> Optional[pd.DataFrame]:
        """
        Score RNA structure(s) using trained scoring tables.

        Evaluates one or more RNA structures by computing a statistical potential
        score based on observed interatomic distances and their frequencies in the
        training set. Lower (more negative) scores indicate structures more similar
        to the training distribution. Equivalent to the CLI command `rna-score score`.

        Parameters
        ----------
        pdb_path : str or Path, optional
            Path to a single PDB/mmCIF file to score. Mutually exclusive with
            `folder` and `list`.
        folder : str or Path, optional
            Directory containing multiple PDB/mmCIF files to score. Mutually
            exclusive with `pdb_path` and `list`.
        list : str or Path, optional
            Path to a text file containing structure file paths (local files or
            PDB IDs) to score, one per line.
            Example lines: `/path/to/file1.pdb`, `/path/to/file2.cif`, `1EHZ`
        tables_dir : str or Path, optional
            Directory containing scoring tables (output from `train_scoring`).
            If not provided, uses `self.scoring_tables_dir` from a previous
            training step.
        file_format : str, default="pdb"
            Input file format. Options: "pdb" or "mmcif". Auto-detected for
            local files.
        cutoff : float, default=20.0
            Maximum distance (Å) to consider when scoring. Should match the
            cutoff used during training.
        seq_sep : int, default=4
            Minimum sequence separation (residues) for intra-chain contacts.
            Should match the seq_sep used during training.
        atom_mode : str or list of str, default="C3'"
            Atom selection mode(s) for scoring. Should match the atom_mode
            used during distance extraction and training.
            Options: "C3'", "centroid", "all", or list of atom names.
        detailed : bool, default=False
            If True, prints detailed per-interaction scores showing individual
            distance contributions to the total score.
        output_csv : str or Path, optional
            Output CSV file path for saving scores. If provided, scores are
            written to this file and loaded into `self.scores`.

        Returns
        -------
        pd.DataFrame or None
            DataFrame containing scores for each structure if `output_csv` is
            specified, otherwise None. Columns typically include structure ID,
            total score, and optionally per-pair-type scores.

        Raises
        ------
        ValueError
            If no scoring tables are provided and model has not been trained,
            or if no structure input (`pdb_path`, `folder`, or `list`) is provided.
        FileNotFoundError
            If the score_structures.py script is not found.

        Examples
        --------
        Score a single structure after training:
        
        >>> scorer = RNAScorer()
        >>> scorer.train_scoring(input_dir="distances", output_dir="tables")
        >>> scorer.score_structure(pdb_path="test_structure.pdb", output_csv="scores.csv")
        
        Score multiple structures from a folder:
        
        >>> scorer.score_structure(
        ...     folder="structures_to_score/",
        ...     tables_dir="training_output",
        ...     file_format="mmcif",
        ...     output_csv="batch_scores.csv"
        ... )
        
        Score with detailed output showing individual contributions:
        
        >>> scorer.score_structure(
        ...     pdb_path="structure.pdb",
        ...     detailed=True,
        ...     output_csv="detailed_scores.csv"
        ... )
        
        Score structures from a list file:
        
        >>> scorer.score_structure(
        ...     list="structures_to_score.txt",
        ...     tables_dir="scoring_tables",
        ...     output_csv="list_scores.csv"
        ... )

        Notes
        -----
        - The method caches the output file path in `self.scores_file` and loads
          scores into `self.scores` for downstream analysis.
        - Ensure scoring parameters (`cutoff`, `seq_sep`, `atom_mode`) match
          those used during distance extraction and training for consistency.
        - The `detailed` flag is useful for debugging but produces verbose output
          for structures with many contacts.
        """
        tables_dir = kwargs.get("tables_dir") or self.scoring_tables_dir
        if tables_dir is None:
            raise ValueError("No scoring tables provided and model not trained yet.")

        pdb_path = kwargs.get("pdb_path")
        if pdb_path is None:
            raise ValueError("pdb_path is required for scoring.")

        output_csv = kwargs.get("output_csv")
        script_path = self.scripts_dir / "score_structures.py"
        if not script_path.exists():
            raise FileNotFoundError(f"score_structures.py not found at {script_path}")

        cmd = [
            sys.executable, str(script_path),
            "--pdb", str(pdb_path),
            "--tables", str(tables_dir),
            "--format", kwargs.get("file_format", "pdb"),
            "--cutoff", str(kwargs.get("cutoff", 20.0)),
            "--seq-sep", str(kwargs.get("seq_sep", 4)),
        ]

        atom_mode = kwargs.get("atom_mode", "C3'")
        if isinstance(atom_mode, list):
            cmd += ["--atom-mode"] + atom_mode
        else:
            cmd += ["--atom-mode", atom_mode]

        if kwargs.get("detailed"):
            cmd.append("--detailed")
        if output_csv:
            cmd += ["--output", str(output_csv)]

        subprocess.run(cmd, check=True)
        if output_csv:
            self.scores_file = Path(output_csv)
            self.scores = pd.read_csv(output_csv)
            return self.scores
        return None

    def run_workflow(self, **kwargs):
        """
        Run the full workflow: extract => train => score
        Memoizes both paths and loaded tables.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments passed to `extract_distances`, `train_scoring`,
            and `score_structure` methods.
        Returns
        -------
        None
            The method runs the full workflow and saves scores to a CSV file
            specified in `score_structure` parameters

        """
        extracted_dir = self.extract_distances(
            pdb_list=kwargs.get("train_list"),
            folder=kwargs.get("train_folder"),
            chains=kwargs.get("chains"),
            atom_mode=kwargs.get("atom_mode", ["C3'"]),
            dist_mode=kwargs.get("dist_mode", "intra"),
            cutoff=kwargs.get("cutoff", 20.0),
            seq_sep=kwargs.get("seq_sep", 4),
            bin_width=kwargs.get("bin_width", 1.0),
            method=kwargs.get("method", "histogram"),
            cores=kwargs.get("cores"),
            out_dir=os.path.join(kwargs.get("output_dir", "workflow_output"), "extracted"),
            save_detailed=kwargs.get("save_detailed", False),
            fmt=kwargs.get("fmt", "pdb")
        )

        training_dir = self.train_scoring(
            input_dir=extracted_dir,
            output_dir=os.path.join(kwargs.get("output_dir", "workflow_output"), "training_output"),
            max_score=kwargs.get("max_score", 10.0),
            pseudocount=kwargs.get("pseudocount", 0.0),
            cutoff=kwargs.get("cutoff", 20.0),
            bin_width=kwargs.get("bin_width", 1.0),
            method=kwargs.get("method", "histogram")
        )

        self.score_structure(
            pdb_path=kwargs.get("score_list") or kwargs.get("score_folder"),
            tables_dir=training_dir,
            cutoff=kwargs.get("cutoff", 20.0),
            seq_sep=kwargs.get("seq_sep", 4),
            atom_mode=kwargs.get("atom_mode", ["C3'"]),
            output_csv=os.path.join(kwargs.get("output_dir", "workflow_output"), "scores.csv")
        )

        print(f"Workflow complete. Scores saved in {self.scores_file}")

    def __call__(self, **kwargs):
        """Callable shortcut to score a given structure ONLY IF already trained

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments passed to `score_structure`.
        Returns
        -------
        pd.DataFrame or None
            DataFrame containing scores for each structure if `output_csv` is
            specified, otherwise None

        
        
        """
        return self.score_structure(**kwargs)