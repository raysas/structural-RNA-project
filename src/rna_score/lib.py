import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Union
import os
import pandas as pd
import glob


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

        # Paths
        self.distances_dir: Optional[Path] = None
        self.training_dir: Optional[Path] = None
        self.kde_dir: Optional[Path] = None
        self.scoring_tables_dir: Optional[Path] = None
        self.scores_file: Optional[Path] = None

        # Loaded tables
        self.histograms: dict[str, pd.DataFrame] = {}
        self.scoring_tables: dict[str, pd.DataFrame] = {}
        self.kde_tables: dict[str, pd.DataFrame] = {}
        self.scores: Optional[pd.DataFrame] = None

    # ----------------------------
    # Distance extraction
    # ----------------------------
    def extract_distances(self, **kwargs) -> Path:
        out_dir_path = Path(kwargs.get("out_dir", "dist_data"))
        script_path = self.base_dir / "src" / "extract_distances.py"
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

        # Load all *_hist.csv as DataFrames
        self.histograms = {}
        for f in glob.glob(str(out_dir_path / "*.csv")):
            self.histograms[Path(f).stem] = pd.read_csv(f)

        return out_dir_path

    # ----------------------------
    # Train scoring function
    # ----------------------------
    def train_scoring(self, **kwargs) -> Path:
        input_dir = kwargs.get("input_dir") or self.distances_dir
        if input_dir is None:
            raise ValueError("No input_dir provided and distances not extracted yet.")

        output_dir_path = Path(kwargs.get("output_dir", "training_output"))
        script_path = self.base_dir / "src" / "train.py"
        if not script_path.exists():
            raise FileNotFoundError(f"train.py not found at {script_path}")

        cmd = [
            sys.executable, str(script_path),
            "--input-dir", str(input_dir),
            "--output-dir", str(output_dir_path),
            "--max-score", str(kwargs.get("max_score", 10.0)),
            "--pseudocount", str(kwargs.get("pseudocount", 0.0)),
            "--cutoff", str(kwargs.get("cutoff", 20.0)),
            "--bin-width", str(kwargs.get("bin_width", 1.0)),
            "--method", kwargs.get("method", "histogram")
        ]
        subprocess.run(cmd, check=True)
        self.training_dir = output_dir_path
        self.scoring_tables_dir = output_dir_path

        # Load scoring tables
        self.scoring_tables = {}
        for f in glob.glob(str(output_dir_path / "score_*.csv")):
            self.scoring_tables[Path(f).stem] = pd.read_csv(f)
        return output_dir_path

    # ----------------------------
    # Score RNA structure
    # ----------------------------
    def score_structure(self, **kwargs) -> Optional[pd.DataFrame]:
        tables_dir = kwargs.get("tables_dir") or self.scoring_tables_dir
        if tables_dir is None:
            raise ValueError("No scoring tables provided and model not trained yet.")

        pdb_path = kwargs.get("pdb_path")
        if pdb_path is None:
            raise ValueError("pdb_path is required for scoring.")

        output_csv = kwargs.get("output_csv")
        script_path = self.base_dir / "src" / "score_structure.py"
        if not script_path.exists():
            raise FileNotFoundError(f"score_structure.py not found at {script_path}")

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

    # ----------------------------
    # Full workflow
    # ----------------------------
    def run_workflow(self, **kwargs):
        """
        Run the full workflow: extract → train → score
        Memoizes both paths and loaded tables.
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
        """Callable shortcut to score a given structure ONLY IF already trained"""
        return self.score_structure(**kwargs)
