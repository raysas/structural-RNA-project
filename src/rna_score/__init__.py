# RNA scoring library - orchestrates src/ modules
from utils.structure_io import FastParser, OnlineFetcher
from . import lib
from .lib import download_rna_structures, RNAScorer

__version__ = "0.1.0"

__all__ = [
    "lib",
    "download_rna_structures",
    "RNAScorer",
    "FastParser",
    "OnlineFetcher",
]