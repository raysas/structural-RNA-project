rna-score CLI Documentation
===========================

rna-score is a Python CLI tool to evaluate predicted RNA tertiary structures
based on interatomic distance distributions.

**Goal:** Creation of an objective function for the RNA folding problem.  
Supervised by Professor Guillaume Postic.  

**Team:** Joelle Assy, Yazid Hoblos, Denys Buryi, Raul Duran De Alba, Rayane Adam

.. button-link:: https://rna-score.onrender.com/
   :color: primary
   :expand:

   Try it online

.. dropdown:: Quick install
   :icon: download
   :color: primary

   .. code-block:: bash

      pip install git+https://github.com/raysas/structural-RNA-project.git

.. dropdown:: Example CLI workflow
   :icon: terminal
   :color: secondary

   .. code-block:: bash

      rna-score access -n 50 --rna-only -f cif -o data/rna_structures --workers 4
      rna-score extract --folder rna_structures/mmcif --format mmcif --out-dir dist_data
      rna-score train --input-dir dist_data --output-dir training_output --method histogram
      rna-score score --folder rna_structures/mmcif --tables training_output --format mmcif --output scores.csv
      rna-score plot --input-dir training_output --output-dir plots --combined

.. grid:: 1 2 2 2
   :gutter: 2
   :margin: 2 0 2 0

   .. grid-item-card:: Background
      :link: background/index
      :link-type: doc
      :shadow: sm
      :class-card: sd-border-1 sd-rounded-2

      Learn the scoring rationale, datasets, and assumptions.


   .. grid-item-card:: Library API
      :link: library
      :link-type: doc
      :shadow: sm
      :class-card: sd-border-1 sd-rounded-2

      Explore Python entry points and helper utilities.

   .. grid-item-card:: CLI Usage
      :link: cli
      :link-type: doc
      :shadow: sm
      :class-card: sd-border-1 sd-rounded-2

      Run the pipeline end-to-end from the command line.

   .. grid-item-card:: Webtool
      :link: webtool/index
      :link-type: doc
      :shadow: sm
      :class-card: sd-border-1 sd-rounded-2

      Score structures from the browser without installing anything.

.. toctree::
   :maxdepth: 2
   :hidden:

   background/index
   library/index
   cli
   webtool/index
