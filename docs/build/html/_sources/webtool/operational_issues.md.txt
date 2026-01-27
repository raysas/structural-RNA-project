# Common operational issues (and fixes)

## “My run seems stuck”
- Check **Status** with the Job ID.
- Large inputs + large cutoff can make extraction/training slow; try a smaller test first.

## “Outputs are empty or extremely small”
- Verify the structure actually contains RNA residues and the selected atom (e.g., `C3'`).
- Reduce **sequence separation** (too large can remove many contacts).
- Ensure chains are valid (in paste mode, wrong chain IDs often yield empty extraction).

## “Scores look inconsistent across runs”
- Ensure you did not change cutoff/bin width/atom mode/separation.
- Ensure the training set is the same (different training structures = different model).

## “I want to score many test sets with the same model”
- Run **Extract + Train** once (or a Pipeline run).
- Then repeatedly use the **Score** tab with new structures, keeping parameters identical.