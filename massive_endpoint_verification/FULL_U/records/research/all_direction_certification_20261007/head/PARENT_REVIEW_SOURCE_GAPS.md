# Independent parent review of finite source-gap caps

The parent agent independently reviewed `finite_gap.py` after the gap
implementation was frozen and the final component jobs had started.
Its findings were:

* The omitted normalized logarithmic coefficient is exactly
  `Q^3 H/(width Rstar)*W`.
* Factoring `W=(1-y^2)P` and using
  `|y+1| sup |x log|x||` plus the symmetric term is a valid uniform bound.
* The required radius less than0.1 makes `-r log r` monotone for intervals
  containing zero, with the continuous zero value at the singularity.
* The Roy source polynomial is identically zero above maximum local spin6,
  and the code asserts this fact.
* Real source terms united with zero correctly enclose the inactive portion
  of a gap crossing a source window endpoint.
* No defect was identified in this new endpoint treatment. The inherited
  overall kernel normalization remains supported by the earlier reviews and
  independent field comparisons; this note is not a separate full derivation.

This review concerns the formula and interval treatment. Full component
records still require all cells to finish, and a complete coupling bound
requires the remaining finite spins and a matching all-spin tail.
