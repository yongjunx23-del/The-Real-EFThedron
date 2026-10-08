# Upper55 source and gap-coarsening review

Internal parent-agent review, 2026-10-07. This reviews the adaptation and exact accounting, not an external or formal proof.

The parent read upper55_head.py, upper55_common.py and the complete assembler. `review_source_diffs.py` reconstructs the assembler from the sealed extra-upper assembler and the declared replacement ledger; equality is exact. Only base paths, general two-vector rational weights, and authenticated gap coarsening change the mathematical inputs.

The weights solve the exact 2x2 objective equation and must be positive. All 220 multipliers use the same exact weights. Every unchanged component is bounded by positive homogeneity and subadditivity of H(U,V)=sqrt(U²+V²)+V. Combined spin0/2 integrals may improve on this convex bound.

The repaired83 spin2 cap ends at beta=1/1024, whereas the85 cap ends at1/256. `align_gaps` moves exactly64 existing, contiguous, nonnegative cell enclosures covering [1/1024,1/256] into the83 gap cost. It adds their entire exact dyadic costs, removes precisely those ordinary cells, checks the remaining total is bounded by the original head component, and makes the two gap lists identical. No new kernel approximation or extrapolation is used. The assembler subtracts the original authenticated selected83 head component, not an invented coarsened cost.

The cell formulas and source activity guard are unchanged from the corrected extra-upper implementation: beta<1/2 universally, or for spin2 only hi²<=79/119 for the normalized shortcut. Ordinary cells and analytic gaps cover [0,1] exactly; both256/384-bit integrations must agree on partition and witness. The retained direct component is max of the pair, optionally capped by the valid original convex cost. Exact remaining costs plus the new components are outward rounded; normal length uses a rigorous rational lower bound.

Source review accepts this adaptation. Successful complete output still depends on paired runs, exact decomposition equality and source/moment authentication. Shared-code replay is not an independent mathematical implementation.
