# Internal review of the continuous boundary-length audit

Reviewed `geometry_continuous_audit.py` on 2026-10-07.
This is exact geometry for the digitized reference inputs, not a statement
that the published curve coordinates themselves are exact or certified.

For each segment p+t*d, the code uses an unnormalized normal in the
axis-range metric. All halfplane intersections along the physical normal
ray are therefore affine rational functions of t. Constraints with zero
normal denominator restrict t directly. Every lower/upper pair restricts
the projected domain by lower(t)<=upper(t), which is precisely the condition
that the normal line has a nonempty polygon intersection.

The outer endpoint is the minimum of finitely many affine upper functions.
On a closed projected domain its maximum occurs at an endpoint or an affine
envelope crossing. The code includes all pairwise upper-function crossings,
and verifies each retained active piece at both endpoints. It thus cannot
miss an interior peak between plotted arclength stations. Multiplying the
normal parameter by sqrt(segment metric length squared) gives the actual
normalized gap; signed squared comparisons order these gaps exactly.

Reference-point exclusion is computed separately by exact clipping against
all halfplanes. Clipped parameter fractions multiply the corresponding
metric segment lengths. Positive square roots and the final length ratios
are enclosed with Arb. Normal lines missing the polygon are not silently
given zero gap. Corners use their two one-sided segment normals and receive
no artificial normal-angle length. The unknown right cap is excluded.

The polygon vertex and recession-cone checks agree with the exact rational
halfplanes. No projection, sign, endpoint, or coverage defect was identified.
These graphical diagnostics measure the geometry of completed physical
inequalities but do not create new ones or certify a digitization error bar.
