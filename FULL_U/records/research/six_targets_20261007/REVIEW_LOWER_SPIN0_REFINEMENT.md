# Exact bisection refinement for lower15/lower45

Parent-agent internal source review, 2026-10-07. Read assemble_spin0_refined.py and the source adaptation description against the reviewed lower mixture formula. This is not external review or proof-assistant verification.

The new producer bisects each ordinary spin0 mixture interval at its exact rational midpoint. Gap intervals, their costs, witness definitions, weights, all220 coefficients and physical field formulas are unchanged. Every child interval lies in the same containing base cell, so the existing full-cell fallback remains a conservative nonnegative bound; it is not divided by2. The normalized shortcut keeps hi<1/2. Both256/384-bit runs enclose the entire refined partition.

The wrapper starts from an authenticated complete original target packet. It reconstructs every required midpoint and demands exact equality of the two new partitions to those children. It checks unchanged exact gaps and definitions, all legal methods and safe domains, exact fallback costs and exact sums. The replacement is max of the two direct upper enclosures and must be smaller than the old retained spin0 component. It subtracts precisely that component, preserves every other head and inherited tail cost, and rechecks the exact total and outward norm conversion. Initial complete packets are preserved unchanged.

This refinement is acceptable after both new runs and the complete original packets pass the asserted checks. It improves numerical enclosure quality for the same witness; it adds no physical constraint and performs no new multiplier optimization.
