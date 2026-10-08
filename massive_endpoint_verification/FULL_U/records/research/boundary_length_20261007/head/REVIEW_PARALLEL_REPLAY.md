# Internal review of exact partition replay in chunks

Reviewed `parallel_replay.py` and `assemble_parallel.py` on 2026-10-07.
The method changes scheduling and interval summation, not the mathematical
integrand or its enclosure. The original complete replay remains a separate
record; no result is overwritten.

Each worker validates the immutable original record, witness, source table
manifest, dependency hashes, exact all-energy partition, and higher working
precision. It re-evaluates a half-open range of saved cell indices with the
unchanged model. Only the range beginning at index zero includes the source
endpoint and threshold gap contributions. All workers independently verify
that the model returns the same exact gap geometry.

The combiner checks common witness/objective/source/record/precision metadata,
sorts the ranges, and requires an exact disjoint cover of every cell index.
It checks that all gaps occur once, then sums the directed dyadic upper
endpoints as exact fractions. A sum of upper bounds remains an upper bound;
the exact sum introduces no new rounding error. It binds every chunk file
and all implementation dependencies into the combined replay record.

The alternate head assembler selects only the explicitly named parallel
record where it exists, retains the maximum of producer and replay upper
bounds, and saves a separate summary. The complete directional assembler
also rechecks the chunk-index cover and exact sum directly. No domain or
summation defect was identified. This is an internal review, not an external
or formal verification, and a single-spin replay still needs all other
spins and tail components before constituting a full directional bound.
