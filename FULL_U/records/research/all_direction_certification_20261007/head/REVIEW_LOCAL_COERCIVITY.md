# Independent review: local coercivity weighting

Reviewed source: `research/all_direction_certification_20261007/real220_local_coercivity.py`.
SHA256 at review: `a72027bf891465c16862b5716962b1930ca53018880a0f11dc34c4cf406e0041`.

This review concerns the new mathematical weighting and coverage logic. The
inherited unweighted real-density formula was reviewed separately in
`REVIEW_REAL220.md`. The new module is presently unused. Its descriptor
assembly is not yet validated, because the semantic proof bindings below
remain to be implemented against the finalized input schemas.

For every spin and energy where V/p <= -kappa < 0,

```
H(U,V) = sqrt(U^2+V^2)+V <= U^2/(-2V) <= U^2/(2*kappa*p).
```

Consequently multiplying each previously reviewed positive real-norm
summand by the exact positive rational weight1/(2*kappa) is valid. The
new `density_weighted` is a textual copy of the reviewed finite density
with only this additional factor in the spin sum. Its unit-weight overlap
checks are useful implementation diagnostics, not a substitute for this
algebraic argument.

With beta^2=1-4/mu and T=ell(ell+1),

```
tau = T/(mu-4) = T*(1-beta^2)/(4*beta^2).
```

This is decreasing on beta>0. A closed beta cell[lo,hi] therefore maps to
[tau(hi),tau(lo)]. The compact boundarymu49 is exactly beta^2=45/49.
On the UV portion, clipping its upper tau endpoint byT/45 is correct.
The unbounded tau portion beyond256 can occur in the UV domain only for
even ell>=108. If the whole cell lies beyond256, no finite-tau lookup is
needed; if it crosses256, the code includes both the large-tau margin and
the finite-cell margins. At beta1 the finite UV lookup includes tau0.

The lookup takes the smallest kappa over every certified closed tau cell
touching the mapped interval, including both sides of an exact boundary.
This is conservative and valid; it interpolates no sample values. The
finite UV loader checks consecutive rational cells from0 through
min(256,T/45), paired replay endpoints, negative directed upper bounds,
and every even spin42..448 exactly once. Taking the worse upper endpoint
from the original proof and its replay preserves validity.

The real integration starts with an exact partition, bisects rational
cells, and verifies exact final coverage of[0,1]. At higher precision it
recomputes each cell with the same exact coercivity weights. Retaining the
maximum of the original and replay upper masses is conservative. Float
heap priorities influence only which interval is split; exact Fractions
and directed Arb endpoints determine the reported integral.

## Unresolved descriptor requirements

Before execution or acceptance, bind the descriptor's
`compact_kappa_exact` and `large_tau_kappa_exact` to the actual negative
bounds proved by their referenced records. The reviewed code only checks
status and witness SHA for these two proofs; a descriptor could otherwise
claim a larger margin or a broader domain than its input proves.

The checks must verify at least:

* The same frozen left witness and exact(-1,0) direction.
* Complete even spin42..448 coverage for each applicable region.
* Compact proof coverage4<mu<=49, with descriptor kappa no greater than
  the proved negative margin.
* Large-tau proof coveragemu>=49 and tau>=256 for every applicable spin,
  with its descriptor kappa no greater than the proved margin.
* The finite UV proof's exact tau convention and mu>=49 domain, as well
  as consistency of record spins with top-level spin lists and matching
  record counts in each original/replay pair.
* Referenced proof and source/code hashes, with finalized schemas tied
  to their actual proving/replay implementations.

The parent explicitly agreed these bindings are required and will be
implemented before this module is executed or accepted. No complete
weighted-tail certificate follows from this review alone.


## Follow-up: semantic bindings resolved

The revised module SHA256
`ff62a3009941936af71f5556ad6f0a3854d571888aad9dd1d77aa6defd896e23`
was reviewed after the parent added the descriptor checks. It now binds
exact objective(-1,0), compact minimum even spin42, compact energy4..49
with the lower endpoint open, and the selected compact kappa to the
completed compact proof's own certified margin. The unbounded-tau proof
must covermu49..infinity, tau256..infinity and a starting spin no larger
than108; its margin is checked against its directed upper endpoint.

For finite UV input it requires the actual `PASS_LEFT_FINITE_UV` schema,
matching objective and witness, empty failure lists, strictly higher replay
precision, a matching reference digest, and exact matching record/spin
lists. Thus the original and replay cannot silently truncate their paired
record lists. The source-moment manifest is read and its file hashes are
also verified.

The finite UV schema does not separately serialize its tau convention.
I inspected the hash-bound producer and field implementation:
`certify_left_finite_uv.py` sets the exact endpointmin(256,T/45), and
`finite_uv_bridge220.py` definesT=ell(ell+1), w=tau/T and accepts only
0<=tau<=min(256,T/45). Herew=1/(mu-4), so this is exactly the convention
andmu>=49 domain used by the weighted integration. This implementation
binding resolves the earlier semantic concern for these specific proof
files; unrelated proof schemas must not be substituted merely because
they have similarly named fields.

No remaining defect was identified in the new weighting/coverage or its
finalized descriptor bindings. Acceptance remains conditional on all
referenced proof files actually being complete, passing their own
arithmetic replays, and retaining the asserted source/code hashes. At the
time of this follow-up, the parent was still waiting for those final files;
this note is not a claim that weighted integration had already finished.
