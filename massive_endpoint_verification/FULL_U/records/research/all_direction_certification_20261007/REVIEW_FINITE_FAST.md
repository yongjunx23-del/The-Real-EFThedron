# Review of the unused contracted finite-head evaluator

The prototype `head/finite_fast.py` was reviewed on 7 October 2026.
It is not part of the completed g2 or b4 certificates and does not replace
any already selected left-head production record.

For q=1/(mu+h), delta=mu-4=[1-(4+h)q]/q. Hence

    sum_j C_j delta^(ell-j) U_j(q) q
      = [q sum_j C_j q^j (1-(4+h)q)^(ell-j) U_j(q)]/q^ell.

This is exactly the `combined_u` polynomial and its normalization. All
exponents ell-j are nonnegative because only j=0..ell enter. The combined
polynomial coefficients are Arb enclosures of the original retained
coefficients. Algebraic rearrangement and interval polynomial arithmetic
therefore enclose the same mathematical polynomial; correlations may
improve or weaken the enclosure but are not treated as independent facts.
The finite-energy normalization branch stays away from q=0.

For x=1/mu and h=x/(1-4x), the UV contraction similarly gives

    sum_j C_j h^j UVa_j(x)
      = [sum_j C_j x^j(1-4x)^(ell-j)UVa_j(x)]/(1-4x)^ell.

The UV denominator stays positive, including at infinite energy x=0.
The extra overall x/16 in the scaled spectral field is unchanged.

The prototype retains the original source geometric error estimates;
it rearranges only the finite retained polynomials. It preserves the
Roy/strip polynomial terms, PV branch selection, logarithmic terms,
active-source real part, transfer derivatives, crossing and g0 subtraction
terms. Threshold caps and endpoint-gap logic are inherited from FiniteGap.
Its exact beta1 endpoint handling agrees with the original evaluator.

No algebraic or error-normalization defect was found. Before selecting this
as a replay implementation, require representative interval and Taylor-jet
comparisons as well as the already reported point checks; in particular,
the normalized PV branch must be selected consistently for the full cell
and midpoint. This document by itself authorizes no completed new bound.
The original selected proof implementations remain the authoritative ones
until any new replay is explicitly assembled and hash-bound.
