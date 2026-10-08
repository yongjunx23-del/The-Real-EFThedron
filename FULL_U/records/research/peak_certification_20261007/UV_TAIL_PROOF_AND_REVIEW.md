# UV tail proof for the two peak witnesses

These records bound the dual imaginary field in the region mu >= 49 for every even spin ell >= 42. They do not by themselves bound the support integral or establish a complete directional certificate. The head, compact-energy coercivity, and real-norm integrations must use the identical exact witness.

The two input NPZs are frozen in the previous weak-region campaign. Every binary64 number is interpreted as its exact dyadic rational. The upper83 objective is `(535119985632953/4503599627370496, 8943390007067311/9007199254740992)`. The lower263 objective is `(-2109641626282615/18014398509481984, -8945221711437375/9007199254740992)`. The adapters change paths and the exact second objective in the reviewed common UV formulas. Prepared source manifests bind the witness, source-generation code, exact moments, norms, and exact Roy/crossing positivity audits.

For finite spins 42,44,...,448, set tau=ell*(ell+1)/(mu-4). Each producer proves strict V/p < -kappa on a partition of `[0,min(256,ell*(ell+1)/45)]`; a separate 480-bit run re-evaluates every accepted 400-bit cell. The saved local margin is minus the larger of the two directed upper endpoints. Producer acceptance uses kappa=1/200 for upper83 and 1/50 for lower263. These acceptance constants need not be the best local margins and are not substituted for the stronger retained bounds when computing energy-dependent real costs.

For all even ell>=450, eta=1/(ell*(ell+1)) lies in `[0,1/(450*451)]`. Exact bivariate source and crossing polynomials plus uniform remainders enclose the entire eta interval, including its zero endpoint. Upper83 uses 16 closed eta panels and lower263 one panel; every panel covers tau in `[0,256]` with a full 400/480-bit producer/replay partition. Point screens and the unselected four-panel upper83 probe are diagnostics, not the selected continuum certificate.

For tau>=256, an analytic Jensen bound applies uniformly for ell>=108 and mu>=49. All shifted crossing powers have the required nonpositive Bernstein coefficient bounds. The finite-spin cases ell<=106 have tau<=106*107/45<256, so the two proofs leave no gap. Both analytic computations are evaluated at 320 and 400 bits and require V/p<-1. The lower263 proof uses the previously reviewed absolute source envelope.

## Sharpened upper83 crossed-source envelope

The absolute crossed-source envelope was too loose for upper83. Its failed diagnostic is retained and does not enter a completed packet. Let z=1/(mu-4/3), so 0<=z<=3/143 for mu>=49, and write the source term as S2(z)=sum m_k*z^(2*k), with exact signed moments from `compact_norms.pkl`. The existing signed-measure remainder uses r=(317/60)*(3/143)<1 and source total-variation norm norm2. Therefore

```text
-S2(z) <= -m_0 + sum_{k=1}^{11} |m_k|*(3/143)^(2*k)
           + norm2*r^24/(1-r^2).
```

The producer checks the witness fingerprint, at least 12 moments, and m_0=Hs from the independent joint-moment representation. Its prepared manifest already binds the exact compact moment packet. The signed constant is not counted twice: the contracted forward numerator contains -Hs, the uniform helper adds +Hs-S2, and the large-tau forward envelope has no Hs term. Thus replacing the absolute bound by the displayed bound for -S2 directly is consistent.

The upper-search subagent independently reviewed this source2 change and its convention against `uniform_source_uv.py` and `ContractedRegion.numerator`, and reported formula review PASS. The reviewed producer `tail_upper83/certify_angleupper83_large_tau_moments.py` has SHA-256 `7a0a03110b849011d5c694d4aa0839da56e0f6a8679e135918f73e96c4e58862`. The first 320-bit resulting upper is approximately -48.1094; this numerical display is not a substitute for the exact directed endpoint retained in the proof JSON.

## Source and arithmetic checks

The exact-rational source cache is rebuilt without mutation and every entry compared. An additional algebra audit independently specializes the uniform source and crossing polynomials at ell=450,452,600 and checks exact equality with the finite-spin polynomials, followed by 18 high-precision enclosure overlaps. These audits do not replace the complete interval cover.

The assembler checks witness/objective consistency, paired cell geometry, exact endpoint continuity, all spin and eta indices, strict margins, analytic overlap, and recursively bound source hashes. The real support contribution uses the already reviewed real-norm integral and the identity bounding the disk support by A/(2*kappa); finite-spin integration may use a separately checked minimum of all local margins touching each energy cell. A mathematically valid positive coercivity constant is required before this division.

This is directed interval arithmetic with an analytic treatment of the infinite domains under the stated dispersive representation. Producer/replay arithmetic and internal code reviews share the underlying formulas; they are not proof-assistant verification or independent external replication.
