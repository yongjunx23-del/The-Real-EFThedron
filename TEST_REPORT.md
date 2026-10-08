# Staging verification report

All three endpoint certificates passed the final default verifier on macOS,
using Python 3.12 and python-flint 0.9.0. The run took about 20 seconds;
that timing is environment-dependent. No optimizer was run.

The checks reproduced all three exact rational support constants and their
outward decimal values:

- C08: `g0 <= 2.705955`.
- C20: `g2 <= 1.052351`.
- C21: `g0 >= -7.751502`.

Checked data included 590 exact multipliers, 1,008 original dependency files,
44 low-spin head partitions comprising 99,965 ordinary cells and 220 endpoint
caps, and the separate C08 refined direction partitions. C08's finite-block
bound was rebuilt by the original directed interval merge at slope zero.

The original C20 sign assemblers checked 43,948 compact bridge cells,
68,030 base UV finite-spin cells, 90,287 stronger UV finite-spin cells,
and the stated analytic domains for the remaining spins and energies.
C21's original coercivity checker validated 191,182 finite-spin UV cells
and their paired precision records for every even spin 42–448.

This report concerns saved-certificate verification and exact assembly.
The optional `--fresh-local` calculation is a separate test, and a complete
fresh replay of all interval integrals was not performed in this staging run.
The enclosing archive's final release manifest is generated after staging;
its validation is performed when the delivered package is checked.
