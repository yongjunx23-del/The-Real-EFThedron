# Actual large Gram matrices underlying the continuum certificate

The head certificate can be expressed using actual large matrix sum-of-squares Gram matrices. It is not necessary to regard it solely as a collection of positive `2x2` coefficients.

On each compact velocity cell, write the transformed kernel matrix as a degree `n=2m` polynomial in `t in [0,1]`:

\[
\widetilde C(t)=\sum_{j=0}^{2m}\binom{2m}{j}
t^j(1-t)^{2m-j}M_j,
\qquad M_j\succeq0.
\]

For even indices `j=2r`, define

\[
q_j(t)=t^r(1-t)^{m-r}.
\]

For odd indices `j=2r+1`, define

\[
q_j(t)=t^r(1-t)^{m-r-1}.
\]

Let `v_j` be the rational monomial coefficient vector of this polynomial. The following are explicit rational Gram matrices:

\[
Q_0=\sum_{j\text{ even}}\binom{2m}{j}
(v_jv_j^T)\otimes M_j\succeq0,
\]

\[
Q_1=\sum_{j\text{ odd}}\binom{2m}{j}
(v_jv_j^T)\otimes M_j\succeq0.
\]

Their dimensions are respectively `2(m+1)` and `2m`. With

\[
Z_k(t)=\begin{pmatrix}I_2\\tI_2\\\vdots\\t^kI_2\end{pmatrix},
\]

the exact polynomial identity is

\[
\boxed{\widetilde C(t)=Z_m(t)^TQ_0Z_m(t)
+t(1-t)Z_{m-1}(t)^TQ_1Z_{m-1}(t).}
\]

This is the interval matrix sum-of-squares representation. The second term imposes the support interval. Because each constant `M_j` has been checked positive semidefinite exactly, both large matrices are positive semidefinite by an explicit rational sum of tensor products; neither floating-point eigenvalues nor an approximate Cholesky factor are needed.

For the twelve-null certificate, the polynomial degrees range from 44 at spin zero to 80 at spin 18. Thus the local certificates use Gram matrices as large as **82 by 82 and 80 by 80**. There is a separate pair on each velocity interval. They are higher-degree matrix-SOS certificates; they should not be described as a single global Hankel moment matrix.

The exporter `certificate_export_sos_gram.py` writes the full dense rational matrices for a selected cell, together with their positive-semidefinite decomposition. It independently checks the polynomial identity coefficient by coefficient in exact rational arithmetic. The complete calculation stores the more compact Bernstein data for every cell; the displayed formula reconstructs every large Gram pair from that data.

This construction is a sufficient subfamily of general interval matrix-SOS certificates. Optimizing unrestricted `Q0,Q1` could improve efficiency or the bound further, but the explicit matrices already establish valid continuum inequalities. The improvement is not a claim that Gram positivity adds new physics beyond unitarity; it comes from retaining the complete unitarity matrix and combining more valid harmonic and transfer-derivative equations.
