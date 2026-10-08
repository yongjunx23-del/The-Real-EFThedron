# Review of the arbitrary-direction finite-head driver

Reviewed `certify_head.py` on 7 October 2026. The driver calls the unchanged,
previously reviewed finite220/centered/quartic/gap/fast evaluators. The exact
direction is taken from the witness model and recorded as rationals; it is
not replaced by a hardcoded upper-g2 or negative-g0 direction.

The initial knots include the complete threshold/source gap endpoints and a
rational beta cutoff whose corresponding energy is checked to be at least
the analytic UV start. Excluding only intervals wholly inside the explicit
gaps, followed by exact partition assertions, covers beta from zero to one.
Adaptive bisection replaces a parent by two children without changing this
coverage. Every checkpoint is written only after checking the complete
partition. Nonfinite ordinary-cell evaluations are subdivided or raise an
error; they are never counted as zero.

Every retained `trapezoid` value is an interval enclosure obtained from the
reviewed support/error machinery. Heap priority uses floating estimates only
to choose subdivision order, not to establish an upper bound. Interval sums
and updates remain directed; dependency loss can widen the result but cannot
justify a smaller invalid upper. The recorded dyadic endpoint is taken from
the upper endpoint of the final total. Quality targets do not affect validity.

Replay verifies the exact same witness, source manifest, implementation hashes,
direction, ordinary partition and analytic gap intervals. It reevaluates every
ordinary cell and every gap. The complete final assembler must additionally
check that replay precision exceeds production precision (the driver records
both but does not enforce that relationship itself), verify all source moment
and remainder metadata, and include all even spins before adding a separate
same-witness omitted-spin certificate. Neither one completed head nor an
adaptive checkpoint is a complete coupling bound.

Checkpoint records can change while a worker is active. Select only final
records or explicitly copy and hash a fixed checkpoint before using it in a
downstream certificate. Record the actual enclosure when a quality target is
missed; do not substitute the requested target. Replays may be slightly wider;
retaining the larger production/replay endpoint is conservative.

No new formula or domain-coverage defect was found in this orchestration
review. The underlying analytic assumptions and shared-kernel implementation
limitations remain those of the previous certificates. This is internal
review, not formal verification or a second independent kernel implementation.
