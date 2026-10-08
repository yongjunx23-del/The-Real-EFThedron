# Direct95 scheduler adoption review

Internal parent-agent review, 2026-10-07. Independently reconstructed all three sources in upper95_direct/parallel_adaptation_manifest.json using review_parallel95.py. They exactly match the declared transformations from the reviewed86 scheduler and reviewed95 assemblers.

Only scheduling width changes from2 to4; the upper95 budgets, exact witness, moments, producer/replay precision and6500/3000 cell limits remain unchanged. The additional ParentProcessId check permits adoption only of children of the identified ordinary head controller. This correctly excludes the separate threshold-refinement producer/replay. The process-handle preservation and explicitly labelled completed-producer exit-code inference follow the reviewed86 protocol. No scientific output folder may be overwritten or silently reclassified.

The new main assembler adds the same authenticated switch/old-controller/new-controller/budget checks as86, with the exact95 paths and worker counts. It retains the complete all21-spin head checker, source rebuild authentication and every compact/UV/tail requirement. The banded assembler changes its complete global packet and source bindings only. The threshold refinement remains a separately checked exact component substitution after a complete original packet exists.

This adaptation may be executed within the current authorized calculation. It establishes no new physical bound by itself; all unchanged producer/replay, coverage, coefficient and sum checks must finish successfully. Internal source review is not independent external or formal verification.
