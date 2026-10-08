# Upper 95% UV assembler metadata repair review

The root reviewer checked the complete replacement source and independently reconstructed both replacement files from every exact substitution in `output/six_targets_20261007/upper95_direct/tail/metadata_reader_fix.json`. Every ledger input SHA-256 matched. The replacement assembler SHA-256 is `0c1eaf8cbe9995bfdf020d2300caba93cf4cae2e403957ae7d49e95beb4ce777`.

The sole reader change applies `sha256_manifest` recursion only to dictionary-root JSON. A list-root metadata document remains authenticated by its entire file hash. Dictionary manifests retain the same recursive authentication, duplicate-hash consistency checks, and expected-hash checks. The additional inputs bind this repair ledger and review. The wrapper selects the new assembler and fresh log paths; it changes no proof commands or numerical parameters.

The failure occurred while reading the adapter-diff metadata list, before a UV completion packet was written. The original source and failed log are preserved. No numerical kernel, coverage partition, coercivity threshold, paired-precision rule, or support-integral formula changed. All finite and infinite UV coverage assertions remain in place. This review approves the metadata repair, not an independent re-execution of the interval kernels; the completed assembly and readback must still pass.
