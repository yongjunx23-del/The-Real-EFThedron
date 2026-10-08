"""Preserve the failed reader and create a list-root-safe metadata reader."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
CODE = ROOT / 'research/six_targets_20261007'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    old = CODE / 'assemble_uv_upper95.py'
    new = CODE / 'assemble_uv_upper95_v2.py'
    substitutions = [
        ("for q,hh in d.get('sha256_manifest',{}).items():ingest(ROOT/q,hh)",
         "for q,hh in (d.get('sha256_manifest',{}).items()if isinstance(d,dict)else []):ingest(ROOT/q,hh)"),
        (" ingest(outdir/'exact_rebuild_and_algebra_audit.json')",
         " ingest(ROOT/'output/six_targets_20261007/upper95_direct/tail/metadata_reader_fix.json')\n"
         " ingest(ROOT/'research/six_targets_20261007/REVIEW_UPPER95_UV_METADATA_FIX.md')\n"
         " ingest(outdir/'exact_rebuild_and_algebra_audit.json')"),
    ]
    text = old.read_text()
    for before, after in substitutions:
        assert text.count(before) == 1
        text = text.replace(before, after)
    assert not new.exists()
    new.write_text(text)

    driver_old = CODE / 'run_tail95_completion.py'
    driver_new = CODE / 'run_tail95_completion_v2.py'
    driver_substitutions = [
        ('assemble_uv_upper95.py', 'assemble_uv_upper95_v2.py'),
        ("label+'.log'", "label+'_v2.log'"),
        ('tail_completion_commands.jsonl', 'tail_completion_commands_v2.jsonl'),
        ('tail_completion_campaign.json', 'tail_completion_campaign_v2.json'),
    ]
    text = driver_old.read_text()
    for before, after in driver_substitutions:
        assert before in text
        text = text.replace(before, after)
    assert not driver_new.exists()
    driver_new.write_text(text)
    ledger = dict(
        status='NARROW_METADATA_ROOT_TYPE_FIX_NO_PROOF_OR_FORMULA_CHANGE',
        changes=[dict(original=str(old.relative_to(ROOT)), replacement=str(new.relative_to(ROOT)),
                      exact_substitutions=substitutions),
                 dict(original=str(driver_old.relative_to(ROOT)), replacement=str(driver_new.relative_to(ROOT)),
                      exact_substitutions=driver_substitutions)],
        reason=('The original reader stopped before creating any UV packet when the independently '
                'reviewed adapter-diff JSON had a list root. List-root metadata is authenticated by '
                'its whole-file hash; dictionary-root sha256_manifest mappings still recurse unchanged.'),
        failed_log='output/six_targets_20261007/upper95_direct/uv_assembly.log',
        sha256_manifest={str(p.relative_to(ROOT)): sha(p) for p in [old, new, driver_old, driver_new,
            Path(__file__), ROOT/'output/six_targets_20261007/upper95_direct/uv_assembly.log',
            ROOT/'output/six_targets_20261007/upper95_direct/tail_completion_commands.jsonl']},
    )
    output = ROOT/'output/six_targets_20261007/upper95_direct/tail/metadata_reader_fix.json'
    assert not output.exists()
    output.write_text(json.dumps(ledger, indent=2) + '\n')
    print(sha(new))


if __name__ == '__main__':
    main()
