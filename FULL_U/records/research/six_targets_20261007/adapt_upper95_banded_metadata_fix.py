"""Create a narrow list-root reader repair without changing saved proofs."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
CODE = ROOT / 'research/six_targets_20261007'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    change = ("for name,expected in d.get('sha256_manifest',{}).items():",
              "for name,expected in (d.get('sha256_manifest',{}).items()if isinstance(d,dict)else []):")
    changes = []
    files = [Path(__file__)]

    def adapt(old_name, new_name, substitutions):
        old, new = CODE / old_name, CODE / new_name
        text = old.read_text()
        for before, after in substitutions:
            assert before in text
            text = text.replace(before, after)
        assert not new.exists()
        new.write_text(text)
        changes.append(dict(original=str(old.relative_to(ROOT)), replacement=str(new.relative_to(ROOT)),
                            exact_substitutions=substitutions))
        files.extend([old, new])

    adapt('banded_upper95.py', 'banded_upper95_v2.py', [change])
    adapt('review_upper95_tail.py', 'review_upper95_tail_v2.py', [
        ('research/six_targets_20261007/banded_upper95.py',
         'research/six_targets_20261007/banded_upper95_v2.py'),
        (".replace('complete_compact_coercivity_stronger.json','complete_compact_coercivity.json')==new.read_text()",
         ".replace('complete_compact_coercivity_stronger.json','complete_compact_coercivity.json')"
         + '.replace(' + repr(change[0]) + ',' + repr(change[1]) + ')==new.read_text()'),
        ("for name in ['assemble_uv_peak.py','review_upper_tail.py']:",
         "read(BASE/'tail/banded_metadata_reader_fix.json')\n"
         "for name in ['assemble_uv_upper95_v2.py','review_upper95_tail_v2.py','REVIEW_UPPER95_BANDED_METADATA_FIX.md']:"),
    ])
    adapt('run_tail95_completion_v2.py', 'run_tail95_completion_v3.py', [
        (" run('uv_assembly',[CODE/'assemble_uv_upper95_v2.py','--angle','upper95'])\n", ''),
        ('banded_upper95.py', 'banded_upper95_v2.py'),
        ('review_upper95_tail.py', 'review_upper95_tail_v2.py'),
        ("label+'_v2.log'", "label+'_v3.log'"),
        ('tail_completion_commands_v2.jsonl', 'tail_completion_commands_v3.jsonl'),
        ('tail_completion_campaign_v2.json', 'tail_completion_campaign_v3.json'),
    ])
    files.extend([ROOT/'output/six_targets_20261007/upper95_direct/banded_cost_v2.log',
                  ROOT/'output/six_targets_20261007/upper95_direct/tail_completion_commands_v2.jsonl'])
    ledger = dict(status='NARROW_BANDED_METADATA_ROOT_TYPE_FIX_NO_MATHEMATICAL_CHANGE',
                  changes=changes,
                  explanation=('The original local-cost reader encountered the same list-root metadata '
                               'already handled in the UV v2 reader. Only its manifest-traversal root-type '
                               'guard changes. The independent coverage checker confirms this exact '
                               'source substitution and otherwise repeats every original exact test. '
                               'The completed UV packets are verified again, never overwritten.'),
                  sha256_manifest={str(p.relative_to(ROOT)): sha(p) for p in files})
    target = ROOT/'output/six_targets_20261007/upper95_direct/tail/banded_metadata_reader_fix.json'
    assert not target.exists()
    target.write_text(json.dumps(ledger, indent=2) + '\n')
    print(sha(CODE/'banded_upper95_v2.py'))


if __name__ == '__main__':
    main()
