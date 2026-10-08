"""Prepare independent fixed-witness heads using the unchanged reviewed kernels."""
from pathlib import Path
from decimal import Decimal
from fractions import Fraction
import argparse, hashlib, json, os, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[3]
CONFIGS = {
    '261p3814': dict(sha='5bd1b99e5fd4e237fe7e98583aa683ab076d4d76319b05b06c023e583f48383e',
        numerical='head384/angle261p3814.json', slack={0: '.0008', 2: '.00022', 4: '.00021'}, other='.000015'),
    '268': dict(sha='1a0754b19846867214d914c595cdc744f433abfcdaff309d645677bbe2602044',
        numerical='head384_268/angle268.json', slack={0: '.0004', 2: '.00014', 4: '.00008'}, other='.00001'),
}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--campaign', choices=CONFIGS, required=True); a = ap.parse_args()
    c = CONFIGS[a.campaign]
    witness = ROOT/f'output/boundary_length_20261007/search/angle{a.campaign}.npz'
    out = ROOT/f'output/boundary_length_20261007/head/angle{a.campaign}'
    numerical_path = ROOT/'output/boundary_length_20261007'/c['numerical']
    assert sha(witness) == c['sha']
    assert not out.exists(), 'Use a fresh campaign directory.'
    out.mkdir(parents=True)
    env = dict(os.environ); env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    numerical = json.loads(numerical_path.read_text()); assert numerical['sha256'] == c['sha']
    rows = []
    for ell in range(0, 41, 2):
        slack = c['slack'].get(ell, c['other']); value = str(numerical['per_spin'][str(ell)])
        rows.append(dict(spin=ell, model='centered' if ell < 4 else 'fast', numerical=value,
                         slack=slack, target=str(Decimal(value) + Decimal(slack))))
    budget = dict(witness_sha256=c['sha'], numerical_record=str(numerical_path.relative_to(ROOT)),
        numerical_record_sha256=sha(numerical_path), rows=rows,
        slack_sum_exact=str(sum((Fraction(r['slack']) for r in rows), Fraction(0))),
        scope='Quality targets only. Certification uses the independently computed directed enclosures.')
    (out/'budgets.json').write_text(json.dumps(budget, indent=2))
    started = time.time(); runs = []
    for label in ['moments128', 'moments128_replay']:
        command = [sys.executable, str(ROOT/'research/all_direction_certification_20261007/head/build_moments.py'),
            '--witness', str(witness), '--out', str(out/label), '--lmax', '40', '--N', '128', '--Nuv', '32', '--bits', '256']
        tic = time.time(); logfile = out/f'{label}.log'
        with logfile.open('w') as fh: subprocess.run(command, cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT, check=True)
        runs.append(dict(label=label, command=command, log=str(logfile.relative_to(ROOT)), elapsed_seconds=time.time()-tic))
        print(label, 'finished', runs[-1]['elapsed_seconds'], flush=True)
    manifest = json.loads((out/'moments128/manifest.json').read_text()); checks = []
    for table in manifest['tables']:
        for suffix, key in [('_moments.json', 'moment_sha256'), ('_uv_moments.json', 'uv_sha256')]:
            name = table['label']+suffix
            assert sha(out/'moments128'/name) == sha(out/'moments128_replay'/name) == table[key]
            checks.append(dict(file=name, sha256=table[key], match=True))
    command = [sys.executable, str(Path(__file__).with_name('preflight.py')), '--witness', str(witness),
               '--moments', str(out/'moments128'), '--out', str(out/'preflight.json')]
    tic = time.time(); logfile = out/'preflight.log'
    with logfile.open('w') as fh: subprocess.run(command, cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT, check=True)
    runs.append(dict(label='preflight', command=command, log=str(logfile.relative_to(ROOT)), elapsed_seconds=time.time()-tic))
    preflight = json.loads((out/'preflight.json').read_text())
    result = dict(status='PASS_PREPARED_HEAD_INPUTS' if preflight['status']=='PASS_PREFLIGHT' else 'PREFLIGHT_FAILED',
        witness_sha256=c['sha'], source_rebuild_checks=checks, runs=runs, elapsed_seconds=time.time()-started,
        preparation_script_sha256=sha(Path(__file__)),
        worker_environment={k:env[k] for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']},
        scope='Source tables and listed finite-width caps only; not an integrated head or coupling bound.')
    (out/'preparation.json').write_text(json.dumps(result, indent=2)); print(json.dumps(result, indent=2), flush=True)

if __name__ == '__main__': main()
