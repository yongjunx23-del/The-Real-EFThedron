"""Exact union of independently proved energy slices, with no kernel calls."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/competitive_20261007/compact'
EXPECTED='933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92'
PAIRS=[('5','5.5'),('5.5','6'),('6','6.5'),('6.5','7'),('7','7.5'),('7.5','8'),
       ('8','12'),('12','24'),('24','36'),('36','48'),('48','49')]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    original=None
    for bits in [512,640]:
        records={l:[]for l in range(42,110,2)};manifest={};parts=[];seconds=0
        for a,b in PAIRS:
            suffix='_parallel' if bits==640 and (a,b)==('5','5.5') else ''
            p=OUT/f'slice_{a}_{b}_{bits}{suffix}.json';j=json.loads(p.read_text())
            assert j['status']=='PASS_FINITE_COMPACT_BRIDGE' and j['witness_sha256']==EXPECTED
            assert j['spins']==list(records) and list(map(F,j['mu_domain']))==[F(a),F(b)]
            assert j['precision_bits']==bits and F(j['kappa_exact'])==1
            if bits==640:
                assert j['every_original_spin_interval_replayed']
                assert j['replay_input_sha256']==sha(OUT/f'slice_{a}_{b}_512.json')
            manifest[str(p.relative_to(ROOT))]=sha(p);parts.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
            for key,h in j['sha256_manifest'].items():
                assert sha(ROOT/key)==h,key
                if key in manifest:assert manifest[key]==h
                manifest[key]=h
            seconds+=j['seconds']
            for l,rows in j['cells_by_spin'].items():records[int(l)].extend(rows)
        partitions={}
        for l,rows in records.items():
            rows.sort(key=lambda c:F(c['a']));intervals=[(F(c['a']),F(c['b']))for c in rows]
            assert intervals[0][0]==5 and intervals[-1][1]==49
            assert all(a<b for a,b in intervals)
            assert all(a[1]==b[0]for a,b in zip(intervals,intervals[1:]))
            assert all(F(c['V_over_p_upper_exact'])<-1 for c in rows)
            partitions[l]=intervals
        if bits==512:original=partitions
        else:assert partitions==original
        manifest[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
        out=dict(status='PASS_FINITE_COMPACT_BRIDGE',witness_sha256=EXPECTED,spins=list(records),
            mu_domain=['5','49'],kappa_exact='1',precision_bits=bits,cells_by_spin={str(l):v for l,v in records.items()},
            worst_upper_exact=str(max(F(c['V_over_p_upper_exact'])for rows in records.values()for c in rows)),
            replay_input_sha256=sha(OUT/'compact_bridge_sliced_512.json')if bits==640 else None,
            every_original_spin_interval_replayed=bits==640,sha256_manifest=manifest,parts=parts,
            summed_worker_seconds=seconds,scope='Exact gapless union of eleven disjoint-energy interval proofs; seam points are included in both adjacent closed intervals.')
        (OUT/f'compact_bridge_sliced_{bits}.json').write_text(json.dumps(out,indent=2))
        print('PASS',bits,sum(map(len,records.values())),float(F(out['worst_upper_exact'])))
if __name__=='__main__':main()
