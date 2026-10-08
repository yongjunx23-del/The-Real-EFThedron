"""Rebuild every source moment at higher precision and prove cache enclosure.

All freshly rebuilt intervals must be contained in the corresponding saved
intervals; overlap alone is not accepted. Thus the actual cached moment
enclosures used by the head verifiers are independently validated in full.
"""
from certified_roy_moments import ROOT,build as build_roy,load
from certified_strip_moments import build as build_strip
from certified_uv_moments import build as build_uv,load_uv
from flint import ctx
import hashlib,json,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def compare(old,new,keys2=('M_u','M_s'),keys1=('T','M_t')):
    count=0
    for key in keys2:
        assert len(old[key])==len(new[key])
        for i,(a,b) in enumerate(zip(old[key],new[key])):
            assert len(a)==len(b)
            for j,(cached,fresh) in enumerate(zip(a,b)):
                assert cached.is_finite() and fresh.is_finite(),(key,i,j,'nonfinite')
                assert cached.contains(fresh),(key,i,j,cached,fresh)
                count+=1
    for key in keys1:
        assert len(old[key])==len(new[key])
        for i,(cached,fresh) in enumerate(zip(old[key],new[key])):
            assert cached.is_finite() and fresh.is_finite(),(key,i,'nonfinite')
            assert cached.contains(fresh),(key,i,cached,fresh)
            count+=1
    return count

def main():
    ctx.prec=340;tic=time.time();wit=ROOT/'output/angles_resummed/angle_000p000.npz'
    records=[];fresh64={}
    for N,folder in [(32,ROOT/'output/certification_head'),(64,ROOT/'output/certification_head/N64')]:
        fresh=build_roy(wit,40,N,40);path=folder/'roy_moments.json'
        records.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),'entries':compare(load(path),fresh)})
        if N==64:fresh64['roy']=fresh
        for fresh in build_strip(wit,40,N):
            label='strip_'+fresh['metadata']['A']+'_'+fresh['metadata']['B'];path=folder/(label+'_moments.json')
            records.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),'entries':compare(load(path),fresh)})
            if N==64:fresh64[label]=fresh
    for label,source in fresh64.items():
        fresh=build_uv(source);path=ROOT/'output/certification_head'/(label+'_uv_moments.json')
        records.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),'entries':compare(load_uv(path),fresh,('A',),('B',))})
    out={'status':'PASS','witness_sha256':sha(wit),'precision_bits':ctx.prec,
         'method':'Rebuild every polynomial moment using certified Gauss nodes at higher precision, require full interval containment in saved tables; then rebuild every UV contraction.',
         'entries_verified':sum(r['entries'] for r in records),'seconds':time.time()-tic,'tables':records}
    (ROOT/'output/certification_head/all_moment_tables_verified.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='tables'},indent=2))

if __name__=='__main__':main()
