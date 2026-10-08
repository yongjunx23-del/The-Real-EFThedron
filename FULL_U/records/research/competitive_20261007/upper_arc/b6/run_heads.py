"""Exact b6 heads: old150 functional plus107/5 X4, objective(1,6)."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse,sys
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'research/full_direction_20261007'))
from head_family_extended import run
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--spin',type=int,required=True)
    p.add_argument('--precision',type=int,default=280);p.add_argument('--subdivide',type=int,default=1);a=p.parse_args()
    run(a.spin,['6'],a.precision,ROOT/f'output/competitive_20261007/upper_arc/b6/head_p{a.precision}_spin{a.spin}.json',a.subdivide,Fraction(107,30))
