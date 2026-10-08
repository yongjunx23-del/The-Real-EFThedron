"""Angleupper95 adapters: only paths and the exact second objective change."""
from angleupper95_campaign import *
import uniform_source_uv as us
import finite_uv_bridge220 as fv
import tail_finite_uv_region as region
from matrix_crossing_series import shifted_certificate as exact_shifted_certificate

def four_panel_shifted_certificate(poly,cutoff,subdivisions):
 # This constructor precondition proves the same polynomial inequality.
 # Four exact adjacent z-panels sharpen the sufficient Bernstein test;
 # no polynomial coefficient, kernel, or remainder is changed or dropped.
 assert cutoff==10 and subdivisions==1
 result=exact_shifted_certificate(poly,cutoff,4)
 assert result is not None
 return result

region.shifted_certificate=four_panel_shifted_certificate

class Uniform(us.UniformSource):
 def __init__(self,L=450):
  us.LEFT=WITNESS;us.FOLDER=OUT;super().__init__(L)
 def field(self,lo,hi,elo,ehi):return super().field(lo,hi,elo,ehi)+us.ar(OBJECTIVE[1])

class Finite(fv.FiniteUV):
 def __init__(self):
  fv.WITNESS=WITNESS;fv.OUT=OUT;super().__init__()
 def field(self,l,lo,hi):return super().field(l,lo,hi)+fv.ar(OBJECTIVE[1]-1)
