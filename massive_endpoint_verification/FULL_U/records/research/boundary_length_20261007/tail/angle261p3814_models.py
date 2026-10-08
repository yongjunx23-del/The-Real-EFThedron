"""Angle261p3814 adapters: only paths and the exact second objective change."""
from angle261p3814_campaign import *
import uniform_source_uv as us
import finite_uv_bridge220 as fv

class Uniform(us.UniformSource):
 def __init__(self,L=450):
  us.LEFT=WITNESS;us.FOLDER=OUT;super().__init__(L)
 def field(self,lo,hi,elo,ehi):return super().field(lo,hi,elo,ehi)+us.ar(OBJECTIVE[1])

class Finite(fv.FiniteUV):
 def __init__(self):
  fv.WITNESS=WITNESS;fv.OUT=OUT;super().__init__()
 def field(self,l,lo,hi):return super().field(l,lo,hi)+fv.ar(OBJECTIVE[1]-1)
