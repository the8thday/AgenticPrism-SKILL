"""Shared solve_ivp primitive for declared surface-binding mechanisms.

Internal concentration nM, response surface units, time seconds. Not a PK model.
"""
import numpy as np
from scipy.integrate import solve_ivp


def integrate(rhs,initial,times,*,rtol=1e-7,atol=1e-9,method='LSODA'):
 times=np.asarray(times,float)
 if len(times)==0:return np.empty((0,len(initial)))
 if np.any(times<0) or np.any(np.diff(times)<=0):raise ValueError('ODE times must be nonnegative and increasing')
 if times[-1]==0:return np.asarray(initial,float)[None,:]
 sol=solve_ivp(rhs,(0,float(times[-1])),np.asarray(initial,float),t_eval=times,rtol=rtol,atol=atol,method=method)
 if not sol.success or not np.isfinite(sol.y).all():raise ArithmeticError('ODE integration failed: '+sol.message)
 return sol.y.T


def response(model,t,c,duration,p,*,gamma=.1,second_rebinding=True,rtol=1e-7,atol=1e-9):
 """Surface response; first-arm bivalent ka has the statistical factor two."""
 t=np.asarray(t,float);out=np.zeros(len(t));active=t>=0
 if not active.any():return out
 tt=t[active];a=tt[tt<duration];b=tt[tt>=duration]-duration
 if model in ('one_to_one','one_to_one_drift','heterogeneous_ligand'):
  from .kinetics_fit import unit_response
  ans=p['rmax']*unit_response(tt,c,duration,p['ka'],p['kd'])
  if model=='heterogeneous_ligand':ans=p['rmax']*p['fraction']*unit_response(tt,c,duration,p['ka'],p['kd'])+p['rmax']*(1-p['fraction'])*unit_response(tt,c,duration,p['ka2'],p['kd2'])
  if model=='one_to_one_drift':ans=ans+p.get('drift',0)*tt
  out[active]=ans;return out
 if model=='bivalent_analyte':
  initial=[0.,0.]
  def rhs(bulk,rebinding):
   def fun(_,state):
    x1,x2=state;free=p['rmax']-x1-2*x2;first=2*p['ka']*bulk*free-p['kd']*x1
    second=(p['ka2']*x1*free if rebinding else 0)-2*p['kd2']*x2
    return [first-second,second]
   return fun
  signal=lambda values:values[:,0]+values[:,1]
 elif model=='mass_transport':
  initial=[0.,0.]
  def rhs(bulk,_):
   def fun(_,state):
    bound,surface_c=state;binding=p['ka']*surface_c*(p['rmax']-bound)-p['kd']*bound
    return [binding,p['km']*(bulk-surface_c)-gamma*binding]
   return fun
  signal=lambda values:values[:,0]
 else:raise ValueError('Unknown ODE mechanism')
 ta=np.unique(np.r_[a,duration]);assoc=integrate(rhs(c,True),initial,ta,rtol=rtol,atol=atol)
 va=signal(assoc[:len(a)])
 vb=signal(integrate(rhs(0,second_rebinding),assoc[-1],b,rtol=rtol,atol=atol)) if len(b) else np.empty(0)
 out[active]=np.r_[va,vb]
 return out
