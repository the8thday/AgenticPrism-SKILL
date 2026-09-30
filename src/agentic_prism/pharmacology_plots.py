"""Saved-result-only pharmacology graphics."""
from pathlib import Path
import numpy as np

def render(run,result,style):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from .report import font_setup
 from . import plot_style as pstyle
 token=pstyle.THEMES[style];rc=pstyle.rc(token,*font_setup())
 out=Path(run)/'figures';out.mkdir(exist_ok=True);names=[]
 if result['analysis_type']=='drug_combination':
  context=result['primary']['assay_context']
  for model in dict.fromkeys(c['model'] for c in result['primary']['cells']):
   cells=[c for c in result['primary']['cells'] if c['model']==model];aa=sorted({c['dose_a'] for c in cells});bb=sorted({c['dose_b'] for c in cells});ix={(c['dose_a'],c['dose_b']):c['mean'] for c in cells}
   x=np.array([[np.nan if ix[a,b] is None else ix[a,b] for b in bb] for a in aa]);finite=x[np.isfinite(x)];lim=max(1,float(np.max(abs(finite)))) if len(finite) else 1
   with plt.rc_context(rc):
    fig,ax=plt.subplots(figsize=(5.6,4.6),layout='constrained')
    cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#dddddd');im=ax.imshow(x,cmap=cm,vmin=-lim,vmax=lim,origin='lower')
    ax.set(xticks=range(len(bb)),yticks=range(len(aa)),xticklabels=[f'{v:g}' for v in bb],yticklabels=[f'{v:g}' for v in aa],xlabel=f"{context['agent_b']} ({context['dose_unit_b']})",ylabel=f"{context['agent_a']} ({context['dose_unit_a']})",title=f'{model}: mean score (percentage points)')
    ax.tick_params(length=0)
    for i in range(len(aa)):
     for j in range(len(bb)):ax.text(j,i,'NA' if np.isnan(x[i,j]) else f'{x[i,j]:.1f}',ha='center',va='center',fontsize=7.5,color='white' if np.isfinite(x[i,j]) and abs(x[i,j])>.6*lim else 'black')
    fig.colorbar(im,ax=ax,label='Score (percentage points)',shrink=.85);name=f'synergy_{model}'
    pstyle.save(fig,out/name,token,300)
   names.append(name)
 else:
  import pandas as pd
  for index,p in enumerate(result['primary']['plates']):
   effects=p['correction']
   with plt.rc_context(rc):
    fig,axes=plt.subplots(1,2,figsize=(7.4,3.2),layout='constrained')
    for k,(axis,which,label) in enumerate(((axes[0],'row_effects','Row'),(axes[1],'column_effects','Column'))):
     marks=pstyle.point_style(token,0,4.5);marks.pop('linestyle')
     axis.plot(list(effects[which]),list(effects[which].values()),'-',color=pstyle.color(token,0),lw=1.2,**marks)
     axis.axhline(0,color=pstyle.MUTED,lw=.8,ls=(0,(3,2)))
     axis.set(xlabel=label,ylabel='Median-polish effect',title=p['plate_id'] if k==0 else f"QC: {p['status']}")
    name=f'plate_{index+1}'
    pstyle.save(fig,out/name,token,300)
   names.append(name)
   data=pd.read_csv(Path(run)/'input.csv',dtype={'plate_id':str});g=data[data.plate_id==p['plate_id']];grid=g.pivot(index='row',columns='column',values='value').sort_index().sort_index(axis=1)
   with plt.rc_context(rc):
    fig,ax=plt.subplots(figsize=(7.2,4.2),layout='constrained');im=ax.imshow(grid.to_numpy(),cmap='viridis')
    ax.set(xticks=range(len(grid.columns)),yticks=range(len(grid.index)),xticklabels=list(grid.columns),yticklabels=list(grid.index),xlabel='Column',ylabel='Row',title=f"{p['plate_id']}: observed plate signal")
    ax.tick_params(length=0)
    fig.colorbar(im,ax=ax,label=result['primary']['assay_context']['response_unit'],shrink=.85);name=f'plate_{index+1}_observed'
    pstyle.save(fig,out/name,token,300)
   names.append(name)
 return names
