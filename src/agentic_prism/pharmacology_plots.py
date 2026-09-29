"""Saved-result-only pharmacology graphics."""
from pathlib import Path
import numpy as np

def render(run,result,style):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 out=Path(run)/'figures';out.mkdir(exist_ok=True);names=[]
 if result['analysis_type']=='drug_combination':
  context=result['primary']['assay_context']
  for model in dict.fromkeys(c['model'] for c in result['primary']['cells']):
   cells=[c for c in result['primary']['cells'] if c['model']==model];aa=sorted({c['dose_a'] for c in cells});bb=sorted({c['dose_b'] for c in cells});ix={(c['dose_a'],c['dose_b']):c['mean'] for c in cells}
   x=np.array([[np.nan if ix[a,b] is None else ix[a,b] for b in bb] for a in aa]);fig,ax=plt.subplots(figsize=(7,6));finite=x[np.isfinite(x)];lim=max(1,float(np.max(abs(finite)))) if len(finite) else 1
   cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#dddddd');im=ax.imshow(x,cmap=cm,vmin=-lim,vmax=lim)
   ax.set(xticks=range(len(bb)),yticks=range(len(aa)),xticklabels=[f'{v:g}' for v in bb],yticklabels=[f'{v:g}' for v in aa],xlabel=f"{context['agent_b']} ({context['dose_unit_b']})",ylabel=f"{context['agent_a']} ({context['dose_unit_a']})",title=f'{model}: mean score (percentage points)')
   for i in range(len(aa)):
    for j in range(len(bb)):ax.text(j,i,'NA' if np.isnan(x[i,j]) else f'{x[i,j]:.1f}',ha='center',va='center',fontsize=8,color='white' if np.isfinite(x[i,j]) and abs(x[i,j])>.6*lim else 'black')
   fig.colorbar(im,ax=ax);fig.tight_layout();name=f'synergy_{model}'
   for ext in ('svg','pdf','png'):fig.savefig(out/f'{name}.{ext}',dpi=180)
   plt.close(fig);names.append(name)
 else:
  for index,p in enumerate(result['primary']['plates']):
   fig,axes=plt.subplots(1,2,figsize=(10,4));effects=p['correction'];axes[0].plot(list(effects['row_effects']),list(effects['row_effects'].values()),'o-');axes[0].set(xlabel='Row',ylabel='Median-polish effect',title=p['plate_id'])
   axes[1].plot(list(effects['column_effects']),list(effects['column_effects'].values()),'o-');axes[1].set(xlabel='Column',ylabel='Median-polish effect',title=f"QC: {p['status']}");fig.tight_layout();name=f'plate_{index+1}'
   for ext in ('svg','pdf','png'):fig.savefig(out/f'{name}.{ext}',dpi=180)
   plt.close(fig);names.append(name)
   import pandas as pd
   data=pd.read_csv(Path(run)/'input.csv',dtype={'plate_id':str});g=data[data.plate_id==p['plate_id']];grid=g.pivot(index='row',columns='column',values='value').sort_index().sort_index(axis=1)
   fig,ax=plt.subplots(figsize=(9,5));im=ax.imshow(grid.to_numpy(),cmap='viridis');ax.set(xticks=range(len(grid.columns)),yticks=range(len(grid.index)),xticklabels=list(grid.columns),yticklabels=list(grid.index),xlabel='Column',ylabel='Row',title=f"{p['plate_id']}: observed plate signal")
   fig.colorbar(im,ax=ax,label=result['primary']['assay_context']['response_unit']);fig.tight_layout();name=f'plate_{index+1}_observed'
   for ext in ('svg','pdf','png'):fig.savefig(out/f'{name}.{ext}',dpi=180)
   plt.close(fig);names.append(name)
 return names
