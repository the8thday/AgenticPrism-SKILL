"""Plots saved directed results; never recomputes blocking or clusters."""
from pathlib import Path
import numpy as np

def render(run,result,style):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from scipy.cluster.hierarchy import dendrogram
 out=Path(run)/'figures';out.mkdir(exist_ok=True)
 labels=result['matrix']['antibodies'];x=np.array([[np.nan if v is None else v for v in row] for row in result['matrix']['blocking']])
 fig,ax=plt.subplots(figsize=(max(6,len(labels)*.55),max(5,len(labels)*.5)))
 cm=plt.get_cmap('viridis').copy();cm.set_bad('#dddddd');im=ax.imshow(x,vmin=0,vmax=1,cmap=cm)
 ax.set(xticks=range(len(labels)),yticks=range(len(labels)),xticklabels=labels,yticklabels=labels,xlabel='Second antibody',ylabel='First antibody',title='Directed normalized blocking')
 for i in range(len(labels)):
  for j in range(len(labels)):
   ax.text(j,i,'NA' if np.isnan(x[i,j]) else f'{x[i,j]:.2f}',ha='center',va='center',fontsize=8,color='black' if np.isnan(x[i,j]) or x[i,j]>.5 else 'white')
 fig.colorbar(im,ax=ax,label='Blocking (values outside scale retained in table)');fig.tight_layout()
 for ext in ('svg','pdf','png'):fig.savefig(out/f'blocking.{ext}',dpi=180)
 plt.close(fig);names=['blocking']
 cl=result['clustering']
 if 'linkage' in cl:
  fig,ax=plt.subplots(figsize=(8,5));dendrogram(np.asarray(cl['linkage']),labels=cl['labels'],ax=ax)
  ax.set(ylabel='Euclidean profile distance',title='Average linkage; branch BP in saved results')
  fig.tight_layout()
  for ext in ('svg','pdf','png'):fig.savefig(out/f'clustering.{ext}',dpi=180)
  plt.close(fig);names.append('clustering')
 return names
