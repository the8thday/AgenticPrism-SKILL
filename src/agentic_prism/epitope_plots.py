"""Plots saved directed results; never recomputes blocking or clusters."""
from pathlib import Path
import numpy as np

def render(run,result,style):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from scipy.cluster.hierarchy import dendrogram
 from .report import font_setup
 from . import plot_style as pstyle
 token=pstyle.THEMES[style];rc=pstyle.rc(token,*font_setup())
 out=Path(run)/'figures';out.mkdir(exist_ok=True)
 labels=result['matrix']['antibodies'];x=np.array([[np.nan if v is None else v for v in row] for row in result['matrix']['blocking']])
 with plt.rc_context(rc):
  fig,ax=plt.subplots(figsize=(max(4.6,len(labels)*.5+1.6),max(4,len(labels)*.45+1.2)),layout='constrained')
  cm=plt.get_cmap('viridis').copy();cm.set_bad('#dddddd');im=ax.imshow(x,vmin=0,vmax=1,cmap=cm)
  ax.set(xticks=range(len(labels)),yticks=range(len(labels)),xlabel='Second antibody',ylabel='First antibody',title='Directed normalized blocking')
  ax.set_xticklabels(labels,rotation=45,ha='right',rotation_mode='anchor');ax.set_yticklabels(labels)
  ax.tick_params(length=0)
  for i in range(len(labels)):
   for j in range(len(labels)):
    ax.text(j,i,'NA' if np.isnan(x[i,j]) else f'{x[i,j]:.2f}',ha='center',va='center',fontsize=7.5,color='black' if np.isnan(x[i,j]) or x[i,j]>.5 else 'white')
  fig.colorbar(im,ax=ax,label='Blocking (values outside scale retained in table)',shrink=.85)
  pstyle.save(fig,out/'blocking',token,300)
 names=['blocking']
 cl=result['clustering']
 if 'linkage' in cl:
  with plt.rc_context(rc):
   fig,ax=plt.subplots(figsize=(6.4,4),layout='constrained')
   dendrogram(np.asarray(cl['linkage']),labels=cl['labels'],ax=ax,color_threshold=0,above_threshold_color='black')
   ax.set(ylabel='Euclidean profile distance',title='Average linkage; branch BP in saved results')
   ax.spines['bottom'].set_visible(False);ax.tick_params(axis='x',length=0)
   pstyle.save(fig,out/'clustering',token,300)
  names.append('clustering')
 return names
