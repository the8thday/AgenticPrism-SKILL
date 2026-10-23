"""Artifact pipeline for independent specialist extensions; saved-result-only facts."""
from pathlib import Path
from datetime import datetime, timezone
import html
import json
import shutil
import platform
import importlib.metadata
from .workflow import dump,sha,implementation_hash,verify_run
from . import __version__
from .routine import TYPES as ROUTINE_TYPES
TYPES = set(ROUTINE_TYPES) | {"epitope_binning", "drug_combination", "hts_qc", "sample_size", "thermal_unfolding", "competing_risks", "nested_comparison", "ancova", "curve_auc", "standard_curve", "qpcr_relative", "nonlinear_fit", "location_test", "generalized_regression", "roc_analysis", "rater_agreement", "outlier_screen", "censored_comparison", "potency_ranking", "growth_rate_inhibition", "panel_comparison", "enzyme_inhibition", "receptor_occupancy", "pca", "clustering", "pk_nca", "population_pk", "fih_dose", "quantal_dose_response", "ordinal_regression"}


def module_for(kind):
 if kind == 'generalized_regression':
  from . import regression
  return regression
 if kind == 'roc_analysis':
  from . import diagnostic_accuracy
  return diagnostic_accuracy
 if kind in ('pca', 'clustering'):
  from . import multivariate
  return multivariate
 if kind == 'ordinal_regression':
  from . import ordinal
  return ordinal
 if kind == 'quantal_dose_response':
  from . import quantal
  return quantal
 if kind == 'fih_dose':
  from . import fih_dose
  return fih_dose
 if kind == 'population_pk':
  from . import population_pk
  return population_pk
 if kind == 'pk_nca':
  from . import pk_nca
  return pk_nca
 if kind == 'receptor_occupancy':
  from . import receptor_occupancy
  return receptor_occupancy
 if kind == 'enzyme_inhibition':
  from . import enzyme_kinetics
  return enzyme_kinetics
 if kind == 'panel_comparison':
  from . import panel
  return panel
 if kind == 'growth_rate_inhibition':
  from . import growth_rate
  return growth_rate
 if kind == 'potency_ranking':
  from . import potency_ranking
  return potency_ranking
 if kind == 'censored_comparison':
  from . import censored
  return censored
 if kind == 'outlier_screen':
  from . import outliers
  return outliers
 if kind == 'rater_agreement':
  from . import agreement
  return agreement
 if kind == 'location_test':
  from . import location_tests
  return location_tests
 from . import routine, epitope, combination, hts, sample_size, thermal, competing, nested, ancova, curve_auc, standard_curve, qpcr, nonlinear
 if kind in nonlinear.TYPES:return nonlinear
 if kind in qpcr.TYPES:return qpcr
 if kind in standard_curve.TYPES:return standard_curve
 if kind in curve_auc.TYPES:return curve_auc
 if kind in ancova.TYPES:return ancova
 if kind in nested.TYPES:return nested
 if kind in competing.TYPES:return competing
 if kind in sample_size.TYPES:return sample_size
 if kind in thermal.TYPES:return thermal
 if kind in combination.TYPES:return combination
 if kind in hts.TYPES:return hts
 if kind in epitope.TYPES:return epitope
 if kind in routine.TYPES:return routine
 raise ValueError('Unsupported extension')


def write_facts(out):
 out=Path(out);r=json.loads((out/'results.json').read_text())
 from .interpretation import disclosures
 dump(out/'interpretation_facts.json',{'schema_version':1,'analysis_type':r['analysis_type'],
  'source_sha256':{n:sha(out/n) for n in ('results.json','config.resolved.json')},
  'primary':{'source':'results.json#/primary','value':r['primary']},
  'states':disclosures(r),'must_mention':r['must_mention'],'failing_items':r['failing_items'],
  'limitations':r['limitations'],'validation_evidence':r.get('validation_evidence',{})})


def analyze_extension(config_path,output,render=True):
 src=Path(config_path).resolve();out=Path(output).resolve()
 if out.exists():raise FileExistsError('Choose a new output directory')
 out.mkdir(parents=True)
 try:
  raw=json.loads(src.read_text());module=module_for(raw['analysis_type']);cfg=module.resolve_config(raw)
  if cfg.get('input') is None:  # design-only analyses (sample size) declare assumptions, not data
   inp=None;data=None;dump(out/'config.resolved.json',cfg)
  else:
   inp=(src.parent/cfg['input']).resolve();data=module.load_data(inp,cfg)
   shutil.copyfile(inp,out/'input.csv');cfg['input']='input.csv';dump(out/'config.resolved.json',cfg)
  r=module.compute(data,cfg)
  try:
   from .evidence_0111 import EVIDENCE
   if r['analysis_type'] in EVIDENCE:
    r['validation_evidence']=EVIDENCE[r['analysis_type']];r['must_mention']+=r['validation_evidence']['must_mention']
  except ImportError:pass
  if r['analysis_type']=='epitope_binning':
   from .evidence_0121 import EVIDENCE
   r['validation_evidence']=EVIDENCE[r['analysis_type']];r['must_mention']+=r['validation_evidence']['must_mention']
  if r['analysis_type'] in ('drug_combination','hts_qc'):
   from .evidence_0130 import EVIDENCE
   r['validation_evidence']=EVIDENCE[r['analysis_type']];r['must_mention']+=r['validation_evidence']['must_mention']
  if r['analysis_type'] in ('sample_size','competing_risks','thermal_unfolding'):
   from .evidence_0131 import EVIDENCE
   if r['analysis_type'] in EVIDENCE:
    r['validation_evidence']=EVIDENCE[r['analysis_type']];r['must_mention']+=r['validation_evidence']['must_mention']
  if r['analysis_type'] in ('nested_comparison','ancova','curve_auc','standard_curve','qpcr_relative') and not (r['analysis_type']=='curve_auc' and cfg['comparison'].get('design')=='paired_t'):
   from .evidence_0132 import EVIDENCE
   r['validation_evidence']=EVIDENCE[r['analysis_type']];r['must_mention']+=r['validation_evidence']['must_mention']
  if r['analysis_type']=='nonlinear_fit':
   from .evidence_0133 import EVIDENCE
   r['validation_evidence']=EVIDENCE['nonlinear_fit'];r['must_mention']+=r['validation_evidence']['must_mention']
  if r['analysis_type']=='hts_qc' and 'hit_reference' in cfg:
   from .evidence_0131 import EVIDENCE
   r['validation_evidence_0131']=EVIDENCE['hts_layout_simulation'];r['must_mention']+=r['validation_evidence_0131']['must_mention']
  if r['analysis_type']=='location_test' or (r['analysis_type']=='curve_auc' and cfg['comparison'].get('design')=='paired_t'):
   from .evidence_0134 import EVIDENCE
   key='location_test' if r['analysis_type']=='location_test' else 'paired_auc'
   r['validation_evidence']=EVIDENCE[key];r['must_mention']+=EVIDENCE[key]['must_mention']
  if r['analysis_type']=='generalized_regression':
   from .evidence_0135 import EVIDENCE
   r['validation_evidence']=EVIDENCE;r['must_mention']+=EVIDENCE['must_mention']
   module.write_artifacts(out,r)
  if r['analysis_type'] in ('roc_analysis','rater_agreement','outlier_screen'):
   from .evidence_0136 import EVIDENCE
   r['validation_evidence']=EVIDENCE[r['analysis_type']];r['must_mention']+=EVIDENCE[r['analysis_type']]['must_mention']
  if r['analysis_type']=='censored_comparison' and cfg.get('paired') is not None:
   from .evidence_0141 import EVIDENCE as E141
   r['validation_evidence']=E141['censored_paired'];r['must_mention']+=E141['censored_paired']['must_mention']
  elif r['analysis_type'] in ('censored_comparison','potency_ranking','growth_rate_inhibition','panel_comparison'):
   from .evidence_0140 import EVIDENCE as E140
   r['validation_evidence']=E140[r['analysis_type']];r['must_mention']+=E140[r['analysis_type']]['must_mention']
   if r['analysis_type']=='growth_rate_inhibition' and cfg['fit'].get('grmax_definition')=='fitted_minimum_top_two':
    from .evidence_0141 import EVIDENCE as E141
    r['validation_evidence_0141']=E141['gr_fitted_grmax'];r['must_mention']+=E141['gr_fitted_grmax']['must_mention']
   if r['analysis_type']=='panel_comparison' and cfg.get('censoring') is not None:
    from .evidence_0141 import EVIDENCE as E141
    r['validation_evidence_0141']=E141['panel_censored'];r['must_mention']+=E141['panel_censored']['must_mention']
  if r['analysis_type']=='receptor_occupancy':
   from .evidence_0142 import EVIDENCE as E142
   r['validation_evidence']=E142['receptor_occupancy'];r['must_mention']+=E142['receptor_occupancy']['must_mention']
   if cfg.get('exposure_response',{}).get('interval')=='cluster_robust_t':
    from .evidence_0143 import EVIDENCE as E143
    r['validation_evidence_0143']=E143['ro_cluster_robust_t'];r['must_mention']+=E143['ro_cluster_robust_t']['must_mention']
  if r['analysis_type']=='enzyme_inhibition' and cfg.get('tight_binding') is not None:
   from .evidence_0143 import EVIDENCE as E143
   r['validation_evidence']=E143['enzyme_morrison'];r['must_mention']+=E143['enzyme_morrison']['must_mention']
  elif r['analysis_type']=='enzyme_inhibition':
   from .evidence_0142 import EVIDENCE as E142
   r['validation_evidence']=E142['enzyme_inhibition'];r['must_mention']+=E142['enzyme_inhibition']['must_mention']
  if r['analysis_type'] in ('quantal_dose_response','ordinal_regression'):
   from .evidence_0153 import EVIDENCE as E153
   r['validation_evidence']=E153[r['analysis_type']];r['must_mention']+=E153[r['analysis_type']]['must_mention']
  if r['analysis_type']=='fih_dose':
   from .evidence_0152 import EVIDENCE as E152
   r['validation_evidence']=E152['fih_dose'];r['must_mention']+=E152['fih_dose']['must_mention']
  if r['analysis_type']=='population_pk':
   from .evidence_0151 import EVIDENCE as E151
   r['validation_evidence']=E151['population_pk'];r['must_mention']+=E151['population_pk']['must_mention']
  if r['analysis_type']=='pk_nca':
   from .evidence_0150 import EVIDENCE as E150
   r['validation_evidence']=E150['pk_nca'];r['must_mention']+=E150['pk_nca']['must_mention']
  if r['analysis_type'] in ('pca','clustering'):
   module.write_artifacts(out,r)
   from .evidence_0144 import EVIDENCE as E144
   r['validation_evidence']=E144[r['analysis_type']];r['must_mention']+=E144[r['analysis_type']]['must_mention']
  if r['analysis_type']=='growth_rate_inhibition':
   import pandas as pd
   pd.DataFrame(r['gr_values']).to_csv(out/'gr_values.csv',index=False)
  if r['analysis_type']=='roc_analysis':
   import pandas as pd
   pd.DataFrame(r['roc_curve']).to_csv(out/'roc_curve.csv',index=False)
  dump(out/'results.json',r);write_facts(out)
  (out/'rerun.txt').write_text('agentic-prism analyze --config config.resolved.json --output ../extension-rerun-new\n')
  dump(out/'manifest.json',{'schema_version':1,'analysis_type':cfg['analysis_type'],'package_version':__version__,
   'created_utc':datetime.now(timezone.utc).isoformat(),'input_sha256':sha(inp) if inp else None,'source':cfg['source'],
   'python':platform.python_version(),'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','scipy','pandas','statsmodels')},
   'implementation_sha256':implementation_hash(),'prism_numerical_equivalence':'not_claimed',
   'scientific_artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}})
 except Exception as e:
  dump(out/'failure.json',{'status':'failed','exception':type(e).__name__,'message':str(e)});raise
 if render:render_extension(out)
 return out


def render_extension(run,style=None,extra_figures=''):
 from . import report_shell as shell
 run=Path(run);verify_run(run);r=json.loads((run/'results.json').read_text());cfg=json.loads((run/'config.resolved.json').read_text());style=style or cfg['report']['plot_style']
 if style not in ('prism_like','standard'):raise ValueError('Unsupported style')
 show=shell.value_html
 names=[]
 if r['analysis_type']=='generalized_regression':
  from .regression import render
  names=render(run,r,style)
 if r['analysis_type']=='roc_analysis':
  from .diagnostic_accuracy import render
  names=render(run,r,style)
 if r['analysis_type'] in ('pca','clustering'):
  from .multivariate import render
  names=render(run,r,style)
 if r['analysis_type']=='ordinal_regression':
  from .ordinal import render
  names=render(run,r,style)
 if r['analysis_type']=='quantal_dose_response':
  from .quantal import render
  names=render(run,r,style)
 if r['analysis_type']=='fih_dose':
  from .fih_dose import render
  names=render(run,r,style)
 if r['analysis_type']=='population_pk':
  from .population_pk import render
  names=render(run,r,style)
 if r['analysis_type']=='pk_nca':
  from .pk_nca import render
  names=render(run,r,style)
 if r['analysis_type']=='receptor_occupancy':
  from .receptor_occupancy import render
  names=render(run,r,style)
 if r['analysis_type']=='enzyme_inhibition':
  from .enzyme_kinetics import render
  names=render(run,r,style)
 if r['analysis_type']=='panel_comparison':
  from .panel import render
  names=render(run,r,style)
 if r['analysis_type']=='growth_rate_inhibition':
  from .growth_rate import render
  names=render(run,r,style)
 if r['analysis_type']=='potency_ranking':
  from .potency_ranking import render
  names=render(run,r,style)
 if r['analysis_type']=='censored_comparison':
  from .censored import render
  names=render(run,r,style)
 if r['analysis_type']=='outlier_screen':
  from .outliers import render
  names=render(run,r,style)
 if r['analysis_type']=='rater_agreement':
  from .agreement import render
  names=render(run,r,style)
 if r['analysis_type']=='location_test':
  from .location_tests import render
  names=render(run,r,style)
 if r['analysis_type']=='epitope_binning':
  from .epitope_plots import render
  names=render(run,r,style)
 if r['analysis_type'] in ('drug_combination','hts_qc'):
  from .pharmacology_plots import render
  names=render(run,r,style)
 if r['analysis_type']=='sample_size':
  from .sample_size import render
  names=render(run,r,style)
 if r['analysis_type']=='thermal_unfolding':
  from .thermal import render
  names=render(run,r,style)
 if r['analysis_type']=='competing_risks':
  from .competing import render
  names=render(run,r,style)
 if r['analysis_type']=='nested_comparison':
  from .nested import render
  names=render(run,r,style)
 if r['analysis_type']=='ancova':
  from .ancova import render
  names=render(run,r,style)
 if r['analysis_type']=='curve_auc':
  from .curve_auc import render
  names=render(run,r,style)
 if r['analysis_type']=='standard_curve':
  from .standard_curve import render
  names=render(run,r,style)
 if r['analysis_type']=='qpcr_relative':
  from .qpcr import render
  names=render(run,r,style)
 if r['analysis_type']=='nonlinear_fit':
  from .nonlinear import render
  names=render(run,r,style)
 figures=''.join(shell.single_figure(run/'figures',name,name.replace('_',' ')) for name in names)+extra_figures
 title=r['analysis_type'].replace('_',' ')
 body=shell.section('primary','主要结果',show(r['primary']))
 body+=shell.section('must','必须说明','<ul class="must">'+''.join('<li>'+html.escape(str(m))+'</li>' for m in r['must_mention'])+'</ul>')
 if figures:body+=shell.section('figures','图形',figures)
 body+=shell.section('details','详细结果、诊断与证据','<details><summary>展开全部字段</summary>'+show({k:v for k,v in r.items() if k not in ('primary','must_mention')})+'</details>')
 body+=shell.section('files','可追溯文件',shell.downloads([(n,n) for n in ('results.json','interpretation_facts.json','input.csv','config.resolved.json') if (run/n).exists()]))
 nav=[('primary','主要结果'),('must','必须说明')]+([('figures','图形')] if figures else [])+[('details','详细结果'),('files','可追溯文件')]
 page=shell.page('分析结果 · '+title,eyebrow='AgenticPrism / '+title,heading=title[:1].upper()+title[1:],
                 lede='保存结果的只读呈现；声明的读数、响应尺度、独立性与设计范围限制同样适用。',nav=nav,body=body,
                 footer=f'AgenticPrism · 图形风格 {html.escape(style)} · 数值分析与图形渲染分离')
 (run/'report.html').write_text(page);dump(run/'render_manifest.json',{'style':style,'scientific_artifacts_changed':False,'report_sha256':sha(run/'report.html')});return run/'report.html'
