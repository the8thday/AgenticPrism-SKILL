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
TYPES = set(ROUTINE_TYPES) | {"epitope_binning", "drug_combination", "hts_qc"}


def module_for(kind):
 from . import routine, epitope, combination, hts
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
  dump(out/'results.json',r);write_facts(out)
  (out/'rerun.txt').write_text('agentic-prism analyze --config config.resolved.json --output ../extension-rerun-new\n')
  dump(out/'manifest.json',{'schema_version':1,'analysis_type':cfg['analysis_type'],'package_version':__version__,
   'created_utc':datetime.now(timezone.utc).isoformat(),'input_sha256':sha(inp),'source':cfg['source'],
   'python':platform.python_version(),'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','scipy','pandas','statsmodels')},
   'implementation_sha256':implementation_hash(),'prism_numerical_equivalence':'not_claimed',
   'scientific_artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}})
 except Exception as e:
  dump(out/'failure.json',{'status':'failed','exception':type(e).__name__,'message':str(e)});raise
 if render:render_extension(out)
 return out


def render_extension(run,style=None):
 run=Path(run);verify_run(run);r=json.loads((run/'results.json').read_text());cfg=json.loads((run/'config.resolved.json').read_text());style=style or cfg['report']['plot_style']
 if style not in ('prism_like','standard'):raise ValueError('Unsupported style')
 def show(v):
  if isinstance(v,dict):return '<table>'+''.join('<tr><th>'+html.escape(str(k))+'</th><td>'+show(x)+'</td></tr>' for k,x in v.items())+'</table>'
  if isinstance(v,list):return '<ol>'+''.join('<li>'+show(x)+'</li>' for x in v)+'</ol>'
  return html.escape(str(v))
 body='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>分析结果</title><style>body{font:16px sans-serif;max-width:1100px;margin:32px auto;padding:16px}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:6px;text-align:left;vertical-align:top}li{margin:6px}h1{color:#24445c}</style>'
 body+='<h1>'+html.escape(r['analysis_type'])+'</h1><h2>主要结果</h2>'+show(r['primary'])+'<h2>必须说明</h2>'+show(r['must_mention'])+'<h2>详细结果、诊断与证据</h2>'+show({k:v for k,v in r.items() if k not in ('primary','must_mention')})
 if r['analysis_type']=='epitope_binning':
  from .epitope_plots import render
  for name in render(run,r,style):
   body+=f'<img src="figures/{name}.svg" style="max-width:100%">'+''.join(f'<a href="figures/{name}.{ext}" download>{ext}</a> ' for ext in ('svg','pdf','png'))
 if r['analysis_type'] in ('drug_combination','hts_qc'):
  from .pharmacology_plots import render
  for name in render(run,r,style):
   body+=f'<img src="figures/{name}.svg" style="max-width:100%">'+''.join(f'<a href="figures/{name}.{ext}" download>{ext}</a> ' for ext in ('svg','pdf','png'))
 body+=''.join(f'<p><a href="{n}" download="{n}">{n}</a></p>' for n in ('results.json','interpretation_facts.json','input.csv','config.resolved.json'))+'</html>'
 (run/'report.html').write_text(body);dump(run/'render_manifest.json',{'style':style,'scientific_artifacts_changed':False,'report_sha256':sha(run/'report.html')});return run/'report.html'
