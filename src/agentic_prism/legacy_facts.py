"""Retrofit facts from saved artifacts only. Audit and withheld values stay separate."""
import json
from pathlib import Path
from copy import deepcopy
from .workflow import dump,sha
from .interpretation import disclosures


def safe_copy(value):
 if isinstance(value,list):return [safe_copy(v) for v in value]
 if not isinstance(value,dict):return value
 if value.get('reportable') is False or value.get('status') in ('failed','withheld'):
  return {k:deepcopy(v) for k,v in value.items() if k in ('status','reportable','diagnostics','reason','curve_id','fit_group_id','plate_id','sample_id')}
 return {k:safe_copy(v) for k,v in value.items() if not k.startswith('audit_')}


def write_facts(out, all_states=False):
 out=Path(out);r=json.loads((out/'results.json').read_text());cfg=json.loads((out/'config.resolved.json').read_text())
 states=disclosures(r)
 evidence=[{'kind':'calibration_evidence','value':m,'source':f'results.json#/validation_evidence/must_mention/{i}'} for i,m in enumerate((r.get('validation_evidence') or {}).get('must_mention',[]))] if r['analysis_type']=='dose_response_4pl' else []
 supplemental={}
 if all_states:supplemental['assay_context']={'source':'config.resolved.json#','value':cfg}
 if r['analysis_type']=='binding_kinetics':supplemental['assay_context']={'source':'config.resolved.json#','value':cfg}
 if r['analysis_type']=='binding_kinetics' and (out/'steady_state.json').exists():
  ss=json.loads((out/'steady_state.json').read_text());supplemental.update({'steady_state':{'source':'steady_state.json#','sha256':sha(out/'steady_state.json'),'value':safe_copy(ss)}})
  states+=[{**state,'source':state['source'].replace('results.json#','steady_state.json#',1)} for state in disclosures(ss)]
 dump(out/'interpretation_facts.json',{'schema_version':1,'analysis_type':r['analysis_type'],
  'source_sha256':{n:sha(out/n) for n in ('results.json','config.resolved.json')},
  'saved_results':{'source':'results.json#','value':safe_copy(r)},'states':states,
  'must_mention':(states if all_states else [s for s in states if s['kind']=='diagnostic' or s['value'] in (False,'failed','withheld','limited')])+evidence,
  **supplemental,
  'limitations':['Saved-result extraction only; no refitting, new inference or promotion of audit estimates.',
   'Independent units and declared preprocessing determine the interpretation; significance does not establish equivalence.']})
