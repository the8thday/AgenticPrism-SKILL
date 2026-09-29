"""Opt-in kinetics workflow; legacy scientific artifacts use their original path."""
from pathlib import Path
from datetime import datetime,timezone
import json,shutil,platform
import numpy as np
import pandas as pd
from . import __version__
from .workflow import dump,sha,implementation_hash,verify_run
from .kinetics_schema import resolve_kinetic_config,load_kinetic_data
from .advanced_kinetics import fit_group,off_rate
from .extension_workflow import write_facts,render_extension


def analyze_advanced(config_path,output,render=True):
 src=Path(config_path).resolve();out=Path(output).resolve()
 if out.exists():raise FileExistsError('Choose a new output directory')
 out.mkdir(parents=True)
 try:
  cfg=resolve_kinetic_config(json.loads(src.read_text()));inp=(src.parent/cfg['input']).resolve();d,pre=load_kinetic_data(inp,cfg);shutil.copy2(inp,out/'input.csv');cfg['input']='input.csv'
  if cfg['provenance']:
   p=(src.parent/cfg['provenance']).resolve();record=json.loads(p.read_text())
   if record.get('data_sha256')!=sha(inp):raise ValueError('Importer hash does not match data')
   shutil.copy2(p,out/'source_provenance.json');cfg['provenance']='source_provenance.json'
  dump(out/'config.resolved.json',cfg);d.to_csv(out/'normalized_data.csv',index=False);dump(out/'preprocessing_log.json',{'configuration':cfg['preprocessing'],'curves':pre,'automatic_smoothing':False,'phase_boundaries':'declared metadata'})
  fits=[(off_rate if cfg['model']=='off_rate_screening' else fit_group)(g,cfg) for _,g in d.groupby('fit_group_id',sort=False)]
  notes=['Complex model predeclared with mechanistic rationale; better fit alone does not establish mechanism.','Time points and sensor curves are not independent experiments.','Readouts must track bound analyte on the declared common response scale.']
  if cfg['model']=='bivalent_analyte':notes+=['First-arm ka includes factor two for free analyte arms; second association is in response^-1 s^-1. A single KD is not defined.']
  if cfg['model']=='heterogeneous_ligand':notes+=['Ordered surface KD1 < KD2 describe heterogeneous surface sites, not two epitopes.']
  if cfg['model']=='off_rate_screening':notes+=['Dissociation-only apparent koff ranking; no affinity or kon is estimated.']
  for f in fits:
   notes+=f['diagnostics']
   for curve in f.get('curves',[]):notes += ['Curve '+curve['curve_id']+': '+reason for reason in curve['diagnostics']]
  r={'schema_version':1,'analysis_type':'binding_kinetics','model':cfg['model'],'fits':fits,'primary':{'status':'estimated' if all(f['status']=='estimated' for f in fits) else 'limited','fits':[{'fit_group_id':f['fit_group_id'],'reportable':f['reportable'],'parameters':f.get('parameters'),'curves':[{k:v for k,v in c.items() if not k.startswith('audit_') and k!='bootstrap'} for c in f.get('curves',[])]} for f in fits]},'must_mention':list(dict.fromkeys(notes)),'failing_items':[{'fit_group_id':f['fit_group_id'],'reasons':f['diagnostics']} for f in fits if not f['reportable']],'limitations':['Bootstrap is conditional on declared concentrations, reference processing and mechanism.','Profile endpoints test support within a declared log span; they are not full confidence intervals.','No instrument-software or Prism equivalence claim.']}
  r['primary']['assay_context']={'assay':cfg['assay'],'advanced_declarations':cfg['advanced'],'weighting':cfg['fit']['weighting'],'rmax_sharing':cfg['fit']['rmax']}
  r['failing_items'] += [{'fit_group_id':f['fit_group_id'],'curve_id':curve['curve_id'],'reasons':curve['diagnostics']} for f in fits for curve in f.get('curves',[]) if not curve['reportable']]
  try:
   from .evidence_0120 import EVIDENCE
   r['validation_evidence']=EVIDENCE;r['must_mention']+=EVIDENCE['must_mention']
  except ImportError:pass
  dump(out/'results.json',r);pd.DataFrame([p for f in fits for p in f['predictions']]).to_csv(out/'predictions.csv',index=False);write_facts(out)
  (out/'rerun.txt').write_text('agentic-prism analyze --config config.resolved.json --output ../advanced-kinetics-rerun-new\n')
  dump(out/'manifest.json',{'schema_version':1,'analysis_type':'binding_kinetics','package_version':__version__,'created_utc':datetime.now(timezone.utc).isoformat(),'input_sha256':sha(inp),'source':cfg['source'],'python':platform.python_version(),'implementation_sha256':implementation_hash(),'scientific_artifacts_sha256':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}})
 except Exception as e:dump(out/'failure.json',{'exception':type(e).__name__,'message':str(e)});raise
 if render:render_advanced(out)
 return out


def render_advanced(run, style=None):
    run = Path(run)
    verify_run(run)
    result = json.loads((run / 'results.json').read_text())
    cfg = json.loads((run / 'config.resolved.json').read_text())
    style = style or cfg['report']['plot_style']
    report = render_extension(run, style)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import data_uri
    width = cfg['report']['figure_width_mm'] / 25.4
    folder = run / 'figures'
    folder.mkdir(exist_ok=True)
    extra = []

    def save_figure(fig, name):
        for extension in ('svg', 'pdf', 'png'):
            fig.savefig(folder / f'{name}.{extension}', dpi=cfg['report']['png_dpi'])
        plt.close(fig)
        extra.append((folder / f'{name}.svg').read_text())
        extra.append('<p>' + ' | '.join(
            f'<a download="{name}.{extension}" href="{data_uri(folder / f"{name}.{extension}")}">{extension.upper()}</a>'
            for extension in ('svg', 'pdf', 'png')) + '</p>')

    fig, (ax, residual_axis) = plt.subplots(
        2, 1, figsize=(width, width * .75), sharex=True, layout='constrained',
        gridspec_kw={'height_ratios': [3, 1]})
    color_index = 0
    for fit in result['fits']:
        if not fit.get('predictions'):
            continue
        data = pd.DataFrame(fit['predictions'])
        for curve, group in data.groupby('curve_id', sort=False):
            color = plt.get_cmap('tab10')(color_index % 10)
            color_index += 1
            ax.plot(group.time_s, group.observed, '.', ms=2, alpha=.5, color=color)
            ax.plot(group.time_s, group.predicted, label=curve, color=color)
            residual_axis.plot(group.time_s, group.residual, '.', ms=2, color=color)
    ax.set(ylabel='Response (declared units)', title='Observed and fitted traces; consult reliability flags')
    xlabel = ('Time since first fitted dissociation point (s)'
              if result['model'] == 'off_rate_screening' else 'Time since association (s)')
    residual_axis.set(xlabel=xlabel, ylabel='Residual')
    residual_axis.axhline(0, color='grey', lw=.7)
    if style == 'prism_like':
        ax.spines[['top', 'right']].set_visible(False)
        residual_axis.spines[['top', 'right']].set_visible(False)
    if 0 < color_index <= 12:
        ax.legend(fontsize=7, ncol=2)
    save_figure(fig, 'kinetics')

    points = []
    for fit in result['fits']:
        p = fit.get('parameters')
        if not p:
            continue
        if 'kd_M' in p:
            points.append((p['ka_M_inv_s_inv'], p['kd'], fit['fit_group_id']))
        elif 'surface_KDs_M' in p:
            points.extend([
                (p['ka_M_inv_s_inv'], p['kd'], fit['fit_group_id'] + ' surface site1'),
                (p['ka2_M_inv_s_inv'], p['kd2'], fit['fit_group_id'] + ' surface site2')])
    if points:
        fig, ax = plt.subplots(figsize=(width, width * .7), layout='constrained')
        xx = np.logspace(2, 8, 100)
        for kd in (1e-12, 1e-10, 1e-8, 1e-6):
            ax.loglog(xx, kd * xx, '--', color='grey', alpha=.5)
            ax.text(xx[-1], kd * xx[-1], f'{kd:g} M', fontsize=8)
        for ka, kd, label in points:
            ax.loglog(ka, kd, 'o', label=label)
        if len(points) <= 12:
            ax.legend(fontsize=7)
        ax.set(xlabel='ka (M^-1 s^-1)', ylabel='kd (s^-1)', title='Reportable kinetic KD only')
        save_figure(fig, 'iso_affinity')
    report.write_text(report.read_text().replace(
        '<h2>详细结果、诊断与证据</h2>', ''.join(extra) + '<h2>详细结果、诊断与证据</h2>'))
    dump(run / 'render_manifest.json', {
        'style': style, 'scientific_artifacts_changed': False, 'report_sha256': sha(report)})
    return report
