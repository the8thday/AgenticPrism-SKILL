"""Saved-result affinity facts, immutable artifacts, and offline reporting."""
from pathlib import Path
from datetime import datetime, timezone
import json
import html
import shutil
import platform
import importlib.metadata
import numpy as np
import pandas as pd
from . import __version__
from .workflow import dump,sha,verify_run,implementation_hash
from .fit import summarize


def write_facts(out):
    """Extract only from saved files; never fit or derive an interval here."""
    out=Path(out); r=json.loads((out/'results.json').read_text()); c=json.loads((out/'config.resolved.json').read_text())
    primary=[]; failing=[]; must=list(r.get('must_mention',[]))
    for i,v in enumerate(r['fits']):
        item={'source':f'results.json#/fits/{i}','curve_id':v['curve_id'],'status':v['status'],
            'reportable':v['reportable'],'interpretation':v.get('interpretation',c['assay']['interpretation']),
            'kd_M':v.get('kd_M') if v['reportable'] else None,'ki_M':v.get('ki_M') if v['reportable'] else None,
            'ci_low_M':v.get('ci_low_M') if v['reportable'] else None,'ci_high_M':v.get('ci_high_M') if v['reportable'] else None,
            'supported_upper_bound_M':v.get('supported_upper_bound_M'), 'ci_status':v['ci_status'],
            'diagnostics':v['diagnostics'],'regime_diagnostics':v.get('regime_diagnostics',[]),
            'design_warnings':v.get('design_warnings',[])}
        primary.append(item)
        if not v['reportable']: failing.append({'source':item['source'],'curve_id':v['curve_id'],'reasons':v['diagnostics']})
        must+=v['diagnostics']
        if v.get('design_warnings'):
            must.append('Design warning: the titrant never reached twice the constant-species concentration, so the '
                        'equivalence region was not spanned and KD information comes only from sub-saturation curvature.')
    must+=['Report only reportable points or explicitly supported bounds; audit optima are not reported affinity.',
           'Profile-F intervals are model-conditional and not exact finite-sample coverage guarantees.',
           'Independent experiments receive equal log-scale weights; wells are not experiments.']
    if c['analysis_type']=='cell_binding': must+=['Always apparent KD; bivalent IgG binding is avidity-influenced, not intrinsic affinity.']
    dump(out/'interpretation_facts.json',{'schema_version':1,'analysis_type':c['analysis_type'],
        'source_sha256':{n:sha(out/n) for n in ('results.json','config.resolved.json')},'primary':primary,
        'summaries':{'source':'results.json#/summaries','value':r['summaries']},'failing_items':failing+r.get('failing_items',[]),
        'must_mention':list(dict.fromkeys(must)), 'validation_evidence':r.get('validation_evidence',{'scope':'Legacy model; consult historical validation; no new calibration claim.'})})


def analyze_affinity(config_path,output,render=True):
    from . import affinity_schema
    from .affinity_fit import fit_affinity
    src,out=Path(config_path).resolve(),Path(output).resolve()
    if out.exists(): raise FileExistsError('Choose a new affinity output directory')
    out.mkdir(parents=True)
    try:
        raw=json.loads(src.read_text())
        if raw.get('analysis_type')=='cell_binding':
            from . import cell_binding as module
        else: module=affinity_schema
        cfg=module.resolve_config(raw); inp=(src.parent/cfg['input']).resolve(); d=module.load_data(inp,cfg)
        shutil.copyfile(inp,out/'input.csv'); cfg['input']='input.csv'; dump(out/'config.resolved.json',cfg)
        d.to_csv(out/'normalized_data.csv',index=False)
        fits=[fit_affinity(g,cfg) for _,g in d.groupby('fit_group_id',sort=False)]
        result={'schema_version':1,'analysis_type':cfg['analysis_type'],'model':cfg['model'],'fits':fits,
                'summaries':summarize(fits,cfg),'must_mention':[],'failing_items':[]}
        if cfg['model']=='competition_exact':
            # Never label Ki summaries as KD.
            sf=[dict(v,kd_M=v['ki_M']) for v in fits]; result['summaries']=summarize(sf,cfg)
            for s in result['summaries']: s['geometric_mean_ki_M']=s.pop('geometric_mean_kd_M')
            if cfg['competition']['cheng_prusoff']:
                from .binding_models import cheng_prusoff
                q=cfg['competition']; result['cheng_prusoff_approximation_M']=cheng_prusoff(q['ic50_M'],q['tracer_total_M'],q['tracer_kd_M'],cfg['active_sites']['concentration_M'])
                result['must_mention'].append('Cheng-Prusoff is a separate approximation diagnostic; exact competition is the primary model.')
        if cfg['model']=='solution_equilibrium_titration' and cfg['valency']['constant_species']=='bivalent':
            result['must_mention'].append('Bivalent constant species analysed as independent binding sites; this assumes the '
                                          'declared readout is proportional to free sites, not to molecules with any free site.')
        try:
            from .evidence_0110 import EVIDENCE
            ev=EVIDENCE.get(cfg['model'],{})
        except ImportError: ev={'status':'pending','must_mention':['0.11.0 validation is pending; no calibration claim.']}
        result['validation_evidence']=ev; result['must_mention']+=ev.get('must_mention',[])
        result['failing_items'] += [{'source':'results.json#/validation_evidence','reason':'calibration_miss','scenario':v['scenario']} for v in ev.get('calibration_rows',[]) if v.get('passed') is False]
        dump(out/'results.json',result)
        dump(out/'diagnostics.json',{f['curve_id']:f['diagnostics'] for f in fits})
        pd.DataFrame([{k:v for k,v in f.items() if k not in ('profile','predictions','parameters')} for f in fits]).to_csv(out/'fit_results.csv',index=False)
        pd.DataFrame(result['summaries']).to_csv(out/'sample_summary.csv',index=False)
        dump(out/'profiles.json',{f['curve_id']:f['profile'] for f in fits})
        pd.DataFrame([p for f in fits for p in f['predictions']]).to_csv(out/'predictions.csv',index=False)
        dump(out/'preprocessing_log.json',{'source':cfg['source'],'exclusions':d.loc[d.exclude,['observation_id','exclusion_reason']].to_dict('records'),
            'automatic_outlier_removal':False,'weighting':cfg['fit']['weighting'],'background':cfg.get('background'),
            'concentration_units':'M internally; original units retained','replicates':cfg['replicates']})
        write_facts(out)
        (out/'rerun.txt').write_text('agentic-prism analyze --config config.resolved.json --output ../affinity-rerun-new\n')
        dump(out/'manifest.json',{'schema_version':1,'analysis_type':cfg['analysis_type'],'package_version':__version__,
            'created_utc':datetime.now(timezone.utc).isoformat(),'input_sha256':sha(inp),'source':cfg['source'],
            'python':platform.python_version(),'platform':platform.platform(),'implementation_sha256':implementation_hash(),
            'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','scipy','pandas','matplotlib')},
            'scientific_artifacts_sha256':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}})
    except Exception as exc:
        dump(out/'failure.json',{'status':'failed','exception':type(exc).__name__,'message':str(exc)}); raise
    if render: render_affinity(out)
    return out


def render_affinity(run,style=None):
    run=Path(run); verify_run(run); cfg=json.loads((run/'config.resolved.json').read_text())
    facts=json.loads((run/'interpretation_facts.json').read_text()); r=json.loads((run/'results.json').read_text())
    style=style or cfg['report']['plot_style']
    if style not in ('standard','prism_like'): raise ValueError('Unsupported style')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    from . import report_shell as shell
    token=pstyle.THEMES[style]
    (run/'figures').mkdir(exist_ok=True); pictures=[]
    for i,f in enumerate(r['fits']):
        with plt.rc_context(pstyle.rc(token,*font_setup())):
            fig,ax=plt.subplots(figsize=(6.2,3.9),layout='constrained')
            for j,(cid,series) in enumerate(dict.fromkeys((p['curve_id'],p.get('series','binding')) for p in f['predictions'])):
                ps=[p for p in f['predictions'] if p['curve_id']==cid and p.get('series','binding')==series]; xx=np.array([p['concentration_M'] for p in ps]); order=np.argsort(xx)
                ax.plot(xx,[p['response'] for p in ps],label=cid+' / '+series,**pstyle.point_style(token,j,4.5))
                # Draw saved predictions only, never refit during rendering.
                ax.plot(xx[order],np.array([p['predicted_response'] for p in ps])[order],color=pstyle.color(token,j),lw=token['line_width'],alpha=.9)
            ax.set_xscale('symlog',linthresh=max(f.get('min_positive_M',1e-15),1e-30))
            ax.set(xlabel='Total ligand (M)',ylabel=f['response_unit'],title=f['curve_id']+' — '+f['interpretation']+' ('+f['status']+')')
            if f['predictions']: ax.legend(fontsize=7)
            pstyle.save(fig,run/'figures'/f'affinity-{i+1:03d}__{style}',token,cfg['report']['png_dpi'])
        pictures.append(shell.single_figure(run/'figures',f'affinity-{i+1:03d}__{style}',f['curve_id']))
    label='细胞结合：表观 KD' if cfg['analysis_type']=='cell_binding' else '平衡结合亲和力'
    page='<section id="must" class="card"><h2>必须披露</h2><ul class="must">'+''.join('<li>'+html.escape(v)+'</li>' for v in facts['must_mention'])+'</ul></section>'
    def molar(v): return '—' if v is None else f'{v:.3g} M'
    head=['曲线','状态','解释','点估计（仅可报告时）','95% 区间','可引用上界','区间状态','Pt/KD','诊断','设计提示']
    body=''
    for v in facts['primary']:
        value=v['ki_M'] if v['interpretation']=='Ki' else v['kd_M']
        ci='—' if v['ci_low_M'] is None and v['ci_high_M'] is None else f"{molar(v['ci_low_M'])} – {molar(v['ci_high_M'])}"
        ratio='; '.join(f"{r['pt_over_kd']:.3g} ({r['regime']})" for r in v['regime_diagnostics']) or '—'
        cells=[v['curve_id'],v['status'],v['interpretation'],molar(value),ci,molar(v['supported_upper_bound_M']),v['ci_status'],ratio,
               ', '.join(v['diagnostics']) or '—',', '.join(w['code'] for w in v.get('design_warnings',[])) or '—']
        body+='<tr>'+''.join('<td>'+html.escape(str(c))+'</td>' for c in cells)+'</tr>'
    page+=('<section id="results" class="card"><h2>可报告结果与未通过项目</h2><p>未达到可报告状态的曲线不显示点估计和区间；优化器最优值仅在 results.json 中作审计用途。</p>'
           '<div class="table-wrap"><table><tr>'+''.join('<th>'+h+'</th>' for h in head)+'</tr>'+body+'</table></div>'
           +shell.json_block(facts['primary'],'interpretation_facts.json 原文')+'</section>'
           +shell.section('figures','图形',''.join(pictures)))
    page+=shell.section('summary','独立实验汇总',shell.json_block(r['summaries']))
    page+=shell.section('config','配置与证据',shell.json_block({'config':cfg,'evidence':r['validation_evidence']},'完整配置与验证证据'))
    page+=shell.section('downloads','下载与复现',shell.downloads([(p.name,shell.data_uri(p)) for p in sorted(run.iterdir()) if p.name in verify_run(run)['scientific_artifacts_sha256']]))
    page=shell.page(label,eyebrow='AgenticPrism / '+cfg['analysis_type'].replace('_',' '),heading=label,
                    lede='结合模型、可报告状态与限制共同呈现；未通过项目只保留审计值。',
                    nav=[('must','必须披露'),('results','结果'),('figures','图形'),('summary','汇总'),('config','配置'),('downloads','下载')],body=page,
                    footer=f'AgenticPrism · 图形风格 {html.escape(style)} · 数值分析与图形渲染分离')
    (run/'report.html').write_text(page)
    dump(run/'render_manifest.json',{'style':style,'scientific_results_unchanged':True,'report_sha256':sha(run/'report.html')})
    return run/'report.html'
