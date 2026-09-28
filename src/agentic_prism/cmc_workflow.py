"""CMC artifacts and facts extracted from saved results; rendering never refits."""
from datetime import datetime, timezone
from pathlib import Path
import html
import importlib.metadata
import json
import platform
import shutil
import numpy as np
from . import __version__, stability, potency_assay, comparability, specifications
from .workflow import dump, sha, implementation_hash, verify_run


def import_dose_runs(d, base, out, cfg):
    """Reuse verified dose_fit/equivalence results without a second curve fitter."""
    rows=[]; copied={}; seen=set(); references=set()
    for row in d.to_dict('records'):
        folder=(base/row['dose_run']).resolve(); manifest=verify_run(folder)
        source_cfg=json.loads((folder/'config.resolved.json').read_text())
        r=json.loads((folder/'results.json').read_text())
        if r['analysis_type']!='dose_response_4pl': raise ValueError('dose_run must be a verified dose-response analysis')
        key=(sha(folder/'results.json'), row['comparison_id'])
        if key in seen: raise ValueError('Duplicate source comparison cannot count as an independent determination')
        seen.add(key)
        comparisons=[c for c in r['comparisons'] if c['comparison_id']==row['comparison_id']]
        if len(comparisons)!=1: raise ValueError('comparison_id must identify exactly one dose comparison')
        c=comparisons[0]
        refkey=(sha(folder/'results.json'),c['reference_curve'])
        if refkey in references: raise ValueError('Shared reference fits induce correlated RP estimates; unsupported')
        references.add(refkey)
        fits={f['curve_id']:f for f in r['fits']}; rf=fits[c['reference_curve']]
        eq=c.get('parallelism_equivalence')
        if source_cfg['comparison_settings']['parallelism_method']!='equivalence':
            raise ValueError('Upstream dose-response must use prespecified equivalence parallelism')
        if folder not in copied:
            dest=out/'upstream'/str(len(copied)); dest.mkdir(parents=True)
            for name in [*manifest['scientific_artifacts_sha256'],'manifest.json']:
                target=dest/name; target.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(folder/name,target)
            copied[folder]=str(dest.relative_to(out))
        rp=c.get('relative_potency')
        # Failed fits have no RP. Keep their diagnostic and withhold the entire
        # analysis; never substitute a numerical potency for an unestimated fit.
        rows.append({**row,'dose_run':copied[folder], 'relative_potency':rp,
            'reference_quality':'pass' if rf['reportable'] and c['reportable'] else 'fail',
            'parallelism':'pass' if eq and eq.get('equivalent') is True else 'fail',
            'fit_source':copied[folder]+'/results.json#/comparisons/'+str(r['comparisons'].index(c)),
            'upstream_diagnostics':c['diagnostics']})
    import pandas as pd
    return pd.DataFrame(rows)


def write_facts(out):
    out=Path(out); r=json.loads((out/'results.json').read_text())
    dump(out/'interpretation_facts.json', {'schema_version':1,'analysis_type':r['analysis_type'],
        'source_sha256':{n:sha(out/n) for n in ('results.json','config.resolved.json')},
        'primary':{'source':'results.json#/primary',**r['primary']},
        'must_mention':r['must_mention'], 'failing_items':{'source':'results.json#/failing_items','items':r['failing_items']},
        'limitations':r['limitations'], 'calibration':{'source':'results.json#/validation_evidence','value':r['validation_evidence']}})


def analyze_cmc(config_path, output, render=True):
    src,out=Path(config_path).resolve(),Path(output).resolve()
    if out.exists(): raise FileExistsError('Choose a new CMC output directory')
    out.mkdir(parents=True)
    try:
        raw=json.loads(src.read_text()); module={'stability':stability,'potency_assay':potency_assay,'comparability':comparability,'specification':specifications}[raw['analysis_type']]
        cfg=module.resolve_config(raw); inp=(src.parent/cfg['input']).resolve(); d=module.load_data(inp,cfg)
        shutil.copyfile(inp,out/'input.csv'); cfg['input']='input.csv'
        if module is potency_assay and cfg['input_mode']=='dose_runs':
            shutil.copyfile(inp,out/'input.original.csv')
            d=import_dose_runs(d,src.parent,out,cfg)
            d[['observation_id','run_id','nominal_rp','dose_run','comparison_id']].to_csv(out/'input.csv',index=False)
            d.to_json(out/'per_run_import.json',orient='records',indent=2)
        dump(out/'config.resolved.json',cfg)
        if module is potency_assay and d.relative_potency.isna().any():
            # No numerical model can be fitted if even one upstream RP is absent.
            result={'analysis_type':'potency_assay','schema_version':1,'mode':cfg['mode'],
                'primary':{'status':'withheld','combined_rp':None,'reason':'upstream_RP_not_estimable'},
                'per_run':d.astype(object).where(d.notna(),None).to_dict('records'),'levels':[],'linearity':None,
                'failing_items':[{'observation_id':v['observation_id'],'reason':'upstream_RP_not_estimable'} for v in d[d.relative_potency.isna()].to_dict('records')],
                'must_mention':['A supplied upstream fit has no RP; all runs retained and combination withheld.'],
                'limitations':['Upstream fit failures must be resolved before combination.']}
        else: result=module.compute(d,cfg)
        from .evidence_010 import EVIDENCE
        from .evidence_0101 import EVIDENCE as EVIDENCE_0101
        evidence={**EVIDENCE,**EVIDENCE_0101}[result['analysis_type']]
        result['validation_evidence']=evidence
        result['must_mention'] += evidence['must_mention']
        result['failing_items'] += [{'reason':'calibration_miss','scenario':r['scenario'],'metric':r['metric']} for r in evidence['calibration_rows'] if not r['passed']]
        dump(out/'results.json',result); write_facts(out)
        (out/'rerun.txt').write_text('agentic-prism analyze --config config.resolved.json --output ../cmc-rerun-new\n')
        dump(out/'manifest.json',{'schema_version':1,'analysis_type':cfg['analysis_type'],'package_version':__version__,
            'created_utc':datetime.now(timezone.utc).isoformat(),'input_sha256':sha(inp),'source':cfg['source'],
            'python':platform.python_version(),'platform':platform.platform(),
            'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','scipy','pandas','matplotlib')},
            'implementation_sha256':implementation_hash(),'prism_numerical_equivalence':'not_claimed',
            'scientific_artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}})
    except Exception as exc:
        dump(out/'failure.json',{'status':'failed','exception':type(exc).__name__,'message':str(exc)}); raise
    if render: render_cmc(out)
    return out


def render_cmc(run, style=None):
    run=Path(run); verify_run(run)
    r=json.loads((run/'results.json').read_text()); cfg=json.loads((run/'config.resolved.json').read_text())
    style=style or cfg['report']['plot_style']
    if style not in ('standard','prism_like'): raise ValueError('Unsupported plot style')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    if r['analysis_type']=='stability':
        d=pd.read_csv(run/'input.csv',dtype={'batch':str}); cap=r['decision_tree']['allowed_horizon_months']; x=np.linspace(0,cap,150)
        for line in r['batch_estimates']:
            g=d[d.batch==line['batch']]; color=ax.scatter(g.time,g.value,label=line['batch']).get_facecolor()[0]
            ax.plot(x,line['coefficients'][0]+line['coefficients'][1]*x,color=color)
            for side in ('lower','upper'):
                if r['specification'][side] is not None:
                    ax.plot(x,[stability.mean_bound(line,t,side,cfg['assay']['direction']=='two_sided') for t in x],ls='--',color=color)
        for side in ('lower','upper'):
            if r['specification'][side] is not None: ax.axhline(r['specification'][side],color='black',ls=':')
        ax.axvline(r['decision_tree']['observed_common_duration_months'],ls=':',color='grey')
        ax.set(xlabel='Time (months)',ylabel=cfg['assay']['attribute']+' ('+cfg['assay']['unit']+')',title='Batch mean and confidence bounds'); ax.legend()
    elif r['analysis_type']=='comparability':
        lots=pd.DataFrame(r['lots']); ref=lots[lots['product']=='reference']; tst=lots[lots['product']=='test']
        ax.scatter(np.zeros(len(ref)),ref.lot_mean,color='grey',label='Reference lots'); ax.scatter(np.ones(len(tst)),tst.lot_mean,color='black',label='Test lots')
        det=r['detail'] or {}
        if 'range' in det:
            for v in det['range']: ax.axhline(v,ls='--',color='grey')
        ax.set_xticks([0,1],['Reference','Test']); ax.set(ylabel=cfg['attribute']['name']+' ('+cfg['attribute']['unit']+')',title='Lot means'); ax.legend()
        if 'interval' in det:
            ax2=ax.inset_axes([.62,.1,.35,.3]); ax2.errorbar([det['difference_test_minus_reference']],[0],xerr=[[det['difference_test_minus_reference']-det['interval'][0]],[det['interval'][1]-det['difference_test_minus_reference']]],fmt='o',color='black')
            for v in (-det['margin'],det['margin']): ax2.axvline(v,ls='--',color='grey')
            ax2.set_yticks([]); ax2.set_title('Difference and margin',fontsize=8)
    elif r['analysis_type']=='specification':
        d=pd.read_csv(run/'input.csv'); d=d[d.exclude.astype(str).str.lower()!='true']
        ax.hist(d.value,bins=min(30,max(5,len(d)//3)),color='lightgrey',edgecolor='grey')
        ti=r['tolerance_interval']
        if ti and ti.get('limits'):
            for v in ti['limits']:
                if v is not None: ax.axvline(v,color='black',ls='--')
        cap=r['capability']
        if cap:
            for k in ('lsl','usl'):
                if cap[k] is not None: ax.axvline(cap[k],color='black',ls=':')
        ax.set(xlabel=cfg['attribute']['name']+' ('+cfg['attribute']['unit']+')',ylabel='Units',title='Tolerance interval (--) and limits (:)')
    elif r['levels']:
        rows=r['levels']; y=np.array([v['combined_rp'] for v in rows]); ci=np.array([v['combined_rp_interval'] for v in rows])
        ax.errorbar([v['nominal_rp'] for v in rows],y,yerr=np.array([y-ci[:,0],ci[:,1]-y]),fmt='o',color='black',capsize=4)
        ax.set(xlabel='Nominal relative potency',ylabel='Combined RP with interval',title='Random-run log RP model')
    else:
        ax.text(.5,.5,'Combined RP withheld\nInspect per-run diagnostics',ha='center',va='center',transform=ax.transAxes); ax.set_axis_off()
    ax.spines[['top','right']].set_visible(style=='standard'); (run/'figures').mkdir(exist_ok=True)
    for ext in ('svg','pdf','png'): fig.savefig(run/'figures'/f'cmc.{ext}',dpi=160)
    plt.close(fig)
    esc=html.escape
    body='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CMC 分析</title><style>body{font:16px sans-serif;max-width:1000px;margin:30px auto;padding:16px}pre{white-space:pre-wrap;overflow-wrap:anywhere}svg{width:100%;height:auto}</style>'
    body+='<h1>'+{'stability':'稳定性与有效期','potency_assay':'跨运行相对效价','comparability':'可比性 / 生物类似性','specification':'容忍区间与过程能力'}[r['analysis_type']]+'</h1>'
    body+='<h2>主要结果</h2><pre>'+esc(json.dumps(r['primary'],indent=2,ensure_ascii=False))+'</pre>'
    body+='<h2>必须说明的事项</h2><ul>'+''.join('<li>'+esc(v)+'</li>' for v in r['must_mention'])+'</ul>'
    body+=(run/'figures/cmc.svg').read_text()+'<h2>详细结果与证据</h2><pre>'+esc(json.dumps(r,indent=2,ensure_ascii=False))+'</pre>'
    body+=''.join(f'<p><a href="{n}" download="{n}">{n}</a></p>' for n in ('results.json','interpretation_facts.json','config.resolved.json','input.csv'))+'</html>'
    (run/'report.html').write_text(body)
    dump(run/'render_manifest.json',{'style':style,'scientific_artifacts_changed':False,'report_sha256':sha(run/'report.html')})
    return run/'report.html'
