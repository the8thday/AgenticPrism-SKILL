"""ADA artifacts and a render-only offline report."""
from datetime import datetime, timezone
from pathlib import Path
import html
import importlib.metadata
import json
import platform
import shutil
import pandas as pd
from . import __version__
from .workflow import dump, sha, implementation_hash, verify_run
from . import ada


def write_facts(out):
    out = Path(out)
    result = json.loads((out/'results.json').read_text())
    fit = result['fits'][0]
    facts = {"schema_version": 1, "analysis_type": "ada_cut_point", "estimand": "Declared upper percentile of drug-naive negative validation-panel responses",
        "source_sha256": {n: sha(out/n) for n in ('results.json', 'config.resolved.json')},
        "primary": {"source": "results.json#/fits/0", "reportable": fit['reportable'], "cut_point": fit['cut_point'],
            "tier": fit['tier'], "type": fit['selected_type'], "scale": fit['scale'], "deployment": fit['deployment'],
            "target_false_positive_rate": fit['false_positive_rate_target']},
        "withholding_reasons": fit['withholding_reasons'], "must_mention": fit['diagnostics'],
        "variance_components": {"source": "results.json#/fits/0/variance_components", "value": fit['variance_components']},
        "limitations": result['limitations'],
        "calibration": {"source": "results.json#/validation_evidence", "value": result["validation_evidence"]}}
    if 'cut_point_bound' in fit:
        facts['cut_point_bound']={'source':'results.json#/fits/0/cut_point_bound','value':fit['cut_point_bound']}
        facts['estimand']='Lower confidence bound on the declared marginal negative-response percentile; FPR at least target'
    dump(out/'interpretation_facts.json', facts)


def analyze_ada(config_path, output, render=True):
    src, out = Path(config_path).resolve(), Path(output).resolve()
    if out.exists():
        raise FileExistsError(f"Output already exists: {out}")
    out.mkdir(parents=True)
    try:
        cfg = ada.resolve_config(json.loads(src.read_text()))
        inp = (src.parent/cfg['input']).resolve()
        d = ada.load_data(inp, cfg)
        shutil.copyfile(inp, out/'input.csv')
        cfg['input'] = 'input.csv'
        dump(out/'config.resolved.json', cfg)
        result = ada.compute(d, cfg)
        pd.DataFrame(result.pop('cells')).to_csv(out/'analysis_cells.csv', index=False)
        pd.DataFrame(result['per_run_audit']).to_csv(out/'per_run_audit.csv', index=False)
        fit = result['fits'][0]
        from .ada_evidence import CALIBRATION
        if cfg['cut_point']['bound']=='lower':
            from .interval_evidence import ADA_BOUNDS as CALIBRATION
            if cfg['cut_point']['method']=='nonparametric' and cfg['cut_point']['nonparametric_bound']=='two_way_bootstrap':
                from .evidence_093 import TWO_WAY_BOOTSTRAP as CALIBRATION
        result['validation_evidence'] = CALIBRATION
        result['limitations'] += CALIBRATION['limitations']
        dump(out/'preprocessing_log.json', {"replicate_aggregation": "arithmetic means before inhibition, normalization and transformation",
            "declared_outlier_policy": cfg['outliers'], "analytical_flags": fit['analytical_flag_observations'],
            "biological_flags": fit['biological_flag_subjects'],
            "explicit_exclusions": d.loc[d.exclude, ['observation_id','exclusion_reason']].to_dict('records')})
        dump(out/'results.json', result)
        dump(out/'diagnostics.json', {k: fit[k] for k in ('reportable','withholding_reasons','diagnostics')})
        write_facts(out)
        (out/'rerun.txt').write_text('agentic-prism analyze --config config.resolved.json --output ../rerun-new\n')
        dump(out/'manifest.json', {"schema_version": 1, "analysis_type": 'ada_cut_point', 'package_version': __version__,
            'created_utc': datetime.now(timezone.utc).isoformat(), 'input_original_path': str(inp), 'input_sha256': sha(inp),
            'source': cfg['source'], 'python': platform.python_version(), 'platform': platform.platform(),
            'dependencies': {p: importlib.metadata.version(p) for p in ('numpy','scipy','pandas','matplotlib')},
            'implementation_sha256': implementation_hash(), 'prism_numerical_equivalence': 'not_claimed',
            'scientific_artifacts_sha256': {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()}})
    except Exception as exc:
        dump(out/'failure.json', {'status':'failed','exception':type(exc).__name__,'message':str(exc)})
        raise
    if render:
        render_ada(out)
    return out


def render_ada(run, style=None):
    run = Path(run)
    verify_run(run)
    cfg = json.loads((run/'config.resolved.json').read_text())
    style = style or cfg['report']['plot_style']
    if style not in ('standard','prism_like'):
        raise ValueError('Unsupported style')
    result = json.loads((run/'results.json').read_text())
    fit = result['fits'][0]
    d = pd.read_csv(run/'analysis_cells.csv', keep_default_na=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    font, cjk = font_setup()
    with plt.rc_context({'font.family':[font]+([cjk] if cjk else [])}):
        fig, ax = plt.subplots(figsize=(8,4), layout='constrained')
        runs = sorted(d.run_id.unique())
        for i, r in enumerate(runs):
            rows = d[(d.run_id == r)&(d.used.astype(str).str.lower() == 'true')]
            ax.scatter([i]*len(rows), rows.normalized_value, alpha=.3, s=10)
        if fit['reportable']:
            ax.axhline(fit['cut_point'], color='black', linestyle='--', label='Cut point'); ax.legend()
        ax.set_xticks(range(len(runs)), runs)
        ax.set(xlabel='Run', ylabel=f"Analysis response ({fit['scale']})", title=f"ADA {fit['tier']}: {fit['status']}")
        ax.spines[['top','right']].set_visible(style == 'standard')
        (run/'figures').mkdir(exist_ok=True)
        for ext in ('svg','png','pdf'):
            fig.savefig(run/'figures'/f'ada.{ext}', dpi=160)
        plt.close(fig)
    esc = html.escape
    svg = (run/'figures/ada.svg').read_text()
    body = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>ADA 切点分析</title><style>body{font:16px sans-serif;max-width:1100px;margin:30px auto}pre{white-space:pre-wrap}svg{width:100%;height:auto}</style><h1>ADA 切点分析</h1>'
    body += '<p>可报告切点：'+esc(str(fit['cut_point']))+'；状态：'+esc(fit['status'])+'</p>'
    body += '<p>切点用于分析判定，不代表 ADA 发生率、浓度或临床风险。审计切点不能覆盖暂停报告的判定。</p>'+svg
    body += '<h2>结果、诊断和适用范围</h2><pre>'+esc(json.dumps(result,ensure_ascii=False,indent=2))+'</pre>'
    body += '<h2>可追溯文件</h2>'+''.join(f'<p><a href="{n}" download="{n}">{n}</a></p>' for n in ('interpretation_facts.json','results.json','analysis_cells.csv','per_run_audit.csv','config.resolved.json','preprocessing_log.json'))+'</html>'
    (run/'report.html').write_text(body)
    dump(run/'render_manifest.json', {'style':style,'scientific_artifacts_changed':False,'report_sha256':sha(run/'report.html')})
    return run/'report.html'
