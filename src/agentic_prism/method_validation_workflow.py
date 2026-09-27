"""Method-validation artifacts, interpretation facts and a render-only report."""
from datetime import datetime, timezone
from pathlib import Path
import html
import importlib.metadata
import json
import platform
import shutil
from . import __version__, method_validation
from .workflow import dump, sha, implementation_hash, verify_run

HEADLINE = {'accuracy_precision': 'all_levels_pass', 'dilution_linearity': 'passes_without_hook',
            'parallelism': 'all_samples_pass_cv', 'selectivity': 'all_groups_pass', 'specificity': 'all_groups_pass',
            'stability': 'all_conditions_pass'}


def write_facts(out):
    out = Path(out)
    r = json.loads((out/'results.json').read_text())
    e = r['experiment']
    must = list(r['design_diagnostics'])
    if e == 'dilution_linearity' and r['hook_effect_suspected']:
        must.insert(0, 'Hook effect suspected: an above-ULOQ QC did not read above the ULOQ.')
    if e == 'accuracy_precision' and r['summary'].get('accuracy_profile', {}).get('status') == 'estimated':
        ap = r['summary']['accuracy_profile']
        if not ap['lower_is_lowest_tested'] or not ap['upper_is_highest_tested']:
            must.append(f"Accuracy profile (supplement) narrows the range to {ap['lower']:.4g}-{ap['upper']:.4g}; the level-based rule and the profile disagree at the ends.")
    if e == 'parallelism' and r['trend'] and r['trend']['p_value_slope_zero'] < .05:
        must.insert(0, 'Dilution trend: the common log slope differs from zero; passing the CV criterion can still hide non-parallelism.')
    if e == 'accuracy_precision':
        failing = [{'level': v['level'], 'failed_checks': [k for k, ok in v['checks'].items() if not ok]} for v in r['levels'] if not v['passes_criteria']]
    elif e == 'dilution_linearity':
        failing = [{'dilution_factor': v['dilution_factor'], 'accuracy_percent': v['accuracy_percent'], 'cv_percent': v['cv_percent']} for v in r['dilutions'] if v['passes_criteria'] is False]
    elif e == 'parallelism':
        failing = [{'sample_id': v['sample_id'], 'cv_percent': v['cv_percent']} for v in r['samples'] if v['passes_criteria'] is False]
    elif e in ('selectivity', 'specificity'):
        failing = [{k: v[k] for k in ('role', 'interferent', 'passing', 'sources', 'required_fraction') if k in v} for v in r['groups'] if not v['passes_criteria']]
    else:
        failing = [{k: v.get(k) for k in ('condition', 'level', 'accuracy_percent', 'note')} for v in r['conditions'] if not v['passes_criteria']]
    facts_failing = {'source': 'results.json', 'items': failing}
    facts = {'schema_version': 1, 'analysis_type': 'method_validation', 'experiment': e, 'failing_items': facts_failing,
             'source_sha256': {n: sha(out/n) for n in ('results.json', 'config.resolved.json')},
             'primary': {'source': f'results.json#/summary/{HEADLINE[e]}', 'meets_declared_criteria': r['summary'][HEADLINE[e]],
                         'criteria_source': r['criteria']['source'], 'summary': r['summary']},
             'must_mention': must, 'limitations': r['limitations'],
             'statistical_supplements': 'Intervals, tolerance intervals, slope trends and pass-rate intervals are supplements; the declared acceptance rule uses point estimates.',
             'calibration': {'source': 'results.json#/validation_evidence', 'value': r.get('validation_evidence')}}
    dump(out/'interpretation_facts.json', facts)


def analyze_method_validation(config_path, output, render=True):
    src, out = Path(config_path).resolve(), Path(output).resolve()
    if out.exists():
        raise FileExistsError(f'Output already exists: {out}')
    out.mkdir(parents=True)
    try:
        cfg = method_validation.resolve_config(json.loads(src.read_text()))
        inp = (src.parent/cfg['input']).resolve()
        d = method_validation.load_data(inp, cfg)
        shutil.copyfile(inp, out/'input.csv'); cfg['input'] = 'input.csv'
        dump(out/'config.resolved.json', cfg)
        result = method_validation.compute(d, cfg)
        from .evidence_093 import METHOD_VALIDATION
        result['validation_evidence'] = METHOD_VALIDATION
        dump(out/'results.json', result)
        write_facts(out)
        (out/'rerun.txt').write_text('agentic-prism analyze --config config.resolved.json --output ../rerun-new\n')
        dump(out/'manifest.json', {'schema_version': 1, 'analysis_type': 'method_validation', 'package_version': __version__,
             'created_utc': datetime.now(timezone.utc).isoformat(), 'input_original_path': str(inp), 'input_sha256': sha(inp),
             'source': cfg['source'], 'python': platform.python_version(), 'platform': platform.platform(),
             'dependencies': {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'pandas', 'matplotlib')},
             'implementation_sha256': implementation_hash(), 'prism_numerical_equivalence': 'not_claimed',
             'scientific_artifacts_sha256': {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()}})
    except Exception as exc:
        dump(out/'failure.json', {'status': 'failed', 'exception': type(exc).__name__, 'message': str(exc)})
        raise
    if render:
        render_method_validation(out)
    return out


def _figure(run, r, style):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    font, cjk = font_setup()
    e = r['experiment']
    with plt.rc_context({'font.family': [font]+([cjk] if cjk else [])}):
        fig, ax = plt.subplots(figsize=(8, 4), layout='constrained')
        if e == 'accuracy_precision':
            rows = r['levels']; x = [v['nominal'] for v in rows]
            ax.plot(x, [v['bias_percent'] for v in rows], 'o-', color='black', label='Mean bias %')
            if 'tolerance_interval' in rows[0]:
                ax.fill_between(x, [v['tolerance_interval']['lower_percent'] for v in rows], [v['tolerance_interval']['upper_percent'] for v in rows], alpha=.2, label='Beta-expectation interval')
            for v in rows:
                ax.hlines([-v['limits']['accuracy'], v['limits']['accuracy']], v['nominal']*.9, v['nominal']*1.1, colors='grey', linestyles='--')
            ax.set_xscale('log'); ax.set(xlabel='Nominal concentration', ylabel='Relative error (%)', title='Accuracy profile')
            ax.legend()
        elif e in ('dilution_linearity', 'parallelism'):
            key = 'dilutions' if e == 'dilution_linearity' else 'samples'
            if e == 'dilution_linearity':
                rows = [v for v in r[key] if v['in_range']]
                ax.plot([v['dilution_factor'] for v in rows], [v['accuracy_percent'] for v in rows], 'o-', color='black')
                ax.set(xlabel='Dilution factor', ylabel='Accuracy after dilution correction (%)', title='Dilution linearity')
                ax.set_xscale('log')
            else:
                rows = [v for v in r[key] if v.get('cv_percent') is not None]
                ax.bar([v['sample_id'] for v in rows], [v['cv_percent'] for v in rows], color='grey')
                ax.axhline(r['criteria']['parallelism_cv_percent'], color='black', linestyle='--')
                ax.set(xlabel='Sample', ylabel='CV of corrected concentrations (%)', title='Parallelism')
        else:
            rows = r['groups'] if 'groups' in r else r['conditions']
            labels = [str(v.get('role', v.get('condition')))+('' if 'level' not in v else ' / '+str(v['level'])) for v in rows]
            vals = [v.get('pass_fraction', v.get('accuracy_percent', 0)) or 0 for v in rows]
            ax.bar(labels, vals, color='grey')
            ax.set(ylabel='Pass fraction' if 'groups' in r else 'Accuracy (%)', title=r['experiment'].replace('_', ' ').title())
        ax.spines[['top', 'right']].set_visible(style == 'standard')
        (run/'figures').mkdir(exist_ok=True)
        for ext in ('svg', 'png', 'pdf'):
            fig.savefig(run/'figures'/f'method_validation.{ext}', dpi=160)
        plt.close(fig)


def render_method_validation(run, style=None):
    run = Path(run)
    verify_run(run)
    cfg = json.loads((run/'config.resolved.json').read_text())
    style = style or cfg['report']['plot_style']
    if style not in ('standard', 'prism_like'):
        raise ValueError('Unsupported style')
    r = json.loads((run/'results.json').read_text())
    _figure(run, r, style)
    esc = html.escape
    facts = json.loads((run/'interpretation_facts.json').read_text())
    body = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>方法学验证</title><style>body{font:16px sans-serif;max-width:1100px;margin:30px auto}pre{white-space:pre-wrap}svg{width:100%;height:auto}</style>'
    body += f'<h1>方法学验证：{esc(r["experiment"])}</h1>'
    body += f'<p>按声明的接受标准（来源：{esc(str(r["criteria"]["source"]))}）判定：<b>{"符合" if facts["primary"]["meets_declared_criteria"] else "不符合"}</b>。这只是验证的一部分，不等于完整的方法学验证或监管认可。</p>'
    if facts['must_mention']:
        body += '<h2>必须说明的事项</h2><ul>'+''.join(f'<li>{esc(m)}</li>' for m in facts['must_mention'])+'</ul>'
    body += (run/'figures/method_validation.svg').read_text()
    body += '<h2>结果</h2><pre>'+esc(json.dumps(r, ensure_ascii=False, indent=2))+'</pre>'
    body += ''.join(f'<p><a href="{n}" download="{n}">{n}</a></p>' for n in ('interpretation_facts.json', 'results.json', 'config.resolved.json', 'input.csv'))+'</html>'
    (run/'report.html').write_text(body)
    dump(run/'render_manifest.json', {'style': style, 'scientific_artifacts_changed': False, 'report_sha256': sha(run/'report.html')})
    return run/'report.html'
