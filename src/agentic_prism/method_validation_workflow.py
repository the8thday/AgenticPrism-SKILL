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
    from . import plot_style as pstyle
    font, cjk = font_setup()
    e = r['experiment']
    token = pstyle.THEMES[style]
    with plt.rc_context(pstyle.rc(token, font, cjk)):
        fig, ax = plt.subplots(figsize=(6.4, 3.8), layout='constrained')
        if e == 'accuracy_precision':
            rows = r['levels']; x = [v['nominal'] for v in rows]
            ax.plot(x, [v['bias_percent'] for v in rows], 'o-', color='black', mfc='black', label='Mean bias %')
            if 'tolerance_interval' in rows[0]:
                ax.fill_between(x, [v['tolerance_interval']['lower_percent'] for v in rows], [v['tolerance_interval']['upper_percent'] for v in rows], color=pstyle.color(token, 1), alpha=.18, lw=0, label='Beta-expectation interval')
            ax.axhline(0, color=pstyle.MUTED, lw=.8)
            for v in rows:
                ax.hlines([-v['limits']['accuracy'], v['limits']['accuracy']], v['nominal']*.9, v['nominal']*1.1, colors=pstyle.EXCLUDED, linestyles=(0, (4, 2)))
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
                ax.bar([v['sample_id'] for v in rows], [v['cv_percent'] for v in rows], color='#bfc5cc', edgecolor='black', linewidth=.8, width=.6)
                ax.axhline(r['criteria']['parallelism_cv_percent'], color=pstyle.EXCLUDED, linestyle=(0, (4, 2)))
                ax.set(xlabel='Sample', ylabel='CV of corrected concentrations (%)', title='Parallelism')
        else:
            rows = r['groups'] if 'groups' in r else r['conditions']
            labels = [str(v.get('role', v.get('condition')))+('' if 'level' not in v else ' / '+str(v['level'])) for v in rows]
            vals = [v.get('pass_fraction', v.get('accuracy_percent', 0)) or 0 for v in rows]
            ax.bar(labels, vals, color='#bfc5cc', edgecolor='black', linewidth=.8, width=.6)
            ax.set(ylabel='Pass fraction' if 'groups' in r else 'Accuracy (%)', title=r['experiment'].replace('_', ' ').title())
        (run/'figures').mkdir(exist_ok=True)
        pstyle.save(fig, run/'figures'/'method_validation', token, 300, ('svg', 'png', 'pdf'))


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
    from . import report_shell as shell
    meets = facts["primary"]["meets_declared_criteria"]
    body = (f'<section id="verdict" class="card"><div class="head"><h2>按声明标准判定</h2>'
            f'<span class="badge{"" if meets else " bad"}">{"符合" if meets else "不符合"}</span></div>'
            f'<p>按声明的接受标准（来源：{esc(str(r["criteria"]["source"]))}）判定：<b>{"符合" if meets else "不符合"}</b>。这只是验证的一部分，不等于完整的方法学验证或监管认可。</p>')
    if facts['must_mention']:
        body += '<h2>必须说明的事项</h2><ul class="must">'+''.join(f'<li>{esc(m)}</li>' for m in facts['must_mention'])+'</ul>'
    body += shell.single_figure(run/'figures', 'method_validation', '方法学验证 ' + r['experiment']) + '</section>'
    body += shell.section('details', '结果', shell.json_block(r, '完整 results.json', open_=True))
    body += shell.section('files', '可追溯文件', shell.downloads([(n, n) for n in ('interpretation_facts.json', 'results.json', 'config.resolved.json', 'input.csv')]))
    body = shell.page('方法学验证', eyebrow='AgenticPrism / method validation', heading=f'方法学验证：{r["experiment"]}',
                      lede='接受标准均由用户声明并注明来源；本报告不代表监管认可。',
                      nav=[('verdict', '判定'), ('details', '结果'), ('files', '可追溯文件')], body=body,
                      footer=f'AgenticPrism · 图形风格 {esc(style)} · 数值分析与图形渲染分离')
    (run/'report.html').write_text(body)
    dump(run/'render_manifest.json', {'style': style, 'scientific_artifacts_changed': False, 'report_sha256': sha(run/'report.html')})
    return run/'report.html'
