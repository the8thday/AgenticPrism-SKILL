"""ADA assay sensitivity and drug tolerance from positive-control experiments.

Both read a cut point that was established separately (for example by an
ada_cut_point run) and never re-estimate it. Decision values use the declared
deployment (raw signal, signal/NC or signal-NC), matching that cut point.

Sensitivity: per run, the PC concentration where the mean decision value first
reaches the cut point for good (all higher concentrations stay at or above it),
by linear interpolation on log concentration between the bracketing levels.
Across runs, a log-normal upper prediction limit for a future run's sensitivity,
exp(m + t(conf, n-1) s sqrt(1+1/n)), is reported as the run-consistent value.
No extrapolation: runs positive at the lowest or negative at the highest tested
concentration are censored and withhold the prediction limit.

Drug tolerance: for each PC level and run, the highest drug concentration at
which the PC still reads at or above the cut point, interpolated on log drug
concentration; bounded results are reported as censored, never extrapolated.
"""
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
import html
import importlib.metadata
import json
import platform
import shutil
import numpy as np
import pandas as pd
from scipy import stats
from . import __version__
from .workflow import dump, sha, implementation_hash, verify_run

TYPES = ('ada_sensitivity', 'ada_drug_tolerance')
DEFAULTS = {'analysis_type': None, 'schema_version': 1, 'input': None, 'source': None,
    'assay': {'signal_unit': None, 'pc_unit': None, 'drug_unit': None, 'positive_control': None,
              'independent_runs': None, 'rationale': None},
    'cut_point': {'value': None, 'deployment': None, 'source': None},
    'statistics': {'prediction_confidence': None},
    'report': {'plot_style': 'prism_like'}}


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported ADA performance configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] not in TYPES or type(c['schema_version']) is not int or c['schema_version'] != 1:
        raise ValueError('Unsupported analysis/schema')
    a, cp, s = c['assay'], c['cut_point'], c['statistics']
    need = ['signal_unit', 'pc_unit', 'positive_control', 'rationale'] + (['drug_unit'] if c['analysis_type'] == 'ada_drug_tolerance' else [])
    if any(not isinstance(c[k], str) or not c[k].strip() for k in ('input', 'source')) or any(not isinstance(a[k], str) or not a[k].strip() for k in need):
        raise ValueError('input, source, units, positive_control and rationale are required strings')
    if c['analysis_type'] == 'ada_sensitivity' and a['drug_unit'] is not None:
        raise ValueError('drug_unit applies only to ada_drug_tolerance')
    if a['independent_runs'] is not True:
        raise ValueError('Declare independent_runs as literal true with evidence in the rationale')
    v = cp['value']
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v):
        raise ValueError('cut_point.value must be a finite number on the deployment scale')
    if cp['deployment'] not in ('raw', 'ratio', 'difference'):
        raise ValueError('cut_point.deployment must be raw, ratio or difference (as established for the cut point)')
    if not isinstance(cp['source'], str) or not cp['source'].strip():
        raise ValueError('cut_point.source must identify where the cut point came from (for example an ada_cut_point run and its hash)')
    if cp['deployment'] == 'ratio' and v <= 0:
        raise ValueError('A ratio cut point must be positive')
    conf = s['prediction_confidence']
    if c['analysis_type'] == 'ada_sensitivity':
        if isinstance(conf, bool) or not isinstance(conf, (int, float)) or not .5 < conf < 1:
            raise ValueError('Declare statistics.prediction_confidence between .5 and 1 (for example .95)')
    elif conf is not None:
        raise ValueError('prediction_confidence applies only to ada_sensitivity')
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Unsupported plot style')
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    drug = cfg['analysis_type'] == 'ada_drug_tolerance'
    required = {'observation_id', 'run_id', 'pc_concentration', 'signal', 'exclude', 'exclusion_reason'} | ({'drug_concentration'} if drug else set()) \
        | ({'nc_signal'} if cfg['cut_point']['deployment'] != 'raw' else set())
    if not required <= set(d.columns) or not len(d):
        raise ValueError('Missing columns or empty data: ' + ', '.join(sorted(required - set(d.columns))))
    if d.observation_id.duplicated().any():
        raise ValueError('Duplicate observation_id')
    for k in required - {'exclusion_reason'}:
        if (d[k].str.strip() == '').any():
            raise ValueError(f'Empty {k}')
    if set(d.exclude) - {'true', 'false'}:
        raise ValueError('exclude must be literal true/false')
    d['exclude'] = d.exclude == 'true'
    if (d.exclude & (d.exclusion_reason.str.strip() == '')).any():
        raise ValueError('Excluded rows require a reason')
    for k in required & {'pc_concentration', 'signal', 'nc_signal', 'drug_concentration'}:
        d[k] = pd.to_numeric(d[k], errors='raise')
        if not np.isfinite(d[k]).all():
            raise ValueError(f'Non-finite {k}')
    if (d.pc_concentration <= 0).any():
        raise ValueError('PC concentrations must be positive (put negative-control wells in nc_signal)')
    if drug and (d.drug_concentration < 0).any():
        raise ValueError('Drug concentrations must be nonnegative')
    if 'nc_signal' in d and ((d.groupby('run_id').nc_signal.nunique() != 1).any() or (d.nc_signal <= 0).any()):
        raise ValueError('nc_signal must be the run NC aggregate: positive and constant within a run')
    return d


def decision(g, cfg):
    cp = cfg['cut_point']
    if cp['deployment'] == 'raw':
        return g.signal
    return g.signal/g.nc_signal if cp['deployment'] == 'ratio' else g.signal-g.nc_signal


def crossing(x, v, cut, increasing=True):
    """Stable crossing on log x. Returns (value, status)."""
    x, v = np.asarray(x, float), np.asarray(v, float)
    pos = v >= cut
    if increasing:
        if pos.all():
            return None, 'positive_at_lowest_tested'
        if not pos[-1]:
            return None, 'negative_at_highest_tested'
        k = len(pos)-1
        while k > 0 and pos[k-1]:
            k -= 1
        i, j = k-1, k
    else:
        if not pos[0]:
            return None, 'negative_at_lowest_tested'
        if pos.all():
            return None, 'positive_at_highest_tested'
        k = 0
        while k < len(pos)-1 and pos[k+1]:
            k += 1
        i, j = k, k+1
    lx = np.log(x[[i, j]])
    w = (cut-v[i])/(v[j]-v[i])
    return float(np.exp(lx[0]+w*(lx[1]-lx[0]))), 'interpolated'


def sensitivity(d, cfg):
    used = d[~d.exclude].copy(); used['decision'] = decision(used, cfg); cut = cfg['cut_point']['value']; rows, diag = [], []
    for run, g in used.groupby('run_id', sort=True):
        m = g.groupby('pc_concentration').decision.mean().sort_index()
        if len(m) < 3:
            raise ValueError(f'Run {run}: at least three PC concentrations are required')
        value, status = crossing(m.index.to_numpy(), m.to_numpy(), cut, increasing=True)
        nonmono = bool(np.count_nonzero(np.diff((m.to_numpy() >= cut).astype(int))) > 1)
        rows.append({'run_id': run, 'sensitivity': value, 'status': status, 'n_levels': len(m), 'non_monotone_response': nonmono,
                     'lowest_tested': float(m.index.min()), 'highest_tested': float(m.index.max())})
        if nonmono:
            diag.append(f'Run {run}: the mean response crosses the cut point more than once; the stable (highest) crossing is used. Check the dilution series.')
    ok = [r['sensitivity'] for r in rows if r['status'] == 'interpolated']
    reasons = []
    if len(rows) < 3:
        reasons.append('At least three independent runs are needed for a run-to-run prediction limit.')
    if len(ok) < len(rows):
        reasons.append('Some runs are censored (outside the tested PC range); the prediction limit is withheld rather than extrapolated. Extend the dilution range.')
    conf = cfg['statistics']['prediction_confidence']; summary = {'n_runs': len(rows), 'n_interpolated': len(ok)}
    if len(ok) >= 2:
        lg = np.log(ok); m, s = float(lg.mean()), float(lg.std(ddof=1)); n = len(ok)
        summary.update(geometric_mean=float(np.exp(m)), log_sd=s, arithmetic_mean=float(np.mean(ok)), minimum=float(min(ok)), maximum=float(max(ok)))
        if not reasons:
            summary['upper_prediction_limit'] = float(np.exp(m+stats.t.ppf(conf, n-1)*s*np.sqrt(1+1/n)))
            summary['prediction_confidence'] = conf
    if len(rows) < 6:
        diag.append('Fewer than 6 runs; the prediction limit is wide and sensitive to the log-normal assumption.')
    return {'runs': rows, 'summary': summary, 'reportable': not reasons, 'withholding_reasons': reasons, 'diagnostics': diag,
            'method': 'Per-run stable crossing, linear in log concentration; log-normal upper prediction limit across runs'}


def drug_tolerance(d, cfg):
    used = d[~d.exclude].copy(); used['decision'] = decision(used, cfg); cut = cfg['cut_point']['value']; rows, diag = [], []
    for (pc, run), g in used.groupby(['pc_concentration', 'run_id'], sort=True):
        m = g.groupby('drug_concentration').decision.mean().sort_index()
        if 0 not in m.index or len(m) < 3:
            raise ValueError(f'PC {pc}, run {run}: need a no-drug reference and at least two nonzero drug levels')
        row = {'pc_concentration': float(pc), 'run_id': run, 'positive_without_drug': bool(m.loc[0] >= cut)}
        nz = m[m.index > 0]
        if not row['positive_without_drug']:
            row.update(tolerated_drug=None, status='not_detected_without_drug')
        else:
            value, status = crossing(nz.index.to_numpy(), nz.to_numpy(), cut, increasing=False)
            if status == 'negative_at_lowest_tested':
                status = 'below_lowest_nonzero_tested'
            row.update(tolerated_drug=value, status=status, lowest_nonzero=float(nz.index.min()), highest_tested=float(nz.index.max()))
            if np.count_nonzero(np.diff((nz.to_numpy() >= cut).astype(int))) > 1:
                diag.append(f'PC {pc}, run {run}: the response crosses the cut point more than once with increasing drug; the first loss is used. Check the plate.')
        rows.append(row)
    table = pd.DataFrame(rows); summary = []
    for pc, g in table.groupby('pc_concentration', sort=True):
        vals = g.tolerated_drug.dropna().to_numpy(float); st = g.status.value_counts().to_dict()
        entry = {'pc_concentration': float(pc), 'runs': len(g), 'status_counts': st}
        if (g.status == 'interpolated').all():
            entry.update(minimum_tolerated_drug=float(vals.min()), median_tolerated_drug=float(np.median(vals)))
        else:
            entry['note'] = 'At least one run is censored or undetected; per-run statuses are the result (no extrapolated summary).'
        summary.append(entry)
    return {'cells': rows, 'summary': summary, 'diagnostics': diag, 'reportable': True, 'withholding_reasons': [],
            'method': 'Per PC level and run: highest drug concentration still at or above the cut point, linear in log drug concentration'}


def analyze_ada_performance(config_path, output, render=True):
    src, out = Path(config_path).resolve(), Path(output).resolve()
    if out.exists():
        raise FileExistsError(f'Output already exists: {out}')
    out.mkdir(parents=True)
    try:
        cfg = resolve_config(json.loads(src.read_text()))
        inp = (src.parent/cfg['input']).resolve()
        d = load_data(inp, cfg)
        shutil.copyfile(inp, out/'input.csv'); cfg['input'] = 'input.csv'
        dump(out/'config.resolved.json', cfg)
        result = sensitivity(d, cfg) if cfg['analysis_type'] == 'ada_sensitivity' else drug_tolerance(d, cfg)
        result.update(schema_version=1, analysis_type=cfg['analysis_type'], cut_point=cfg['cut_point'],
                      excluded=d.loc[d.exclude, ['observation_id', 'exclusion_reason']].to_dict('records'),
                      limitations=['The cut point is taken as given from its declared source; its own uncertainty is not propagated.',
                                   'Positive-control sensitivity and drug tolerance depend on the PC reagent; they are not the sensitivity for patient antibodies.',
                                   'Interpolation is linear on log concentration between tested levels; no extrapolation beyond the tested range.'])
        if cfg['analysis_type'] == 'ada_sensitivity':
            from .evidence_093 import SENSITIVITY
            result['validation_evidence'] = SENSITIVITY
        dump(out/'results.json', result)
        facts = {'schema_version': 1, 'analysis_type': cfg['analysis_type'],
                 'source_sha256': {n: sha(out/n) for n in ('results.json', 'config.resolved.json')},
                 'primary': {'source': 'results.json#/summary', 'reportable': result['reportable'], 'value': result['summary'],
                             'cut_point': cfg['cut_point']},
                 'withholding_reasons': result['withholding_reasons'], 'must_mention': result['diagnostics'],
                 'limitations': result['limitations'],
                 'calibration': {'source': 'results.json#/validation_evidence', 'value': result.get('validation_evidence')}}
        dump(out/'interpretation_facts.json', facts)
        (out/'rerun.txt').write_text('agentic-prism analyze --config config.resolved.json --output ../rerun-new\n')
        dump(out/'manifest.json', {'schema_version': 1, 'analysis_type': cfg['analysis_type'], 'package_version': __version__,
             'created_utc': datetime.now(timezone.utc).isoformat(), 'input_original_path': str(inp), 'input_sha256': sha(inp),
             'source': cfg['source'], 'python': platform.python_version(), 'platform': platform.platform(),
             'dependencies': {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'pandas', 'matplotlib')},
             'implementation_sha256': implementation_hash(), 'prism_numerical_equivalence': 'not_claimed',
             'scientific_artifacts_sha256': {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()}})
    except Exception as exc:
        dump(out/'failure.json', {'status': 'failed', 'exception': type(exc).__name__, 'message': str(exc)})
        raise
    if render:
        render_ada_performance(out)
    return out


def render_ada_performance(run, style=None):
    run = Path(run)
    verify_run(run)
    cfg = json.loads((run/'config.resolved.json').read_text())
    style = style or cfg['report']['plot_style']
    if style not in ('standard', 'prism_like'):
        raise ValueError('Unsupported style')
    r = json.loads((run/'results.json').read_text())
    d = pd.read_csv(run/'input.csv', keep_default_na=False)
    d = d[d.exclude.astype(str).str.lower() != 'true'].copy()
    for k in ('signal', 'pc_concentration', 'nc_signal', 'drug_concentration'):
        if k in d:
            d[k] = pd.to_numeric(d[k])
    d['decision'] = decision(d, cfg)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    font, cjk = font_setup()
    with plt.rc_context({'font.family': [font]+([cjk] if cjk else [])}):
        fig, ax = plt.subplots(figsize=(8, 4), layout='constrained')
        if cfg['analysis_type'] == 'ada_sensitivity':
            for run_id, g in d.groupby('run_id'):
                m = g.groupby('pc_concentration').decision.mean()
                ax.plot(m.index, m.values, 'o-', alpha=.6, label=str(run_id))
            ax.set(xlabel=f"PC concentration ({cfg['assay']['pc_unit']})", title='ADA sensitivity')
        else:
            for (pc, run_id), g in d.groupby(['pc_concentration', 'run_id']):
                m = g[g.drug_concentration > 0].groupby('drug_concentration').decision.mean()
                ax.plot(m.index, m.values, 'o-', alpha=.6, label=f'PC {pc} / {run_id}')
            ax.set(xlabel=f"Drug concentration ({cfg['assay']['drug_unit']})", title='ADA drug tolerance')
        ax.axhline(cfg['cut_point']['value'], color='black', linestyle='--', label='Cut point')
        ax.set_xscale('log'); ax.set_ylabel(f"Decision value ({cfg['cut_point']['deployment']})")
        ax.legend(fontsize=7, ncol=2)
        ax.spines[['top', 'right']].set_visible(style == 'standard')
        (run/'figures').mkdir(exist_ok=True)
        for ext in ('svg', 'png', 'pdf'):
            fig.savefig(run/'figures'/f'ada_performance.{ext}', dpi=160)
        plt.close(fig)
    esc = html.escape
    title = 'ADA 灵敏度' if cfg['analysis_type'] == 'ada_sensitivity' else 'ADA 药物耐受'
    body = f'<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>{title}</title><style>body{{font:16px sans-serif;max-width:1100px;margin:30px auto}}pre{{white-space:pre-wrap}}svg{{width:100%;height:auto}}</style><h1>{title}</h1>'
    body += f'<p>切点 {esc(str(cfg["cut_point"]["value"]))}（{esc(cfg["cut_point"]["deployment"])}），来源：{esc(cfg["cut_point"]["source"])}。阳性对照结果不代表患者体内抗体的灵敏度。</p>'
    if r['withholding_reasons']:
        body += '<h2>暂不报告的原因</h2><ul>'+''.join(f'<li>{esc(x)}</li>' for x in r['withholding_reasons'])+'</ul>'
    body += (run/'figures/ada_performance.svg').read_text()
    body += '<h2>结果</h2><pre>'+esc(json.dumps(r, ensure_ascii=False, indent=2))+'</pre>'
    body += ''.join(f'<p><a href="{n}" download="{n}">{n}</a></p>' for n in ('interpretation_facts.json', 'results.json', 'config.resolved.json', 'input.csv'))+'</html>'
    (run/'report.html').write_text(body)
    dump(run/'render_manifest.json', {'style': style, 'scientific_artifacts_changed': False, 'report_sha256': sha(run/'report.html')})
    return run/'report.html'
