"""Exploratory principal component analysis and clustering of independent units (0.14.4).

Both analyses describe structure; neither tests a hypothesis. Rows are declared independent units;
columns are declared numeric variables. Rows with a missing value in any declared variable are
excluded and listed (no imputation). Optional declared log10 transform, then centering and, when
declared, scaling to unit variance (standard deviation with n - 1).

PCA (`pca`): singular value decomposition of the processed matrix; eigenvalues = s^2/(n - 1);
loadings are unit vectors (as R prcomp rotation) with a fixed sign convention (the largest absolute
loading of each component is positive); scores = U S. Components are retained by a declared integer
or by Horn's parallel analysis with column permutation (Buja and Eyuboglu 1992): component j is kept
while its eigenvalue exceeds the declared quantile of permuted eigenvalue j, stopping at the first
failure. Eigenvalues closer than North's rule (SE = lambda sqrt(2/n)) are flagged as possibly mixed.
Annotation columns (batch, group) are summarized against retained scores descriptively (eta squared
or correlation), never tested.

Clustering (`clustering`): k-means (k-means++ starts, Lloyd iterations, best of declared starts) or
agglomerative hierarchical clustering (scipy linkage: ward on Euclidean distance, average or complete
on Euclidean or 1 - Pearson correlation distance), on the processed variables or on declared leading
PC scores (k-means fixed points are refined by Hartigan single-point transfers). k is declared or chosen by the gap statistic (Tibshirani, Walther and Hastie 2001;
W = within-cluster sum of squares, reference uniform over the PC-aligned box as cluster::clusGap
"scaledPCA", s_k = sd sqrt(1 + 1/B), smallest k with gap(k) >= gap(k + 1) - s(k + 1)); k = 1 is a
possible answer. Cluster stability by Hennig's (2007) bootstrap Jaccard; silhouette widths
(Rousseeuw 1987) are descriptive.
"""
from copy import deepcopy
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform

TYPES = {'pca', 'clustering'}
COMMON = {'schema_version': 1, 'analysis_type': None, 'input': None, 'source': None,
          'design': {'unit': None, 'unit_rationale': None, 'purpose_rationale': None},
          'data': {'id_column': None, 'variables': None, 'annotations': [], 'log10': [], 'log_rationale': None},
          'scaling': {'method': None, 'rationale': None},
          'report': {'plot_style': 'prism_like'}}
PCA_DEFAULTS = {**deepcopy(COMMON), 'components': {'retain': 'parallel_analysis', 'replicates': 1000, 'seed': 14401, 'quantile': .95}}
CL_DEFAULTS = {**deepcopy(COMMON), 'clustering': {'method': None, 'linkage': None, 'distance': 'euclidean', 'space': 'variables', 'pca_components': None,
                                                  'rationale': None, 'starts': 100, 'seed': 14402},
               'k': {'fixed': None, 'select': None, 'max': 8, 'reference_replicates': 100},
               'stability': {'bootstrap_replicates': 100}}
LINKAGES = ('ward', 'average', 'complete')


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def _int(v, lo):
    return type(v) is int and v >= lo


def resolve_config(raw):
    if not isinstance(raw, dict) or raw.get('analysis_type') not in TYPES: raise ValueError('analysis_type must be pca or clustering')
    base = PCA_DEFAULTS if raw['analysis_type'] == 'pca' else CL_DEFAULTS
    if set(raw) - set(base): raise ValueError('Unsupported configuration keys: ' + ', '.join(sorted(set(raw) - set(base))))
    c = deepcopy(base)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1: raise ValueError('Unsupported schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + name)
    d = c['data']
    if not _text(d['id_column']): raise ValueError('Declare data.id_column (one row per independent unit)')
    if not isinstance(d['variables'], list) or len(d['variables']) < 2 or not all(_text(x) for x in d['variables']) or len(set(d['variables'])) != len(d['variables']):
        raise ValueError('Declare at least two distinct data.variables before looking at the results')
    for key in ('annotations', 'log10'):
        if not isinstance(d[key], list) or not all(_text(x) for x in d[key]): raise ValueError(f'data.{key} must be a list of column names')
    if set(d['annotations']) & set(d['variables']) or d['id_column'] in d['variables'] + d['annotations']:
        raise ValueError('Annotation columns describe units and are not analysed; they cannot also be variables or the id column')
    if set(d['log10']) - set(d['variables']): raise ValueError('data.log10 lists columns that are not variables')
    if d['log10'] and not _text(d['log_rationale']): raise ValueError('Declare data.log_rationale for the log10 transform')
    s = c['scaling']
    if s['method'] not in ('standardize', 'none') or not _text(s['rationale']):
        raise ValueError('Declare scaling.method (standardize when variables have different units or ranges; none only for commensurate variables) with a rationale')
    if c['analysis_type'] == 'pca':
        p = c['components']
        if not (p['retain'] == 'parallel_analysis' or _int(p['retain'], 1)): raise ValueError('components.retain must be parallel_analysis or a positive integer')
        if not _int(p['replicates'], 200) or type(p['seed']) is not int or type(p['quantile']) not in (int, float) or not .5 < p['quantile'] < 1:
            raise ValueError('components needs replicates >= 200, an integer seed and a quantile between 0.5 and 1')
    else:
        m = c['clustering']; k = c['k']
        if m['method'] not in ('kmeans', 'hierarchical') or not _text(m['rationale']): raise ValueError('Declare clustering.method (kmeans or hierarchical) with a rationale, before the results')
        if m['method'] == 'hierarchical':
            if m['linkage'] not in LINKAGES: raise ValueError('Declare clustering.linkage: ward, average or complete')
            if m['distance'] not in ('euclidean', 'correlation'): raise ValueError('distance must be euclidean or correlation')
            if m['linkage'] == 'ward' and m['distance'] != 'euclidean': raise ValueError('Ward linkage requires Euclidean distance')
        else:
            if m['linkage'] is not None or m['distance'] != 'euclidean': raise ValueError('k-means uses Euclidean distance and no linkage')
            if not _int(m['starts'], 10): raise ValueError('starts must be an integer of at least 10')
        if m['space'] not in ('variables', 'pca_scores'): raise ValueError('clustering.space must be variables or pca_scores')
        if m['space'] == 'pca_scores':
            if not _int(m['pca_components'], 2) or m['pca_components'] > len(d['variables']): raise ValueError('pca_scores needs pca_components (2 to the number of variables), declared before the results')
            if m['distance'] == 'correlation': raise ValueError('Correlation distance is not used on PC scores')
        elif m['pca_components'] is not None: raise ValueError('pca_components applies only to space pca_scores')
        if type(m['seed']) is not int: raise ValueError('clustering.seed must be an integer')
        if (k['fixed'] is None) == (k['select'] is None): raise ValueError('Declare exactly one of k.fixed (an integer >= 2) or k.select: gap')
        if k['fixed'] is not None and not _int(k['fixed'], 2): raise ValueError('k.fixed must be an integer of at least 2')
        if k['select'] is not None and (k['select'] != 'gap' or not _int(k['max'], 2) or not _int(k['reference_replicates'], 50)):
            raise ValueError('k.select must be gap with max >= 2 and reference_replicates >= 50')
        if not _int(c['stability']['bootstrap_replicates'], 50): raise ValueError('stability.bootstrap_replicates must be at least 50')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={c['data']['id_column']: str})
    dd = c['data']; need = [dd['id_column'], *dd['variables'], *dd['annotations']]
    miss = [x for x in need if x not in d.columns]
    if miss: raise ValueError('Missing columns: ' + ', '.join(miss))
    extra = set(d.columns) - set(need) - {'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns (declare them as variables or annotations, or remove them): ' + ', '.join(sorted(extra)))
    if d[dd['id_column']].isna().any() or d[dd['id_column']].duplicated().any(): raise ValueError('id_column must be unique and present: one row per independent unit; average technical replicates first by a declared rule')
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()): raise ValueError('Every excluded row needs an exclusion_reason')
    else: d['exclude'] = False
    for v in dd['variables']:
        x = pd.to_numeric(d[v], errors='coerce')
        if (x.isna() & d[v].notna()).any(): raise ValueError(f'{v} has non-numeric entries; text such as "<LLOQ" needs a censored-data analysis, not PCA or clustering')
        if not np.isfinite(x.dropna()).all(): raise ValueError(f'{v} has non-finite values')
        d[v] = x
    for v in dd['log10']:
        if (d.loc[~d.exclude, v].dropna() <= 0).any(): raise ValueError(f'{v} has non-positive values; log10 is undefined (declare no log, or handle values below a limit as censored elsewhere)')
    return d


def preprocess(d, c):
    dd = c['data']; u = d[~d.exclude]; complete = u[dd['variables']].notna().all(axis=1)
    dropped = u.loc[~complete, dd['id_column']].tolist(); u = u[complete]
    X = u[dd['variables']].to_numpy(float).copy()
    for j, v in enumerate(dd['variables']):
        if v in dd['log10']: X[:, j] = np.log10(X[:, j])
    center = X.mean(0); sd = X.std(0, ddof=1) if len(X) > 1 else np.zeros(X.shape[1])
    zero = [v for v, s in zip(dd['variables'], sd) if not s > 0]
    if zero: raise ValueError('Zero variance after exclusions (remove or do not declare): ' + ', '.join(zero))
    Z = (X - center)/(sd if c['scaling']['method'] == 'standardize' else 1.)
    return u, X, Z, center, sd, dropped


def svd_pca(Z):
    U, S, Vt = np.linalg.svd(Z, full_matrices=False); V = Vt.T
    for j in range(V.shape[1]):  # sign convention: largest absolute loading positive
        i = int(np.argmax(np.abs(V[:, j])))
        if V[i, j] < 0: V[:, j] *= -1; U[:, j] *= -1
    return S**2/(len(Z) - 1), V, U*S


def parallel_analysis(Z, replicates, seed, q):
    rng = np.random.default_rng(seed); n, p = Z.shape; ev = np.empty((replicates, min(n - 1, p)))
    for b in range(replicates):
        P = np.column_stack([rng.permutation(Z[:, j]) for j in range(p)]); P = P - P.mean(0)
        ev[b] = (np.linalg.svd(P, compute_uv=False)**2/(n - 1))[:ev.shape[1]]
    return np.quantile(ev, q, axis=0)


def _annotation_summary(u, scores, names, k):
    out = []
    for a in names:
        col = u[a]; num = pd.to_numeric(col, errors='coerce')
        if num.notna().all() and col.nunique() > 6:
            for j in range(k):
                r = float(np.corrcoef(num, scores[:, j])[0, 1]) if num.std() > 0 else None
                out.append({'annotation': a, 'component': f'PC{j+1}', 'type': 'numeric', 'pearson_r': r})
        else:
            g = col.astype(str).to_numpy()
            for j in range(k):
                s = scores[:, j]; tot = ((s - s.mean())**2).sum()
                between = sum((g == lv).sum()*(s[g == lv].mean() - s.mean())**2 for lv in np.unique(g))
                out.append({'annotation': a, 'component': f'PC{j+1}', 'type': 'categorical', 'eta_squared': float(between/tot) if tot > 0 else None,
                            'level_mean_scores': {lv: float(s[g == lv].mean()) for lv in np.unique(g)}, 'level_n': {lv: int((g == lv).sum()) for lv in np.unique(g)}})
    return out


def compute_pca(d, c):
    u, X, Z, center, sd, dropped = preprocess(d, c); n, p = Z.shape; dd = c['data']; failing = []; diag = []
    if n < 3: raise ValueError('PCA needs at least 3 complete units')
    ev, V, scores = svd_pca(Z); m = len(ev); tot = ev.sum()
    pc = c['components']
    if pc['retain'] == 'parallel_analysis':
        thr = parallel_analysis(Z, pc['replicates'], pc['seed'], pc['quantile']); k = 0
        while k < m and ev[k] > thr[k]: k += 1
        retain = {'rule': 'parallel_analysis', 'retained': k, 'thresholds': thr.tolist(), 'replicates': pc['replicates'], 'seed': pc['seed'], 'quantile': pc['quantile']}
        if k == 0: diag.append('no_component_exceeds_parallel_analysis')
    else:
        k = min(pc['retain'], m); retain = {'rule': 'declared', 'retained': k}
    se = ev*np.sqrt(2/n); close = [f'PC{j+1}-PC{j+2}' for j in range(min(max(k, 1), m - 1)) if ev[j] - ev[j + 1] < se[j]]
    if close: diag.append('eigenvalues_within_north_rule:' + ','.join(close))
    if n < 5*p: diag.append('fewer_than_5_units_per_variable')
    if c['scaling']['method'] == 'none':
        share = sd**2/np.sum(sd**2)
        if share.max() > .5: diag.append('one_variable_holds_over_half_the_variance_unscaled:' + dd['variables'][int(np.argmax(share))])
    names = dd['variables']; show = max(k, min(2, m))
    comps = [{'component': f'PC{j+1}', 'eigenvalue': float(ev[j]), 'sd': float(np.sqrt(ev[j])), 'proportion': float(ev[j]/tot), 'cumulative': float(ev[:j+1].sum()/tot),
              'north_se': float(se[j]), 'retained': j < k, 'loadings': {v: float(V[i, j]) for i, v in enumerate(names)},
              'correlations': {v: float(V[i, j]*np.sqrt(ev[j])/(sd[i] if c['scaling']['method'] == 'none' else 1.)) for i, v in enumerate(names)}} for j in range(m)]
    primary = {'status': 'estimated', 'reportable': True, 'n_units': n, 'n_variables': p, 'excluded_incomplete_units': dropped, 'scaling': c['scaling']['method'],
               'log10_variables': dd['log10'], 'center': dict(zip(names, center.tolist())), 'scale_sd': dict(zip(names, sd.tolist())),
               'components': comps, 'retention': retain, 'diagnostics': diag, 'annotation_summary': _annotation_summary(u, scores, dd['annotations'], show),
               'scores': [{'id': i, **{f'PC{j+1}': float(scores[r, j]) for j in range(m)}, **{a: (None if pd.isna(u[a].iloc[r]) else u[a].iloc[r]) for a in dd['annotations']}}
                          for r, i in enumerate(u[dd['id_column']])]}
    notes = [f"Exploratory PCA of {p} declared variables over {n} independent units ({c['design']['unit']}; {c['design']['unit_rationale']}); PCA describes correlation structure and tests nothing.",
             ('Variables were standardized (correlation-matrix PCA): ' if c['scaling']['method'] == 'standardize' else 'Variables were centred but not scaled (covariance PCA); variables with larger variance dominate: ') + c['scaling']['rationale'] + '.',
             'Component signs are arbitrary; here the largest absolute loading of each component is positive. Interpret a component from its loadings, not from its sign.']
    if dd['log10']: notes.append('log10 was applied to ' + ', '.join(dd['log10']) + ' (' + dd['log_rationale'] + ').')
    if dropped: notes.append(f'{len(dropped)} unit(s) with a missing value were excluded (no imputation): ' + ', '.join(map(str, dropped[:10])) + ('...' if len(dropped) > 10 else '') + '.')
    notes.append((f"Parallel analysis ({pc['replicates']} column permutations, {pc['quantile']:.0%} quantile, seed {pc['seed']}) retained {k} component(s)." if pc['retain'] == 'parallel_analysis'
                  else f'{k} component(s) retained by declaration.') + f' The first {max(k, 1)} explain {ev[:max(k, 1)].sum()/tot:.1%} of the total variance.')
    if close: notes.append('Neighbouring eigenvalues are within North\'s rule of thumb (' + ', '.join(close) + '); those components may mix or swap between samples, so do not interpret them separately.')
    if 'fewer_than_5_units_per_variable' in diag: notes.append(f'Only {n} units for {p} variables (fewer than 5 per variable); loadings are unstable.')
    if any(x.startswith('one_variable_holds') for x in diag): notes.append('Unscaled, one variable holds over half of the variance and will dominate PC1; standardize unless the variables share a unit and scale for a reason.')
    if dd['annotations']:
        notes.append('Annotations (' + ', '.join(dd['annotations']) + ') were not used to compute the components; their association with scores is descriptive (eta squared or correlation). A component aligned with batch, run or date reflects that factor, not biology.')
    return {'analysis_type': 'pca', 'primary': primary, 'must_mention': notes, 'failing_items': failing,
            'limitations': ['Linear PCA on complete cases; no imputation, robust PCA, factor rotation or sparse loadings.',
                            'Exploratory: no p-values, no confirmation of groups or mechanisms; components found in these data need independent confirmation.',
                            'No Prism or software-equivalence claim beyond the recorded benchmarks.']}


# ---- clustering ----

def kmeans(X, k, starts, rng, max_iter=300):
    """Best of `starts` k-means++ starts, each run by Lloyd iterations and then Hartigan single-point transfers."""
    best = None; n = len(X)
    for _ in range(starts):
        C = X[[rng.integers(n)]]
        for _ in range(1, k):  # k-means++
            d2 = ((X[:, None, :] - C[None])**2).sum(-1).min(1)
            C = np.vstack([C, X[rng.choice(n, p=d2/d2.sum()) if d2.sum() > 0 else rng.integers(n)]])
        lab = None
        for _ in range(max_iter):
            new = ((X[:, None, :] - C[None])**2).sum(-1).argmin(1)
            if lab is not None and (new == lab).all(): break
            lab = new
            for j in range(k):
                if (lab == j).any(): C[j] = X[lab == j].mean(0)
        lab = hartigan_refine(X, lab, k); w = within_ss(X, lab)
        if best is None or w < best[0] - 1e-12*max(1., abs(w)): best = (w, lab)
    return relabel(best[1])


def hartigan_refine(X, lab, k, max_moves=10000):
    """Hartigan (1975) single-point transfers: moving x from cluster i to j changes the within-cluster SS by
    n_j/(n_j + 1)|x - c_j|^2 - n_i/(n_i - 1)|x - c_i|^2. The best improving transfer is made until none remains.
    Lloyd's fixed points are not always transfer-stable, which left k-means short of R's Hartigan-Wong optimum."""
    lab = lab.copy(); n = np.bincount(lab, minlength=k).astype(float); C = np.array([X[lab == j].mean(0) if n[j] else X[0] for j in range(k)]); r = np.arange(len(X))
    for _ in range(max_moves):
        d2 = ((X[:, None, :] - C[None])**2).sum(-1); own = d2[r, lab]; no = n[lab]
        rem = np.where(no > 1, no/np.maximum(no - 1, 1)*own, -np.inf); add = n/(n + 1)*d2; add[r, lab] = np.inf
        j = add.argmin(1); gain = rem - add[r, j]; t = int(np.argmax(gain))
        if not gain[t] > 1e-12*max(rem[t], 1e-300): break
        i, jj = lab[t], j[t]; C[i] = (C[i]*n[i] - X[t])/(n[i] - 1); C[jj] = (C[jj]*n[jj] + X[t])/(n[jj] + 1); n[i] -= 1; n[jj] += 1; lab[t] = jj
    return lab


def relabel(lab):
    """Clusters numbered 1..k by decreasing size, ties by first appearance in row order."""
    ids = list(dict.fromkeys(lab.tolist())); size = {g: int((lab == g).sum()) for g in ids}
    order = sorted(ids, key=lambda g: (-size[g], ids.index(g))); m = {g: i + 1 for i, g in enumerate(order)}
    return np.array([m[g] for g in lab])


def within_ss(X, lab):
    return float(sum(((X[lab == g] - X[lab == g].mean(0))**2).sum() for g in np.unique(lab)))


def cluster_once(X, c, k, rng):
    m = c['clustering']
    if k == 1: return np.ones(len(X), int)
    if m['method'] == 'kmeans': return kmeans(X, k, m['starts'], rng)
    if m['distance'] == 'correlation': L = hierarchy.linkage(pdist(X, 'correlation'), method=m['linkage'])
    else: L = hierarchy.linkage(X, method=m['linkage'], metric='euclidean')
    return relabel(hierarchy.fcluster(L, k, criterion='maxclust'))


def gap_statistic(X, c, kmax, B, rng, return_references=False):
    n = len(X); Xc = X - X.mean(0); _, _, Vt = np.linalg.svd(Xc, full_matrices=False); R = Xc@Vt.T
    lo, hi = R.min(0), R.max(0); logW = []; refs = []
    for k in range(1, kmax + 1): logW.append(np.log(within_ss(X, cluster_once(X, c, k, rng))))
    sim = np.empty((B, kmax))
    for b in range(B):
        Y = rng.uniform(lo, hi, size=R.shape)@Vt + X.mean(0)
        if return_references: refs.append(Y)
        for k in range(1, kmax + 1): sim[b, k - 1] = np.log(within_ss(Y, cluster_once(Y, c, k, rng)))
    gap = sim.mean(0) - np.array(logW); s = sim.std(0, ddof=1)*np.sqrt(1 + 1/B)
    k_hat = next((k for k in range(1, kmax) if gap[k - 1] >= gap[k] - s[k]), kmax)
    out = {'k': list(range(1, kmax + 1)), 'logW': logW, 'E_logW_reference': sim.mean(0).tolist(), 'gap': gap.tolist(), 's': s.tolist(), 'selected_k': int(k_hat),
           'rule': 'smallest k with gap(k) >= gap(k+1) - s(k+1) (Tibshirani 2001)', 'reference': 'uniform over the PC-aligned bounding box (scaledPCA)', 'replicates': B}
    return (out, refs) if return_references else out


def silhouette(D, lab):
    n = len(lab); s = np.zeros(n); ks = np.unique(lab)
    if len(ks) < 2: return None
    for i in range(n):
        same = lab == lab[i]
        if same.sum() == 1: s[i] = 0.; continue
        a = D[i, same].sum()/(same.sum() - 1); b = min(D[i, lab == g].mean() for g in ks if g != lab[i])
        s[i] = (b - a)/max(a, b) if max(a, b) > 0 else 0.
    return s


def jaccard_stability(X, c, k, lab, B, rng):
    n = len(X); scores = {g: [] for g in np.unique(lab)}
    for _ in range(B):
        idx = rng.integers(0, n, n); bl = cluster_once(X[idx], c, k, rng); uniq = np.unique(idx)
        first = {i: bl[np.flatnonzero(idx == i)[0]] for i in uniq}
        bsets = {g: {i for i in uniq if first[i] == g} for g in np.unique(bl)}
        for g in scores:
            orig = {i for i in uniq if lab[i] == g}
            if not orig: continue
            scores[g].append(max(len(orig & s)/len(orig | s) for s in bsets.values()))
    return {int(g): {'mean_jaccard': float(np.mean(v)) if v else None, 'n_bootstrap_with_members': len(v)} for g, v in scores.items()}


def compute_clustering(d, c):
    u, X, Z, center, sd, dropped = preprocess(d, c); n, p = Z.shape; dd = c['data']; m = c['clustering']; kk = c['k']; diag = []; failing = []
    if n < 10: raise ValueError('Clustering needs at least 10 complete units')
    space = Z
    if m['space'] == 'pca_scores': _, _, sc = svd_pca(Z); space = sc[:, :m['pca_components']]
    rng = np.random.default_rng(m['seed']); gap = None
    if kk['select'] == 'gap':
        kmax = min(kk['max'], n - 1); gap = gap_statistic(space, c, kmax, kk['reference_replicates'], rng); k = gap['selected_k']
        if k == kmax: diag.append('gap_selected_the_largest_k_considered')
    else: k = kk['fixed']
    lab = cluster_once(space, c, k, rng)
    D = squareform(pdist(space, 'correlation' if m['distance'] == 'correlation' else 'euclidean'))
    sil = silhouette(D, lab) if k > 1 else None
    stab = jaccard_stability(space, c, k, lab, c['stability']['bootstrap_replicates'], rng) if k > 1 else {}
    sizes = {int(g): int((lab == g).sum()) for g in np.unique(lab)}
    if k > 1 and min(sizes.values()) < 5: diag.append('a_cluster_has_fewer_than_5_units')
    unstable = [g for g, v in stab.items() if v['mean_jaccard'] is not None and v['mean_jaccard'] < .6]
    if unstable: diag.append('unstable_clusters_jaccard_below_0.6:' + ','.join(map(str, unstable)))
    if n < 5*p and m['space'] == 'variables': diag.append('fewer_than_5_units_per_variable')
    names = dd['variables']
    profiles = [{'cluster': int(g), 'n': sizes[int(g)], 'mean': dict(zip(names, X[lab == g].mean(0).tolist())), 'sd': dict(zip(names, (X[lab == g].std(0, ddof=1) if sizes[int(g)] > 1 else np.full(p, np.nan)).tolist())),
                 'mean_silhouette': None if sil is None else float(sil[lab == g].mean()), **({'stability': stab[int(g)]} if stab else {})} for g in np.unique(lab)]
    crosstab = [{'annotation': a, 'counts': {str(g): u.loc[lab == g, a].astype(str).value_counts().to_dict() for g in np.unique(lab)}} for a in dd['annotations']]
    primary = {'status': 'estimated', 'reportable': True, 'n_units': n, 'n_variables': p, 'excluded_incomplete_units': dropped, 'method': m['method'], 'linkage': m['linkage'],
               'distance': m['distance'], 'space': m['space'], 'pca_components': m['pca_components'], 'k': int(k), 'k_rule': 'gap' if gap else 'declared', 'gap': gap,
               'mean_silhouette': None if sil is None else float(sil.mean()), 'cluster_sizes': sizes, 'profiles': profiles, 'annotation_crosstab': crosstab, 'diagnostics': diag,
               'assignments': [{'id': i, 'cluster': int(g), 'silhouette': None if sil is None else float(sil[r]), **{a: (None if pd.isna(u[a].iloc[r]) else u[a].iloc[r]) for a in dd['annotations']}}
                               for r, (i, g) in enumerate(zip(u[dd['id_column']], lab))]}
    desc = (f"k-means (best of {m['starts']} k-means++ starts, seed {m['seed']})" if m['method'] == 'kmeans' else f"hierarchical clustering ({m['linkage']} linkage, {m['distance']} distance)")
    notes = [f"Exploratory {desc} of {n} independent units ({c['design']['unit']}) on " + (f"the first {m['pca_components']} PC scores" if m['space'] == 'pca_scores' else f'{p} declared variables')
             + f" ({c['scaling']['method']}; {c['scaling']['rationale']}). Method declared before the results: {m['rationale']}.",
             'Clustering always returns groups, even from data without groups; the clusters are a description of these data, not discovered subtypes, until confirmed in independent data.',
             'Do not test the clustering variables between clusters (t tests, ANOVA): the clusters were built to differ on them, so such p-values are invalid. Cluster profiles here are descriptive.']
    if gap:
        notes.append(f"Gap statistic ({gap['replicates']} uniform references, {gap['rule']}) chose k = {k}" + (' - no cluster structure beyond a single group is supported; no clusters are reported.' if k == 1 else '.'))
        if 'gap_selected_the_largest_k_considered' in diag: notes.append('The gap statistic chose the largest k considered; the structure is not resolved within the declared range.')
    else: notes.append(f'k = {k} was declared, not estimated; the data were not asked whether {k} groups exist.')
    if stab:
        notes.append("Bootstrap stability (Hennig's Jaccard, " + str(c['stability']['bootstrap_replicates']) + ' resamples): ' + ', '.join(f"cluster {g} {v['mean_jaccard']:.2f}" for g, v in stab.items() if v['mean_jaccard'] is not None)
                     + '. Below 0.6 a cluster is not stable; above 0.75 it is reasonably stable; above 0.85 highly stable.')
    if sil is not None: notes.append(f'Mean silhouette width {sil.mean():.2f} (descriptive; above 0.5 suggests reasonable separation, below 0.25 little structure).')
    if 'a_cluster_has_fewer_than_5_units' in diag: notes.append('A cluster has fewer than 5 units; treat it as possible outliers rather than a group.')
    if dropped: notes.append(f'{len(dropped)} unit(s) with a missing value were excluded (no imputation).')
    if dd['log10']: notes.append('log10 was applied to ' + ', '.join(dd['log10']) + ' (' + dd['log_rationale'] + ').')
    if dd['annotations']: notes.append('Annotations (' + ', '.join(dd['annotations']) + ') were not used for clustering; cross-tabulations are descriptive. Clusters that follow batch or run reflect that factor.')
    return {'analysis_type': 'clustering', 'primary': primary, 'must_mention': notes, 'failing_items': failing,
            'limitations': ['k-means and agglomerative hierarchical clustering only; no model-based (Gaussian mixture), density-based or consensus clustering.',
                            'Complete cases only; Euclidean or correlation distance; results depend on the declared variables, transform and scaling.',
                            'No Prism or software-equivalence claim beyond the recorded benchmarks.']}


def compute(d, c):
    return compute_pca(d, c) if c['analysis_type'] == 'pca' else compute_clustering(d, c)


def write_artifacts(out, result):
    p = result['primary']
    if result['analysis_type'] == 'pca':
        pd.DataFrame(p['scores']).to_csv(Path(out)/'pca_scores.csv', index=False)
        pd.DataFrame([{'variable': v, **{c['component']: c['loadings'][v] for c in p['components']}} for v in p['center']]).to_csv(Path(out)/'pca_loadings.csv', index=False)
    else:
        pd.DataFrame(p['assignments']).to_csv(Path(out)/'cluster_assignments.csv', index=False)


def render(run, result, style):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        if result['analysis_type'] == 'pca':
            fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), layout='constrained'); ax = axes[0]; comps = p['components']
            ax.plot(range(1, len(comps) + 1), [c['eigenvalue'] for c in comps], 'o-', color=pstyle.color(token, 0), label='Observed')
            if p['retention'].get('thresholds'): ax.plot(range(1, len(comps) + 1), p['retention']['thresholds'], 's--', color=pstyle.MUTED, label='Parallel analysis')
            ax.set(xlabel='Component', ylabel='Eigenvalue', title='Scree'); ax.legend(frameon=False, fontsize=7)
            ax = axes[1]; sc = pd.DataFrame(p['scores']); ann = next(iter(p['annotation_summary']), None)
            key = ann['annotation'] if ann and ann['type'] == 'categorical' else None
            groups = sc.groupby(key, sort=False) if key else [('units', sc)]
            for i, (g, part) in enumerate(groups): ax.scatter(part.PC1, part.PC2, s=14, color=pstyle.color(token, i % 8), label=str(g))
            ax.axhline(0, color=pstyle.MUTED, lw=.6); ax.axvline(0, color=pstyle.MUTED, lw=.6)
            ax.set(xlabel=f"PC1 ({comps[0]['proportion']:.0%})", ylabel=f"PC2 ({comps[1]['proportion']:.0%})", title='Scores (exploratory)')
            if key: ax.legend(frameon=False, fontsize=7, title=key, title_fontsize=7)
            pstyle.save(fig, out/'pca', token, 300); return ['pca']
        a = pd.DataFrame(p['assignments']); n = 2 if p['gap'] else 1
        fig, axes = plt.subplots(1, n, figsize=(4.3*n, 3.4), layout='constrained', squeeze=False)
        if p['gap']:
            ax = axes[0, 0]; g = p['gap']; ax.errorbar(g['k'], g['gap'], yerr=g['s'], fmt='o-', color=pstyle.color(token, 0), capsize=2)
            ax.axvline(p['k'], ls=':', color=pstyle.MUTED); ax.set(xlabel='k', ylabel='Gap', title=f"Gap statistic (k = {p['k']})")
        ax = axes[0, -1]
        if p['mean_silhouette'] is not None:
            a = a.sort_values(['cluster', 'silhouette'], ascending=[True, False]).reset_index(drop=True)
            for g, part in a.groupby('cluster'): ax.barh(part.index, part.silhouette, color=pstyle.color(token, (g - 1) % 8), height=1.)
            ax.axvline(p['mean_silhouette'], ls=':', color=pstyle.MUTED); ax.set(xlabel='Silhouette width', yticks=[], title=f"k = {p['k']}, {p['method']}")
        else:
            ax.text(.5, .5, 'k = 1: no clusters reported', ha='center', va='center', transform=ax.transAxes); ax.set_axis_off()
        pstyle.save(fig, out/'clustering', token, 300); return ['clustering']
