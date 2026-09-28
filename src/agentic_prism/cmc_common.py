"""Strict configuration helpers shared by the scoped CMC modules."""
from copy import deepcopy
import numpy as np
import pandas as pd


def merge(raw, defaults):
    if not isinstance(raw, dict) or set(raw) - set(defaults):
        raise ValueError('Unsupported configuration keys')
    out = deepcopy(defaults)
    for k, v in raw.items():
        if isinstance(out[k], dict):
            out[k] = merge(v, out[k])
        else:
            out[k] = v
    return out


def text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{name} requires a nonempty source or declaration')


def number(value, name, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or (positive and value <= 0):
        raise ValueError(f'{name} requires a finite' + (' positive' if positive else '') + ' number')


def common(cfg, kind):
    if cfg['analysis_type'] != kind or type(cfg['schema_version']) is not int or cfg['schema_version'] != 1:
        raise ValueError('Unsupported analysis_type/schema_version')
    for k in ('input', 'source'):
        text(cfg[k], k)
    if cfg['report']['plot_style'] not in ('standard', 'prism_like'):
        raise ValueError('Unsupported plot style')


def table(path, required, numeric):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    if not len(d) or not set(required) <= set(d):
        raise ValueError('Empty data or missing required columns: ' + ', '.join(required))
    if d.observation_id.duplicated().any():
        raise ValueError('Duplicate observation_id')
    for k in required:
        if (d[k].str.strip() == '').any():
            raise ValueError('Missing ' + k)
    for k in numeric:
        d[k] = pd.to_numeric(d[k], errors='raise')
        if not np.isfinite(d[k]).all():
            raise ValueError('Non-finite ' + k)
    return d
