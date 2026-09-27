"""Optional R methods: fixed packaged scripts; JSON data, never generated R code."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
SCRIPTS=Path(__file__).with_name('r_scripts')
import csv
PINS={r['package']:r['version'] for r in csv.DictReader((SCRIPTS/'packages.csv').open())}


def library_path():
    configured=os.environ.get('AGENTIC_PRISM_R_LIB')
    if configured:return Path(configured).expanduser().resolve()
    # Editable collections and installed wheels inside <collection>/.venv.
    for root in (Path(__file__).resolve().parents[2],Path(__import__('sys').prefix).parent):
        if (root/'.r-lib').is_dir():return root/'.r-lib'
    return None


def call(method,payload):
    if method not in ('mmrm','doctor'):raise ValueError('Unsupported R bridge method')
    executable=shutil.which('Rscript')
    if executable is None:raise RuntimeError('This method requires optional R. Install R and run install.py --with-r; Python methods remain available.')
    with tempfile.TemporaryDirectory(prefix='prism-r-') as tmp:
        inp,out=Path(tmp)/'input.json',Path(tmp)/'output.json'
        inp.write_text(json.dumps({'method':method,'pins':PINS,'data':payload},allow_nan=False))
        args=[executable,'--vanilla',str(SCRIPTS/'runner.R'),str(inp),str(out),str(library_path() or '')]
        proc=subprocess.run(args,capture_output=True,text=True,timeout=180)
        if proc.returncode or not out.exists():
            raise RuntimeError('Optional R method unavailable or failed: '+proc.stderr[-1500:]+'; check agentic-prism doctor and install.py --with-r')
        result=json.loads(out.read_text())
        if 'r_environment' in result:
            result['r_environment']['scripts_sha256']={p.name:__import__('hashlib').sha256(p.read_bytes()).hexdigest() for p in sorted(SCRIPTS.iterdir()) if p.is_file()}
        if 'error' in result:
            error=ArithmeticError(result['error']);error.r_environment=result.get('r_environment');raise error
        return result


def availability():
    try:return {'available':True,**call('doctor',{})}
    except (RuntimeError,ArithmeticError,subprocess.TimeoutExpired) as exc:return {'available':False,'reason':str(exc)}
