"""Immutable variance-component artifacts and render-only report."""
from pathlib import Path
from datetime import datetime,timezone
import json,shutil,platform,html,importlib.metadata
from . import precision,__version__
from .workflow import dump,sha,implementation_hash,verify_run


def analyze_precision(config_path,output,render=True):
    src,out=Path(config_path).resolve(),Path(output).resolve()
    if out.exists():raise FileExistsError(f'Output already exists: {out}')
    out.mkdir(parents=True)
    try:
        c=precision.resolve_config(json.loads(src.read_text()));inp=(src.parent/c['input']).resolve()
        d=precision.load_data(inp,c);shutil.copyfile(inp,out/'input.csv');c['input']='input.csv';dump(out/'config.resolved.json',c)
        result=precision.compute(d,c)
        from .precision_evidence import CALIBRATION
        if c['inference']['sum_interval']=='mls':
            from .interval_evidence import PRECISION as CALIBRATION
        result['validation_evidence']=CALIBRATION
        result['limitations']+=CALIBRATION['limitations']
        dump(out/'results.json',result)
        dump(out/'interpretation_facts.json',{'schema_version':1,'analysis_type':'variance_components',
             'source_sha256':{'results.json':sha(out/'results.json'),'config.resolved.json':sha(out/'config.resolved.json')},
             'primary':{'source':'results.json#/intermediate_precision','value':result['intermediate_precision']},
             'components':{'source':'results.json#/components','value':result['components']},
             'repeatability':{'source':'results.json#/repeatability','value':result['repeatability']},
             **({'positive_sums':result['positive_sums'],'mean_square_inference':result['mean_square_inference']} if 'mean_square_inference' in result else {}),
             'boundary_components':result['boundary_components'],'limitations':result['limitations'],
             'calibration':{'source':'results.json#/validation_evidence','value':CALIBRATION}})
        dump(out/'manifest.json',{'schema_version':1,'analysis_type':'variance_components','package_version':__version__,
             'created_utc':datetime.now(timezone.utc).isoformat(),'input_original_path':str(inp),'source':c['source'],
             'python':platform.python_version(),'platform':platform.platform(),'input_sha256':sha(inp),
             'dependencies':{k:importlib.metadata.version(k) for k in ('numpy','scipy','pandas')},
             'implementation_sha256':implementation_hash(),
             'scientific_artifacts_sha256':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}})
    except Exception as exc:
        dump(out/'failure.json',{'exception':type(exc).__name__,'message':str(exc)});raise
    if render:render_precision(out)
    return out


def render_precision(run,style=None):
    run=Path(run);verify_run(run)
    if style is not None and style not in ('standard','prism_like'):raise ValueError('Unsupported style')
    r=json.loads((run/'results.json').read_text());esc=html.escape
    from . import report_shell as shell
    body=shell.section('details','结果',shell.json_block(r,'完整 results.json',open_=True))
    body+=shell.section('files','可追溯文件',shell.downloads([('interpretation_facts.json','interpretation_facts.json')]))
    body=shell.page('精密度与方差组件',eyebrow='AgenticPrism / variance components',heading='精密度与方差组件',
                    lede='方差组件描述所声明设计的变异；未声明接受限，不作方法验证通过判定。区间是近似的，零边界和校准偏差需一并解释。',
                    nav=[('details','结果'),('files','可追溯文件')],body=body,footer='AgenticPrism · 数值分析与报告渲染分离')
    (run/'report.html').write_text(body);dump(run/'render_manifest.json',{'scientific_artifacts_changed':False,'report_sha256':sha(run/'report.html')})
    return run/'report.html'
