"""Only verified Duke T200 XY text and Carterra XY workbook layouts."""
from pathlib import Path
import csv,json,zipfile,re,shutil
from xml.etree import ElementTree as ET
import pandas as pd
import numpy as np
from .workflow import dump,sha


def workbook_rows(path):
    # Both transitional and strict OOXML are present in the pinned public exports.
    with zipfile.ZipFile(path) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for si in ET.fromstring(z.read('xl/sharedStrings.xml')):
                strings.append(''.join(n.text or '' for n in si.iter() if n.tag.rsplit('}',1)[-1]=='t'))
        root=ET.fromstring(z.read('xl/worksheets/sheet1.xml'));rows=[]
        for row in root.iter():
            if row.tag.rsplit('}',1)[-1]!='row':continue
            cells={}
            for c in row:
                ref=c.attrib.get('r','');letters=re.match('[A-Z]+',ref)
                if not letters:continue
                col=0
                for ch in letters[0]:col=col*26+ord(ch)-64
                value=next((v.text for v in c if v.tag.rsplit('}',1)[-1]=='v'),None)
                if c.attrib.get('t')=='s' and value is not None:value=strings[int(value)]
                elif c.attrib.get('t')=='inlineStr':value=''.join(n.text or '' for n in c.iter() if n.tag.rsplit('}',1)[-1]=='t')
                cells[col-1]=value
            rows.append([cells.get(i) for i in range(max(cells,default=-1)+1)])
        return rows


def import_surface(manifest,output):
    src=Path(manifest).resolve();m=json.loads(src.read_text());out=Path(output).resolve()
    allowed={'layout','file','source','time_source','curves'}
    if set(m)!=allowed:raise ValueError('Declare layout, file, source, time_source and curves only')
    if m['layout'] not in ('duke_t200_xy_text','duke_carterra_xy_workbook'):raise ValueError('Unsupported export layout; provide a real sample before adding another parser')
    if not all(isinstance(m[k],str) and m[k].strip() for k in ('source','time_source')):raise ValueError('Sources for export and injection times required')
    path=(src.parent/m['file']).resolve()
    if m['layout']=='duke_t200_xy_text':
        with path.open(encoding='latin1',newline='') as f:rows=list(csv.reader(f,delimiter='\t'))
        header=rows[0];data=rows[1:]
    else:
        rows=workbook_rows(path)
        if len(rows)<4 or rows[2][:2]!=['X','Y']:raise ValueError('Expected verified Carterra X/Y row 3')
        header=rows[1];data=rows[3:]
    if not isinstance(m['curves'],list) or not m['curves']:raise ValueError('Declare every selected column pair')
    output_rows=[];ids=set();pairs=set()
    for c in m['curves']:
        keys={'pair_index','curve_id','sample_id','experiment_id','fit_group_id','concentration_M','association_start_s','dissociation_start_s','response_unit'}
        if set(c)!=keys:raise ValueError('Each trace needs explicit identity, pair index, molarity, response unit and injection times')
        if c['curve_id'] in ids:raise ValueError('Duplicate curve identity')
        ids.add(c['curve_id']);idx=c['pair_index']
        if type(idx)is not int or idx<0 or 2*idx+1>=len(header):raise ValueError('Invalid zero-based column pair')
        if idx in pairs:raise ValueError('A source trace cannot be duplicated as another curve')
        pairs.add(idx)
        if m['layout']=='duke_carterra_xy_workbook' and rows[2][2*idx:2*idx+2]!=['X','Y']:raise ValueError('Selected Carterra columns must be an X/Y pair')
        if m['layout']=='duke_t200_xy_text' and not (str(header[2*idx]).endswith('_X') and str(header[2*idx+1]).endswith('_Y') and str(header[2*idx])[:-2]==str(header[2*idx+1])[:-2]):raise ValueError('Expected matched T200 X/Y headers')
        a,b,conc=[float(c[k]) for k in ('association_start_s','dissociation_start_s','concentration_M')]
        if not np.isfinite([a,b,conc]).all() or b<=a or conc<=0:raise ValueError('Invalid declared times/concentration')
        for row in data:
            xvalue=row[2*idx] if len(row)>2*idx else None
            yvalue=row[2*idx+1] if len(row)>2*idx+1 else None
            if xvalue in (None,'') and yvalue in (None,''):continue
            if xvalue in (None,'') or yvalue in (None,''):raise ValueError('Unpaired missing time/response; no silent observation exclusion')
            t,y=float(row[2*idx]),float(row[2*idx+1])
            if not np.isfinite([t,y]).all():raise ValueError('Nonfinite export values')
            output_rows.append({k:c[k] for k in ('curve_id','sample_id','experiment_id','fit_group_id','response_unit')}|dict(time=t,time_unit='s',response=y,concentration=conc,concentration_unit='M',association_start=a,dissociation_start=b,phase='baseline' if t<a else 'association' if t<b else 'dissociation'))
    if out.exists():raise FileExistsError('Choose a new import directory')
    out.mkdir(parents=True);(out/'original').mkdir();shutil.copy2(path,out/'original'/path.name);pd.DataFrame(output_rows).to_csv(out/'data.csv',index=False)
    dump(out/'import_manifest.json',dict(m,source_copy='original/'+path.name,source_sha256=sha(path),data_sha256=sha(out/'data.csv'),selected_column_labels={c['curve_id']:header[2*c['pair_index']:2*c['pair_index']+2] for c in m['curves']},automatic_phase_detection=False))
    dump(out/'config.template.json',{'analysis_type':'binding_kinetics','input':'data.csv','source':m['source'],'provenance':'import_manifest.json','assay':{'one_to_one_supported':None,'independent_cycles':None,'concentration_known':None,'rationale':''}})
    return out
