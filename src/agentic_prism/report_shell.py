"""Shared HTML shell for every report: one stylesheet, one page scaffold, one figure block."""
import base64
import html
import json
import re
from pathlib import Path

CSS = """
:root{--ink:#1d2329;--muted:#5b6670;--faint:#8a939c;--line:#e2e6ea;--bg:#f5f6f8;--card:#fff;--accent:#0b5cad;
--ok:#1f6e43;--ok-bg:#e6f4ec;--warn:#8a5712;--warn-bg:#fff3de;--bad:#9e2a1f;--bad-bg:#fdeceb;color-scheme:light}
*{box-sizing:border-box}html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.7 -apple-system,BlinkMacSystemFont,"Helvetica Neue",Arial,"PingFang SC","Hiragino Sans GB","Microsoft YaHei",sans-serif}
p,li,h1,h2,h3,footer,summary{overflow-wrap:anywhere}
header{background:var(--card);border-top:4px solid var(--accent);border-bottom:1px solid var(--line);padding:34px max(4vw,20px) 22px}
header>*{max-width:1120px;margin-left:auto;margin-right:auto}
.eyebrow{font-size:11.5px;font-weight:600;letter-spacing:.14em;text-transform:uppercase;color:var(--accent)}
header h1{font-size:30px;line-height:1.25;font-weight:650;letter-spacing:-.01em;margin-top:8px;margin-bottom:10px}
header p{color:var(--muted);margin-top:0}
header nav{display:flex;flex-wrap:wrap;gap:8px;margin-top:16px}
header nav a{padding:4px 12px;border:1px solid var(--line);border-radius:999px;font-size:13px;color:var(--ink);background:#fafbfc}
header nav a:hover{border-color:var(--accent);color:var(--accent);text-decoration:none}
main{max-width:1160px;margin:auto;padding:22px 20px 32px}
section{margin:24px 0}
h2{font-size:20px;line-height:1.35;font-weight:650;margin:4px 0 12px}
h3{font-size:17px;line-height:1.35;font-weight:650;margin:0}
.panel,.curve,.dose-curve,.kinetic-group,.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:24px;margin:18px 0;box-shadow:0 1px 2px rgba(16,24,40,.04);min-width:0}
.card>h2:first-child,.panel>h2:first-child{margin-top:0}
.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}.grid>.curve{margin:0;scroll-margin-top:16px}
.head,.curve-head{display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap;margin-bottom:6px}
.head h2{margin:0}
.badge{display:inline-block;font-size:12px;font-weight:600;line-height:1;padding:6px 11px;border-radius:999px;background:var(--ok-bg);color:var(--ok);white-space:nowrap}
.badge.limited,.badge.estimated_with_diagnostics,.badge.warn{background:var(--warn-bg);color:var(--warn)}
.badge.failed,.badge.bad{background:var(--bad-bg);color:var(--bad)}
.metric{margin:14px 0 2px;color:var(--muted)}.metric strong{font-size:28px;font-weight:650;color:var(--ink);font-variant-numeric:tabular-nums}
.notice{border-left:3px solid #d49a3a;background:#fff8ec;padding:9px 14px;margin:10px 0;font-size:13.5px;border-radius:0 6px 6px 0}
.muted,.caption{color:var(--muted);font-size:12.5px;line-height:1.65}
.toolbar{display:flex;flex-wrap:wrap;gap:10px 18px;align-items:center;font-size:13px;color:var(--muted)}
.toolbar label{color:var(--ink);font-weight:600}
figure{margin:16px 0 8px;padding:14px 14px 10px;background:#fff;border:1px solid var(--line);border-radius:8px}
figure img{display:block;max-width:100%;height:auto;margin:0 auto}
figure svg{display:block;width:100%;max-width:820px;height:auto;margin:0 auto}
figcaption{margin-top:8px}
.downloads{display:flex;flex-wrap:wrap;gap:8px;font-size:12.5px}
.downloads a{padding:2px 10px;border:1px solid var(--line);border-radius:6px;color:var(--ink);background:#fafbfc}
.downloads a:hover{border-color:var(--accent);color:var(--accent);text-decoration:none}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}
.table-wrap{overflow-x:auto;border:1px solid var(--line);border-radius:8px;margin:10px 0}
.table-wrap table{margin:0}
table{width:100%;border-collapse:collapse;font-size:13.5px;font-variant-numeric:tabular-nums}
th,td{padding:8px 12px;text-align:left;border-bottom:1px solid var(--line);white-space:nowrap;vertical-align:top}
th{background:#f3f5f7;font-weight:600;color:#39424b;font-size:12.5px}
tr:nth-child(even)>td{background:#fafbfc}tr:last-child>td{border-bottom:0}
td table{font-size:12.5px}td td,td th{white-space:normal}
select{font:inherit;font-size:13px;padding:5px 10px;border:1px solid #c9d0d6;border-radius:6px;background:#fff;color:var(--ink)}
pre{background:#f6f8fa;border:1px solid var(--line);border-radius:8px;padding:12px 14px;white-space:pre-wrap;overflow-wrap:anywhere;font:12px/1.55 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;max-height:560px;overflow:auto}
details{margin:14px 0;border:1px solid var(--line);border-radius:8px;background:#fcfcfd}
details>summary{cursor:pointer;padding:9px 14px;font-weight:600;font-size:14px}
details[open]>summary{border-bottom:1px solid var(--line)}
details>:not(summary){margin-left:14px;margin-right:14px}
details>pre{border:0;background:transparent;padding:0}
ul,ol{padding-left:22px}li{margin:4px 0}
.must li{margin:6px 0}
.results>h2~h2,.card>h2~h2{font-size:16px;margin:24px 0 8px}
footer{color:var(--faint);font-size:12px;padding:16px 0 8px;border-top:1px solid var(--line);margin-top:28px}
[hidden]{display:none!important}
@media(max-width:760px){.grid{grid-template-columns:1fr}main{padding:14px 12px 24px}.panel,.curve,.dose-curve,.kinetic-group,.card{padding:16px}header{padding:24px 16px 16px}header h1{font-size:25px}}
@media print{body{background:#fff}header{border-top:0}header nav,select,.toolbar,.downloads{display:none}
.panel,.curve,.dose-curve,.kinetic-group,.card{box-shadow:none;break-inside:avoid}figure{break-inside:avoid}pre{max-height:none}}
"""

STYLE_LABELS = {"prism_like": "Prism-like（出版风格）", "standard": "标准"}
_SWITCH_SCRIPT = ("<script>(function(){const s=document.getElementById('style');"
                  "function set(v){document.querySelectorAll('[data-theme]').forEach(e=>e.hidden=e.dataset.theme!==v);}"
                  "if(s){s.value=%s;s.addEventListener('change',()=>set(s.value));}set(%s);})();</script>")


def data_uri(path):
    mime = {".svg": "image/svg+xml", ".png": "image/png", ".pdf": "application/pdf", ".csv": "text/csv",
            ".json": "application/json", ".txt": "text/plain"}.get(Path(path).suffix, "application/octet-stream")
    return "data:" + mime + ";base64," + base64.b64encode(Path(path).read_bytes()).decode()


def _display_style(svg_path, scale=1.4):
    """Show a figure at a fixed multiple of its physical size (never wider than the card)."""
    head = Path(svg_path).read_text(errors="ignore")[:600]
    match = re.search(r'width="([0-9.]+)pt"', head)
    if not match:
        return ""
    return f' style="width:min(100%,{float(match.group(1)) * 4 / 3 * scale:.0f}px)"'


def style_select(themes, style, note="切换风格只改变图形外观；数据、坐标范围与拟合结果不变。"):
    options = "".join(f'<option value="{t}"{" selected" if t == style else ""}>{STYLE_LABELS.get(t, t)}</option>' for t in themes)
    return (f'<div class="toolbar"><label for="style">图形风格</label><select id="style">{options}</select>'
            f'<span>{html.escape(note)}</span></div>')


def switch_script(style):
    value = json.dumps(style)
    return _SWITCH_SCRIPT % (value, value)


def theme_figures(folder, key, style, alt, themes, formats=("svg", "pdf", "png")):
    """One <figure> per theme; only the selected one is visible. Images and downloads are embedded."""
    blocks = []
    for theme in themes:
        links = "".join(f'<a download="{key}__{theme}.{ext}" href="{data_uri(folder / f"{key}__{theme}.{ext}")}">{ext.upper()}</a>'
                        for ext in formats)
        blocks.append(f'<figure class="theme-figure" data-theme="{theme}"{"" if theme == style else " hidden"}>'
                      f'<img alt="{html.escape(alt)} · {theme}"{_display_style(folder / f"{key}__{theme}.svg")} '
                      f'src="{data_uri(folder / f"{key}__{theme}.svg")}">'
                      f'<figcaption class="downloads">{links}</figcaption></figure>')
    return "".join(blocks)


def single_figure(folder, name, alt, formats=("svg", "pdf", "png")):
    """A figure rendered in the requested style only, with embedded downloads."""
    links = "".join(f'<a download="{name}.{ext}" href="{data_uri(folder / f"{name}.{ext}")}">{ext.upper()}</a>'
                    for ext in formats if (folder / f"{name}.{ext}").exists())
    return (f'<figure><img alt="{html.escape(alt)}"{_display_style(folder / f"{name}.svg")} src="{data_uri(folder / f"{name}.svg")}">'
            f'<figcaption class="downloads">{links}</figcaption></figure>')


def table_wrap(table_html):
    return f'<div class="table-wrap">{table_html}</div>'


def page(title, *, eyebrow, heading, lede="", nav=(), body, footer="", style=None, themes=None):
    """Complete standalone HTML document. ``body`` is trusted HTML; text arguments are escaped here."""
    nav_html = ("<nav>" + "".join(f'<a href="#{html.escape(anchor)}">{html.escape(label)}</a>' for anchor, label in nav)
                + "</nav>") if nav else ""
    toolbar = f'<div class="panel">{style_select(themes, style)}</div>' if style and themes and len(themes) > 1 else ""
    script = switch_script(style) if toolbar else ""
    return ('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1"><link rel="icon" href="data:,">'
            f'<title>{html.escape(title)}</title><style>{CSS}</style></head><body>'
            f'<header><div class="eyebrow">{html.escape(eyebrow)}</div><h1>{html.escape(heading)}</h1>'
            + (f'<p>{html.escape(lede)}</p>' if lede else "") + nav_html + '</header>'
            f'<main>{toolbar}{body}<footer>{footer}</footer></main>{script}</body></html>')


def section(anchor, heading, content, card=True):
    cls = ' class="card"' if card else ""
    ident = f' id="{html.escape(anchor)}"' if anchor else ""
    return f'<section{ident}{cls}><h2>{html.escape(heading)}</h2>{content}</section>'


def downloads(pairs):
    """``pairs`` of (label, href) → a row of download chips; href may be a data URI or a relative path."""
    return '<div class="downloads">' + "".join(f'<a download="{html.escape(label)}" href="{href}">{html.escape(label)}</a>'
                                               for label, href in pairs) + '</div>'


def fmt(value, digits=4):
    """Display text for one saved value. Floats are shown to ``digits`` significant figures;
    the downloadable CSV/JSON artifacts keep full precision."""
    if isinstance(value, bool) or value is None:
        return "—" if value is None else str(value)
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def value_html(v):
    """Nested dict/list rendering as tables and lists, with rounded floats for display."""
    if isinstance(v, dict):
        return ('<div class="table-wrap"><table>' + "".join(f"<tr><th>{html.escape(str(k))}</th><td>{value_html(x)}</td></tr>"
                                                            for k, x in v.items()) + "</table></div>")
    if isinstance(v, list):
        if v and all(not isinstance(x, (dict, list)) for x in v):
            return html.escape(", ".join(fmt(x) for x in v))
        return "<ol>" + "".join(f"<li>{value_html(x)}</li>" for x in v) + "</ol>"
    return html.escape(fmt(v))


def json_block(obj, summary=None, open_=False):
    text = html.escape(json.dumps(obj, ensure_ascii=False, indent=2))
    if summary is None:
        return f"<pre>{text}</pre>"
    return f'<details{" open" if open_ else ""}><summary>{html.escape(summary)}</summary><pre>{text}</pre></details>'
