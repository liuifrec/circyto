"""Prepare/verify a local Chrome test page when CDP sockets are unavailable.

python tests/explorer_headless.py prepare REPORT_DIR TEST.html
google-chrome --headless --no-sandbox --disable-gpu --virtual-time-budget=30000 \
  --run-all-compositor-stages-before-draw --dump-dom file:///absolute/TEST.html > /tmp/explorer-dom.html
python tests/explorer_headless.py verify /tmp/explorer-dom.html RESULTS.json
"""
from pathlib import Path
import argparse
import base64
import hashlib
from html.parser import HTMLParser
import json


PROBE = r"""
(async()=>{
const result={status:'running',checks:[],resources:[],exports:{}},out=document.createElement('pre');out.id='browser-test-result';document.body.append(out);
const started=performance.now(),E=id=>document.getElementById('ex-'+id),D=JSON.parse(E('data').textContent),blobs=new Map();
const originalCreate=URL.createObjectURL;URL.createObjectURL=blob=>{const url=originalCreate(blob);blobs.set(url,blob);return url;};
HTMLAnchorElement.prototype.click=function(){if(this.download)result.exports[this.download]=this.href;};
const errors=[];window.addEventListener('error',e=>errors.push(e.message));
const check=(ok,label)=>{if(!ok)throw Error(label);result.checks.push(label);};
const change=(id,v)=>{E(id).value=v;E(id).dispatchEvent(new Event('change'));};
const search=(id,v)=>{E(id+'-query').value=v;E(id+'-query').dispatchEvent(new Event('input'));};
try{
check(document.getElementById('explorer').dataset.ready==='true','Interface initializes under strict CSP');
check(document.querySelectorAll('#ex-left [data-cell]').length===D.embedding.plotted_cells,'Every finite coordinate is plotted');
if(D.embedding.plotted_cells){
 const left=[...document.querySelectorAll('#ex-left [data-cell]')],right=[...document.querySelectorAll('#ex-right [data-cell]')];
 check(left.every(n=>{const r=right.find(o=>o.dataset.cell===n.dataset.cell);if(!r)return false;const a=n.getBBox(),b=r.getBBox();return Math.abs(a.x+a.width/2-b.x-b.width/2)<.6&&Math.abs(a.y+a.height/2-b.y-b.height/2)<.6;}),'Gene and circRNA points share exact RNA positions');
 const cell=left[0].dataset.cell;left[0].dispatchEvent(new Event('pointerenter'));
 check(E('hover').textContent.includes(D.cells[cell]),'Linked point inspection exposes exact biological identity');
}
const c=D.circs[0],g=D.genes[0];
search('circ',c.id);change('circ','0');check(E('candidate-title').textContent===c.id,'Exact BSJ search selects the candidate');
search('circ',c.aliases.find(a=>a!==c.id)||c.id);check(E('circ').options.length>0,'Catalogue or host alias is searchable');
search('gene',g.id);change('gene','0');check(E('left-unit').textContent===g.unit||!D.embedding.plotted_cells,'Stable gene ID search and explicit expression unit');
if(D.genes.some(g=>g.metadata.gene_name==='DUPLICATE')){search('gene','DUPLICATE');check(E('gene').options.length===2,'Duplicate symbols retain separate stable IDs');}
search('gene','');change('gene','0');search('circ','');change('circ','0');
search('gene','__NOT_EMBEDDED__');check(E('gene').disabled&&E('gene').textContent.includes('No embedded match'),'Search makes bounded feature coverage explicit');search('gene','');change('gene','0');
const group=D.metadata.find(m=>m.kind==='category'&&/condition/.test(m.key));if(group)change('group',group.key);
change('mode','detection');change('threshold','2');
E('export-groups').click();E('export-cells').click();
document.querySelector('[data-export="left"][data-format="svg"]').click();document.querySelector('[data-export="right"][data-format="png"]').click();
for(let i=0;i<50&&!result.exports['circyto_right.png']&&D.embedding.plotted_cells;i++)await new Promise(r=>setTimeout(r,50));
const cells=await blobs.get(result.exports['circyto_selected_cells.tsv']).text(),summary=await blobs.get(result.exports['circyto_detection_summary.tsv']).text();
const rows=summary.trim().split('\n').slice(1).map(r=>r.split('\t'));
const complete=c.states.filter(s=>s.startsWith('completed_')).length,detected=c.values.filter((v,i)=>c.states[i].startsWith('completed_')&&v>=2).length;
check(rows.reduce((n,r)=>n+Number(r[3]),0)===complete&&rows.reduce((n,r)=>n+Number(r[4]),0)===detected,'TSV frequency denominator excludes failed/unknown and retains completed zeros');
const cellRows=cells.trim().split('\n').slice(1).map(r=>r.split('\t'));
check(cellRows.length===D.cells.length,'TSV includes selected cells without UMAP coordinates');
check(cellRows.every((r,i)=>r[0]===D.cells[i]&&r[11]===c.states[i]&&r[9]===(c.values[i]===null?'':String(c.values[i]))),'TSV cell identities, native values and processing states are exact');
if(group){const desired=String(group.values.find(v=>v!==null));change('filter',JSON.stringify(desired));E('export-cells').click();const filtered=(await blobs.get(result.exports['circyto_selected_cells.tsv']).text()).trim().split('\n').length-1;check(filtered===group.values.filter(v=>String(v)===desired).length,'Group filter and exported cell scope agree');change('filter','*');}
if(D.embedding.plotted_cells){
 const svg=await blobs.get(result.exports['circyto_left.svg']).text();check(svg.includes(g.unit)&&svg.includes(g.id),'SVG export includes selected feature and measurement units');
 const png=blobs.get(result.exports['circyto_right.png']);check(png&&png.type==='image/png'&&png.size>1000,'PNG export produces a nonempty image');
 const numeric=D.metadata.find(m=>m.kind==='numeric');if(numeric){change('color',numeric.key);check(E('left-unit').textContent.includes(numeric.key),'QC selector uses stored values');}
 const cluster=D.metadata.find(m=>/cluster|leiden/.test(m.key));if(cluster){change('color',cluster.key);check(E('left-unit').textContent==='Existing cell annotation','Existing RNA annotations selectable');}
}
change('color','gene');change('threshold','1');change('mode','support');
check(errors.length===0,'No JavaScript runtime exceptions');
check(document.documentElement.scrollWidth<=window.innerWidth,'Responsive layout has no horizontal page overflow');
const resources=performance.getEntriesByType('resource').map(e=>e.name);check(resources.filter(u=>/^https?:/.test(u)).length===0,'No HTTP or CDN requests');
result.resources=resources;result.status='pass';result.browser=navigator.userAgent;result.viewport=[innerWidth,innerHeight];result.script_ms=performance.now()-started;result.js_heap_bytes=performance.memory?.usedJSHeapSize;
result.downloads=Object.fromEntries([...Object.entries(result.exports)].map(([k,v])=>[k,{bytes:blobs.get(v).size,type:blobs.get(v).type}]));delete result.exports;
}catch(error){result.status='fail';result.error=String(error);result.stack=error.stack;}
out.textContent=JSON.stringify(result);document.documentElement.dataset.testStatus=result.status;
document.body.replaceChildren(out);document.head.replaceChildren();
})();
"""


class ResultParser(HTMLParser):
    collecting = False
    result = ''
    def handle_starttag(self, tag, attrs):
        if tag == 'pre' and dict(attrs).get('id') == 'browser-test-result':
            self.collecting = True
    def handle_endtag(self, tag):
        if tag == 'pre':
            self.collecting = False
    def handle_data(self, data):
        if self.collecting:
            self.result += data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'verify'])
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.command == 'prepare':
        document = (args.source/'report.html').read_text()
        digest = base64.b64encode(hashlib.sha256(PROBE.encode()).digest()).decode()
        document = document.replace('script-src ', f"script-src 'sha256-{digest}' ", 1)
        document = document.replace('</body>', f'<script>{PROBE}</script></body>')
        args.output.write_text(document)
    else:
        p = ResultParser(); p.feed(args.source.read_text())
        result = json.loads(p.result)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result, indent=2))
        if result['status'] != 'pass':
            raise SystemExit(1)


if __name__ == '__main__':
    main()
