"""Self-contained explorer assets. No runtime JavaScript or plotting dependency."""
from __future__ import annotations

import base64
import hashlib
import html
import json


STYLE = r"""
#explorer{--ex-teal:#087f83;--ex-purple:#7056ba;padding:0;overflow:hidden}
#explorer .ex-intro{padding:28px 30px 22px;background:linear-gradient(110deg,#eaf5f3,#f3f1fa);border-bottom:1px solid #dce4e9}
#explorer h1{font-size:32px;margin:6px 0 12px}#explorer h2{font-size:20px}#explorer h3{font-size:16px}
.ex-kicker{font-size:11px;letter-spacing:1.8px;font-weight:800;color:#087f83}.ex-sub{font-size:14px;color:#506375;max-width:940px}
.ex-badges{display:flex;gap:8px;flex-wrap:wrap;margin-top:18px}.ex-badge{font-size:12px;background:white;border:1px solid #d5e4e5;border-radius:20px;padding:4px 12px}
.ex-body{padding:24px 30px}.ex-toolbar{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px;margin-bottom:16px}
#explorer label{display:block;font-size:12px;font-weight:650;color:#40576b}#explorer input,#explorer select{display:block;width:100%;font:13px system-ui;color:#20384b;border:1px solid #c7d5de;background:#fff;border-radius:6px;padding:9px 10px;margin:5px 0;min-height:38px}
#explorer input:focus,#explorer select:focus,#explorer button:focus-visible{outline:2px solid #087f83;outline-offset:2px}
#explorer button{cursor:pointer;background:white;color:#234558;border:1px solid #bdcfd9;border-radius:6px;font:600 12px system-ui;padding:9px 12px}#explorer button:hover{background:#edf6f6}#explorer button:disabled{opacity:.5;cursor:default}
.ex-plot-grid{display:grid;grid-template-columns:1fr 1fr;gap:20px}.ex-plot{border:1px solid #d7e1e8;border-radius:9px;min-width:0;overflow:hidden}.ex-plot-head{padding:16px 18px 10px;border-top:3px solid #087f83}.ex-plot:nth-child(2) .ex-plot-head{border-top-color:#7056ba}
.ex-plot svg{display:block;width:100%;height:auto;touch-action:manipulation;background:#fff}.ex-plot-tools{display:flex;gap:7px;padding:0 18px 14px;flex-wrap:wrap}.ex-unit{color:#526677;font-size:12px;margin:3px 0}.ex-status{font-size:12px;color:#506375;background:#f3f7f9;border-radius:6px;padding:9px 13px;margin:12px 0;min-height:38px;overflow-wrap:anywhere}
.ex-detail{margin-top:22px;padding-top:24px;border-top:1px solid #dce4e9}.ex-detail-grid{display:grid;grid-template-columns:1fr 1.2fr;gap:25px}.ex-counts{display:flex;gap:25px;margin:14px 0}.ex-counts strong{font-size:26px;display:block;color:#294856}.ex-counts span{font-size:11px;color:#526677;display:block}
.ex-detail-grid>div{min-width:0}.ex-meta{display:grid;grid-template-columns:130px minmax(0,1fr);gap:6px 16px;font-size:12px;margin:12px 0}.ex-meta dt{color:#526677;overflow-wrap:anywhere}.ex-meta dd{margin:0;overflow-wrap:anywhere}.ex-small{font-size:12px;color:#526677;overflow-wrap:anywhere}.ex-note{background:#f4f7f9;border-left:3px solid #92aabb;padding:10px 14px;font-size:12px;margin:12px 0}.ex-actions{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:12px 0}
#ex-candidate-title{overflow-wrap:anywhere}#ex-cell-table td,#ex-group-table td{font-variant-numeric:tabular-nums}.ex-hidden{display:none!important}#explorer a{color:#086b77}.ex-no-data{padding:26px;color:#526677;background:#f5f7f9;border-radius:8px}
#ex-group-table th{white-space:normal;letter-spacing:0;font-size:10px}#ex-group-table th,#ex-group-table td{padding:8px 7px}
@media(max-width:900px){.ex-toolbar{grid-template-columns:1fr 1fr}.ex-detail-grid{grid-template-columns:1fr}.ex-body{padding:20px}.ex-plot-grid{gap:12px}}
@media(max-width:620px){.ex-plot-grid,.ex-toolbar{grid-template-columns:1fr}#explorer .ex-intro{padding:22px 20px}#explorer h1{font-size:26px}.ex-counts{gap:16px}.ex-body{padding:16px}}
@media print{.ex-toolbar,.ex-plot-tools,.ex-actions{display:none}.ex-plot-grid{grid-template-columns:1fr 1fr}#explorer{overflow:visible}}
"""

FRAGMENT = r"""
<section id="explorer" aria-label="Single-cell circRNA explorer">
<div class="ex-intro"><div class="ex-kicker">CIRCYTO / SINGLE-CELL EXPLORER</div><h1>RNA context. Junction evidence.</h1>
<p id="ex-dataset"></p><p class="ex-sub">Explore gene expression and circRNA candidates on the same saved RNA embedding. Each point is a biological cell. Native measurements and processing states stay distinct.</p>
<div class="ex-badges" id="ex-badges"></div></div>
<div class="ex-body"><noscript>This optional explorer requires JavaScript enabled in your browser. No network connection is needed. The exported explorer.json retains all selected values.</noscript>
<p class="ex-small" id="ex-scope"></p><div id="ex-warnings"></div>
<div class="ex-toolbar">
<label>Search embedded genes<input id="ex-gene-query" type="search" placeholder="Symbol or stable gene ID" autocomplete="off"><select id="ex-gene" aria-label="Gene"></select></label>
<label>Search embedded circRNAs<input id="ex-circ-query" type="search" placeholder="BSJ, catalogue alias or host gene" autocomplete="off"><select id="ex-circ" aria-label="circRNA candidate"></select></label>
<label>RNA plot colour<select id="ex-color" aria-label="RNA plot colour"></select><span class="ex-small">Expression, annotations, exposure or QC</span></label>
<label>Junction view<select id="ex-mode"><option value="support">Native junction support</option><option value="detection">Detection at threshold</option></select><span class="ex-small">Measured zeros remain visible</span></label>
<label>Group detection frequencies by<select id="ex-group"></select></label>
<label>Show cells in group<select id="ex-filter"></select></label>
<label>Detection threshold (native support ≥)<input id="ex-threshold" type="number" value="1" min="1" max="9007199254740991" step="1"></label>
<label>Inspect biological cell<input id="ex-cell-query" type="search" placeholder="Exact cell ID; Enter to inspect" autocomplete="off"></label>
</div>
<div id="ex-no-umap" class="ex-no-data ex-hidden">RNA UMAP unavailable. Candidate details and exact-cell tables remain available; no embedding was computed.</div>
<div class="ex-plot-grid" id="ex-plots">
<article class="ex-plot" id="ex-rna-panel"><div class="ex-plot-head"><h3 id="ex-left-heading">RNA expression</h3><p class="ex-unit" id="ex-left-unit"></p></div><div id="ex-left"></div><div class="ex-plot-tools"><button data-export="left" data-format="svg">Export SVG</button><button data-export="left" data-format="png">Export PNG</button></div></article>
<article class="ex-plot" id="ex-circ-panel"><div class="ex-plot-head"><h3>circRNA junction evidence</h3><p class="ex-unit" id="ex-right-unit"></p></div><div id="ex-right"></div><div class="ex-plot-tools"><button data-export="right" data-format="svg">Export SVG</button><button data-export="right" data-format="png">Export PNG</button></div></article>
</div>
<div class="ex-status" id="ex-hover" role="status" aria-live="polite">Hover or click a point to inspect its exact cell identity in both views.</div>
<div class="ex-detail"><div class="ex-kicker">CANDIDATE DETAIL</div><h2 id="ex-candidate-title"></h2><p class="ex-small" id="ex-candidate-aliases"></p>
<div class="ex-detail-grid"><div><div class="ex-counts" id="ex-counts"></div><dl class="ex-meta" id="ex-meta"></dl><div id="ex-evidence"></div></div>
<div><h3>Detection frequencies</h3><p class="ex-small" id="ex-denominator"></p><div class="table-scroll" id="ex-group-table"></div><div class="ex-actions"><button id="ex-export-groups">Export summary TSV</button><button id="ex-export-cells">Export selected cells TSV</button></div></div></div>
<p class="ex-note">Completed-zero means measured native support is zero. Missing, failed, unprocessed and unknown observations are excluded from detection denominators. Completion does not establish locus-level callability; junction support and catalogue matches do not confirm circularity or biological function. Group frequencies are descriptive.</p></div>
<details id="ex-cells"><summary>Underlying selected measurements and processing states</summary><p class="ex-small" id="ex-table-note"></p><div class="table-scroll" id="ex-cell-table"></div></details>
<details><summary>Selection, identity joins and provenance</summary><p class="ex-small" id="ex-policy"></p><div class="table-scroll" id="ex-sources"></div></details>
</div></section>
"""

SCRIPT = r"""
(() => {
'use strict';
const D=JSON.parse(document.getElementById('ex-data').textContent), E=id=>document.getElementById('ex-'+id);
const NS='http://www.w3.org/2000/svg', palette=['#087f83','#9566b8','#dc8b38','#4d78af','#c65570','#6e8f41','#7e6387','#92714d'];
const stateColor={missing:'#8c9aa6',failed:'#ba414c',unprocessed:'#ae7b1b',unknown:'#795996'};
let gene=0,circ=0,selected=-1,active=[],groups=[],cellRows=[],threshold=1;
const fmt=v=>v===null?'Not available':typeof v==='number'?(Number.isInteger(v)?String(v):v.toPrecision(5)):String(v);
const label=f=>f?((f.metadata.gene_name||f.metadata.gene_symbol||f.metadata.symbol||f.metadata.host_gene||'')+' · '+f.id).replace(/^ · /,''):'';
const node=(tag,text,attrs={})=>{const n=document.createElement(tag);if(text!==null)n.textContent=text;for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);return n;};
const svgNode=(tag,attrs={},text=null)=>{const n=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);if(text!==null)n.textContent=text;return n;};
const option=(el,text,value)=>el.append(node('option',text,{value}));
const meta=key=>D.metadata.find(m=>m.key===key);
function table(el,heads,rows){el.replaceChildren();const t=node('table',null),h=node('tr',null);for(const x of heads)h.append(node('th',x,{scope:'col'}));const th=node('thead',null);th.append(h);t.append(th);const body=node('tbody',null);for(const row of rows){const tr=node('tr',null);for(const v of row)tr.append(node('td',fmt(v)));body.append(tr);}t.append(body);el.append(t);}
function search(kind){const list=kind==='gene'?D.genes:D.circs,q=E(kind+'-query').value.toLocaleLowerCase(),sel=E(kind),current=kind==='gene'?gene:circ;sel.replaceChildren();list.forEach((f,i)=>{if(f.aliases.some(a=>a.toLocaleLowerCase().includes(q)))option(sel,label(f),i);});if([...sel.options].some(o=>Number(o.value)===current))sel.value=String(current);else if(sel.options.length){if(kind==='gene')gene=Number(sel.value);else circ=Number(sel.value);}else option(sel,'No embedded match — regenerate with --'+kind,'');sel.disabled=!list.length||!sel.options[0].value&&sel.options[0].text.startsWith('No embedded');render();}
function groupName(i){const m=meta(E('group').value);return m?(m.values[i]===null?'(missing annotation)':String(m.values[i])):'All cells';}
function refreshFilter(){E('filter').replaceChildren();option(E('filter'),'All groups','*');const names=[...new Set(D.cells.map((_,i)=>groupName(i)))].sort();for(const n of names)option(E('filter'),n,JSON.stringify(n));}
function isCompleted(f,i){return f&&['completed_zero','completed_detected'].includes(f.states[i])&&f.values[i]!==null;}
function summarize(f){const by=new Map();for(const i of active){const name=groupName(i);if(!by.has(name))by.set(name,{group:name,cells:0,completed:0,detected:0,zero:0,missing:0,failed:0,unprocessed:0,unknown:0,support:0});const g=by.get(name);g.cells++;if(isCompleted(f,i)){g.completed++;g.detected+=Number(f.values[i]>=threshold);g.zero+=Number(f.values[i]===0);g.support+=f.values[i];}else{const state=f?f.states[i]:'missing';g[state in stateColor?state:'unknown']++;}}return [...by.values()].sort((a,b)=>a.group.localeCompare(b.group));}
function inspect(i){selected=i;const g=D.genes[gene],c=D.circs[circ];E('hover').textContent=D.cells[i]+' | RNA: '+(g?fmt(g.values[i])+' '+g.unit+' ['+g.states[i]+']':'unavailable')+' | Junction: '+(c?fmt(c.values[i])+' '+c.unit+' ['+c.states[i]+']':'unavailable')+(D.coordinates[i]?'':' | No RNA coordinates');document.querySelectorAll('#ex-plots [data-cell]').forEach(n=>{n.setAttribute('stroke-width',Number(n.dataset.cell)===i?'2.5':'0.7');n.setAttribute('stroke',Number(n.dataset.cell)===i?'#142e3f':n.dataset.stroke);});}
const validXY=D.coordinates.filter(Boolean);
const bounds=validXY.length?[0,1].map(k=>{let low=Infinity,high=-Infinity;for(const p of validXY){low=Math.min(low,p[k]);high=Math.max(high,p[k]);}const pad=(high-low||1)*.1;return [low-pad,high+pad];}):[[0,1],[0,1]];
const scale=(v,k)=>k===0?60+(v-bounds[0][0])/(bounds[0][1]-bounds[0][0])*420:330-(v-bounds[1][0])/(bounds[1][1]-bounds[1][0])*260;
function gradient(t,right){const a=right?[235,229,247]:[225,242,239],b=right?[101,65,166]:[0,114,118];return '#'+a.map((v,i)=>Math.round(v+(b[i]-v)*t).toString(16).padStart(2,'0')).join('');}
function draw(side){const right=side==='right',f=right?D.circs[circ]:D.genes[gene],m=!right&&E('color').value!=='gene'?meta(E('color').value):null,category=m?.kind==='category',detection=right&&E('mode').value==='detection';
 const values=m?m.values:f?f.values:D.cells.map(()=>null),states=m?D.cells.map(()=> 'available'):f?f.states:D.cells.map(()=> 'missing');
 const eligible=i=>values[i]!==null&&!Object.hasOwn(stateColor,states[i]);
 const numeric=active.filter(eligible).map(i=>values[i]).filter(v=>typeof v==='number');let lo=0,hi=0;for(const v of numeric){lo=Math.min(lo,v);hi=Math.max(hi,v);}const cats=category?[...new Set(active.map(i=>values[i]).filter(v=>v!==null).map(String))].sort():[];
 const heading=m?m.key:f?label(f):'No feature embedded',unit=m?(category?'Existing cell annotation':'Stored observation value · '+m.key):f?f.unit:'Not available';
 E(right?'right-unit':'left-unit').textContent=unit; if(!right)E('left-heading').textContent=m?'RNA annotation / QC':'RNA gene expression';
 const s=svgNode('svg',{viewBox:'0 0 540 480',role:'img','aria-label':heading+'; '+unit});s.append(svgNode('title',{},heading+' | '+unit+' | '+D.embedding.key));s.append(svgNode('rect',{width:540,height:480,fill:'white'}));
 const text=(x,y,t,size=11,fill='#506375')=>s.append(svgNode('text',{x,y,'font-family':'system-ui,sans-serif','font-size':size,fill},t));
 text(24,24,heading.length>68?heading.slice(0,65)+'…':heading,12,'#18354a');text(24,43,(unit.length>84?unit.slice(0,81)+'…':unit),10);
 s.append(svgNode('path',{d:'M60 65 V340 H490',stroke:'#becbd4',fill:'none'}));text(240,364,'RNA UMAP 1');s.append(svgNode('text',{x:14,y:204,transform:'rotate(-90 14 204)','text-anchor':'middle','font-family':'system-ui,sans-serif','font-size':11,fill:'#506375'},'RNA UMAP 2'));
 for(let k=0;k<=2;k++){const vx=bounds[0][0]+(bounds[0][1]-bounds[0][0])*k/2,vy=bounds[1][0]+(bounds[1][1]-bounds[1][0])*k/2;text(scale(vx,0)-10,352,vx.toFixed(1),9);s.append(svgNode('text',{x:53,y:scale(vy,1)+3,'text-anchor':'end','font-family':'system-ui,sans-serif','font-size':9,fill:'#506375'},vy.toFixed(1)));}
 const order=active.filter(i=>D.coordinates[i]).sort((a,b)=>Number(eligible(a))-Number(eligible(b))||(Number(values[a])||0)-(Number(values[b])||0));
 for(const i of order){const x=scale(D.coordinates[i][0],0),y=scale(D.coordinates[i][1],1),state=values[i]===null&&eligible(i)?'missing':states[i],ok=eligible(i);let fill=ok?(category?palette[cats.indexOf(String(values[i]))%palette.length]:detection?(values[i]>=threshold?'#7056ba':'#dce2e8'):gradient((values[i]-lo)/(hi-lo||1),right)):'white';if(values[i]===null)fill='white';const stroke=stateColor[state]|| (values[i]===null?stateColor.missing:'#fff');
 const p=(state==='failed')?svgNode('path',{d:`M${x-3} ${y-3}l6 6m0 -6l-6 6`,fill:'none'}):(state==='unprocessed')?svgNode('path',{d:`M${x} ${y-4}l4 7h-8z`,fill}):(state==='unknown')?svgNode('path',{d:`M${x} ${y-4}l4 4l-4 4l-4 -4z`,fill}):svgNode('circle',{cx:x,cy:y,r:3.5,fill});
 p.dataset.cell=i;p.dataset.stroke=stroke;p.setAttribute('stroke',stroke);p.setAttribute('stroke-width',stateColor[state]?1.4:.7);p.append(svgNode('title',{},D.cells[i]+' | '+fmt(values[i])+' | '+states[i]));p.addEventListener('pointerenter',()=>inspect(i));p.addEventListener('click',()=>{inspect(i);E('cell-query').value=D.cells[i];});s.append(p);}
 if(category){const shown=cats.slice(0,12);shown.forEach((c,i)=>{const x=24+(i%3)*170,y=389+Math.floor(i/3)*17;s.append(svgNode('circle',{cx:x,cy:y-3,r:4,fill:palette[i%palette.length]}));text(x+9,y,c.length>22?c.slice(0,20)+'…':c,10);});if(cats.length>12)text(24,463,`${cats.length} categories; first 12 shown. Colours repeat after 8.`);}
 else if(detection){text(24,390,`● Detected ≥ ${threshold}    ○ Below threshold (completed)`,11,'#7056ba');}
 else{for(let i=0;i<100;i++)s.append(svgNode('rect',{x:24+i*2,y:378,width:2,height:9,fill:gradient(i/99,right)}));text(24,405,fmt(lo));text(198,405,fmt(hi));text(245,388,'Linear colour scale',10);}
 if(!category){text(24,433,'○ Missing    × Failed    △ Unprocessed    ◇ Unknown',10);text(24,451,`${order.length} plotted / ${active.length} selected cells · ${D.embedding.key}`,10);text(24,468,'Zero ≠ biological absence. States retained in exported tables.',9);}
 E(side).replaceChildren(s);
}
function detail(){const c=D.circs[circ];E('candidate-title').textContent=c?c.id:'No circRNA feature embedded';E('candidate-aliases').textContent=c?'Search aliases: '+c.aliases.join(' · '):'RNA exploration remains available.';groups=summarize(c);let completed=0,detected=0,zero=0;for(const g of groups){completed+=g.completed;detected+=g.detected;zero+=g.zero;}
 E('counts').replaceChildren();for(const [value,title] of [[detected,'supporting cells ≥ '+threshold],[completed,'completed measured cells'],[completed?(100*detected/completed).toFixed(1)+'%':'N/A','detection frequency']]){const n=node('div',null);n.append(node('strong',String(value)),node('span',title));E('counts').append(n);}
 E('meta').replaceChildren();if(c)for(const [k,v] of Object.entries(c.metadata)){E('meta').append(node('dt',k.replaceAll('_',' ')),node('dd',String(v)));}
 E('evidence').replaceChildren(node('a','Inspect underlying cell evidence',{href:'#ex-cells'}));E('evidence').firstChild.addEventListener('click',()=>E('cells').open=true);
 if(c?.evidence_links.length){E('evidence').append(node('p','Source evidence links (open externally only on request)',{class:'ex-small'}));for(const l of c.evidence_links){const p=node('p',null);p.append(node('a',l.label+': '+l.url,{href:l.url,target:'_blank',rel:'noopener noreferrer',referrerpolicy:'no-referrer'}));E('evidence').append(p);}}
 E('denominator').textContent=`Support ≥ ${threshold}. Denominator: completed cells with a finite value for this candidate, including ${zero} measured zeros. All selected cells count, including those without coordinates. Native unit: ${c?c.unit:'unavailable'}.`;
 table(E('group-table'),['Group','Detected / completed','Frequency','Missing','Failed','Unprocessed','Unknown'],groups.map(g=>[g.group,g.detected+' / '+g.completed,g.completed?(100*g.detected/g.completed).toFixed(1)+'%':null,g.missing,g.failed,g.unprocessed,g.unknown]));
 E('export-groups').disabled=!c;
}
function measurements(){const g=D.genes[gene],c=D.circs[circ],m=meta(E('color').value);cellRows=active.map(i=>[D.cells[i],D.coordinates[i]?.[0]??null,D.coordinates[i]?.[1]??null,groupName(i),g?.id??null,g?.values[i]??null,g?.unit??null,g?.states[i]??'missing',c?.id??null,c?.values[i]??null,c?.unit??null,c?.states[i]??'missing',m?.key??g?.id??null,m?m.values[i]:g?.values[i]??null,m?(m.kind==='category'?'annotation':'stored observation value'):g?.unit??null]);table(E('cell-table'),['Cell','UMAP 1','UMAP 2','Group','Gene','Expression','Expression unit','RNA state','circRNA','Native support','Support unit','circ state','RNA plot field','RNA plot value','RNA plot unit'],cellRows.slice(0,100));E('table-note').textContent=`Showing ${Math.min(100,active.length)} of ${active.length} selected cells. TSV exports all selected rows, including supplied values from failed/unknown cells; those values are excluded from detection rates.`;}
function render(){const v=Number(E('threshold').value);if(!Number.isSafeInteger(v)||v<1){E('threshold').setCustomValidity('Use a positive safe integer');E('threshold').reportValidity();return;}E('threshold').setCustomValidity('');threshold=v;active=D.cells.map((_,i)=>i).filter(i=>E('filter').value==='*'||JSON.stringify(groupName(i))===E('filter').value);if(validXY.length){draw('left');draw('right');}detail();measurements();if(selected>=0)inspect(selected);E('explorer-ready')?.remove();document.getElementById('explorer').dataset.ready='true';}
function save(blob,name){const a=node('a',null,{href:URL.createObjectURL(blob),download:name});document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(a.href),1000);}
function tsv(rows){return rows.map(r=>r.map(v=>{let s=v===null?'':String(v);if(/^[=+@-]/.test(s)&&typeof v!=='number')s="'"+s;return s.replaceAll('\t',' ').replaceAll('\r',' ').replaceAll('\n',' ');}).join('\t')).join('\n')+'\n';}
E('export-cells').addEventListener('click',()=>save(new Blob([tsv([['cell_id','umap_1','umap_2','group','gene_id','expression','expression_unit','rna_state','circ_id','native_support','support_unit','circ_state','rna_plot_field','rna_plot_value','rna_plot_unit'],...cellRows])],{type:'text/tab-separated-values'}),'circyto_selected_cells.tsv'));
E('export-groups').addEventListener('click',()=>{const c=D.circs[circ];save(new Blob([tsv([['circ_id','group','selected_cells','completed_measured_cells','detected_cells','completed_zero_cells','detection_frequency','native_support_completed','support_unit','threshold','missing','failed','unprocessed','unknown','group_by'],...groups.map(g=>[c.id,g.group,g.cells,g.completed,g.detected,g.zero,g.completed?g.detected/g.completed:null,g.support,c.unit,threshold,g.missing,g.failed,g.unprocessed,g.unknown,E('group').value||'all_cells'])])],{type:'text/tab-separated-values'}),'circyto_detection_summary.tsv');});
document.querySelectorAll('#explorer [data-export]').forEach(b=>b.addEventListener('click',()=>{const svg=E(b.dataset.export).querySelector('svg');if(!svg)return;const clone=svg.cloneNode(true);clone.setAttribute('xmlns',NS);const raw=new XMLSerializer().serializeToString(clone);if(b.dataset.format==='svg')save(new Blob([raw],{type:'image/svg+xml'}),'circyto_'+b.dataset.export+'.svg');else{const img=new Image();img.onload=()=>{const canvas=document.createElement('canvas');canvas.width=1080;canvas.height=960;canvas.getContext('2d').drawImage(img,0,0,1080,960);canvas.toBlob(blob=>{if(blob)save(blob,'circyto_'+b.dataset.export+'.png');});};img.src='data:image/svg+xml;charset=utf-8,'+encodeURIComponent(raw);}}));
for(const kind of ['gene','circ']){E(kind+'-query').addEventListener('input',()=>search(kind));E(kind).addEventListener('change',()=>{if(kind==='gene'){gene=Number(E(kind).value);E('color').value='gene';}else circ=Number(E(kind).value);render();});}
option(E('color'),'Selected gene expression','gene');option(E('group'),'All cells','');
for(const m of D.metadata){const unique=new Set(m.values.filter(v=>v!==null).map(String));if(m.kind==='numeric'||unique.size<=128)option(E('color'),m.key,m.key);if(m.kind==='category'&&unique.size<=128)option(E('group'),m.key,m.key);}
const preferred=D.metadata.find(m=>m.kind==='category'&&/condition|exposure/.test(m.key));if(preferred)E('group').value=preferred.key;
E('group').addEventListener('change',()=>{refreshFilter();render();});for(const k of ['filter','color','mode','threshold'])E(k).addEventListener('change',render);
E('cell-query').addEventListener('keydown',e=>{if(e.key==='Enter'){const i=D.cells.indexOf(E('cell-query').value);if(i>=0)inspect(i);else E('hover').textContent='No exact cell identity matches this query.';}});
for(const text of [D.cells.length+' biological cells',D.join.overlap+' paired cells',D.embedding.plotted_cells+' saved RNA coordinates',D.genes.length+' genes / '+D.circs.length+' circRNAs','Offline · no network'])E('badges').append(node('span',text,{class:'ex-badge'}));
E('dataset').textContent=D.title;
E('scope').textContent=`Embedded feature subset: ${D.genes.length} / ${D.selection.genes_total} RNA genes and ${D.circs.length} / ${D.selection.circs_total} circRNA candidates. Search covers this subset only. Regenerate with --gene or --circ to include other exact IDs or aliases.`;
for(const warning of D.warnings)E('warnings').append(node('p',warning,{class:'ex-note'}));
E('policy').textContent=D.join.policy+'. '+D.selection.policy+' RNA values: '+D.selection.rna_layer+'. Embedding: '+D.embedding.key+'. No normalization, clustering or embedding fitting is performed.';
table(E('sources'),['Input basename','Bytes','SHA-256'],D.sources.map(s=>[s.name,s.bytes,s.sha256]));
E('plots').classList.toggle('ex-hidden',!validXY.length);E('no-umap').classList.toggle('ex-hidden',!!validXY.length);E('rna-panel').classList.toggle('ex-hidden',!D.join.rna_cells);E('circ-panel').classList.toggle('ex-hidden',!D.circs.length);
if(!D.genes.length){E('gene-query').disabled=true;E('gene').disabled=true;if(D.metadata.length)E('color').value=D.metadata[0].key;}
refreshFilter();search('gene');search('circ');render();
})();
"""


def script_policy():
    digest = base64.b64encode(hashlib.sha256(SCRIPT.encode()).digest()).decode()
    return f"default-src 'none'; script-src 'sha256-{digest}'; style-src 'unsafe-inline'; img-src data: blob:; connect-src 'none'; base-uri 'none'; form-action 'none'"


def render_explorer(payload):
    data = json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    for char, escaped in (("&", "\\u0026"), ("<", "\\u003c"), (">", "\\u003e"), ("\u2028", "\\u2028"), ("\u2029", "\\u2029")):
        data = data.replace(char, escaped)
    return f'<style>{STYLE}</style>{FRAGMENT}<script type="application/json" id="ex-data">{data}</script><script>{SCRIPT}</script>'


def render_explorer_document(payload):
    from circyto.pipeline.qc_report_style import STYLE as BASE_STYLE

    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="{html.escape(script_policy(), quote=True)}"><meta name="referrer" content="no-referrer">
<title>CIRCYTO · Single-cell explorer</title><style>{BASE_STYLE}</style></head><body>
<header><div class="brand">CIRCYTO <span>SINGLE-CELL REPORT</span></div><nav aria-label="Report sections"><a href="#explorer">Explore</a><a href="#ex-cells">Cell evidence</a></nav></header>
<main>{render_explorer(payload)}</main><footer>CIRCYTO · Optional exploration report · Saved RNA embedding · Native evidence · No network required</footer></body></html>\n'''
