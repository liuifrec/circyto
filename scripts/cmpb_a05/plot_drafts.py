"""Deterministic SVG plots from small A05 TSVs; stdlib only, no fitted models."""
import argparse
import csv
import html
import hashlib
import json
import math
from pathlib import Path
import subprocess


def read(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class SVG:
    def __init__(self, title, subtitle, width, height, metadata):
        self.w, self.h = width, height
        self.parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            '<style>text{font-family:DejaVu Sans,Arial,sans-serif;fill:#23313d} .axis{stroke:#7c8a96;stroke-width:1} .grid{stroke:#e6ebef;stroke-width:1}</style>',
            '<rect width="100%" height="100%" fill="white"/>',
            '<metadata>'+html.escape(json.dumps(metadata, sort_keys=True))+'</metadata>']
        self.text(35, 31, title, 21, weight="bold")
        self.text(35, 55, subtitle, 12)

    def text(self, x, y, text, size=12, anchor="start", weight="normal"):
        self.parts.append(f'<text x="{x:.4f}" y="{y:.4f}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}">{html.escape(str(text))}</text>')

    def line(self, x1, y1, x2, y2, color="#7c8a96"):
        self.parts.append(f'<line x1="{x1:.4f}" y1="{y1:.4f}" x2="{x2:.4f}" y2="{y2:.4f}" stroke="{color}"/>')

    def dot(self, x, y, label, color="#176c89", radius=3, identifier="", xr=0, yr=0):
        # Self-closing circle also works in the installed ImageMagick SVG reader,
        # which silently drops circles containing a nested SVG title element.
        self.parts.append(f'<circle cx="{x:.6f}" cy="{y:.6f}" r="{radius}" fill="{color}" fill-opacity="0.68" aria-label="{html.escape(label, quote=True)}" data-id="{html.escape(identifier, quote=True)}" data-x="{xr}" data-y="{yr}"/>')

    def save(self, path, footer):
        for i, line in enumerate(footer): self.text(35, self.h-40+i*15, line, 10)
        path.write_text("\n".join(self.parts+["</svg>"])+"\n", encoding="utf-8")


def ticks(top):
    if top <= 0: return 1
    power = 10**math.floor(math.log10(top))
    return math.ceil(top/power*2)*power/2


def axes(svg, x, y, w, h, xmax, ymax, xlabel, ylabel, title):
    svg.text(x, y-15, title, 14, weight="bold")
    for i in range(5):
        fx, fy = x+w*i/4, y+h-h*i/4
        svg.line(x, fy, x+w, fy, "#e6ebef")
        svg.text(fx, y+h+19, f"{xmax*i/4:,.3g}", 10, "middle")
        svg.text(x-8, fy+4, f"{ymax*i/4:,.3g}", 10, "end")
    svg.line(x,y+h,x+w,y+h); svg.line(x,y,x,y+h)
    svg.text(x+w/2,y+h+40,xlabel,11,"middle")
    svg.text(x,y-34,ylabel,11)
    return lambda a,b:(x+w*a/xmax,y+h-h*b/ymax)


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--tables",type=Path,required=True);ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--renderer-commit",help="Pin the recorded renderer commit when recreating from a ZIP without Git metadata")
    args=ap.parse_args();args.out.mkdir(exist_ok=True,parents=True)
    cells=read(args.tables/"descriptive_cell_metrics.tsv");candidates=read(args.tables/"descriptive_candidate_metrics.tsv");thresholds=read(args.tables/"threshold_sensitivity.tsv")
    provenance=json.loads((args.tables/"analysis_provenance.json").read_text())
    observations=[]; manifest=[]
    # Renderer commit binds plotting code separately from analysis helper commit.
    repo=Path(__file__).resolve().parents[2]
    renderer_commit=args.renderer_commit or (subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip() if (repo/".git").exists() else provenance['helper_commit'])
    renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    assert json.loads((args.tables/'cnv_gate.json').read_text())['decision']=='NOT VERIFIED', 'Review gate decision before choosing this RNA/circ-only draft'
    def record(fig,panel,dataset,oid,x,y,source,order):
        observations.append(dict(figure=fig,panel=panel,dataset=dataset,observation_id=oid,display_order=order,
            x=x,y=y,source_table=source,source_object_sha256=provenance['source_hashes'][dataset],analysis_helper_commit=provenance['helper_commit'],renderer_commit=renderer_commit))
    def save(svg, name, dataset, source_tables, note):
        sha=provenance['source_hashes'][dataset]
        footer=[note, f"Object SHA256 {sha} | helper {provenance['helper_commit'][:12]}; renderer {renderer_commit[:12]}"]
        svg.save(args.out/name,footer)
        manifest.append(dict(figure=name,dataset=dataset,source_tables=source_tables,source_object_sha256=sha,
            analysis_helper_commit=provenance['helper_commit'],renderer_commit=renderer_commit,renderer_sha256=renderer_sha256,scope_note=note))
    for dataset,slug in [("Smart-seq3","smartseq3"),("IMR90","imr90")]:
        cc=[r for r in cells if r['dataset']==dataset];cv=[r for r in candidates if r['dataset']==dataset]
        matched=[r for r in cc if r['matched_rna_circ']=='True']
        meta=dict(dataset=dataset,**provenance,renderer_commit=renderer_commit)
        name=slug+'_rna_circ_scatter_draft.svg'
        s=SVG(dataset+" | stored RNA and circ burden",f"All {len(matched)} recorded shared identities; historical circ completion unknown",1140,480,meta)
        for panel,(xm,ym,xl,yl) in enumerate([
            ('rna_overlap_signal','circ_support','RNA overlap-profile total (stored counts)','circ junction-read support (stored units)'),
            ('rna_positive_features','circ_positive_candidates','RNA positive features (count)','positive circ candidates (count)')]):
            x=[int(r[xm]) for r in matched]; y=[int(r[ym]) for r in matched]
            coord=axes(s,85+panel*555,110,440,245,ticks(max(x)),ticks(max(y)),xl,yl,'A' if panel==0 else 'B')
            for i,r in enumerate(matched):
                xp,yp=coord(x[i],y[i]);s.dot(xp,yp,f"{r['cell_id']}: {x[i]}, {y[i]}",identifier=r['cell_id'],xr=x[i],yr=y[i])
                record(name,str(panel),dataset,r['cell_id'],x[i],y[i],'descriptive_cell_metrics.tsv',i)
        note='Descriptive stored measurements; no correlation test, causal inference or cross-assay comparison.'
        if dataset=='Smart-seq3': note+=' RNA-only DIYHEK_192 is retained in tables/distributions.'
        save(s,name,dataset,['descriptive_cell_metrics.tsv'],note)
        name=slug+'_distributions_draft.svg'
        s=SVG(dataset+' | stored-signal distributions',f"All stored rows (N={len(cc)}) and candidates (N={len(cv)}); zeros retained; completion unknown",1140,830,meta)
        for panel,(rows,field,label) in enumerate([(cc,'circ_positive_candidates','positive circ candidates per stored row'),(cc,'circ_support','circ support per stored row'),(cv,'positive_stored_cells','positive stored rows per candidate'),(cv,'pooled_support','pooled support per candidate')]):
            vals=sorted([(int(r[field]),i,r) for i,r in enumerate(rows)],key=lambda z:(z[0],z[1]))
            coord=axes(s,85+(panel%2)*555,110+(panel//2)*330,440,225,ticks(max(v[0] for v in vals)),1,label,'fraction of stored '+('rows' if panel<2 else 'candidates'),chr(65+panel))
            path=[]; prev=0
            for i,(value,_,r) in enumerate(vals):
                fraction=(i+1)/len(vals); xp,yp=coord(value,fraction)
                xb,yb=coord(value,prev)
                if not path:
                    x0,y0=coord(0,0);path.append(f'M{x0:.6f},{y0:.6f}')
                path.extend([f'L{xb:.6f},{yb:.6f}',f'L{xp:.6f},{yp:.6f}']);prev=fraction
                oid=r['cell_id'] if panel<2 else r['candidate_id']
                record(name,str(panel),dataset,oid,value,fraction,'descriptive_cell_metrics.tsv' if panel<2 else 'descriptive_candidate_metrics.tsv',i)
            s.parts.append('<path d="'+' '.join(path)+'" fill="none" stroke="#176c89" stroke-width="1.6"/>')
        save(s,name,dataset,['descriptive_cell_metrics.tsv','descriptive_candidate_metrics.tsv'],'Empirical cumulative distributions; stored-row fractions are not completed-cohort or biological prevalence.')
        name=slug+'_threshold_sensitivity_draft.svg'
        s=SVG(dataset+' | prespecified support thresholds',f"Each threshold starts from original X; {len(cc)} stored rows and {len(cv)} original candidates",1140,490,meta)
        rr=[r for r in thresholds if r['dataset']==dataset]
        for panel,field,label in [(0,'retained_candidates','retained candidate features (count)'),(1,'retained_support','retained support (stored junction-read units)')]:
            coord=axes(s,85+panel*555,110,440,245,5,ticks(max(int(r[field]) for r in rr)),'support threshold (units)',label,chr(65+panel))
            for mode,color in [('individual_entry','#176c89'),('pooled_candidate','#bb6d23')]:
                points=[]
                for i,r in enumerate([x for x in rr if x['threshold_mode']==mode]):
                    x,y=int(r['threshold']),int(r[field]);xp,yp=coord(x,y);points.append((xp,yp))
                    s.dot(xp,yp,mode+f': {x}, {y}',color,4,mode+':'+str(x),x,y)
                    record(name,str(panel),dataset,mode+':'+str(x),x,y,'threshold_sensitivity.tsv',i)
                s.parts.append('<polyline points="'+' '.join(f'{x:.6f},{y:.6f}' for x,y in points)+'" fill="none" stroke="'+color+'"/>')
        s.text(85,417,'Blue: individual entry >= threshold. Orange: pooled candidate sum >= threshold; all its original entries retained.',12)
        save(s,name,dataset,['threshold_sensitivity.tsv'],'No biological cutoff is recommended; zero rows remain in every threshold denominator.')
    dataset='IMR90';cc=[r for r in cells if r['dataset']==dataset];name='figure2c_rna_circ_only_alternative.svg'
    s=SVG('Figure 2C alternative | IMR90 stored RNA and circ signals','23 shared identities in original source order. CNV callability NOT VERIFIED; RNA/circ-only review draft.',1250,880,provenance)
    fields=[('circ_positive_candidates','Positive circ candidates','count'),('circ_support','circ support','junction-read units'),('rna_positive_features','RNA positive features','count'),('rna_overlap_signal','RNA overlap-profile total','stored counts')]
    for panel,(field,title,units) in enumerate(fields):
        x0=190+panel*260;w=210;top=ticks(max(int(r[field]) for r in cc))
        s.text(x0,90,title,13,weight='bold');s.text(x0,109,units,11)
        for tick in range(3):
            x=x0+w*tick/2;s.line(x,126,x,775,'#e6ebef');s.text(x,798,f'{top*tick/2:,.3g}',10,'middle')
        for i,r in enumerate(cc):
            y=140+i*28
            if panel==0:s.text(172,y+4,r['cell_id'],11,'end')
            v=int(r[field]);x=x0+w*v/top
            s.dot(x,y,r['cell_id']+f': {v}',radius=4,identifier=r['cell_id'],xr=v,yr=i)
            record(name,str(panel),dataset,r['cell_id'],v,i,'descriptive_cell_metrics.tsv',i)
    save(s,name,dataset,['descriptive_cell_metrics.tsv'],'Historical simple-overlap RNA; circ assay completion and original reference digests unknown. No CNV/association inference.')
    with (args.out/'figure_observations.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(observations[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(observations)
    (args.out/'figure_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(figures=len(manifest),plotted_observations=len(observations))))


if __name__=='__main__':main()
