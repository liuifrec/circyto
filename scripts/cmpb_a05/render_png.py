"""Rasterize authored scientific SVG drafts with existing ImageMagick; no installs."""
import json
from pathlib import Path
import struct
import subprocess
import sys

from run_stage import digest


def main():
    folder=Path(sys.argv[1]);manifest=json.loads((folder/'figure_manifest.json').read_text());records=[]
    version=subprocess.check_output(['/usr/bin/convert','-version'],text=True).splitlines()[0]
    for item in manifest:
        source=folder/item['figure'];target=source.with_suffix('.png')
        assert not target.exists(),target
        argv=['/usr/bin/convert','-background','white','-density','120',str(source),'-strip',str(target)]
        subprocess.run(argv,check=True)
        data=target.read_bytes();assert data[:8]==b'\x89PNG\r\n\x1a\n'
        width,height=struct.unpack('>II',data[16:24])
        records.append(dict(svg=source.name,png=target.name,svg_sha256=digest(source),png_sha256=digest(target),width=width,height=height,renderer=version,argv=argv))
    (folder/'png_render_record.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(dict(png_previews=len(records),renderer=version)))


if __name__=='__main__':main()
