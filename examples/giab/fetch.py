"""Download the pinned public files for the GIAB example."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('directory', type=Path)
parser.add_argument('--inputs', action='store_true', help='Also download both gVCFs and the reference (about 5.2 GB unpacked)')
args = parser.parse_args()
args.directory.mkdir(parents=True, exist_ok=True)
source_dir = Path(__file__).resolve().parent
sources = json.loads((source_dir / 'sources.json').read_text())
for source in sources:
    if source['kind'] != 'regions' and not args.inputs:
        continue
    name = source['file']
    if Path(name).name != name:
        raise ValueError('Download filenames must be plain basenames')
    path = args.directory / name
    partial = path.with_name(name + '.partial')
    if not path.exists():
        print('Downloading', name, flush=True)
        with urllib.request.urlopen(source['url'], timeout=120) as response, partial.open('wb') as handle:
            stream = gzip.GzipFile(fileobj=response) if source.get('gunzip') else response
            shutil.copyfileobj(stream, handle, 1024 * 1024)
        candidate = partial
    else:
        candidate = path
    if candidate.stat().st_size != source['bytes']:
        raise ValueError(f'Incorrect file size: {candidate}; remove it and download again')
    if 'sha256' in source:
        digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
        if digest != source['sha256']:
            raise ValueError(f'Checksum mismatch: {candidate}')
    if candidate == partial:
        partial.rename(path)
shutil.copyfile(source_dir / 'autosomes.bed', args.directory / 'autosomes.bed')
print('Ready. Region files have checksums; large inputs were checked by size only.')
