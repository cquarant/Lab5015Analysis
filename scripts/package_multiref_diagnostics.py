"""Package original step2 PNGs for the static viewer, one download per DUT bar."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import tarfile


def package(sample, dest):
    receipt = json.loads((sample/'step2_diagnostics.json').read_text())
    archive = sample/'step2_diagnostics.tar.gz'
    pattern = re.compile(r'c_(.+)_bar(\d+)_Vov3\.00_th(\d+)_energyBin01\.png')
    bars = defaultdict(dict)
    with tarfile.open(archive, 'r|gz') as tar:
        for member in tar:
            if not member.isfile():
                continue
            quantity, bar, threshold = pattern.fullmatch(Path(member.name).name).groups()
            key = f'{int(threshold)}/{quantity}'
            image = tar.extractfile(member).read()
            if key in bars[bar] or not image.startswith(b'\x89PNG\r\n\x1a\n'):
                raise ValueError(member.name)
            bars[bar][key] = image
    if sum(map(len, bars.values())) != receipt['plots']:
        raise ValueError(str(archive))
    dest.mkdir(parents=True, exist_ok=True)
    quantities = set()
    for bar, images in bars.items():
        offsets = {}
        size = 0
        for key, image in images.items():
            offsets[key] = [size, len(image)]
            size += len(image)
            quantities.add(key.split('/')[1])
        header = json.dumps(offsets, separators=(',', ':')).encode()
        path = dest/f'bar{bar}.bin'
        with path.with_suffix('.partial').open('wb') as out:
            out.write(len(header).to_bytes(4, 'little'))
            out.write(header)
            out.writelines(images.values())
        path.with_suffix('.partial').replace(path)
    return quantities


def package_all(production, output, summaries=None):
    inventory = set()
    for run in ('4194', '4257', '4281', '4291'):
        for ref in ('04','05','07','08','09','10','11','12'):
            folder = Path(run)/f'ref{ref}'
            if summaries is not None:
                summary = json.loads((summaries/f'{run}_ref{ref}_summary.json').read_text())
                receipt = json.loads((production/folder/'step2_diagnostics.json').read_text())
                if summary['input_uuids']['corrections'] != receipt['step2_uuid']:
                    raise ValueError(f'Wrong diagnostic archive for {folder}')
            inventory.update(package(production/folder, output/folder))
    (output/'quantities.json').write_text(json.dumps(sorted(inventory)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    package_all(args.production, args.output)
