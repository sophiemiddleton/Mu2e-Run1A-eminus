"""Dump the muon stops of each stopping-target geometry to parquet, with the stops per generated event.

For each geometry, runs dump_stops.fcl (Offline StoppedParticlesDumper) on the files of
sim.mu2e.MuminusStopsCat.MDC2025ba_<geometry>.art (the nominal target: MDC2025ac, the SimEfficiencies2
sample) and writes

    <outdir>/stops/stops_<geometry>.parquet   x, y, z (detector coordinates, mm), t (ns), one row per stop
    <outdir>/stops/stops_summary.yaml         per geometry: files, stops, gen.count, prescale, MuBeamCat events,
                                              stops per MuBeamCat event

The detector frame is the Mu2e frame shifted by DETECTOR_ORIGIN, the frame of the EventNtuple truth
(trkmcsim pos), so the stops can be compared with the CE/DIO start positions.

Every stops stage resampled the same MuBeamCat.MDC2025ab. The nominal one kept one stop in
NOMINAL_STOP_PRESCALE (1000), the geometries every stop, so stops per MuBeamCat event
(stops x prescale / gen.count) compare between all of them. stop_rate_rel in the registry is the
same ratio from rse.nevent, which counts events, not stops.

    python dump_stops_st.py --outdir /exp/mu2e/data/users/$USER/stoppingtarget --max-files 3 --jobs 3

Each file is ~1.3 GB and ~585k stops; a file takes a few minutes and ~1 GB of memory.
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import uproot
import yaml

PYFITTER_DIR = Path(os.environ.get('PYFITTER_DIR', Path(__file__).resolve().parents[1] / 'pyfitter'))
sys.path.insert(0, str(PYFITTER_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_registry import SETUP_MU2E, SETUP_OPS, run_mu2e  # noqa: E402
from make_registry_st import (CAMPAIGN, GEOMETRIES, GEOMETRY_STOP_PRESCALE, NOMINAL, NOMINAL_STOP_CAMPAIGN,  # noqa: E402
                              NOMINAL_STOP_PRESCALE)

FCL = Path(__file__).resolve().parent / 'dump_stops.fcl'
DETECTOR_ORIGIN = np.array([-3904.0, 0.0, 10171.0])  # Mu2e frame position of the detector-frame origin (mm)


def sim_path(filename, location):
    """/pnfs path of a sim-tier file (as make_registry.pnfs_path, for phy-sim)."""
    tier, owner, desc, config, _seq, fmt = filename.split('.')
    digest = hashlib.sha256(filename.encode()).hexdigest()
    area = {'tape': 'tape/phy-sim', 'disk': 'persistent/datasets/phy-sim', 'scratch': 'scratch/datasets/phy-sim'}[location]
    return f'/pnfs/mu2e/{area}/{tier}/{owner}/{desc}/{config}/{fmt}/{digest[:2]}/{digest[2:4]}/{filename}'


def stop_source(geometry):
    """(campaign, MuminusStopsCat dataset, prescale) of a geometry's stops."""
    if geometry == NOMINAL:
        return NOMINAL_STOP_CAMPAIGN, f'mu2e:sim.mu2e.MuminusStopsCat.{NOMINAL_STOP_CAMPAIGN}.art', NOMINAL_STOP_PRESCALE
    return CAMPAIGN, f'mu2e:sim.mu2e.MuminusStopsCat.{CAMPAIGN}_{geometry}.art', GEOMETRY_STOP_PRESCALE


def stop_files(geometry, location):
    """[(path, gen.count)] of the geometry's MuminusStopsCat files, sorted by name."""
    dataset = stop_source(geometry)[1]
    files = json.loads(run_mu2e(SETUP_OPS, f'metacat query -m gen.count -j "files from {dataset}"'))
    return sorted((sim_path(f['name'], location), f['metadata']['gen.count']) for f in files)


def dump(path, out, campaign):
    """Run the dumper on one file with the campaign's SimJob; return the stops as a DataFrame in detector coordinates."""
    if not out.exists():
        setup = f'{SETUP_MU2E} && muse setup SimJob {campaign} >/dev/null 2>&1'
        run_mu2e(setup, f'mu2e -c {FCL} -s {path} -T {out}.tmp > {out}.log 2>&1 && mv {out}.tmp {out}')
    stops = uproot.open(out)['stops/stops']['stops'].array(library='np')
    return pd.DataFrame({
        'x': stops['x'] - DETECTOR_ORIGIN[0],
        'y': stops['y'] - DETECTOR_ORIGIN[1],
        'z': stops['z'] - DETECTOR_ORIGIN[2],
        't': stops['time'],
    }).astype(np.float32)


def main(args):
    out_dir = Path(args.outdir).expanduser().resolve() / 'stops'
    work = out_dir / 'root'
    work.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / 'stops_summary.yaml'
    summary = yaml.safe_load(summary_path.read_text()) if summary_path.exists() else {}

    for geometry in args.geometries:
        campaign, _, prescale = stop_source(geometry)
        files = stop_files(geometry, args.loc)[:args.max_files]
        print(f'[{geometry}] dumping {len(files)} files')
        with ThreadPoolExecutor(args.jobs) as pool:
            frames = list(pool.map(lambda f: dump(f[0], work / f'{Path(f[0]).stem}.root', campaign), files))
        df = pd.concat(frames, ignore_index=True)
        df.to_parquet(out_dir / f'stops_{geometry}.parquet', index=False)
        gen = int(sum(n for _, n in files))
        mubeam = gen / prescale
        summary[geometry] = {'campaign': campaign, 'files': len(files), 'stops': len(df), 'gen_count': gen,
                             'prescale': prescale, 'mubeam_events': mubeam, 'stops_per_mubeam': len(df) / mubeam}
        print(f'[{geometry}] {len(df)} stops, gen.count {gen}, prescale {prescale}: {len(df) / mubeam:.4g} per MuBeamCat event')

    summary_path.write_text(yaml.safe_dump(summary, sort_keys=False))
    print(f'Wrote {summary_path}')
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Muon stop positions of the stopping-target geometry samples')
    parser.add_argument('--outdir', required=True, help='Study directory (the make_registry_st.py --outdir)')
    parser.add_argument('--geometries', nargs='+', default=GEOMETRIES, help='Geometries to dump')
    parser.add_argument('--max-files', type=int, default=None, help='Files per geometry (default: all)')
    parser.add_argument('--jobs', type=int, default=3, help='Files dumped at once (~1 GB of memory each)')
    parser.add_argument('--loc', default='disk', choices=['tape', 'disk', 'scratch'], help='Location of the stop files')
    raise SystemExit(main(parser.parse_args()))
