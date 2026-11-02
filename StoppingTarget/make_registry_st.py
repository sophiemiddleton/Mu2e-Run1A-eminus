"""Write a pyfitter dataset registry for the MDC2025ba stopping-target geometry samples.

For each geometry in GEOMETRIES and each process in PROCESSES (CE and DIO, unmixed OnSpill):

    filelist   = the nts dataset's /pnfs paths, written to <outdir>/filelists/
    acceptance = sum(rse.nevent) / sum(gen.count) over the mcs dataset's files (metacat)
    expected   = generated-phase-space yield for the run, from pyfitter's make_registry.expected_yield

The run is --days of real time in beam mode --beam-mode at --beam-power kW (default 28 days,
1BB, 3.84 kW). POT follows pyfitter's normalization.py: the intensity is scaled from the mode's
quoted power to --beam-power at fixed timing.

The run's stopped muons are the nominal-target SimEfficiencies2 value (Production normalizations.py):
MuBeamCat x MuminusStopsCat x 1000, the 1000 undoing the stops stage's MuminusTargetStopPrescale.
Each geometry's stops stage resampled the same MuBeamCat.MDC2025ab with nPrescale 1
(cnf.mu2e.TargetStops.MDC2025ba_<geometry>.tar), so its stops per POT are
MuBeamCat x MuminusStopsCat with no prescale. stop_rate_rel is that over the nominal value
(nominal: 1), and plot_st.py scales each geometry's yields by it.

    python make_registry_st.py --outdir /exp/mu2e/data/users/$USER/stoppingtarget
"""

import argparse
import json
import os
import sys
from pathlib import Path

import yaml

PYFITTER_DIR = Path(os.environ.get('PYFITTER_DIR', Path(__file__).resolve().parents[1] / 'pyfitter'))
sys.path.insert(0, str(PYFITTER_DIR))

from make_registry import (PRODUCTION_DIR, SETUP_OPS, expected_yield, load_production_constants,  # noqa: E402
                           mcs_acceptance, run_mu2e, sim_efficiencies, write_filelist)
from normalization import beam_power_kW, get_pot, on_spill_s, protons_per_spill_for_power  # noqa: E402

CAMPAIGN = 'MDC2025ba'
RECO = 'best_v1_5'
NOMINAL = 'nominal'
MUBEAM_DATASET = 'mu2e:sim.mu2e.MuBeamCat.MDC2025ab.art'  # input of every stops stage
GEOMETRY_STOP_PRESCALE = 1  # TargetStopPrescaleFilter.nPrescale of the geometry stops stage
NOMINAL_STOP_CAMPAIGN = 'MDC2025ac'  # the SimEfficiencies2 stops sample (Sim_best v1_1)
NOMINAL_STOP_PRESCALE = 1000  # its MuminusTargetStopPrescale (Production pileup/prolog.fcl)
STOP_GEOMETRIES = ['st37plus6up6dn', 'st49hole1to5', 'st49hole9to41']  # with a MuminusStopsCat dataset
GEOMETRIES = [NOMINAL] + STOP_GEOMETRIES
PROCESSES = {
    'CE':  {'desc': 'CeEndpointOnSpill', 'process': 'CE'},
    'DIO': {'desc': 'DIOtail95OnSpill',  'process': 'DIO', 'e_min': 95.0},
}


def dataset_rate(dataset):
    """Sum of rse.nevent over sum of gen.count, as mu2eGenFilterEff."""
    files = json.loads(run_mu2e(SETUP_OPS, f'metacat query -m rse.nevent,gen.count -j "files from {dataset}"'))
    return sum(f['metadata']['rse.nevent'] for f in files) / sum(f['metadata']['gen.count'] for f in files)


def stops_per_pot(geometry):
    """Stopped mu- per POT for a geometry: MuBeamCat x MuminusStopsCat x prescale."""
    stops = dataset_rate(f'mu2e:sim.mu2e.MuminusStopsCat.{CAMPAIGN}_{geometry}.art')
    return dataset_rate(MUBEAM_DATASET) * stops * GEOMETRY_STOP_PRESCALE


def sample_config(geometry):
    """Dataset config of a geometry's mcs and nts samples (the nominal target has no geometry tag)."""
    return f'{CAMPAIGN}_{RECO}' if geometry == NOMINAL else f'{CAMPAIGN}_{geometry}_{RECO}'


def main(args):
    outdir = Path(args.outdir).expanduser().resolve()
    c = load_production_constants()
    sim = sim_efficiencies()

    pps = protons_per_spill_for_power(args.beam_power, args.beam_mode)
    onspill = on_spill_s(args.days, args.beam_mode)
    pot = get_pot(onspill, args.beam_mode, pps)
    print(f'Run: {args.days} days {args.beam_mode} at {beam_power_kW(args.beam_mode, pps):.3g} kW, '
          f'on-spill {onspill:.4g} s, POT {pot:.4g}, stopped mu {pot * sim["target_stopped_muons_per_pot"]:.4g}')

    nominal_stops = sim['target_stopped_muons_per_pot']
    stop_rate_rel = {g: stops_per_pot(g) / nominal_stops for g in STOP_GEOMETRIES}
    stop_rate_rel[NOMINAL] = 1.0
    for g, r in stop_rate_rel.items():
        print(f'[{g}] stopped mu- per POT {r * nominal_stops:.4g}, relative to nominal {r:.4f}')
    datasets, groups = {}, {}
    for geometry in GEOMETRIES:
        groups[geometry] = []
        for label, sample in PROCESSES.items():
            name = f'{label}_{geometry}'
            config = sample_config(geometry)
            nevent, ngen, nfiles = mcs_acceptance(sample['desc'], config)
            path, nts, nlist = write_filelist(name, sample['desc'], config, outdir / 'filelists', args.loc)
            datasets[name] = {
                'filelist': str(path),
                'process': sample['process'],
                'acceptance': nevent / ngen,
                'expected': expected_yield(sample, pot, sim, c, args.rmue),
                'geometry': geometry,
                'stop_rate_rel': stop_rate_rel[geometry],
                'mcs_dataset': f'mcs.mu2e.{sample["desc"]}.{config}.art',
                'mcs_events': nevent,
                'generated_events': ngen,
                'nts_dataset': nts,
            }
            groups[geometry].append(name)
            print(f'[{name}] {nlist} nts files, acceptance {nevent / ngen:.4g}, expected {datasets[name]["expected"]:.4g}')
    groups['all'] = list(datasets)

    run = {
        'beam_mode': args.beam_mode,
        'realtime_days': args.days,
        'beam_power_kW': args.beam_power,
        'on_spill_s': float(onspill),
        'pot': float(pot),
        'stopped_muons': float(pot * sim['target_stopped_muons_per_pot']),
        'captures_per_stopped_muon': float(c.CAPTURES_PER_STOPPED_MUON),
        'rmue': args.rmue,
    }
    header = (f'# Stopping-target geometry registry, GENERATED by StoppingTarget/make_registry_st.py\n'
              f'#   samples: {CAMPAIGN}[_<geometry>]_{RECO}, normalizations: {PRODUCTION_DIR}\n'
              f'#   run: {args.days} days {args.beam_mode} at {args.beam_power} kW, POT {pot:.4g}\n\n')
    registry = outdir / 'datasets_st.yaml'
    registry.write_text(header + yaml.safe_dump({'run': run, 'datasets': datasets, 'groups': groups}, sort_keys=False))
    print(f'Wrote {registry}')
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Registry for the MDC2025ba stopping-target geometry samples')
    parser.add_argument('--outdir', required=True, help='Directory for the registry and file lists')
    parser.add_argument('--days', type=float, default=28, help='Real (on+off-spill) running time in days')
    parser.add_argument('--beam-mode', default='1BB', help='Beam mode in pyfitter normalization.py')
    parser.add_argument('--beam-power', type=float, default=3.84, help='Beam power in kW')
    parser.add_argument('--rmue', type=float, default=1e-15, help='R_mue for the CE expected yield')
    parser.add_argument('--loc', default='tape', choices=['tape', 'disk', 'scratch'], help='Location of the nts files')
    raise SystemExit(main(parser.parse_args()))
