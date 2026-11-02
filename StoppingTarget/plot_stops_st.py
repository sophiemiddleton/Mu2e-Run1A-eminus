"""Muon stops per stopping-target geometry: where they stop, how many, and how much CE signal each foil gives.

Reads the dump_stops_st.py output (<outdir>/stops) and, with --tag, the run_chain.py CE output. The stop
rates are per MuBeamCat event: every stops stage (the nominal MDC2025ac one, prescale 1000, and the
geometries', prescale 1) resampled the same MuBeamCat, so they compare directly.

Foils are found from the stop z positions (a gap of more than --foil-gap mm starts a new foil),
numbered 1, 2, ... from upstream.

Writes to <outdir>/stops/plots/:
    stops_z.{png,pdf}        stops per foil against foil z, ratio to the first geometry
    stops_r.{png,pdf}        stop radius, ratio to the first geometry
    stops_rz.{png,pdf}       stop r against z for each geometry (the holes show up as empty bands)
    stops_foils.csv          per foil: z, inner and outer stop radius, stops per MuBeamCat event

and with --tag, from the CE start position (the stop position, pos_z_true in the post-cut parquet,
written by process.py since it records the truth position):
    ce_foils.{png,pdf}       per foil: CE signal efficiency, mean momentum response, and CE signal for the run
    ce_foils.csv

    python plot_stops_st.py --outdir /exp/mu2e/data/users/$USER/stoppingtarget --tag /exp/mu2e/data/users/$USER/stoppingtarget/run1a
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

PYFITTER_DIR = Path(os.environ.get('PYFITTER_DIR', Path(__file__).resolve().parents[1] / 'pyfitter'))
sys.path.insert(0, str(PYFITTER_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from config import get_version_cuts  # noqa: E402
from plot_st import (GEOMETRY_COLORS, GEOMETRY_LABELS, GEOMETRY_MARKERS, GEOMETRY_STYLES,  # noqa: E402
                     generated_in_sample, load_columns, save)
from style import COLORS, apply_root_text, format_run_info  # noqa: E402

INFO = 'Same MuBeamCat input; nominal from MDC2025ac'


def mubeam_events(s):
    """MuBeamCat events behind a geometry's stops: gen.count / prescale (summaries before the prescale had 1)."""
    return s.get('mubeam_events', s['gen_count'] / s.get('prescale', 1))


def stops_per_mubeam(s):
    return s['stops'] / mubeam_events(s)


def mu2e_label(ax):
    """Bold 'Mu2e' and italic 'Simulation' above the top-left of the frame, as style.apply_root_text."""
    mu2e = ax.text(0.0, 1.03, 'Mu2e', transform=ax.transAxes, ha='left', va='bottom', fontsize=16, fontweight='bold')
    ax.annotate(' Simulation', xy=(1, 0), xycoords=mu2e, ha='left', va='bottom', fontsize=16, style='italic')


def find_foils(z, gap):
    """Foil edges in z: the stops sorted in z, split where consecutive stops are more than gap mm apart."""
    zs = np.sort(z)
    breaks = np.nonzero(np.diff(zs) > gap)[0]
    starts = np.concatenate([[zs[0]], zs[breaks + 1]])
    ends = np.concatenate([zs[breaks], [zs[-1]]])
    return np.column_stack([starts, ends])


def foil_index(z, foils):
    """1-based foil number of each z (0 outside every foil)."""
    i = np.searchsorted(foils[:, 0], z, side='right') - 1
    inside = (i >= 0) & (z <= foils[np.clip(i, 0, None), 1])
    return np.where(inside, i + 1, 0)


def per_foil(stops, summary, gap):
    """{geometry: DataFrame per foil} of the stops."""
    out = {}
    for geometry, df in stops.items():
        foils = find_foils(df['z'].to_numpy(), gap)
        index = foil_index(df['z'].to_numpy(), foils)
        r = np.hypot(df['x'], df['y']).to_numpy()
        per_gen = 1.0 / mubeam_events(summary[geometry])
        rows = []
        for i, (z0, z1) in enumerate(foils, start=1):
            sel = index == i
            rows.append({'foil': i, 'z': 0.5 * (z0 + z1), 'stops': int(sel.sum()), 'stops_per_gen': sel.sum() * per_gen,
                         'r_min': np.percentile(r[sel], 0.1), 'r_max': np.percentile(r[sel], 99.9)})
        out[geometry] = pd.DataFrame(rows)
    return out


def ratio_points(rax, x, y, ref, color, marker, label):
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = y / ref
    rax.plot(x, ratio, marker=marker, ls='none', ms=4, color=color, label=label)


def ratio_legend(rax, geometries):
    rax.legend(loc='upper left', frameon=False, fontsize=9, ncol=2,
               title=f'Ratio to {GEOMETRY_LABELS.get(geometries[0], geometries[0])}', title_fontsize=9)


def plot_z(foils, geometries, summary, out_dir):
    """Stops per foil against foil z, with the totals in the legend and the ratio to the first geometry."""
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    ref = stops_per_mubeam(summary[geometries[0]])
    for geometry, color, marker in zip(geometries, GEOMETRY_COLORS, GEOMETRY_MARKERS):
        f = foils[geometry]
        total = stops_per_mubeam(summary[geometry])
        name = GEOMETRY_LABELS.get(geometry, geometry)
        name = name if f'{len(f)} foils' in name else f'{name} ({len(f)} foils)'
        label = f'{name}: {total * 1e3:.1f} per $10^3$ ({total / ref - 1:+.1%})'
        ax.plot(f['z'], f['stops_per_gen'] * 1e3, marker=marker, ms=5, lw=1, color=color, label=label)
    # Ratio per foil only for foils within 2 mm in z of one of the first geometry's foils
    ref_f = foils[geometries[0]]
    ratios = []
    for geometry, color, marker in zip(geometries[1:], GEOMETRY_COLORS[1:], GEOMETRY_MARKERS[1:]):
        f = foils[geometry]
        nearest = np.abs(f['z'].to_numpy()[:, None] - ref_f['z'].to_numpy()[None, :]).argmin(axis=1)
        matched = np.abs(f['z'].to_numpy() - ref_f['z'].to_numpy()[nearest]) < 2
        ref_y = np.where(matched, ref_f['stops_per_gen'].to_numpy()[nearest], np.nan)
        ratio_points(rax, f['z'], f['stops_per_gen'], ref_y, color, marker, GEOMETRY_LABELS.get(geometry, geometry))
        ratios.append(f['stops_per_gen'].to_numpy() / ref_y)
    rax.axhline(1.0, color='black', lw=0.8)
    finite = np.concatenate(ratios)[np.isfinite(np.concatenate(ratios))] if ratios else np.array([1.0])
    rax.set_ylim(0.0, max(2.0, 1.15 * finite.max()))
    ax.set_ylim(bottom=0, top=ax.get_ylim()[1] * 1.35)
    apply_root_text(fig, ax, rax, 'Foil z, detector frame (mm)', r'Stops per foil per $10^3$ MuBeamCat',
                    'Ratio, same z', info_text=INFO)
    rax.yaxis.label.set_size(14)
    ratio_legend(rax, geometries)
    ax.legend(loc='upper right', frameon=False, fontsize=11, title=r'$\mu^-$ stops; total per $10^3$ MuBeamCat', title_fontsize=11)
    save(fig, out_dir, 'stops_z')


def plot_r(stops, geometries, summary, out_dir, edges):
    """Stop radius per generated event, ratio to the first geometry."""
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    ref = None
    centres = 0.5 * (edges[1:] + edges[:-1])
    for geometry, color, ls in zip(geometries, GEOMETRY_COLORS, GEOMETRY_STYLES):
        df = stops[geometry]
        y = np.histogram(np.hypot(df['x'], df['y']), edges)[0] * 1e3 / mubeam_events(summary[geometry])
        ax.stairs(y, edges, color=color, ls=ls, lw=1.5, label=GEOMETRY_LABELS.get(geometry, geometry))
        if ref is None:
            ref = y
        else:
            ratio_points(rax, centres, y, ref, color, 'o', GEOMETRY_LABELS.get(geometry, geometry))
    rax.axhline(1.0, color='black', lw=0.8)
    rax.set_ylim(0.0, 2.5)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(bottom=0, top=ax.get_ylim()[1] * 1.25)
    apply_root_text(fig, ax, rax, 'Stop radius (mm)', fr'Stops / {edges[1] - edges[0]:.2g} mm per $10^3$ MuBeamCat',
                    'Ratio', info_text=INFO)
    ratio_legend(rax, geometries)
    ax.legend(loc='upper left', frameon=False, fontsize=12, title=r'$\mu^-$ stops', title_fontsize=12)
    save(fig, out_dir, 'stops_r')


def plot_rz(stops, geometries, summary, out_dir, z_edges, r_edges):
    """Stop r against z for each geometry per 10^3 MuBeamCat events, on a common colour scale."""
    fig, axes = plt.subplots(len(geometries), 1, figsize=(11, 3.2 * len(geometries)), sharex=True,
                             gridspec_kw={'hspace': 0.08})
    axes = np.atleast_1d(axes)
    hists = [np.histogram2d(stops[g]['z'], np.hypot(stops[g]['x'], stops[g]['y']), [z_edges, r_edges])[0]
             * 1e3 / mubeam_events(summary[g]) for g in geometries]
    vmax = max(h.max() for h in hists)
    for ax, geometry, h in zip(axes, geometries, hists):
        mesh = ax.pcolormesh(z_edges, r_edges, np.where(h > 0, h, np.nan).T, cmap='viridis', vmin=0, vmax=vmax)
        ax.text(0.01, 0.92, GEOMETRY_LABELS.get(geometry, geometry), transform=ax.transAxes, ha='left', va='top',
                fontsize=12, color='white', bbox={'facecolor': 'black', 'alpha': 0.5, 'lw': 0})
        ax.set_ylabel('r (mm)', fontsize=14)
        ax.set_facecolor('#f0f0f0')
    mu2e_label(axes[0])
    axes[-1].set_xlabel('z, detector frame (mm)', loc='right', fontsize=14)
    fig.colorbar(mesh, ax=axes, label=r'Stops per bin per $10^3$ MuBeamCat', pad=0.01)
    save(fig, out_dir, 'stops_rz')


def ce_per_foil(ce, manifest, registry, foils, geometries, p_window, t_window):
    """Per foil: CE efficiency into the signal window, mean momentum response, and CE signal for the run.

    The CE are generated at the stops, so the CE generated in a foil is expected x (stops in the foil /
    all stops). The efficiency is the signal-window yield from the foil over that.
    """
    entries = {registry['datasets'][e['name']]['geometry']: e for e in manifest['datasets']
               if e.get('status') == 'ok' and e['process'] == 'CE'}
    out = {}
    for geometry in geometries:
        if geometry not in ce:
            continue
        entry, f = entries[geometry], foils[geometry].copy()
        stop_rate_rel = registry['datasets'][entry['name']].get('stop_rate_rel', 1.0)
        expected = float(entry['expected']) * stop_rate_rel
        weight = expected * float(entry['acceptance']) / generated_in_sample(entry['cutflow'])
        df = ce[geometry]
        df = df[(df['time'] > t_window[0]) & (df['time'] < t_window[1])]
        index = foil_index(df['pos_z_true'].to_numpy(), find_foils_from(f))
        dp = (df['momentum'] - df['mom_true']).to_numpy()
        in_window = ((df['momentum'] > p_window[0]) & (df['momentum'] < p_window[1])).to_numpy()
        generated = expected * f['stops'] / f['stops'].sum()
        n_window = np.array([np.sum(in_window & (index == i)) for i in f['foil']])
        f['ce_signal'] = n_window * weight
        f['ce_signal_err'] = np.sqrt(n_window) * weight
        f['efficiency'] = f['ce_signal'] / generated
        f['efficiency_err'] = f['ce_signal_err'] / generated
        f['dp_mean'] = [dp[index == i].mean() if np.any(index == i) else np.nan for i in f['foil']]
        f['dp_median'] = [np.median(dp[index == i]) if np.any(index == i) else np.nan for i in f['foil']]
        outside = np.mean(index == 0)
        if outside > 0.01:
            print(f'[{geometry}] {outside:.1%} of CE tracks start outside the foils found from the stops')
        out[geometry] = f
    return out


def find_foils_from(f):
    """Foil edges back from the per-foil table: midpoints between foil centres, with half a pitch at the ends."""
    z = f['z'].to_numpy()
    mid = 0.5 * (z[1:] + z[:-1])
    half = 0.5 * np.min(np.diff(z)) if len(z) > 1 else 5.0
    return np.column_stack([np.concatenate([[z[0] - half], mid]), np.concatenate([mid, [z[-1] + half]])])


def plot_ce_foils(per, geometries, info, out_dir):
    """Three panels against foil z: CE efficiency into the signal window, mean response, CE signal for the run."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True, gridspec_kw={'hspace': 0.07})
    for geometry, color, marker in zip(geometries, GEOMETRY_COLORS, GEOMETRY_MARKERS):
        if geometry not in per:
            continue
        f = per[geometry]
        label = f'{GEOMETRY_LABELS.get(geometry, geometry)}: {f["ce_signal"].sum():.3f} CE'
        axes[0].errorbar(f['z'], f['efficiency'] * 100, yerr=f['efficiency_err'] * 100, marker=marker, ms=4, lw=1, color=color, label=label)
        axes[1].plot(f['z'], f['dp_median'], marker=marker, ms=4, lw=1, color=color)
        axes[2].errorbar(f['z'], f['ce_signal'] * 1e3, yerr=f['ce_signal_err'] * 1e3, marker=marker, ms=4, lw=1, color=color)
    axes[0].set_ylabel('Signal-window\nefficiency (%)', fontsize=14)
    axes[1].set_ylabel('Median\n$p_{\\mathrm{reco}} - p_{\\mathrm{gen}}$ (MeV/c)', fontsize=14)
    axes[2].set_ylabel('CE in window\nper foil ($10^{-3}$)', fontsize=14)
    for ax in axes:
        ax.set_ylim(bottom=0 if ax is not axes[1] else None)
        ax.grid(alpha=0.3)
    axes[0].set_ylim(top=axes[0].get_ylim()[1] * 1.7)
    axes[0].legend(loc='upper left', frameon=False, fontsize=11, title=r'CE, $R_{\mu e} = 10^{-15}$; total in window', title_fontsize=11)
    mu2e_label(axes[0])
    axes[0].text(1.0, 1.03, info, transform=axes[0].transAxes, ha='right', va='bottom', fontsize=11)
    axes[-1].set_xlabel('Stop foil z, detector frame (mm)', loc='right', fontsize=14)
    save(fig, out_dir, 'ce_foils')


def main(args):
    stops_dir = Path(args.outdir).expanduser().resolve() / 'stops'
    summary = yaml.safe_load((stops_dir / 'stops_summary.yaml').read_text())
    stops = {g: pd.read_parquet(stops_dir / f'stops_{g}.parquet') for g in summary}
    geometries = [g for g in GEOMETRY_LABELS if g in stops] + [g for g in stops if g not in GEOMETRY_LABELS]
    out_dir = stops_dir / 'plots'
    out_dir.mkdir(parents=True, exist_ok=True)

    foils = per_foil(stops, summary, args.foil_gap)
    table = pd.concat([f.assign(geometry=g) for g, f in foils.items()], ignore_index=True)
    table.to_csv(out_dir / 'stops_foils.csv', index=False)
    print('Stops per MuBeamCat event:')
    for g in geometries:
        s = summary[g]
        print(f'  {g:16s} {stops_per_mubeam(s):.4g} ({s["stops"]} stops, {s["files"]} files), {len(foils[g])} foils, '
              f'ratio {stops_per_mubeam(s) / stops_per_mubeam(summary[geometries[0]]):.4f}')
    print(f'Wrote {out_dir / "stops_foils.csv"}')

    z_all = np.concatenate([stops[g]['z'] for g in geometries])
    z_edges = np.arange(np.floor(z_all.min()) - 10, np.ceil(z_all.max()) + 10 + 1e-9, args.z_bin)
    r_edges = np.arange(0, args.r_max + 1e-9, args.r_bin)
    plot_z(foils, geometries, summary, out_dir)
    plot_r(stops, geometries, summary, out_dir, r_edges)
    plot_rz(stops, geometries, summary, out_dir, z_edges, r_edges)

    if args.tag:
        tag_dir = Path(args.tag).expanduser().resolve()
        manifest = yaml.safe_load((tag_dir / 'manifest.yaml').read_text())
        registry = yaml.safe_load(Path(manifest['registry']).read_text())
        cuts = get_version_cuts(manifest.get('version') or 'run1a')
        p_window = tuple(args.p_window or (cuts['signal_region_p_lower'], cuts['signal_region_p_upper']))
        t_window = tuple(args.t_window or (cuts['signal_region_t_lower'], cuts['signal_region_t_upper']))
        ce = load_columns(manifest, registry, 'CE', ['momentum', 'time', 'mom_true', 'pos_z_true'])
        if ce:
            per = ce_per_foil(ce, manifest, registry, foils, geometries, p_window, t_window)
            run = registry['run']
            info = format_run_info(run['beam_power_kW'], run['realtime_days'], run['pot']) + '; stops scaled'
            plot_ce_foils(per, geometries, info, out_dir)
            pd.concat([f.assign(geometry=g) for g, f in per.items()], ignore_index=True).to_csv(out_dir / 'ce_foils.csv', index=False)
            print(f'Wrote {out_dir / "ce_foils.csv"}')
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Muon stops and CE signal per foil for the stopping-target geometries')
    parser.add_argument('--outdir', required=True, help='Study directory (dump_stops_st.py --outdir)')
    parser.add_argument('--tag', default=None, help='run_chain.py tag directory, for the CE per foil (needs pos_z_true in the parquet)')
    parser.add_argument('--p-window', type=float, nargs=2, default=None, metavar=('LO', 'HI'),
                        help='Signal momentum window for the CE per foil (default: signal_region_p of the cut set)')
    parser.add_argument('--t-window', type=float, nargs=2, default=None, metavar=('LO', 'HI'),
                        help='Signal time window for the CE per foil (default: signal_region_t of the cut set)')
    parser.add_argument('--foil-gap', type=float, default=3.0, help='z gap (mm) between stops that starts a new foil')
    parser.add_argument('--z-bin', type=float, default=2.0, help='z bin width (mm) of the r-z plot')
    parser.add_argument('--r-bin', type=float, default=2.0, help='Radius bin width (mm)')
    parser.add_argument('--r-max', type=float, default=120.0, help='Largest radius plotted (mm)')
    raise SystemExit(main(parser.parse_args()))
