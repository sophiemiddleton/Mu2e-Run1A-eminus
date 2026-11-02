"""CE and DIO after the analysis cuts for each stopping-target geometry, scaled to the run.

Reads <tag>/manifest.yaml from run_chain.py --skip-fit on the make_registry_st.py registry.
Each track is weighted as in run_chain.py and make_page.py:

    weight = expected * acceptance / N_generated_in_sample

so a histogram is in events for the run in the registry "run" section (default 28 days at
3.84 kW, CE at R_mue = 1e-15).

Writes to <tag>/plots/:
    mom_<geometry>.{png,pdf}   momentum in the signal time window, DIO and nominal cosmics stacked, CE line
    time_<geometry>.{png,pdf}  time in the signal momentum window
    mom_compare.{png,pdf}      the three geometries overlaid, ratio to the first
    sensitivity_st.csv         per geometry: stops, CE and DIO in the signal window (scaled by the stop rate),
                               CE reco and selection efficiency, SES and the expected R_mue limit
    sensitivity_st_opt.csv     the same in each geometry's optimized momentum window (--optimize-p)
    mom_/time_<geometry>_opt   the per-geometry plots in that window
    window_scan_plo, window_scan_2d, window_scan.csv
                               limit against the window edges, per geometry (--optimize-p)
    window_yields              CE and backgrounds against the lower edge, at the best shared upper edge
    ce_eff_compare.{png,pdf,csv}
                               CE efficiency per generated event after each cut and in the windows

    python plot_st.py --tag /exp/mu2e/data/users/$USER/stoppingtarget/run1a
"""

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

PYFITTER_DIR = Path(os.environ.get('PYFITTER_DIR', Path(__file__).resolve().parents[1] / 'pyfitter'))
sys.path.insert(0, str(PYFITTER_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'sensitivity' / 'cutcount'))

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from config import get_version_cuts  # noqa: E402
from style import COLORS, apply_root_text, format_run_info, line_color  # noqa: E402

GEOMETRY_LABELS = {
    'nominal': 'Nominal (37 foils)',
    'st37plus6up6dn': '37 foils + 6 up + 6 down',
    'st49hole1to5': '49 foils, hole r 1-5 mm',
    'st49hole9to41': '49 foils, hole r 9-41 mm',
}
GEOMETRY_STYLES = ['-', '--', ':', '-.']
GEOMETRY_COLORS = ['#000000', '#e66100', '#5d3a9b', '#1a85ff']
GEOMETRY_MARKERS = ['o', 's', '^', 'D']


def generated_in_sample(cutflow_path):
    """Events before any cut ('No cuts' row): the denominator of the run_chain cut efficiency."""
    with open(cutflow_path, newline='') as f:
        return float(next(csv.DictReader(f))['events_passing'])


def load(manifest, registry, scale_stops):
    """{geometry: {process: (momentum, time, weight per track)}}."""
    out = {}
    for entry in manifest['datasets']:
        if entry.get('status') != 'ok':
            print(f'[{entry["name"]}] skipped, status {entry.get("status")}')
            continue
        reg = registry['datasets'][entry['name']]
        weight = float(entry['expected']) * float(entry['acceptance']) / generated_in_sample(entry['cutflow'])
        if scale_stops:
            weight *= reg['stop_rate_rel']
        df = pd.read_parquet(entry['parquet'], columns=['momentum', 'time'])
        out.setdefault(reg['geometry'], {})[entry['process']] = (
            df['momentum'].to_numpy(np.float64), df['time'].to_numpy(np.float64), weight)
    return out


def load_cosmic(tag_dir):
    """(momentum, time, weight per track) of the main-analysis cosmics; weight = acceptance (on-spill time / livetime)."""
    with open(Path(tag_dir) / 'manifest.yaml') as f:
        manifest = yaml.safe_load(f)
    entries = [e for e in manifest['datasets'] if e.get('status') == 'ok' and e.get('process') == 'Cosmic']
    if len(entries) != 1:
        raise ValueError(f'Expected one Cosmic dataset in {tag_dir}/manifest.yaml, found {len(entries)}')
    df = pd.read_parquet(entries[0]['parquet'], columns=['momentum', 'time'])
    return df['momentum'].to_numpy(np.float64), df['time'].to_numpy(np.float64), float(entries[0]['acceptance'])


def load_columns(manifest, registry, process, columns):
    """{geometry: DataFrame of the post-cut parquet columns} for one process; geometries without them are left out."""
    out = {}
    for entry in manifest['datasets']:
        reg = registry['datasets'].get(entry['name'], {})
        if entry.get('status') != 'ok' or entry['process'] != process:
            continue
        try:
            out[reg['geometry']] = pd.read_parquet(entry['parquet'], columns=columns).dropna()
        except (KeyError, ValueError) as e:
            print(f'[{entry["name"]}] parquet lacks some of {columns}: rerun run_chain.py with the current process.py')
    return out


def hist(values, weight, edges):
    """Weighted counts and their MC statistical errors (constant weight per sample)."""
    n = np.histogram(values, edges)[0]
    return n * weight, np.sqrt(n) * weight


def save(fig, out_dir, name):
    for ext in ('png', 'pdf'):
        fig.savefig(out_dir / f'{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'Wrote {out_dir / name}.png')


def plot_geometry(geometry, data, variable, edges, window, other_window, info, out_dir, suffix=''):
    """DIO and cosmics stacked and CE line for one geometry, in the other variable's signal window."""
    index = 0 if variable == 'mom' else 1
    fig, ax = plt.subplots(figsize=(10, 6))
    centres = 0.5 * (edges[1:] + edges[:-1])
    stack = np.zeros(len(edges) - 1)
    for process, label in (('DIO', 'DIO'), ('Cosmic', 'Cosmic (nominal target)'), ('CE', r'CE, $R_{\mu e} = 10^{-15}$')):
        if process not in data:
            continue
        values = data[process][index]
        other = data[process][1 - index]
        keep = (other > other_window[0]) & (other < other_window[1])
        y, err = hist(values[keep], data[process][2], edges)
        if process == 'CE':
            ax.stairs(y, edges, color=COLORS['Signal'], lw=1.5, label=label, zorder=5)
            ax.errorbar(centres, y, yerr=err, fmt='none', ecolor=COLORS['Signal'], lw=0.8, zorder=6)
            continue
        ax.stairs(stack + y, edges, baseline=stack.copy(), fill=True, color=COLORS[process], edgecolor='black', lw=0.5, label=label)
        stack += y
        ax.errorbar(centres, stack, yerr=err, fmt='none', ecolor=line_color(COLORS[process]), lw=0.8, zorder=6)
    for edge in window:
        ax.axvline(edge, color='gray', ls='--', lw=1)
    ax.set_yscale('log')
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(bottom=1e-4)
    width = edges[1] - edges[0]
    if variable == 'mom':
        xlabel, ylabel = 'Momentum (MeV/c)', f'Events / {width:.2g} MeV/c'
        title = f'{GEOMETRY_LABELS.get(geometry, geometry)}; {other_window[0]:.0f} < t < {other_window[1]:.0f} ns'
    else:
        xlabel, ylabel = 'Time (ns)', f'Events / {width:.0f} ns'
        title = f'{GEOMETRY_LABELS.get(geometry, geometry)}; {other_window[0]:.2f} < p < {other_window[1]:.2f} MeV/c'
    apply_root_text(fig, ax, None, xlabel, ylabel, None, info_text=info)
    ax.legend(loc='upper right', frameon=False, title=title, title_fontsize=12, fontsize=14)
    ax.set_ylim(top=ax.get_ylim()[1] * 100)  # room for the legend above the histograms
    save(fig, out_dir, f'{variable}_{geometry}{suffix}')


def plot_compare(data, geometries, edges, window, t_window, info, out_dir):
    """CE and DIO momentum for every geometry, with the ratio to the first geometry."""
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    ratios = []
    for process, color in (('DIO', line_color(COLORS['DIO'])), ('CE', COLORS['Signal'])):
        ref = None
        for geometry, ls in zip(geometries, GEOMETRY_STYLES):
            if process not in data[geometry]:
                continue
            p, t, w = data[geometry][process]
            keep = (t > t_window[0]) & (t < t_window[1])
            y, err = hist(p[keep], w, edges)
            ax.stairs(y, edges, color=color, ls=ls, lw=1.5, label=f'{process}: {GEOMETRY_LABELS.get(geometry, geometry)}')
            if ref is None:
                ref = (y, err)
                continue
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = y / ref[0]
                ratio_err = ratio * np.sqrt((err / y) ** 2 + (ref[1] / ref[0]) ** 2)
            centres = 0.5 * (edges[1:] + edges[:-1])
            ratios.append(ratio[np.isfinite(ratio) & (ratio > 0)])
            rax.errorbar(centres, ratio, yerr=ratio_err, fmt='o' if process == 'CE' else 's', ms=3, color=color,
                         mfc='none' if ls == ':' else color, lw=0.8, label=f'{process}: {GEOMETRY_LABELS.get(geometry, geometry)}')
    cosmic = data[geometries[0]].get('Cosmic')
    if cosmic is not None:
        p, t, w = cosmic
        keep = (t > t_window[0]) & (t < t_window[1])
        ax.stairs(hist(p[keep], w, edges)[0], edges, color=line_color(COLORS['Cosmic']), lw=1.5, label='Cosmic: nominal target')
    for a in (ax, rax):
        for edge in window:
            a.axvline(edge, color='gray', ls='--', lw=1)
    rax.axhline(1.0, color='black', lw=0.8)
    # Log ratio axis wide enough for every point, symmetric about 1
    values = np.concatenate(ratios) if ratios else np.array([1.0])
    span = max(2.0, values.max(), 1.0 / values.min()) * 1.3
    rax.set_yscale('log')
    rax.set_ylim(1.0 / span, span)
    ax.set_yscale('log')
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(bottom=1e-4)
    apply_root_text(fig, ax, rax, 'Momentum (MeV/c)', f'Events / {edges[1] - edges[0]:.2g} MeV/c',
                    f'Ratio to {GEOMETRY_LABELS.get(geometries[0], geometries[0])}', info_text=info)
    rax.yaxis.label.set_size(11)
    rax.legend(loc='upper left', frameon=False, fontsize=8, ncol=4, columnspacing=1.0, handletextpad=0.3)
    ax.legend(loc='upper left', frameon=False, fontsize=12, ncol=2,
              title=f'{t_window[0]:.0f} < t < {t_window[1]:.0f} ns', title_fontsize=12)
    ax.set_ylim(top=ax.get_ylim()[1] * 1e3)
    save(fig, out_dir, 'mom_compare')


def fwhm(counts, edges):
    """Full width at half maximum of a histogram, from linear interpolation of the half-maximum crossings."""
    centres = 0.5 * (edges[1:] + edges[:-1])
    peak = np.argmax(counts)
    half = counts[peak] / 2
    lo = np.nonzero(counts[:peak] < half)[0][-1]
    hi = peak + np.nonzero(counts[peak:] < half)[0][0]
    left = np.interp(half, counts[lo:lo + 2], centres[lo:lo + 2])
    right = np.interp(half, counts[hi - 1:hi + 1][::-1], centres[hi - 1:hi + 1][::-1])
    return right - left


def plot_resolution(ce, geometries, p_window, t_window, out_dir, edges):
    """CE momentum response (reco - generated momentum) per geometry, unit area, ratio to the first geometry.

    mom_true is the production momentum, so the response includes the energy loss in the target and
    upstream material as well as the tracker resolution. The signal window is drawn at its offset from
    the generated momentum.
    """
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    centres = 0.5 * (edges[1:] + edges[:-1])
    rows, ref, p_gen = [], None, None
    for geometry, ls in zip(geometries, GEOMETRY_STYLES):
        df = ce[geometry]
        df = df[(df['time'] > t_window[0]) & (df['time'] < t_window[1])]
        dp = (df['momentum'] - df['mom_true']).to_numpy()
        p_gen = float(np.median(df['mom_true']))
        n = np.histogram(dp, edges)[0].astype(float)
        density, err = n / n.sum(), np.sqrt(n) / n.sum()
        ax.stairs(density, edges, color=COLORS['Signal'], ls=ls, lw=1.5, label=GEOMETRY_LABELS.get(geometry, geometry))
        core = dp[(dp > -1.5) & (dp < 1.0)]
        rows.append({'geometry': geometry, 'tracks': len(dp), 'peak': centres[np.argmax(n)], 'median': np.median(dp),
                     'fwhm': fwhm(n, edges), 'core_rms': core.std(),
                     'frac_below_-1': np.mean(dp < -1.0), 'frac_in_window': np.mean((dp > p_window[0] - p_gen) & (dp < p_window[1] - p_gen))})
        if ref is None:
            ref = (density, err)
            continue
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = density / ref[0]
            ratio_err = ratio * np.sqrt((err / density) ** 2 + (ref[1] / ref[0]) ** 2)
        rax.errorbar(centres, ratio, yerr=ratio_err, fmt='o', ms=3, color=COLORS['Signal'], mfc='none' if ls == ':' else COLORS['Signal'], lw=0.8,
                     label=GEOMETRY_LABELS.get(geometry, geometry))
    for a in (ax, rax):
        for edge in p_window:
            a.axvline(edge - p_gen, color='gray', ls='--', lw=1)
    rax.axhline(1.0, color='black', lw=0.8)
    rax.set_ylim(0.5, 1.5)
    ax.set_yscale('log')
    ax.set_xlim(edges[0], edges[-1])
    apply_root_text(fig, ax, rax, r'$p_{\mathrm{reco}} - p_{\mathrm{gen}}$ (MeV/c)', f'Fraction / {edges[1] - edges[0]:.2g} MeV/c',
                    'Ratio', info_text='Unit area')
    rax.legend(loc='upper left', frameon=False, fontsize=9, ncol=2,
               title=f'Ratio to {GEOMETRY_LABELS.get(geometries[0], geometries[0])}', title_fontsize=9)
    ax.legend(loc='upper left', frameon=False, fontsize=13,
              title=f'CE, $p_{{\\mathrm{{gen}}}}$ = {p_gen:.2f} MeV/c; {t_window[0]:.0f} < t < {t_window[1]:.0f} ns', title_fontsize=12)
    ax.set_ylim(top=ax.get_ylim()[1] * 5)
    save(fig, out_dir, 'resolution_compare')

    table = pd.DataFrame(rows)
    table.to_csv(out_dir / 'resolution_st.csv', index=False)
    print('\nCE momentum response (MeV/c)')
    with pd.option_context('display.width', 200, 'display.float_format', '{:.4g}'.format):
        print(table.to_string(index=False))
    print(f'Wrote {out_dir / "resolution_st.csv"}')


def main_backgrounds(main_tag, processes):
    """{process: (momentum, time, weight per track)} of the main-analysis tag's other backgrounds (RPC, RMC),
    weighted as run_cutcount.count_window: expected * acceptance / N_generated_in_sample."""
    from run_cutcount import canonical_process_name

    with open(Path(main_tag) / 'manifest.yaml') as f:
        manifest = yaml.safe_load(f)
    out = {}
    for entry in manifest['datasets']:
        process = canonical_process_name(entry.get('process'))
        if entry.get('status') != 'ok' or process not in processes:
            continue
        df = pd.read_parquet(entry['parquet'], columns=['momentum', 'time']).astype('float64')
        weight = float(entry['expected']) * float(entry['acceptance']) / generated_in_sample(entry['cutflow'])
        out.setdefault(process, []).append((df['momentum'].to_numpy(), df['time'].to_numpy(), weight))
    return out


def in_window(p, t, p_window, t_window):
    return (p > p_window[0]) & (p < p_window[1]) & (t > t_window[0]) & (t < t_window[1])


class WindowModel:
    """Yields and the expected limit of one geometry in any momentum window, for the table and the window scan.

    CE and DIO are the geometry's samples (scaled by stop_rate_rel), RPC and RMC the main-analysis tag's.
    The cosmics are counted in the window, or with cosmic_flat taken as flat in momentum: the rate per
    MeV/c over 100-110 MeV/c in the time window, so a scanned window does not follow MC fluctuations.
    """

    def __init__(self, card_path, main_tag, t_window, cosmic, cosmic_flat=False):
        from counting_stats import load_card

        self.card, self.card_yields, self.signal = load_card(card_path)
        self.cl = float(self.card.get('confidence_level', 0.90))
        self.others = [k for k in self.card_yields if k not in (self.signal, 'DIO', 'Cosmic')]
        self.backgrounds = main_backgrounds(main_tag, self.others)
        self.t_window = t_window
        self.cosmic = cosmic
        p, t, w = cosmic
        n_flat = int(in_window(p, t, (100.0, 110.0), t_window).sum())
        self.cosmic_flat = (n_flat * w / 10.0, np.sqrt(n_flat) * w / 10.0) if cosmic_flat else None

    def yields(self, processes, p_window):
        """{process: (yield, error, MC tracks)} in the window."""
        out = {}
        for process in ('CE', 'DIO'):
            if process in processes:
                p, t, w = processes[process]
                n = int(in_window(p, t, p_window, self.t_window).sum())
                out[process] = (n * w, np.sqrt(n) * w, n)
        p, t, w = self.cosmic
        n = int(in_window(p, t, p_window, self.t_window).sum())
        if self.cosmic_flat:
            width = p_window[1] - p_window[0]
            out['Cosmic'] = (self.cosmic_flat[0] * width, self.cosmic_flat[1] * width, n)
        else:
            out['Cosmic'] = (n * w, np.sqrt(n) * w, n)
        for process in self.others:
            parts = [(int(in_window(p, t, p_window, self.t_window).sum()), w) for p, t, w in self.backgrounds.get(process, [])]
            out[process] = (sum(n * w for n, w in parts), np.sqrt(sum(n * w * w for n, w in parts)), sum(n for n, _ in parts))
        return out

    def limit_events(self, yields):
        """Median expected limit in signal events (Combine AsymptoticLimits, counting_stats.py)."""
        from counting_stats import CountingModel, asimov_expected_limits

        counts = dict(self.card_yields, **{k: v[0] for k, v in yields.items()})
        counts.setdefault('DIO', 0.0)
        counts[self.signal] = yields['CE'][0]
        mu = asimov_expected_limits(CountingModel(counts, self.signal, self.card.get('systematics')), self.cl)[0]
        return mu * yields['CE'][0]

    def scan(self, processes, lows, highs):
        """Limit in signal events per CE event (R_mue / 1e-15 at R_mue = 1e-15) for every window, as a
        DataFrame with p_lo, p_hi, the yields and limit_per_ce."""
        rows = []
        for lo in lows:
            for hi in highs:
                if hi <= lo:
                    continue
                y = self.yields(processes, (lo, hi))
                row = {'p_lo': round(lo, 3), 'p_hi': round(hi, 3), **{f'N_{k}': v[0] for k, v in y.items()}}
                row['limit_per_ce'] = self.limit_events(y) / y['CE'][0] if y['CE'][0] > 0 else np.inf
                rows.append(row)
        return pd.DataFrame(rows)


def plot_ce_efficiency(manifest, registry, geometries, windows, model, data, out_dir):
    """CE efficiency per generated event after each cut (acceptance x cut flow), then in the signal windows,
    for each geometry, with the ratio to the first geometry. windows: {label: {geometry: (lo, hi)}}."""
    curves = {}
    for entry in manifest['datasets']:
        reg = registry['datasets'][entry['name']]
        if entry.get('status') != 'ok' or entry['process'] != 'CE' or reg['geometry'] not in geometries:
            continue
        cutflow = pd.read_csv(entry['cutflow']).set_index('name')['events_passing']
        acceptance = float(entry['acceptance'])
        steps = {'mcs acceptance': acceptance}
        steps.update({name: acceptance * n / cutflow['No cuts'] for name, n in cutflow.items() if name != 'No cuts'})
        for label, per_geometry in windows.items():
            window = per_geometry[reg['geometry']]
            n = model.yields(data[reg['geometry']], window)['CE'][2]
            steps[f'{label} window'] = acceptance * n / cutflow['No cuts']
        curves[reg['geometry']] = steps
    geometries = [g for g in geometries if g in curves]
    names = list(curves[geometries[0]])
    x = np.arange(len(names))
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(12, 9), sharex=True, gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.05})
    ref = np.array([curves[geometries[0]].get(n, np.nan) for n in names])
    for geometry, color, marker in zip(geometries, GEOMETRY_COLORS, GEOMETRY_MARKERS):
        y = np.array([curves[geometry].get(n, np.nan) for n in names])
        label = GEOMETRY_LABELS.get(geometry, geometry)
        ax.plot(x, y, marker=marker, color=color, lw=1, ms=5, label=label)
        if geometry != geometries[0]:
            rax.plot(x, y / ref, marker=marker, color=color, lw=1, ms=5, label=label)
    rax.axhline(1.0, color='black', lw=0.8)
    ax.set_ylim(0, 0.45)
    ax.set_ylabel('CE efficiency per generated event')
    rax.set_ylabel(f'Ratio to {GEOMETRY_LABELS.get(geometries[0], geometries[0])}', fontsize=11)
    rax.set_xticks(x)
    rax.set_xticklabels(names, rotation=60, ha='right', fontsize=10)
    ax.grid(axis='x', color='0.9')
    rax.grid(axis='x', color='0.9')
    ax.legend(loc='upper right', frameon=False, fontsize=12)
    ax.text(0.01, 0.97, 'Mu2e $\\it{Simulation}$', transform=ax.transAxes, va='top', fontsize=16, fontweight='bold')
    save(fig, out_dir, 'ce_eff_compare')
    table = pd.DataFrame({g: curves[g] for g in geometries})
    table.index.name = 'step'
    table.to_csv(out_dir / 'ce_eff_compare.csv')


def plot_window_scan(scans, best, info, out_dir):
    """Limit against the window edges: each geometry's limit (x 1e-15) against p_lo at its best p_hi, and
    the 2D map over (p_lo, p_hi) per geometry with its optimum marked."""
    geometries = list(scans)
    fig, ax = plt.subplots(figsize=(10, 6))
    for geometry, color, marker in zip(geometries, GEOMETRY_COLORS, GEOMETRY_MARKERS):
        df = scans[geometry]
        hi = best[geometry][1]
        d = df[np.isclose(df['p_hi'], hi)].sort_values('p_lo')
        ax.plot(d['p_lo'], d['limit_per_ce'], marker=marker, color=color, ms=4, lw=1.2,
                label=f'{GEOMETRY_LABELS.get(geometry, geometry)} ($p_{{hi}}$ = {hi:.1f})')
    ax.set_xlabel('Lower momentum edge (MeV/c)')
    ax.set_ylabel(r'Expected 90% CL $R_{\mu e}$ ($10^{-15}$)')
    apply_root_text(fig, ax, None, ax.get_xlabel(), ax.get_ylabel(), None, info_text=info)
    ax.legend(loc='upper left', frameon=False, fontsize=12)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.4 * (hi - lo))
    save(fig, out_dir, 'window_scan_plo')

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True, constrained_layout=True)
    vmin = min(scans[g]['limit_per_ce'].min() for g in geometries)
    for ax, geometry in zip(axes.flat, geometries):
        df = scans[geometry]
        grid = df.pivot(index='p_hi', columns='p_lo', values='limit_per_ce')
        mesh = ax.pcolormesh(grid.columns, grid.index, grid.values, shading='nearest', cmap='viridis_r',
                             vmin=vmin, vmax=vmin * 1.3)
        ax.plot(*best[geometry], marker='*', color='red', ms=14, mec='white')
        ax.set_title(f'{GEOMETRY_LABELS.get(geometry, geometry)}: {df["limit_per_ce"].min():.2f} at '
                     f'{best[geometry][0]:.1f}-{best[geometry][1]:.1f}', fontsize=12)
    for ax in axes[1]:
        ax.set_xlabel('Lower edge (MeV/c)')
    for ax in axes[:, 0]:
        ax.set_ylabel('Upper edge (MeV/c)')
    fig.colorbar(mesh, ax=axes, label=r'Expected 90% CL $R_{\mu e}$ ($10^{-15}$)')
    save(fig, out_dir, 'window_scan_2d')

    # Best single window for all the geometries: lowest mean limit relative to each geometry's own best
    merged = None
    for geometry in geometries:
        d = scans[geometry][['p_lo', 'p_hi', 'limit_per_ce']].rename(columns={'limit_per_ce': geometry})
        d[geometry] = d[geometry] / d[geometry].min()
        merged = d if merged is None else merged.merge(d, on=['p_lo', 'p_hi'])
    merged['mean'] = merged[geometries].mean(axis=1)
    shared = merged.loc[merged['mean'].idxmin()]
    print(f'Best shared window {shared["p_lo"]:.1f} < p < {shared["p_hi"]:.1f} MeV/c; limit / own best: ' +
          ', '.join(f'{g} {shared[g]:.3f}' for g in geometries))
    out = pd.concat([scans[g].assign(geometry=g) for g in geometries], ignore_index=True)
    out.to_csv(out_dir / 'window_scan.csv', index=False)
    return float(shared['p_lo']), float(shared['p_hi'])


def plot_window_yields(scans, p_hi, info, out_dir):
    """CE and backgrounds in the window against its lower edge, at a common upper edge."""
    geometries = list(scans)
    fig, ax = plt.subplots(figsize=(10, 6))
    for geometry, ls in zip(geometries, GEOMETRY_STYLES):
        d = scans[geometry][np.isclose(scans[geometry]['p_hi'], p_hi)].sort_values('p_lo')
        label = GEOMETRY_LABELS.get(geometry, geometry)
        ax.plot(d['p_lo'], d['N_CE'], color=COLORS['Signal'], ls=ls, lw=1.5, label=f'CE: {label}')
        ax.plot(d['p_lo'], d['N_DIO'], color=line_color(COLORS['DIO']), ls=ls, lw=1.5, label=f'DIO: {label}')
    d = scans[geometries[0]][np.isclose(scans[geometries[0]]['p_hi'], p_hi)].sort_values('p_lo')
    ax.plot(d['p_lo'], d['N_Cosmic'], color=line_color(COLORS['Cosmic']), lw=1.5, label='Cosmic (flat), all')
    if 'N_RPC' in d:
        ax.plot(d['p_lo'], d['N_RPC'], color=line_color(COLORS['RPC']), lw=1.5, label='RPC, all')
    ax.set_yscale('log')
    ax.set_xlabel('Lower momentum edge (MeV/c)')
    ax.set_ylabel(f'Events in the window ($p$ < {p_hi:.1f} MeV/c)')
    apply_root_text(fig, ax, None, ax.get_xlabel(), ax.get_ylabel(), None, info_text=info)
    ax.legend(loc='upper center', frameon=False, fontsize=10, ncol=3, title=r'CE at $R_{\mu e} = 10^{-15}$', title_fontsize=11)
    ax.set_ylim(top=ax.get_ylim()[1] * 30)
    save(fig, out_dir, 'window_yields')


def sensitivity_table(manifest, registry, data, p_windows, t_window, model, out_dir, name='sensitivity_st'):
    """Stops, signal and backgrounds in each geometry's signal window, CE efficiencies and R_mue per geometry.

    data is load(..., scale_stops=True): the yields are always scaled by each geometry's stop rate
    (stop_rate_rel). The CE efficiency splits into reco (a reconstructed track per generated CE:
    acceptance * has_a_track) and selection (in the signal window after the cuts, per reconstructed track).

    R_mue is the Combine AsymptoticLimits expected limit (counting_stats.py) on the count card with CE,
    DIO, cosmics, RPC and RMC replaced by their yields in the window: the card's processes and lnN
    systematics.
    """
    run = registry['run']
    rows = {}
    for entry in manifest['datasets']:
        if entry.get('status') != 'ok' or entry['process'] != 'CE':
            continue
        reg = registry['datasets'][entry['name']]
        geometry = reg['geometry']
        cutflow = pd.read_csv(entry['cutflow']).set_index('name')['events_passing']
        p_window = p_windows[geometry]
        yields = model.yields(data[geometry], p_window)
        row = {'geometry': geometry, 'p_lo': p_window[0], 'p_hi': p_window[1],
               'N_stops': run['stopped_muons'] * reg['stop_rate_rel']}
        for process, (y, err, n) in yields.items():
            label = process.lower() if process == 'Cosmic' else process
            row[f'N_{label}'], row[f'N_{label}_err'], row[f'{process}_tracks'] = y, err, n
        row['CE_reco_eff'] = float(entry['acceptance']) * cutflow['has_a_track'] / cutflow['No cuts']
        row['CE_sel_eff'] = yields['CE'][2] / cutflow['has_a_track']
        row['SES'] = 1.0 / (row['N_stops'] * run['captures_per_stopped_muon'] * row['CE_reco_eff'] * row['CE_sel_eff'])
        row['limit_events'] = model.limit_events(yields)
        row[f'Rmue_{model.cl:.0%}CL'] = row['limit_events'] * row['SES']
        rows[geometry] = row
    others = [c for p in model.others for c in (f'N_{p}', f'N_{p}_err')]
    columns = (['geometry', 'p_lo', 'p_hi', 'N_stops', 'N_CE', 'N_CE_err', 'N_DIO', 'N_DIO_err', 'N_cosmic', 'N_cosmic_err']
               + others + ['CE_reco_eff', 'CE_sel_eff', 'SES', 'limit_events', f'Rmue_{model.cl:.0%}CL', 'CE_tracks',
                           'DIO_tracks', 'Cosmic_tracks'])
    geometries = [g for g in GEOMETRY_LABELS if g in rows] + [g for g in rows if g not in GEOMETRY_LABELS]
    table = pd.DataFrame([rows[g] for g in geometries]).reindex(columns=columns)
    table.to_csv(out_dir / f'{name}.csv', index=False)
    cosmic = 'flat in momentum' if model.cosmic_flat else 'counted in the window'
    print(f'\n{name}: {t_window[0]:.0f} < t < {t_window[1]:.0f} ns, cosmics {cosmic}; '
          f'{run["realtime_days"]} days at {run["beam_power_kW"]} kW, CE at R_mue = {run["rmue"]:g}')
    with pd.option_context('display.width', 300, 'display.max_columns', 40, 'display.float_format', '{:.4g}'.format):
        print(table.to_string(index=False))
    print(f'Wrote {out_dir / name}.csv')


def main(args):
    tag_dir = Path(args.tag).expanduser().resolve()
    with open(tag_dir / 'manifest.yaml') as f:
        manifest = yaml.safe_load(f)
    with open(args.registry or manifest['registry']) as f:
        registry = yaml.safe_load(f)
    cuts = get_version_cuts(args.version or manifest.get('version') or 'run1a')
    p_window = tuple(args.p_window or (cuts['signal_region_p_lower'], cuts['signal_region_p_upper']))
    t_window = tuple(args.t_window or (cuts['signal_region_t_lower'], cuts['signal_region_t_upper']))

    data = load(manifest, registry, args.scale_stops)
    cosmic = load_cosmic(Path(args.main_tag).expanduser())
    for processes in data.values():
        processes['Cosmic'] = cosmic
    geometries = [g for g in GEOMETRY_LABELS if g in data] + [g for g in data if g not in GEOMETRY_LABELS]
    run = registry['run']
    info = format_run_info(run['beam_power_kW'], run['realtime_days'], run['pot'])
    if args.scale_stops:
        info += '; stops scaled'

    out_dir = tag_dir / 'plots'
    out_dir.mkdir(parents=True, exist_ok=True)
    mom_lo, mom_hi = manifest['fitrange_low'][0], manifest['fitrange_hi'][0]
    mom_edges = np.arange(mom_lo, mom_hi + 1e-9, args.mom_bin)
    time_edges = np.arange(manifest['fitrange_low'][1], manifest['fitrange_hi'][1] + 1e-9, args.time_bin)
    for geometry in geometries:
        plot_geometry(geometry, data[geometry], 'mom', mom_edges, p_window, t_window, info, out_dir)
        plot_geometry(geometry, data[geometry], 'time', time_edges, t_window, p_window, info, out_dir)
    plot_compare(data, geometries, mom_edges, p_window, t_window, info, out_dir)

    card = Path(args.card or Path(args.main_tag).expanduser() / 'cutcount' / 'count_card.yaml').expanduser()
    scaled = load(manifest, registry, True) if not args.scale_stops else data
    model = WindowModel(card, Path(args.main_tag).expanduser(), t_window, cosmic)
    sensitivity_table(manifest, registry, scaled, {g: p_window for g in geometries}, t_window, model, out_dir)
    if args.optimize_p:
        # Each geometry's momentum window with the lowest expected limit, cosmics flat in momentum
        flat = WindowModel(card, Path(args.main_tag).expanduser(), t_window, cosmic, cosmic_flat=True)
        lows = np.arange(args.optimize_p[0], args.optimize_p[1] + 1e-9, args.optimize_step)
        highs = np.arange(args.optimize_p[2], args.optimize_p[3] + 1e-9, args.optimize_step)
        scans = {g: flat.scan(scaled[g], lows, highs) for g in geometries if 'CE' in scaled[g]}
        windows = {g: tuple(df.loc[df['limit_per_ce'].idxmin(), ['p_lo', 'p_hi']]) for g, df in scans.items()}
        for g, w in windows.items():
            print(f'[{g}] optimized window {w[0]:.2f} < p < {w[1]:.2f} MeV/c')
        shared = plot_window_scan(scans, windows, info, out_dir)
        plot_window_yields(scans, round(shared[1], 3), info, out_dir)
        plot_ce_efficiency(manifest, registry, geometries, {'fixed': {g: p_window for g in geometries}, 'own': windows},
                           flat, scaled, out_dir)
        sensitivity_table(manifest, registry, scaled, windows, t_window, flat, out_dir, name='sensitivity_st_opt')
        for geometry, window in windows.items():
            plot_geometry(geometry, data[geometry], 'mom', mom_edges, window, t_window, info + '; optimized window',
                          out_dir, suffix='_opt')
            plot_geometry(geometry, data[geometry], 'time', time_edges, t_window, window, info + '; optimized window',
                          out_dir, suffix='_opt')
    ce = load_columns(manifest, registry, 'CE', ['momentum', 'time', 'mom_true'])
    if ce:
        plot_resolution(ce, [g for g in geometries if g in ce], p_window, t_window, out_dir,
                        np.arange(args.res_range[0], args.res_range[1] + 1e-9, args.res_bin))
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='CE and DIO per stopping-target geometry, scaled to the run')
    parser.add_argument('--tag', required=True, help='run_chain.py tag directory made from the make_registry_st.py registry')
    parser.add_argument('--registry', default=None, help='Registry (default: the one recorded in the manifest)')
    parser.add_argument('--version', default=None, help='Cut set for the signal window (default: the manifest version)')
    parser.add_argument('--mom-bin', type=float, default=0.1, help='Momentum bin width in MeV/c')
    parser.add_argument('--time-bin', type=float, default=25.0, help='Time bin width in ns')
    parser.add_argument('--res-bin', type=float, default=0.05, help='CE momentum response bin width in MeV/c')
    parser.add_argument('--res-range', type=float, nargs=2, default=(-4.0, 1.0), help='CE momentum response range in MeV/c')
    parser.add_argument('--main-tag', default=str(PYFITTER_DIR / 'au_run1a'),
                        help='Main-analysis tag: cosmics from its manifest, systematics and other yields from cutcount/count_card.yaml')
    parser.add_argument('--p-window', type=float, nargs=2, default=None, metavar=('LO', 'HI'),
                        help='Signal momentum window in MeV/c (default: signal_region_p of the cut set)')
    parser.add_argument('--t-window', type=float, nargs=2, default=None, metavar=('LO', 'HI'),
                        help='Signal time window in ns (default: signal_region_t of the cut set)')
    parser.add_argument('--optimize-p', type=float, nargs=4, default=None, metavar=('LO_MIN', 'LO_MAX', 'HI_MIN', 'HI_MAX'),
                        help='Also optimize each geometry\'s momentum window over these lower and upper edge ranges (MeV/c): '
                             'sensitivity_st_opt.csv and mom_/time_<geometry>_opt plots')
    parser.add_argument('--optimize-step', type=float, default=0.1, help='Step of the --optimize-p scan in MeV/c')
    parser.add_argument('--card', default=None,
                        help='Count card for the systematics and the other yields (default: <main-tag>/cutcount/count_card.yaml)')
    parser.add_argument('--scale-stops', action='store_true',
                        help='Also scale the plots by each geometry\'s muon-stop rate relative to the nominal target (the table is always scaled)')
    raise SystemExit(main(parser.parse_args()))
