"""Build a self-contained web page with the momentum and time distributions and the
cut-and-count and shape-fit sensitivities of one tag.

Inputs (all from the tag made by run_chain.py and cutcount/run_cutcount.py):
    <tag>/manifest.yaml                     datasets, acceptances, expected yields, parquet paths
    <tag>/cutcount/count_card.yaml          counting window
    <tag>/cutcount/count_limits_asymptotic.yaml   cut-and-count limit
    <tag>/limits/card_{2d,1d}_cls_limits.npz      shape-fit limit
Optional:
    --runtime-summary  studies/runtime/summary.csv of pyfitter's run_studies.py: limit vs run time
    --beam LABEL=TAG   shape-fit limit of other tags (e.g. 1.5 kW, 1BB, 2BB) for the beam-mode comparison

Each track is weighted to the run exactly as in the counting (run_cutcount.py):
    MC:     expected * acceptance / N_generated_in_sample
    Cosmic: acceptance (on-spill time / cosmic sample livetime)

The page has no external dependencies, so it opens from disk or any web server.
The two projection plots are also saved next to it as PNG and PDF, for talks and notes:
    limit_vs_runtime.{png,pdf}   limit_by_beam.{png,pdf}

Usage:
    python make_page.py --tag ~/pyfitter/run1a
    python make_page.py --tag ~/pyfitter/run1a --out run1a.html --cutcount-dir run1a_count/cutcount
    python make_page.py --tag ~/pyfitter/run1a --runtime-summary ~/pyfitter/studies/runtime/summary.csv \
        --beam '1.5 kW=~/pyfitter/studies/beampower/1BB1500W' '1BB=~/pyfitter/run1a' '2BB=~/pyfitter/run1a_2BB'
"""

import argparse
import csv
import datetime
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# Fixed order: each process keeps its color whatever is present
PROCESS_ORDER = ['CE', 'DIO', 'Cosmic', 'RPC', 'RMC']
MOM_BIN = 0.1    # MeV/c
TIME_BIN = 25.0  # ns


def canonical_process_name(process):
    """Same mapping as run_chain.py and run_cutcount.py."""
    name = str(process).strip().lower()
    if name.startswith('rpc'):
        return 'RPC'
    if name.startswith('rmc'):
        return 'RMC'
    if name in ('cosmic', 'cosmics'):
        return 'Cosmic'
    if name in ('ce', 'cem'):
        return 'CE'
    return str(process)


def generated_in_sample(cutflow_path):
    with open(cutflow_path, newline='') as f:
        return float(next(csv.DictReader(f))['events_passing'])


def histograms(manifest, p_window, t_window, mom_edges, time_edges):
    """Weighted histograms per process: all tracks, and tracks inside the other variable's window."""
    out = {}
    for entry in manifest['datasets']:
        if entry.get('status') != 'ok' or not entry.get('process'):
            continue
        process = canonical_process_name(entry['process'])
        df = pd.read_parquet(entry['parquet'], columns=['momentum', 'time'])
        if process == 'Cosmic':
            weight = float(entry['acceptance'])
        elif entry.get('expected') is None:
            continue
        else:
            weight = float(entry['expected']) * float(entry['acceptance']) / generated_in_sample(entry['cutflow'])
        # float64, as in run_cutcount.py: a float32 compare rounds the window edges and moves tracks
        p, t = df['momentum'].to_numpy(np.float64), df['time'].to_numpy(np.float64)
        in_p = (p > p_window[0]) & (p < p_window[1])
        in_t = (t > t_window[0]) & (t < t_window[1])
        h = out.setdefault(process, {'mom_all': 0, 'mom_win': 0, 'time_all': 0, 'time_win': 0,
                                     'fit_yield': 0.0, 'window_yield': 0.0, 'window_tracks': 0, 'tracks': 0})
        h['mom_all'] = h['mom_all'] + np.histogram(p, mom_edges)[0] * weight
        h['mom_win'] = h['mom_win'] + np.histogram(p[in_t], mom_edges)[0] * weight
        h['time_all'] = h['time_all'] + np.histogram(t, time_edges)[0] * weight
        h['time_win'] = h['time_win'] + np.histogram(t[in_p], time_edges)[0] * weight
        h['fit_yield'] += len(df) * weight
        h['window_yield'] += int(np.sum(in_p & in_t)) * weight
        h['window_tracks'] += int(np.sum(in_p & in_t))
        h['tracks'] += len(df)
    for h in out.values():
        for key in ('mom_all', 'mom_win', 'time_all', 'time_win'):
            h[key] = [float(v) for v in np.asarray(h[key], dtype=float)]
    return {p: out[p] for p in PROCESS_ORDER if p in out} | {p: v for p, v in out.items() if p not in PROCESS_ORDER}


def count_limit(cutcount_dir):
    path = cutcount_dir / 'count_limits_asymptotic.yaml'
    if not path.exists():
        return None
    with open(path) as f:
        res = yaml.safe_load(f)
    asimov = res['asimov']
    return {
        'method': 'Asymptotic CLs (Asimov), one bin',
        'events': asimov['band_+0']['events'], 'rmue': asimov['band_+0']['rmue'],
        'rmue_lo': asimov['band_-1']['rmue'], 'rmue_hi': asimov['band_+1']['rmue'],
        'efficiency': res['signal_efficiency'], 'background': res['background_yield'],
        'source': str(path),
    }


def shape_limit(tag_dir):
    for dim in ('2d', '1d'):
        path = tag_dir / 'limits' / f'card_{dim}_cls_limits.npz'
        if path.exists():
            d = np.load(path)
            events, rmue = d['asimov_events'], d['asimov_rmue']
            return {
                'method': f'Asymptotic CLs (Asimov), {dim.upper()} shape fit',
                'events': float(events[2]), 'rmue': float(rmue[2]),
                'rmue_lo': float(rmue[1]), 'rmue_hi': float(rmue[3]),
                'efficiency': float(d['signal_efficiency']), 'background': None,
                'source': str(path),
            }
    return None


def runtime_study(path):
    """Points of the run_studies.py runtime study that have a limit, in run-time order."""
    with open(path, newline='') as f:
        rows = [r for r in csv.DictReader(f) if r.get('status') == 'ok' and r.get('rmue_median')]
    rows.sort(key=lambda r: float(r['realtime_days']))
    return [{'point': r['point'], 'days': float(r['realtime_days']), 'beam_mode': r['beam_mode'],
             'power': float(r['beam_power_kW']), 'pot': float(r['pot']),
             'rmue': float(r['rmue_median']), 'rmue_lo': float(r['rmue_low']), 'rmue_hi': float(r['rmue_high'])}
            for r in rows]


def beam_comparison(specs):
    """Shape-fit limit of each LABEL=TAG, in the order given. A tag without a limit is kept as not run."""
    out = []
    for spec in specs:
        label, sep, tag = spec.partition('=')
        if not sep:
            raise SystemExit(f'--beam expects LABEL=TAG, got {spec!r}')
        tag_dir = Path(tag).expanduser().resolve()
        res = shape_limit(tag_dir)
        if res:
            res['source'] = relative_to_tag(res['source'], tag_dir)
        else:
            print(f'No limit in {tag_dir}/limits: {label} shown as not run')
        out.append({'label': label.strip(), 'limit': res})
    return out


def _limit_axes():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.set_yscale('log')
    ax.set_ylabel(r'90% CL upper limit on $R_{\mu e}$', fontsize=13)
    return plt, fig, ax


def _save(plt, fig, ax, out_dir, name):
    ax.tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out_dir / f'{name}.{ext}', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Wrote {out_dir / name}.png and .pdf')


def plot_runtime(points, out_dir):
    """Median expected limit with its 68% band against real run time, as run_studies.py draws it."""
    plt, fig, ax = _limit_axes()
    days = np.array([p['days'] for p in points])
    med, lo, hi = (np.array([p[k] for p in points]) for k in ('rmue', 'rmue_lo', 'rmue_hi'))
    ax.fill_between(days, lo, hi, color='#66cc66', alpha=0.5, lw=0, label=r'Expected $\pm 1\sigma$')
    ax.plot(days, med, 'k--', marker='o', ms=5, label='Expected median')
    ax.plot(days, med[0] * days[0] / days, color='gray', lw=1, ls=':', label=r'$\propto 1/T$ from first point')
    ax.minorticks_on()
    ax.set_xlim(0, days[-1] * 1.04)
    ax.set_xlabel('Real run time [days]', fontsize=13)
    ax.set_title(f'{points[0]["beam_mode"]}, {points[0]["power"]:.2g} kW', loc='right', fontsize=12)
    ax.legend(frameon=False)
    _save(plt, fig, ax, out_dir, 'limit_vs_runtime')


def plot_beams(beams, out_dir):
    """Median expected limit and 68% band of each beam mode that has a limit."""
    beams = [b for b in beams if b['limit']]
    plt, fig, ax = _limit_axes()
    x = np.arange(len(beams))
    med = np.array([b['limit']['rmue'] for b in beams])
    lo = np.array([b['limit']['rmue_lo'] for b in beams])
    hi = np.array([b['limit']['rmue_hi'] for b in beams])
    ax.errorbar(x, med, yerr=[med - lo, hi - med], fmt='o', color='k', ms=7, capsize=6, lw=1.5,
                label=r'Expected median, $\pm 1\sigma$')
    for xi, m in zip(x, med):
        e = int(np.floor(np.log10(m)))
        ax.annotate(rf'${m / 10 ** e:.1f}\times10^{{{e}}}$', (xi, m), xytext=(10, 0),
                    textcoords='offset points', va='center', fontsize=11)
    ax.set_xticks(x, [b['label'] for b in beams], fontsize=12)
    ax.set_xlim(-0.5, len(beams) - 0.2)  # room for the last value label
    ax.set_xlabel('Beam mode', fontsize=13)
    ax.legend(frameon=False)
    ax.tick_params(axis='x', which='minor', bottom=False, top=False)  # categories: no minor ticks
    _save(plt, fig, ax, out_dir, 'limit_by_beam')


def relative_to_tag(path, tag_dir):
    """<tag name>/<path inside the tag>, or just the file name for a file outside the tag."""
    path = Path(path)
    try:
        return str(Path(tag_dir.name) / path.resolve().relative_to(tag_dir))
    except ValueError:
        return path.name


def main():
    parser = argparse.ArgumentParser(description='Web page with the distributions and sensitivities of a tag')
    parser.add_argument('--tag', required=True, help='Tag directory (run_chain.py output)')
    parser.add_argument('--cutcount-dir', default=None, help='Default: <tag>/cutcount')
    # The page is published from the display submodule (github.com/sophiemiddleton/senisitivity_display)
    parser.add_argument('--out', default=str(Path(__file__).resolve().parent / 'display' / 'index.html'))
    parser.add_argument('--title', default='Run-1A sensitivity')
    parser.add_argument('--runtime-summary', default=None,
                        help='studies/runtime/summary.csv from run_studies.py, for the limit vs run time plot')
    parser.add_argument('--beam', nargs='+', default=[], metavar='LABEL=TAG',
                        help='Tags to compare by beam mode, e.g. "1.5 kW=studies/beampower/1BB1500W" "2BB=run1a_2BB"')
    args = parser.parse_args()

    tag_dir = Path(args.tag).expanduser().resolve()
    cutcount_dir = Path(args.cutcount_dir).expanduser().resolve() if args.cutcount_dir else tag_dir / 'cutcount'
    if not (tag_dir / 'manifest.yaml').exists():
        raise SystemExit(f'No manifest.yaml in {tag_dir}: run the selection for this tag first '
                         '(is PYFITTER_DIR set? source setup.sh)')
    with open(tag_dir / 'manifest.yaml') as f:
        manifest = yaml.safe_load(f)
    card_path = cutcount_dir / 'count_card.yaml'
    if not card_path.exists():
        raise SystemExit(f'No counting card at {card_path}: run cutcount/run_cutcount.py first')
    with open(card_path) as f:
        card = yaml.safe_load(f)
    p_window, t_window = card['window']['mom'], card['window']['time']

    mom_range = (float(manifest['fitrange_low'][0]), float(manifest['fitrange_hi'][0]))
    time_range = (float(manifest['fitrange_low'][1]), float(manifest['fitrange_hi'][1]))
    mom_edges = np.arange(mom_range[0], mom_range[1] + MOM_BIN / 2, MOM_BIN)
    time_edges = np.arange(time_range[0], time_range[1] + TIME_BIN / 2, TIME_BIN)

    data = {
        'title': args.title,
        'tag': manifest.get('tag', tag_dir.name),
        'version': card.get('metadata', {}).get('version'),
        'generated': datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
        'window': {'mom': p_window, 'time': t_window},
        'range': {'mom': mom_range, 'time': time_range},
        'edges': {'mom': [float(v) for v in mom_edges], 'time': [float(v) for v in time_edges]},
        'processes': histograms(manifest, p_window, t_window, mom_edges, time_edges),
        'count': count_limit(cutcount_dir),
        'shape': shape_limit(tag_dir),
        'runtime': runtime_study(Path(args.runtime_summary).expanduser()) if args.runtime_summary else [],
        'beams': beam_comparison(args.beam),
        'sources': {'manifest': str(tag_dir / 'manifest.yaml'), 'card': str(card_path)},
    }
    # The page is public: show inputs relative to the tag, not as paths in someone's area
    for res in [data['sources'], data['count'] or {}, data['shape'] or {}]:
        for key in ('manifest', 'card', 'source'):
            if key in res:
                res[key] = relative_to_tag(res[key], tag_dir)
    payload = json.dumps(data).replace('</', '<\\/')
    html = TEMPLATE.replace('__TITLE__', args.title).replace('__DATA__', payload)
    out = Path(args.out)
    if not out.parent.exists():
        raise SystemExit(f'{out.parent} does not exist: run git submodule update --init, or pass --out')
    out.write_text(html)
    print(f'Wrote {args.out}')
    if len(data['runtime']) >= 2:
        plot_runtime(data['runtime'], out.parent)
    if any(b['limit'] for b in data['beams']):
        plot_beams(data['beams'], out.parent)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7;
  --surface-1: #fcfcfb;
  --text-primary: #0b0b0b;
  --text-secondary: #52514e;
  --text-muted: #898781;
  --grid: #e1e0d9;
  --axis: #c3c2b7;
  --border: rgba(11,11,11,0.10);
  --window: rgba(11,11,11,0.05);
  --s-CE: #2a78d6; --s-DIO: #eb6834; --s-Cosmic: #1baf7a; --s-RPC: #eda100; --s-RMC: #e87ba4;
  --band: rgba(42,120,214,0.18);
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface-1: #1a1a19;
    --text-primary: #f0efec; --text-secondary: #c3c2b7; --text-muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10); --window: rgba(255,255,255,0.06);
    --s-CE: #3987e5; --s-DIO: #d95926; --s-Cosmic: #199e70; --s-RPC: #c98500; --s-RMC: #d55181;
    --band: rgba(57,135,229,0.24);
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface-1: #1a1a19;
  --text-primary: #f0efec; --text-secondary: #c3c2b7; --text-muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10); --window: rgba(255,255,255,0.06);
  --s-CE: #3987e5; --s-DIO: #d95926; --s-Cosmic: #199e70; --s-RPC: #c98500; --s-RMC: #d55181;
  --band: rgba(57,135,229,0.24);
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--text-primary);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1080px; margin: 0 auto; padding: 32px 16px 48px; }
h1 { font-size: 24px; font-weight: 600; margin: 0 0 4px; }
h2 { font-size: 16px; font-weight: 600; margin: 0 0 2px; }
.sub { color: var(--text-secondary); margin: 0 0 24px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 16px; margin-bottom: 12px; }
.tile, .panel { background: var(--surface-1); border: 1px solid var(--border); border-radius: 12px; padding: 20px; }
.tile .label { color: var(--text-secondary); font-size: 14px; }
.tile .value { font-size: 36px; font-weight: 600; margin: 4px 0; }
.tile .detail { color: var(--text-secondary); font-size: 14px; }
.tile .method { color: var(--text-muted); font-size: 13px; margin-top: 8px; }
.compare { color: var(--text-secondary); margin: 0 0 24px; }
.controls { display: flex; flex-wrap: wrap; gap: 12px 24px; align-items: center; margin-bottom: 16px; }
.seg { display: inline-flex; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }
.seg button { font: inherit; font-size: 14px; background: var(--surface-1); color: var(--text-secondary);
  border: 0; padding: 6px 12px; cursor: pointer; }
.seg button[aria-pressed="true"] { background: var(--text-primary); color: var(--surface-1); }
.legend { display: flex; flex-wrap: wrap; gap: 4px 16px; font-size: 14px; color: var(--text-secondary); }
.legend span { display: inline-flex; align-items: center; gap: 6px; }
.key { width: 16px; height: 2px; border-radius: 1px; display: inline-block; }
.windowkey { width: 14px; height: 12px; background: var(--window); border: 1px solid var(--axis); display: inline-block; }
.charts { display: grid; gap: 16px; margin-bottom: 24px; }
.projections { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; margin-bottom: 24px; }
.projections .legend { margin-bottom: 8px; }
.bandkey { width: 16px; height: 10px; background: var(--band); display: inline-block; }
.files { font-size: 13px; color: var(--text-muted); margin: 8px 0 0; }
.files a { color: var(--s-CE); }
.dashkey { width: 16px; border-top: 2px dotted var(--text-muted); display: inline-block; }
.panel .note { color: var(--text-secondary); font-size: 14px; margin: 0 0 8px; }
.chart { position: relative; }
.chart svg { display: block; width: 100%; height: 320px; }
.tooltip { position: absolute; pointer-events: none; background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 8px; padding: 8px 10px; font-size: 13px; box-shadow: 0 2px 8px rgba(0,0,0,0.12);
  display: none; min-width: 160px; z-index: 2; }
.tooltip .row { display: flex; justify-content: space-between; gap: 12px; font-variant-numeric: tabular-nums; }
.tooltip .row span:first-child { display: inline-flex; align-items: center; gap: 6px; color: var(--text-secondary); }
.tooltip .head { color: var(--text-primary); font-weight: 600; margin-bottom: 4px; }
table { width: 100%; border-collapse: collapse; font-size: 14px; }
th, td { text-align: right; padding: 8px 10px; border-bottom: 1px solid var(--grid); font-variant-numeric: tabular-nums; }
th:first-child, td:first-child { text-align: left; }
th { color: var(--text-secondary); font-weight: 600; }
tfoot td { font-weight: 600; border-bottom: 0; }
.tablewrap { overflow-x: auto; }
footer { color: var(--text-muted); font-size: 13px; margin-top: 24px; word-break: break-all; }
</style>
</head>
<body>
<main>
  <h1 id="title"></h1>
  <p class="sub" id="subtitle"></p>

  <div class="tiles" id="tiles"></div>
  <p class="compare" id="compare"></p>

  <div class="controls">
    <div class="seg" role="group" aria-label="Y scale">
      <button data-scale="log" aria-pressed="true">Log</button>
      <button data-scale="linear" aria-pressed="false">Linear</button>
    </div>
    <div class="seg" role="group" aria-label="Tracks shown">
      <button data-sel="all" aria-pressed="true">All selected tracks</button>
      <button data-sel="win" aria-pressed="false">Inside the other window</button>
    </div>
    <div class="legend" id="legend"></div>
  </div>

  <div class="charts">
    <section class="panel">
      <h2>Reconstructed momentum</h2>
      <p class="note" id="mom-note"></p>
      <div class="chart" id="chart-mom"></div>
    </section>
    <section class="panel">
      <h2>Reconstructed time</h2>
      <p class="note" id="time-note"></p>
      <div class="chart" id="chart-time"></div>
    </section>
  </div>

  <div class="projections" id="projections" hidden>
    <section class="panel" id="panel-runtime" hidden>
      <h2>Limit vs run time</h2>
      <p class="note" id="runtime-note"></p>
      <div class="legend">
        <span><i class="key" style="background:var(--s-CE)"></i>Expected median</span>
        <span><i class="bandkey"></i>68% band</span>
        <span><i class="dashkey"></i>∝ 1/T from the first point</span>
      </div>
      <div class="chart" id="chart-runtime"></div>
      <p class="files">Figure: <a href="limit_vs_runtime.png">PNG</a> · <a href="limit_vs_runtime.pdf">PDF</a></p>
    </section>
    <section class="panel" id="panel-beams" hidden>
      <h2>Limit by beam mode</h2>
      <p class="note">Expected 90% CL upper limit from the shape fit of each tag. Bars: 68% band.</p>
      <div class="chart" id="chart-beams"></div>
      <p class="files">Figure: <a href="limit_by_beam.png">PNG</a> · <a href="limit_by_beam.pdf">PDF</a></p>
    </section>
  </div>

  <section class="panel">
    <h2>Expected yields</h2>
    <p class="note">Events expected in the run. MC tracks is the raw number of simulated tracks in the counting window.</p>
    <div class="tablewrap"><table id="yields"></table></div>
  </section>

  <footer id="footer"></footer>
</main>

<script>
const DATA = __DATA__;
const NS = 'http://www.w3.org/2000/svg';
const state = { scale: 'log', sel: 'all' };
const procs = Object.keys(DATA.processes);

function sci(x, digits = 2) {
  if (x == null || !isFinite(x)) return '–';
  const e = Math.floor(Math.log10(Math.abs(x)));
  return `${(x / 10 ** e).toFixed(digits)} × 10<sup>${e}</sup>`;
}
function num(x) {
  if (x === 0) return '0';
  if (Math.abs(x) >= 100) return x.toFixed(0);
  if (Math.abs(x) >= 1) return x.toFixed(2);
  if (Math.abs(x) >= 0.001) return x.toPrecision(2);
  return x.toExponential(1);
}
function decade(x) {
  const e = Math.round(Math.log10(x));
  if (e >= 0 && e <= 3) return String(10 ** e);
  const sup = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  return '10' + String(e).split('').map(c => sup[c]).join('');
}
function el(tag, attrs = {}, parent) {
  const node = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (parent) parent.appendChild(node);
  return node;
}
const color = p => `var(--s-${p}, var(--text-muted))`;

// ---------- header, tiles, legend, table ----------
document.getElementById('title').textContent = DATA.title;
const w = DATA.window;
document.getElementById('subtitle').textContent =
  `Tag ${DATA.tag}${DATA.version ? ', cut set ' + DATA.version : ''}. ` +
  `Counting window ${w.mom[0]}–${w.mom[1]} MeV/c and ${w.time[0]}–${w.time[1]} ns. Expected 90% CL upper limits.`;

function tile(name, res) {
  const div = document.createElement('div');
  div.className = 'tile';
  if (!res) {
    div.innerHTML = `<div class="label">${name}</div><div class="value">Not run</div>`;
    return div;
  }
  div.innerHTML = `<div class="label">${name}: expected R<sub>μe</sub> limit</div>
    <div class="value">${sci(res.rmue)}</div>
    <div class="detail">${res.events.toFixed(2)} signal events · 68% band ${sci(res.rmue_lo)} to ${sci(res.rmue_hi)}</div>
    <div class="detail">CE efficiency ${(100 * res.efficiency).toFixed(1)}%${res.background != null ? ` · background ${res.background.toFixed(2)} events in the window` : ''}</div>
    <div class="method">${res.method}</div>`;
  return div;
}
const tiles = document.getElementById('tiles');
tiles.appendChild(tile('Cut and count', DATA.count));
tiles.appendChild(tile('Shape fit', DATA.shape));
if (DATA.count && DATA.shape) {
  const r = DATA.count.rmue / DATA.shape.rmue;
  document.getElementById('compare').textContent = r >= 1
    ? `The shape fit limit is ${r.toFixed(2)}× stronger than the cut and count.`
    : `The cut and count limit is ${(1 / r).toFixed(2)}× stronger than the shape fit.`;
}

const legend = document.getElementById('legend');
for (const p of procs) legend.insertAdjacentHTML('beforeend', `<span><i class="key" style="background:${color(p)}"></i>${p}</span>`);
legend.insertAdjacentHTML('beforeend', `<span><i class="windowkey"></i>Counting window</span>`);

const table = document.getElementById('yields');
let rows = '<thead><tr><th>Process</th><th>Counting window</th><th>Fit range</th><th>MC tracks in window</th></tr></thead><tbody>';
let bW = 0, bF = 0;
for (const p of procs) {
  const h = DATA.processes[p];
  if (p !== 'CE') { bW += h.window_yield; bF += h.fit_yield; }
  rows += `<tr><td><span class="legend"><span><i class="key" style="background:${color(p)}"></i>${p}</span></span></td>
    <td>${num(h.window_yield)}</td><td>${num(h.fit_yield)}</td><td>${h.window_tracks.toLocaleString()}</td></tr>`;
}
rows += `</tbody><tfoot><tr><td>Total background</td><td>${num(bW)}</td><td>${num(bF)}</td><td></td></tr></tfoot>`;
table.innerHTML = rows;

const r = DATA.range;
document.getElementById('footer').innerHTML =
  `Generated ${DATA.generated} by make_page.py. Selection range ${r.mom[0]}–${r.mom[1]} MeV/c, ${r.time[0]}–${r.time[1]} ns.<br>` +
  `Manifest: ${DATA.sources.manifest}<br>Counting card: ${DATA.sources.card}` +
  (DATA.count ? `<br>Cut-and-count limit: ${DATA.count.source}` : '') +
  (DATA.shape ? `<br>Shape-fit limit: ${DATA.shape.source}` : '') +
  (DATA.runtime.length ? '<br>Limit vs run time: studies/runtime/summary.csv from run_studies.py' : '') +
  DATA.beams.filter(b => b.limit).map(b => `<br>${b.label} limit: ${b.limit.source}`).join('');

// ---------- charts ----------
function niceStep(range, target) {
  const raw = range / target, mag = 10 ** Math.floor(Math.log10(raw)), f = raw / mag;
  return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag;
}
function linearTicks(lo, hi, target) {
  const step = niceStep(hi - lo, target), out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}

function draw(key) {
  const host = document.getElementById(`chart-${key}`);
  host.innerHTML = '';
  const edges = DATA.edges[key], win = DATA.window[key];
  const field = `${key}_${state.sel}`;
  const series = procs.map(p => ({ p, v: DATA.processes[p][field] }));
  const W = Math.max(host.clientWidth, 280), H = 320;
  const m = { l: 64, r: 12, t: 12, b: 40 };
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img',
    'aria-label': `${key === 'mom' ? 'Momentum' : 'Time'} distribution per process` }, host);

  const x0 = edges[0], x1 = edges[edges.length - 1];
  const X = v => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r);
  const all = series.flatMap(s => s.v), vmax = Math.max(...all, 1e-9);
  let Y, yTicks, yLo;
  if (state.scale === 'log') {
    const pos = all.filter(v => v > 0);
    const lo = Math.max(Math.min(...pos), vmax * 1e-5);
    const dLo = Math.floor(Math.log10(lo)), dHi = Math.ceil(Math.log10(vmax * 1.5));
    yLo = 10 ** dLo;
    Y = v => H - m.b - (Math.log10(v) - dLo) / (dHi - dLo) * (H - m.t - m.b);
    yTicks = [];
    const every = Math.ceil((dHi - dLo + 1) / 7);
    for (let d = dLo; d <= dHi; d += every) yTicks.push(10 ** d);
  } else {
    yTicks = linearTicks(0, vmax * 1.1, 5);
    const top = Math.max(yTicks[yTicks.length - 1], vmax * 1.1);
    yLo = 0;
    Y = v => H - m.b - v / top * (H - m.t - m.b);
  }

  // window band, behind everything
  el('rect', { x: X(win[0]), y: m.t, width: X(win[1]) - X(win[0]), height: H - m.t - m.b, fill: 'var(--window)' }, svg);

  // grid and y axis
  for (const t of yTicks) {
    const y = Y(t);
    el('line', { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: 'var(--grid)', 'stroke-width': 1 }, svg);
    const lab = el('text', { x: m.l - 8, y: y + 4, 'text-anchor': 'end', fill: 'var(--text-muted)', 'font-size': 12 }, svg);
    lab.textContent = state.scale === 'log' ? decade(t) : num(t);
  }
  el('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: 'var(--axis)', 'stroke-width': 1 }, svg);
  for (const t of linearTicks(x0, x1, W < 500 ? 5 : 10)) {
    const lab = el('text', { x: X(t), y: H - m.b + 18, 'text-anchor': 'middle', fill: 'var(--text-muted)', 'font-size': 12 }, svg);
    lab.textContent = t;
  }
  const xl = el('text', { x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle', fill: 'var(--text-secondary)', 'font-size': 13 }, svg);
  xl.textContent = key === 'mom' ? 'Momentum at tracker front [MeV/c]' : 'Time at tracker front [ns]';
  const binw = edges[1] - edges[0];
  const yl = el('text', { transform: `translate(14 ${(m.t + H - m.b) / 2}) rotate(-90)`, 'text-anchor': 'middle',
    fill: 'var(--text-secondary)', 'font-size': 13 }, svg);
  yl.textContent = `Events / ${+binw.toFixed(3)} ${key === 'mom' ? 'MeV/c' : 'ns'}`;

  // step lines; on log scale bins below the axis floor break the line
  for (const s of series) {
    let d = '', open = false;
    s.v.forEach((v, i) => {
      const ok = state.scale === 'linear' || v >= yLo;
      if (!ok) { open = false; return; }
      const y = Y(v);
      d += open ? `V${y}H${X(edges[i + 1])}` : `M${X(edges[i])} ${y}H${X(edges[i + 1])}`;
      open = true;
    });
    el('path', { d, fill: 'none', stroke: color(s.p), 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }, svg);
  }

  // hover: bin highlight and tooltip
  const hl = el('rect', { y: m.t, height: H - m.t - m.b, fill: 'var(--text-primary)', opacity: 0.06, visibility: 'hidden' }, svg);
  const tip = document.createElement('div');
  tip.className = 'tooltip';
  host.appendChild(tip);
  const hit = el('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: 'transparent' }, svg);
  const unit = key === 'mom' ? 'MeV/c' : 'ns';
  function show(evt) {
    const box = svg.getBoundingClientRect();
    const px = (evt.clientX - box.left) * W / box.width;
    const i = Math.min(edges.length - 2, Math.max(0, Math.floor((px - m.l) / (W - m.l - m.r) * (edges.length - 1))));
    hl.setAttribute('x', X(edges[i])); hl.setAttribute('width', Math.max(1, X(edges[i + 1]) - X(edges[i])));
    hl.setAttribute('visibility', 'visible');
    let html = `<div class="head">${+edges[i].toFixed(2)}–${+edges[i + 1].toFixed(2)} ${unit}</div>`;
    let bkg = 0;
    for (const s of series) {
      if (s.p !== 'CE') bkg += s.v[i];
      html += `<div class="row"><span><i class="key" style="background:${color(s.p)}"></i>${s.p}</span><span>${num(s.v[i])}</span></div>`;
    }
    html += `<div class="row"><span>Background</span><span>${num(bkg)}</span></div>`;
    tip.innerHTML = html;
    tip.style.display = 'block';
    const left = (evt.clientX - box.left) + 16;
    tip.style.left = `${left + tip.offsetWidth > box.width ? left - tip.offsetWidth - 32 : left}px`;
    tip.style.top = '16px';
  }
  hit.addEventListener('pointermove', show);
  hit.addEventListener('pointerleave', () => { tip.style.display = 'none'; hl.setAttribute('visibility', 'hidden'); });
}

function notes() {
  const inner = state.sel === 'win';
  document.getElementById('mom-note').textContent = inner
    ? `Tracks with ${w.time[0]} < t < ${w.time[1]} ns. Shaded: momentum window.`
    : 'All tracks passing the selection. Shaded: momentum window.';
  document.getElementById('time-note').textContent = inner
    ? `Tracks with ${w.mom[0]} < p < ${w.mom[1]} MeV/c. Shaded: time window.`
    : 'All tracks passing the selection. Shaded: time window.';
}
// ---------- projections: limit vs run time and by beam mode ----------
function logScale(lo, hi, y0, y1) {
  const dLo = Math.floor(Math.log10(lo)), dHi = Math.ceil(Math.log10(hi));
  const Y = v => y0 - (Math.log10(v) - dLo) / (dHi - dLo) * (y0 - y1);
  // decades, plus 2 and 5 when the axis spans few decades
  const ticks = [];
  for (let d = dLo; d <= dHi; d++)
    for (const f of dHi - dLo <= 2 ? [1, 2, 5] : [1]) if (f * 10 ** d <= 10 ** dHi) ticks.push({ v: f * 10 ** d, f, d });
  return { Y, ticks };
}
function tickLabel(t) {
  const p = decade(10 ** t.d);
  return t.f === 1 ? p : `${t.f}×${p}`;
}
function frame(host, label, H = 300) {
  host.innerHTML = '';
  const W = Math.max(host.clientWidth, 280);
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': label }, host);
  svg.style.height = `${H}px`;
  const tip = document.createElement('div');
  tip.className = 'tooltip';
  host.appendChild(tip);
  return { svg, tip, W, H, m: { l: 88, r: 16, t: 12, b: 44 } };
}
function yAxis(svg, sc, m, W) {
  for (const t of sc.ticks) {
    const y = sc.Y(t.v);
    el('line', { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: 'var(--grid)', 'stroke-width': 1 }, svg);
    const lab = el('text', { x: m.l - 8, y: y + 4, 'text-anchor': 'end', fill: 'var(--text-muted)', 'font-size': 12 }, svg);
    lab.textContent = tickLabel(t);
  }
}
function axisTitles(svg, W, H, m, xText) {
  const xl = el('text', { x: (m.l + W - m.r) / 2, y: H - 6, 'text-anchor': 'middle', fill: 'var(--text-secondary)', 'font-size': 13 }, svg);
  xl.textContent = xText;
  const yl = el('text', { transform: `translate(14 ${(m.t + H - m.b) / 2}) rotate(-90)`, 'text-anchor': 'middle',
    fill: 'var(--text-secondary)', 'font-size': 13 }, svg);
  yl.innerHTML = '90% CL limit on R<tspan font-size="10" dy="3">μe</tspan>';
}
function placeTip(tip, host, px, py) {
  tip.style.display = 'block';
  const left = px + 16;
  tip.style.left = `${left + tip.offsetWidth > host.clientWidth ? px - tip.offsetWidth - 16 : left}px`;
  tip.style.top = `${Math.max(0, py - tip.offsetHeight / 2)}px`;
}

function drawRuntime() {
  const pts = DATA.runtime;
  const host = document.getElementById('chart-runtime');
  const { svg, tip, W, H, m } = frame(host, 'Expected R_mue limit against real run time');
  const xMax = pts[pts.length - 1].days * 1.04;
  const X = v => m.l + v / xMax * (W - m.l - m.r);
  const ref = pts.map(p => pts[0].rmue * pts[0].days / p.days);
  const sc = logScale(Math.min(...pts.map(p => p.rmue_lo), ...ref), Math.max(...pts.map(p => p.rmue_hi)), H - m.b, m.t);
  yAxis(svg, sc, m, W);
  el('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: 'var(--axis)', 'stroke-width': 1 }, svg);
  for (const t of linearTicks(0, xMax, W < 500 ? 4 : 8)) {
    const lab = el('text', { x: X(t), y: H - m.b + 18, 'text-anchor': 'middle', fill: 'var(--text-muted)', 'font-size': 12 }, svg);
    lab.textContent = t;
  }
  axisTitles(svg, W, H, m, 'Real run time [days]');

  const band = pts.map(p => `${X(p.days)} ${sc.Y(p.rmue_hi)}`).join('L') + 'L' +
    pts.slice().reverse().map(p => `${X(p.days)} ${sc.Y(p.rmue_lo)}`).join('L');
  el('path', { d: `M${band}Z`, fill: 'var(--band)' }, svg);
  el('path', { d: 'M' + pts.map((p, i) => `${X(p.days)} ${sc.Y(ref[i])}`).join('L'), fill: 'none',
    stroke: 'var(--text-muted)', 'stroke-width': 1.5, 'stroke-dasharray': '2 4', 'stroke-linecap': 'round' }, svg);
  el('path', { d: 'M' + pts.map(p => `${X(p.days)} ${sc.Y(p.rmue)}`).join('L'), fill: 'none',
    stroke: 'var(--s-CE)', 'stroke-width': 2, 'stroke-linejoin': 'round' }, svg);
  const dots = pts.map(p => el('circle', { cx: X(p.days), cy: sc.Y(p.rmue), r: 4, fill: 'var(--s-CE)',
    stroke: 'var(--surface-1)', 'stroke-width': 2 }, svg));

  const hit = el('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: 'transparent' }, svg);
  hit.addEventListener('pointermove', evt => {
    const box = svg.getBoundingClientRect();
    const px = (evt.clientX - box.left) * W / box.width;
    let i = 0;
    pts.forEach((p, j) => { if (Math.abs(X(p.days) - px) < Math.abs(X(pts[i].days) - px)) i = j; });
    dots.forEach((d, j) => d.setAttribute('r', j === i ? 6 : 4));
    const p = pts[i];
    tip.innerHTML = `<div class="head">${+p.days.toFixed(1)} days (${(p.days / (365 / 12)).toFixed(1)} months)</div>
      <div class="row"><span>Median</span><span>${sci(p.rmue)}</span></div>
      <div class="row"><span>68% band</span><span>${sci(p.rmue_lo, 1)} – ${sci(p.rmue_hi, 1)}</span></div>
      <div class="row"><span>POT</span><span>${sci(p.pot)}</span></div>`;
    placeTip(tip, host, X(p.days) * box.width / W, sc.Y(p.rmue) * box.height / H);
  });
  hit.addEventListener('pointerleave', () => { tip.style.display = 'none'; dots.forEach(d => d.setAttribute('r', 4)); });
}

function drawBeams() {
  const beams = DATA.beams;
  const host = document.getElementById('chart-beams');
  const { svg, tip, W, H, m } = frame(host, 'Expected R_mue limit by beam mode');
  const done = beams.filter(b => b.limit).map(b => b.limit);
  const sc = logScale(Math.min(...done.map(r => r.rmue_lo)), Math.max(...done.map(r => r.rmue_hi)), H - m.b, m.t);
  yAxis(svg, sc, m, W);
  el('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: 'var(--axis)', 'stroke-width': 1 }, svg);
  axisTitles(svg, W, H, m, 'Beam mode');
  const step = (W - m.l - m.r) / beams.length;
  beams.forEach((b, i) => {
    const x = m.l + step * (i + 0.5);
    const lab = el('text', { x, y: H - m.b + 18, 'text-anchor': 'middle', fill: 'var(--text-secondary)', 'font-size': 13 }, svg);
    lab.textContent = b.label;
    const r = b.limit;
    if (!r) {
      const na = el('text', { x, y: (m.t + H - m.b) / 2, 'text-anchor': 'middle', fill: 'var(--text-muted)', 'font-size': 13 }, svg);
      na.textContent = 'Not run';
      return;
    }
    el('line', { x1: x, x2: x, y1: sc.Y(r.rmue_lo), y2: sc.Y(r.rmue_hi), stroke: 'var(--s-CE)', 'stroke-width': 2 }, svg);
    for (const v of [r.rmue_lo, r.rmue_hi])
      el('line', { x1: x - 6, x2: x + 6, y1: sc.Y(v), y2: sc.Y(v), stroke: 'var(--s-CE)', 'stroke-width': 2 }, svg);
    el('circle', { cx: x, cy: sc.Y(r.rmue), r: 5, fill: 'var(--s-CE)', stroke: 'var(--surface-1)', 'stroke-width': 2 }, svg);
    const val = el('text', { x: x + 12, y: sc.Y(r.rmue) + 4, fill: 'var(--text-primary)', 'font-size': 12 }, svg);
    const e = Math.floor(Math.log10(r.rmue));
    val.textContent = `${(r.rmue / 10 ** e).toFixed(1)}×${decade(10 ** e)}`;
    const hit = el('rect', { x: x - step / 2, y: m.t, width: step, height: H - m.t - m.b, fill: 'transparent' }, svg);
    hit.addEventListener('pointermove', () => {
      const box = svg.getBoundingClientRect();
      tip.innerHTML = `<div class="head">${b.label}</div>
        <div class="row"><span>Median</span><span>${sci(r.rmue)}</span></div>
        <div class="row"><span>68% band</span><span>${sci(r.rmue_lo, 1)} – ${sci(r.rmue_hi, 1)}</span></div>
        <div class="row"><span>Signal events</span><span>${r.events.toFixed(2)}</span></div>
        <div class="row"><span>CE efficiency</span><span>${(100 * r.efficiency).toFixed(1)}%</span></div>`;
      placeTip(tip, host, x * box.width / W, sc.Y(r.rmue) * box.height / H);
    });
    hit.addEventListener('pointerleave', () => { tip.style.display = 'none'; });
  });
}

const hasRuntime = DATA.runtime.length >= 2, hasBeams = DATA.beams.some(b => b.limit);
document.getElementById('projections').hidden = !(hasRuntime || hasBeams);
document.getElementById('panel-runtime').hidden = !hasRuntime;
document.getElementById('panel-beams').hidden = !hasBeams;
if (hasRuntime) {
  const p0 = DATA.runtime[0];
  document.getElementById('runtime-note').textContent =
    `Expected 90% CL upper limit from the shape fit, ${p0.beam_mode} at ${+p0.power.toFixed(2)} kW. ` +
    'Only the normalization changes between points: the selection is the same.';
}

function render() {
  notes(); draw('mom'); draw('time');
  if (hasRuntime) drawRuntime();
  if (hasBeams) drawBeams();
}

for (const b of document.querySelectorAll('[data-scale]')) b.addEventListener('click', () => {
  state.scale = b.dataset.scale;
  document.querySelectorAll('[data-scale]').forEach(o => o.setAttribute('aria-pressed', o === b));
  render();
});
for (const b of document.querySelectorAll('[data-sel]')) b.addEventListener('click', () => {
  state.sel = b.dataset.sel;
  document.querySelectorAll('[data-sel]').forEach(o => o.setAttribute('aria-pressed', o === b));
  render();
});
let resizeTimer;
window.addEventListener('resize', () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(render, 100); });
render();
</script>
</body>
</html>
"""


if __name__ == '__main__':
    main()
