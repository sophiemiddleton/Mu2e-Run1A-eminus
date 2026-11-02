# Set up the environment for the Run-1A scripts. Source it, from any directory:
#
#     source setup.sh
#
# It
#   1. sets up the Mu2e environment (mu2einit) if that has not been done,
#   2. activates the Mu2e python environment (pyenv ana <version>),
#   3. checks out the pyfitter submodule at the commit this repository pins, and exports
#      PYFITTER_DIR and RUN1A_DIR,
#   4. checks that the python packages and pyfitter pieces the scripts use are there.
#
# The pyfitter version is the submodule commit recorded in this repository (git submodule status).
# To move to a new pyfitter tag:
#     git -C pyfitter fetch --tags && git -C pyfitter checkout <tag> && git add pyfitter && git commit
#
# Override any of these before sourcing:
#   MU2E_ANA_VERSION  pyenv ana version (default 2.8.0, the version the scripts are tested with)
#   PYFITTER_DIR      another pyfitter checkout to use instead of the submodule (for development;
#                     results made this way are not tied to the pinned version)

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "setup.sh must be sourced: source ${0}" >&2
    exit 1
fi

_run1a_setup() {
    local repo_dir
    repo_dir="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
    export RUN1A_DIR="$repo_dir"

    # 1. Mu2e environment
    if ! declare -F pyenv > /dev/null; then
        echo "[setup] Setting up the Mu2e environment (mu2einit)"
        source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh || return 1
    fi

    # 2. Python environment
    local version="${MU2E_ANA_VERSION:-2.8.0}"
    echo "[setup] Activating pyenv ana ${version}"
    pyenv ana "$version" || return 1

    # 3. pyfitter: the pinned submodule unless PYFITTER_DIR points elsewhere
    if [ -z "$PYFITTER_DIR" ]; then
        export PYFITTER_DIR="$repo_dir/pyfitter"
        if [ ! -f "$PYFITTER_DIR/run_chain.py" ]; then
            echo "[setup] Checking out the pyfitter submodule"
            git -C "$repo_dir" submodule update --init pyfitter || {
                echo "[setup] Checkout failed. Mu2e/pyfitter is private: check your GitHub ssh key (ssh -T git@github.com)" >&2
                return 1
            }
        fi
        local pinned actual
        pinned="$(git -C "$repo_dir" ls-tree HEAD pyfitter | awk '{print $3}')"
        actual="$(git -C "$PYFITTER_DIR" rev-parse HEAD 2>/dev/null)"
        if [ -n "$pinned" ] && [ "$pinned" != "$actual" ]; then
            echo "[setup] WARNING: pyfitter is at ${actual:0:7}, but this repository pins ${pinned:0:7}." >&2
            echo "[setup]          Run 'git submodule update pyfitter' to use the pinned version." >&2
        fi
    else
        echo "[setup] Using PYFITTER_DIR from the environment (not the pinned submodule)"
    fi
    echo "[setup] pyfitter: $PYFITTER_DIR ($(git -C "$PYFITTER_DIR" describe --tags --always --dirty 2>/dev/null))"

    # 4. Checks
    python - "$PYFITTER_DIR" << 'EOF' || return 1
import importlib, sys
from pathlib import Path
missing = []
for module in ('pyutils', 'zfit', 'hepstats', 'uproot', 'awkward', 'numpy', 'scipy', 'pandas', 'pyarrow', 'yaml', 'matplotlib'):
    try:
        importlib.import_module(module)
    except ImportError:
        missing.append(module)
pyfitter = Path(sys.argv[1])
for name in ('run_chain.py', 'process.py', 'analyze.py', 'config.py', 'normalization.py',
             'uncertainties_config.py', 'datasets.yaml', 'sensitivity_scan/simple_limit_combine_card.py'):
    if not (pyfitter / name).exists():
        missing.append(f'pyfitter/{name}')
if missing:
    print('[setup] MISSING: ' + ', '.join(missing), file=sys.stderr)
    sys.exit(1)
print(f'[setup] python {sys.version.split()[0]}: all packages and pyfitter files found')
EOF

    echo "[setup] Ready. For example:"
    echo "    python $repo_dir/sensitivity/cutcount/run_cutcount.py --tag run1a --version run1a --group mix --loc tape --jobs 16"
    echo "    python $repo_dir/sensitivity/make_page.py --tag run1a"
    echo "[setup] Reading files from tape or disk (--loc tape/disk) also needs a grid token: getToken"
}

_run1a_setup
_run1a_status=$?
unset -f _run1a_setup
return $_run1a_status
