# Kept so older instructions still work: the setup now lives at the top of the repository.
#
#     source setup.sh

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "setup.sh must be sourced: source ${0}" >&2
    exit 1
fi

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/setup.sh"
