#!/bin/sh
set -eu
srcdir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
: "${OFFLINE_MAIN:?Source the sPHENIX environment first}"
(cd "$srcdir" && autoreconf --install --force -I "$OFFLINE_MAIN/share")
"$srcdir/configure" "$@"
