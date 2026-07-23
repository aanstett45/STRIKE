#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
mkdir -p data

curl -L --fail --retry 3 \
  https://snap.stanford.edu/data/email-Enron.txt.gz \
  -o data/email-Enron.txt.gz
curl -L --fail --retry 3 \
  https://snap.stanford.edu/data/roadNet-CA.txt.gz \
  -o data/roadNet-CA.txt.gz

gzip -df data/email-Enron.txt.gz
gzip -df data/roadNet-CA.txt.gz

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum data/email-Enron.txt data/roadNet-CA.txt > data/checksums.sha256
else
  shasum -a 256 data/email-Enron.txt data/roadNet-CA.txt > data/checksums.sha256
fi

date -u +'%Y-%m-%dT%H:%M:%SZ' > data/download_timestamp_utc.txt
printf 'Downloaded datasets and wrote data/checksums.sha256\n'
