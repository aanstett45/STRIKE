# Graph datasets

The benchmark commands expect the following local files:

```text
data/email-Enron.txt
data/roadNet-CA.txt
```

They are SNAP graph edge lists. Download them with:

```bash
bash scripts/download_snap_data.sh
```

Do not silently replace, filter, deduplicate, reorder, or otherwise preprocess these files. The C++ loader ignores comment lines, remaps vertex identifiers, removes self-loops, merges duplicate edges, and symmetrizes the graph internally.
