#!/usr/bin/env python3
import random
import subprocess
from collections import defaultdict
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR    = Path('/home/senekowitsch/Thesis/Functional/01_prokka')
OUTPUT_BASE = BASE_DIR / 'output'
GENOMES_DIR = Path('/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/genomes')
SWEEPS_FILE = Path('/home/senekowitsch/Thesis/Sweeps/05_check_distance/output/sweeps_bottomup_clonal_5x.txt')
SWEEP7_FILE = BASE_DIR / 'sweep_7.txt'
LABELS_FILE = BASE_DIR / 'genome_sweep_labels.txt'
REPS_FILE   = BASE_DIR / 'prokka_representatives.txt'
THREADS     = 30
RANDOM_SEED = 42

BASE_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

# ── Step 1: Build sweep labels from sweeps_bottomup_clonal_5x.txt ──────────
print("=== Step 1: Building sweep labels ===")

sweep_map = {}
with open(SWEEPS_FILE) as f:
    next(f)  # skip header
    for line in f:
        sweep_id, genome, *_ = line.strip().split('\t')
        genome_id = genome.removeprefix('out_')
        sweep_map[genome_id] = f'sweep_{sweep_id}'

# Assign labels to all genomes in the genomes directory
labels = {}
for fna in sorted(GENOMES_DIR.glob('*.fna')):
    genome_id = fna.stem
    labels[genome_id] = sweep_map.get(genome_id, 'no_sweep')

# ── Step 2: Add sweep_7 ────────────────────────────────────────────────────
print("=== Step 2: Adding sweep_7 labels ===")

sweep7 = set()
with open(SWEEP7_FILE) as f:
    for line in f:
        g = line.strip()
        if g:
            sweep7.add(g)

updated = 0
for genome_id in sweep7:
    if genome_id in labels:
        labels[genome_id] = 'sweep_7'
        updated += 1
    else:
        print(f"  WARNING: {genome_id} from sweep_7 not found in genomes directory")

print(f"  Updated {updated} genomes to sweep_7")

# Write labels file
with open(LABELS_FILE, 'w') as f:
    f.write('genome\tsweep\n')
    for genome_id, label in sorted(labels.items()):
        f.write(f'{genome_id}\t{label}\n')

# Print summary
from collections import Counter
counts = Counter(labels.values())
for label, count in sorted(counts.items()):
    print(f"  {label}: {count}")

# ── Step 3: Select representatives ────────────────────────────────────────
print("\n=== Step 3: Selecting representatives ===")

random.seed(RANDOM_SEED)

groups = defaultdict(list)
for genome_id, label in labels.items():
    groups[label].append(genome_id)

selected = []
for label, genomes in sorted(groups.items()):
    chosen = random.sample(genomes, min(3, len(genomes)))
    selected.extend((g, label) for g in chosen)
    print(f"  {label}: {chosen}")

with open(REPS_FILE, 'w') as f:
    f.write('genome\tsweep\n')
    for genome_id, label in selected:
        f.write(f'{genome_id}\t{label}\n')

print(f"\n  Total representatives: {len(selected)}")

# ── Step 4: Run Prokka ─────────────────────────────────────────────────────
print("\n=== Step 4: Running Prokka ===")

for genome_id, label in selected:
    fna = GENOMES_DIR / f'{genome_id}.fna'
    outdir = OUTPUT_BASE / label / genome_id

    if not fna.exists():
        print(f"  WARNING: {fna} not found, skipping.")
        continue

    if outdir.exists():
        print(f"  Skipping {genome_id} ({label}) — already done.")
        continue

    outdir.mkdir(parents=True, exist_ok=True)
    print(f"  Annotating {genome_id} ({label})...")

    cmd = [
        'prokka',
        '--outdir', str(outdir),
        '--prefix', genome_id,
        '--genus', 'Salmonella',
        '--species', 'Infantis',
        '--kingdom', 'Bacteria',
        '--cpus', str(THREADS),
        '--force',
        str(fna)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR on {genome_id}:\n{result.stderr}")
    else:
        print(f"  Done: {genome_id}")

print("\n=== All done! ===")
print(f"Output in: {OUTPUT_BASE}")