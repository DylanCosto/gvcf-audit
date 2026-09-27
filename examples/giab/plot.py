"""Plot a checked GIAB overlap result (requires Matplotlib)."""

import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('overlap_json', type=Path)
parser.add_argument('independent_check_json', type=Path)
parser.add_argument('output_prefix', type=Path)
args = parser.parse_args()
data = json.loads(args.overlap_json.read_text())
check = json.loads(args.independent_check_json.read_text())
r = data['results']
if not check['match'] or check['bases_by_state'] != r['refseq_cds']['bases_by_state']:
    raise ValueError('The independent CDS check does not match this result')
groups = [
    ('refseq_cds', 'All coding regions'),
    ('cds_giab_high_confidence', 'Coding ∩ GIAB high confidence'),
    ('cds_difficult_regions', 'Coding ∩ difficult regions'),
    ('cds_low_mappability', 'Coding ∩ low mappability'),
    ('cds_segmental_duplications', 'Coding ∩ segmental duplications'),
]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(11, 5.8))
fig.subplots_adjust(left=.34, right=.97, top=.72, bottom=.22)
for y, (key, label) in enumerate(groups):
    value = r[key]['callable_percent']
    ax.barh(y, 100, color='#e5e9ed', height=.57)
    ax.barh(y, value, color='#237f88', height=.57)
    ax.text(value-1.3, y, f'{value:.2f}%', ha='right', va='center', color='white', weight='bold')
ax.set_yticks(range(len(groups)), [f"{label}\n{r[key]['requested_bases']:,} bases" for key,label in groups])
ax.invert_yaxis()
ax.set_xlim(0, 100)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel('Bases meeting the audit policy (%)')
ax.tick_params(axis='y', length=0, pad=12)
ax.spines[['top', 'right', 'left']].set_visible(False)
fig.text(.04, .94, 'Where does a gVCF lack usable genotype evidence?', fontsize=18, weight='bold')
fig.text(.04, .875, f"{data['sample']} • {data['caller']} • GRCh38 chromosomes 1–22", fontsize=12)
fig.legend(handles=[Patch(color='#237f88', label='Callable'), Patch(color='#e5e9ed', label='Unresolved')],
           loc='upper left', bbox_to_anchor=(.035, .84), frameon=False, ncol=2)
unresolved = r['refseq_cds']['requested_bases'] - r['refseq_cds']['callable_bases']
difficult = r['cds_difficult_regions']['requested_bases'] - r['cds_difficult_regions']['callable_bases']
fig.text(.04, .08, f'{unresolved:,} coding bases were unresolved; {100*difficult/unresolved:.2f}% fall in GIAB difficult regions.', weight='bold')
policy = data['policy']
fig.text(.04, .025, f"DP/MIN_DP ≥ {policy['min_dp']:g}; GQ ≥ {policy['min_gq']:g}. Categories overlap. Base-level callability is not variant accuracy.\n"
         f"gvcf-audit {data['audit_version']}; GIAB benchmark v4.2.1 and stratifications v3.6.", fontsize=9, color='#444444')
for ext in ['png', 'svg', 'pdf']:
    fig.savefig(str(args.output_prefix) + '.' + ext, dpi=180)
