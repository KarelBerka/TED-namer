import json
from collections import Counter

with open('examples/random_1000_results.json', encoding='utf-8') as f:
    data = json.load(f)

low_conf = sorted(data, key=lambda x: x['confidence_score'])

print(f"Total results: {len(data)}")
print(f"Confidence < 0.55: {sum(1 for x in data if x['confidence_score'] < 0.55)}")
print(f"Confidence < 0.50: {sum(1 for x in data if x['confidence_score'] < 0.50)}")
print(f"Confidence < 0.40: {sum(1 for x in data if x['confidence_score'] < 0.40)}")

print("\n--- 15 NEJNIŽŠÍCH CONFIDENCE SKÓRE ---")
for x in low_conf[:15]:
    mem = x['members'][0]
    top_hit = x['top_candidates'][0] if x['top_candidates'] else None
    top_name = top_hit['canonical_name'] if top_hit else 'None'
    top_iou = top_hit['mean_iou'] if top_hit else 0
    dbs = top_hit['source_databases'] if top_hit else []
    print(f"{x['cluster_id']}: conf={x['confidence_score']:.3f} | name='{x['consensus_name']}' | tier={x['specificity_tier']} | CATH={x['dominant_cath_label']} | chopping={mem['domain']['chopping']} ({mem['domain']['nres_domain']} res, plddt={mem['domain']['plddt']:.1f}) | dbs={dbs} | top_hit='{top_name}' (IoU={top_iou})")

print("\n--- DOMÉNY S NÁZVEM OBSAHUJÍCÍM 'Uncharacterized' NEBO 'unknown' ---")
unchar = [x for x in data if 'uncharacterized' in x['consensus_name'].lower() or 'unknown' in x['consensus_name'].lower() or 'duf' in x['consensus_name'].lower()]
print(f"Počet: {len(unchar)}")
for x in unchar[:10]:
    mem = x['members'][0]
    top_hit = x['top_candidates'][0] if x['top_candidates'] else None
    top_name = top_hit['canonical_name'] if top_hit else 'None'
    top_iou = top_hit['mean_iou'] if top_hit else 0
    print(f"{x['cluster_id']}: conf={x['confidence_score']:.3f} | name='{x['consensus_name']}' | tier={x['specificity_tier']} | CATH={x['dominant_cath_label']} | chopping={mem['domain']['chopping']} (plddt={mem['domain']['plddt']:.1f}) | top_hit='{top_name}' (IoU={top_iou})")

print("\n--- PŘÍPADY S VELMI NÍZKÝM PkLDDT (< 60) ---")
low_plddt = [x for x in data if x['members'][0]['domain']['plddt'] < 60]
print(f"Počet s pLDDT < 60: {len(low_plddt)}")
for x in low_plddt[:5]:
    mem = x['members'][0]
    print(f"{x['cluster_id']}: plddt={mem['domain']['plddt']:.1f} | name='{x['consensus_name']}' | CATH={x['dominant_cath_label']} | chopping={mem['domain']['chopping']}")
