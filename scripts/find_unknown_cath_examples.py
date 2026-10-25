import json

data = json.load(open('examples/random_10000_results.json', encoding='utf-8'))
unknown_cath = []
for r in data:
    for m in r['members']:
        dom = m['domain']
        raw = dom.get('cath_label')
        src = dom.get('cath_estimation_source')
        if raw == '-' or (src and src != 'ted_official'):
            unknown_cath.append((r, m))

print(f"Total unknown/inferred CATH domains: {len(unknown_cath)}")
by_source = {}
for r, m in unknown_cath:
    src = m['domain'].get('cath_estimation_source') or 'unknown'
    by_source.setdefault(src, []).append((r, m))

for src, items in sorted(by_source.items()):
    print(f"\n=======================================================")
    print(f"Source: {src} ({len(items)} instances)")
    print(f"=======================================================")
    for i in range(min(3, len(items))):
        r, m = items[i]
        dom = m['domain']
        p_ctx = m.get('protein_context') or {}
        p_name = p_ctx.get('protein_name') or p_ctx.get('gene_name') or 'Uncharacterized'
        print(f"  Example {i+1}: {dom['ted_id']}")
        print(f"    UniProt:       {dom['uniprot_acc']} ({p_name})")
        print(f"    Chopping:      {dom['chopping']} (Length: {dom['nres_domain']} res, pLDDT: {dom['plddt']:.1f})")
        print(f"    Raw CATH:      {dom.get('cath_label')} -> Estimated CATH: {dom.get('inferred_cath_label')}")
        print(f"    Proposed Name: {r['consensus_name']} [Tier: {r['specificity_tier']}]")
        print(f"    Confidence:    {r['confidence_score']:.3f}")
        print(f"    Helices/Strands: {dom.get('num_helix')} / {dom.get('num_strand')}")
