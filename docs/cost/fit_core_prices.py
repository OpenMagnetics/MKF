# Fit the COST filter's core price laws, price = c * mass^k (US$, kg), one per MAS material type and
# materialComposition ("ferrite/MnZn", "powder/FeSiAl", ...). "proprietary" names no alloy, so its laws also
# carry the manufacturer they were fitted on ("powder/proprietary/Micrometals"): they never price another
# maker's proprietary grade. Run from this folder with the MAS submodule checked out.
#
# Prices, all at the tier a 1000-piece order pays (the highest break <= 1000) and only where that tier is a
# break of at least 10 pieces; single retail prices (qty 1-9) are left out:
#   core_prices.csv        DigiKey list prices read through findchips.com, 2026-09-30.
#   core_prices_extra.csv  NiZn, carbonyl iron, iron and nanocrystalline cores from DigiKey, Farnell, TME,
#                          rf-microwave.com and others, 2026-09-30; EUR/GBP at the ECB reference rate of
#                          2026-09-30 (unit_price_usd). One price per part: the median over its distributors.
#   core_prices_gap.csv    Micrometals iron and proprietary-grade toroids from Power Magnetics (UK, GBP; its
#                          highest break is 50 pieces) and Feryster (PL, PLN; its breaks are a fixed discount
#                          on the 1-piece price, not negotiated volume prices), Proterial AMCC amorphous cut
#                          cores from eu.mouser.com (EUR; per core pair, unverified), 2026-10-01.
#   core_prices_eu_check.csv  EU prices of parts already priced at DigiKey, used only as a check (printed).
# Mass is what MKF's Core::get_mass() returns: for cores with a datasheet Ve, Ve times the MAS density of the
# material (Ve is always per set / per toroid); otherwise (the toroids of the extra file) the datasheet mass,
# which for a toroid is the geometric volume times the density that get_mass() uses.
import csv, collections, statistics as st, math, json

MIN_PARTS = 10
CLASS = {'MnZn ferrite': 'ferrite/MnZn', 'powder: Kool Mu (Sendust)': 'powder/FeSiAl', 'powder: MPP': 'powder/FeNiMo',
         'powder: High Flux': 'powder/FeNi', 'powder: Edge': 'powder/FeNi', 'powder: XFlux': 'powder/FeSi',
         'powder: iron (Micrometals mix)': 'powder/iron', 'powder: Micrometals MS (Sendust-type)': 'powder/FeSiAl',
         'powder: Micrometals MP (MPP-type)': 'powder/FeNiMo'}
EXTRA_TYPE = {'NiZn': 'ferrite', 'MgZn': 'ferrite', 'carbonylIron': 'powder', 'iron': 'powder', 'proprietary': 'nanocrystalline'}
GAP_TYPE = {'iron': 'powder', 'carbonylIron': 'powder', 'proprietary': 'powder', 'FeSi': 'amorphous'}

def law_key(material_type, composition, manufacturer):
    return f"{material_type}/{composition}" + (f"/{manufacturer}" if composition == 'proprietary' else '')
# Price per set: a per-half quote doubles. Ferroxcube's shape datasheets that state it print "Selling unit: PCS"
# for E and EQ (E14/3.5/5, E18/4/10, EQ13, EQ38/8/25), so Ferroxcube E/EQ quotes with no stated basis are per
# half. Other unstated TDK/Ferroxcube quotes follow the datasheet mass unit (g/pcs -> half, g/set -> set).
BASIS = {'per toroid': 1, 'per set': 1, 'per half (piece)': 2, 'likely per set (unverified)': 1,
         'likely per half (unverified)': 2, 'unverified': 1, 'per piece (single core, not a set)': 1}

def per_set(x):
    if x['price_basis'] == 'unverified' and x['manufacturer'] == 'Ferroxcube' and x['shape_family'] in ('E', 'EQ'):
        return 2
    return BASIS[x['price_basis']]

density = {}
for line in open('../../MAS/data/core_materials.ndjson'):
    d = json.loads(line)
    if 'density' in d:
        density[d['name']] = d['density']

offers = collections.defaultdict(list)   # (law, mpn) -> [(mass, price)]
skipped = collections.Counter()
for x in csv.DictReader(open('core_prices.csv')):
    if x['applies_to_1000pc_order'] != '1' or x['material_class'] not in CLASS:
        continue
    if float(x['qty_break']) < 10:
        skipped['main: break below 10 pieces'] += 1
        continue
    if not x['ve_mm3']:
        skipped['main: no Ve'] += 1
        continue
    mass = float(x['ve_mm3']) * 1e-9 * density[x['material']]
    offers[(CLASS[x['material_class']], x['mpn'])].append((mass, float(x['unit_price']) * per_set(x)))
for x in csv.DictReader(open('core_prices_extra.csv')):
    comp = x['material_composition']
    if x['applies_to_1000pc_order'] != '1' or comp not in EXTRA_TYPE or not x['unit_price_usd']:
        continue
    if float(x['qty_break'] or 1) < 10:
        skipped['extra: break below 10 pieces'] += 1
        continue
    if x['shape_family'] != 'T' or not x['mass_g'] or x['price_basis'] not in BASIS:
        skipped['extra: not a toroid with a datasheet mass'] += 1
        continue
    offers[(law_key(EXTRA_TYPE[comp], comp, x['manufacturer']), x['mpn'])].append((float(x['mass_g']) / 1000, float(x['unit_price_usd']) * per_set(x)))
for x in csv.DictReader(open('core_prices_gap.csv')):
    comp = x['material_composition']
    if x['applies_to_1000pc_order'] != '1' or comp not in GAP_TYPE or not x['unit_price_usd'] or not x['mass_g']:
        continue
    if float(x['qty_break'] or 1) < 10:
        skipped['gap: break below 10 pieces'] += 1
        continue
    if float(x['unit_price_usd']) <= 0:
        skipped['gap: listed at 0 (no offer)'] += 1
        continue
    # Toroids (Micrometals OC/OD/GX grades are toroids too) and the amorphous C-core pairs; E/U cores
    # are left out because it is not stated whether a piece is a half or a set.
    if comp == 'FeSi':
        factor = 1   # "per core set (C-core pair)"
    elif x['shape_family'] in ('T', 'O', 'G'):
        factor = 1
    else:
        skipped['gap: E/U core, half or set not stated'] += 1
        continue
    offers[(law_key(GAP_TYPE[comp], comp, x['manufacturer']), x['mpn'])].append((float(x['mass_g']) / 1000, float(x['unit_price_usd']) * factor))
print('skipped', dict(skipped))

points = collections.defaultdict(list)
for (law, mpn), o in offers.items():
    points[law].append((st.median(m for m, _ in o), st.median(p for _, p in o)))

def ols(d):
    n = len(d)
    mx = sum(a for a, _ in d) / n
    my = sum(b for _, b in d) / n
    b = sum((a - mx) * (c - my) for a, c in d) / sum((a - mx) ** 2 for a, _ in d)
    return my - b * mx, b

for law, d in sorted(points.items()):
    if len(d) < MIN_PARTS:
        print(f"{law:28s} n={len(d)}: too few parts, no law")
        continue
    la, k = ols([(math.log(m), math.log(p)) for m, p in d])
    c = math.exp(la)
    err = sorted(abs(math.log(p / (c * m ** k))) for m, p in d)
    print(f"{law:28s} c={c:8.3f} k={k:.4f} n={len(d):3d} mass {min(m for m, _ in d)*1e3:.2f}-{max(m for m, _ in d)*1e3:.0f} g "
          f"spread median x{math.exp(st.median(err)):.2f} p90 x{math.exp(err[int(0.9 * len(err))]):.2f}")

# Check: EU volume prices over DigiKey's for the same part.
dk = {x['mpn']: float(x['unit_price']) for x in csv.DictReader(open('core_prices.csv'))
      if x['applies_to_1000pc_order'] == '1' and x['distributor'] == 'DigiKey'}
ratio = collections.defaultdict(list)
for x in csv.DictReader(open('core_prices_eu_check.csv')):
    if x['applies_to_1000pc_order'] == '1' and x['mpn'] in dk and x['unit_price_usd'] and float(x['qty_break'] or 1) >= 10:
        ratio[x['distributor'].split(' (')[0]].append(float(x['unit_price_usd']) / dk[x['mpn']])
for shop, r in sorted(ratio.items()):
    print(f"check: {shop:32s} n={len(r):3d} EU/DigiKey geometric mean x{math.exp(st.mean(map(math.log, r))):.2f}")
