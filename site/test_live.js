const fs = require('fs'), path = require('path');
const site = process.argv[2];
global.window = global;
const stub = () => ({ after() {}, addEventListener() {}, setAttribute() {}, style: {}, querySelectorAll: () => [], append() {} });
global.$ = stub; global.esc = s => s;
global.document = { createElement: stub, head: { append() {} } };
eval(fs.readFileSync(path.join(site, 'explorer/explorer.js'), 'utf8'));
global.E = window.EXPLORER;
eval(fs.readFileSync(path.join(site, 'explorer-live.js'), 'utf8'));
const { align, kabsch, superpose, caMap } = window.__live;
for (const [gene, tax] of [['APP', '10090'], ['PSEN1', '7955'], ['PSEN1', '6239'], ['APP', '7227']]) {
  const G = E.genes[gene], o = G.organisms[tax];
  const a = align(G.seq, o.seq);
  const keys = Object.keys(a.map), py = Object.keys(o.map);
  const sameMap = keys.length === py.length && keys.every(k => o.map[k] === a.map[k]);
  const ident = keys.filter(h => G.seq[h - 1] === o.seq[a.map[h] - 1]).length / keys.length * 100;
  const hca = caMap(fs.readFileSync(path.join(site, G.human_model), 'utf8')), oca = caMap(fs.readFileSync(path.join(site, o.model), 'utf8'));
  const regions = G.features.filter(f => f.type === 'Domain' && f.end - f.start >= 30).map(f => [f.description || 'Domain', f.start, f.end]);
  if (gene === 'APP') regions.push(['amyloid-beta and membrane helix', 672, 723]);
  const sp = superpose(hca, oca, a.map, regions);
  console.log(gene, o.name.padEnd(12), 'JS id', ident.toFixed(1), 'py id', o.identity, 'same map as python:', sameMap, '| JS rmsd', sp.rmsd, sp.n, sp.region, '| py', o.superpose.rmsd, o.superpose.n, o.superpose.region);
}
const G = E.genes.PSEN1, s = align(G.seq, G.seq); console.log('self-align identity map ok:', Object.entries(s.map).every(([h, o]) => +h === o) && Object.keys(s.map).length === G.seq.length);

// live path against the real services (needs network): an organism that is not a preset
(async () => {
  if (!process.argv.includes('--live')) return;
  const realFetch = global.fetch;
  global.fetch = (u, o = {}) => realFetch(String(u).startsWith('http') ? u : 'http://127.0.0.1:8765/' + u, { ...o, headers: { ...(o.headers || {}), 'User-Agent': 'Mozilla/5.0' } });
  eval(fs.readFileSync(path.join(site, 'explorer/orthologs.js'), 'utf8'));
  for (const [g, name] of [['PSEN1', 'Gorilla gorilla'], ['APP', 'Felis catus']]) {
    global.gene = g;
    const entry = window.ORTHOLOGS[g].find(e => e[1].startsWith(name));
    const rec = await window.__live.liveCompute(entry);
    console.log('LIVE', g, rec.name, rec.uniprot, 'id', rec.identity, 'disease', rec.disease_same + '/' + rec.disease_total, 'core', rec.core_same + '/' + rec.core_total, 'rmsd', rec.superpose && rec.superpose.rmsd, rec.alphafold, '|', rec.note);
  }
})();
