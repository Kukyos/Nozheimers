// Any-organism mode: NCBI ortholog list -> UniProt sequence -> AlphaFold model, aligned and superposed in the browser
// with the same rules as site/build_explorer.py (global BLOSUM62 alignment, gap -11/-1; Kabsch on Ca with pLDDT >= 70).
(() => {
  const B62 = (() => {
    const aa = 'ARNDCQEGHILKMFPSTWYV';
    const rows = `4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4`.split('\n').map(r => r.trim().split(/\s+/).map(Number));
    const idx = {}; [...aa].forEach((a, i) => idx[a] = i);
    return (a, b) => (a in idx && b in idx) ? rows[idx[a]][idx[b]] : -1;
  })();

  // Gotoh global alignment, BLOSUM62, gap open -11 (first residue), extend -1. Returns {score, map: human pos -> org pos}.
  function align(h, o) {
    const n = h.length, m = o.length, NEG = -1e9, W = m + 1;
    const M = new Float64Array((n + 1) * W).fill(NEG), X = new Float64Array((n + 1) * W).fill(NEG), Y = new Float64Array((n + 1) * W).fill(NEG);
    M[0] = 0;
    for (let i = 1; i <= n; i++) X[i * W] = -11 - (i - 1);
    for (let j = 1; j <= m; j++) Y[j] = -11 - (j - 1);
    for (let i = 1; i <= n; i++) for (let j = 1; j <= m; j++) {
      const k = i * W + j, d = (i - 1) * W + (j - 1), u = (i - 1) * W + j, l = i * W + j - 1;
      M[k] = Math.max(M[d], X[d], Y[d]) + B62(h[i - 1], o[j - 1]);
      X[k] = Math.max(M[u] - 11, X[u] - 1, Y[u] - 11);
      Y[k] = Math.max(M[l] - 11, Y[l] - 1, X[l] - 11);
    }
    let i = n, j = m, k = n * W + m, st = M[k] >= X[k] && M[k] >= Y[k] ? 0 : X[k] >= Y[k] ? 1 : 2;
    const score = Math.max(M[k], X[k], Y[k]), map = {};
    while (i > 0 || j > 0) {
      k = i * W + j;
      if (st === 0 && i > 0 && j > 0) {
        map[i] = j;
        const d = (i - 1) * W + (j - 1), s = M[k] - B62(h[i - 1], o[j - 1]);
        st = s === M[d] ? 0 : s === X[d] ? 1 : 2; i--; j--;
      } else if (st === 1 || j === 0) {
        const u = (i - 1) * W + j; st = i > 0 && X[k] === M[u] - 11 ? 0 : X[k] === X[u] - 1 ? 1 : 2; i--;
        if (i === 0 && j > 0) st = 2;
      } else {
        const l = i * W + j - 1; st = j > 0 && Y[k] === M[l] - 11 ? 0 : Y[k] === Y[l] - 1 ? 2 : 1; j--;
        if (j === 0 && i > 0) st = 1;
      }
    }
    return { score, map };
  }

  // Kabsch via Horn's quaternion method. Returns {rot, tran, rmsd} with y_on_x = y . rot + tran (Biopython convention).
  function kabsch(X, Y) {
    const n = X.length, cx = [0, 0, 0], cy = [0, 0, 0];
    X.forEach(p => p.forEach((v, i) => cx[i] += v / n)); Y.forEach(p => p.forEach((v, i) => cy[i] += v / n));
    const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let k = 0; k < n; k++) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) S[a][b] += (Y[k][a] - cy[a]) * (X[k][b] - cx[b]);
    const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = S;
    const N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
               [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]];
    const q = topEigenvector(N);
    const [w, a, b, c] = q;
    const Rm = [[w*w + a*a - b*b - c*c, 2 * (a*b - w*c), 2 * (a*c + w*b)],
                [2 * (a*b + w*c), w*w - a*a + b*b - c*c, 2 * (b*c - w*a)],
                [2 * (a*c - w*b), 2 * (b*c + w*a), w*w - a*a - b*b + c*c]];
    const rot = [0, 1, 2].map(r => [0, 1, 2].map(cc => Rm[cc][r]));   // row-vector convention
    const tran = [0, 1, 2].map(cc => cx[cc] - (cy[0] * rot[0][cc] + cy[1] * rot[1][cc] + cy[2] * rot[2][cc]));
    let ss = 0;
    for (let k = 0; k < n; k++) for (let cc = 0; cc < 3; cc++) {
      const v = Y[k][0] * rot[0][cc] + Y[k][1] * rot[1][cc] + Y[k][2] * rot[2][cc] + tran[cc] - X[k][cc]; ss += v * v;
    }
    return { rot, tran, rmsd: Math.sqrt(ss / n) };
  }
  function topEigenvector(A) {  // Jacobi rotations on a symmetric 4x4
    const a = A.map(r => r.slice()), v = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    for (let sweep = 0; sweep < 60; sweep++) {
      let off = 0; for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) off += a[p][q] ** 2;
      if (off < 1e-14) break;
      for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) {
        if (Math.abs(a[p][q]) < 1e-18) continue;
        const th = (a[q][q] - a[p][p]) / (2 * a[p][q]), t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 4; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
        for (let k = 0; k < 4; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < 4; k++) { const vkp = v[k][p], vkq = v[k][q]; v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq; }
      }
    }
    let best = 0; for (let i = 1; i < 4; i++) if (a[i][i] > a[best][best]) best = i;
    return [0, 1, 2, 3].map(k => v[k][best]);
  }

  const caMap = pdb => {
    const out = {};
    pdb.split('\n').forEach(l => { if (l.startsWith('ATOM') && l.slice(12, 16).trim() === 'CA') out[+l.slice(22, 26)] = [[+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)], +l.slice(60, 66)]; });
    return out;
  };
  function superpose(hca, oca, map, regions) {
    const fit = span => {
      const pairs = Object.entries(map).filter(([hp, op]) => hca[hp] && oca[op] && hca[hp][1] >= 70 && oca[op][1] >= 70 && (!span || (hp >= span[0] && hp <= span[1])));
      if (pairs.length < 10) return null;
      const r = kabsch(pairs.map(([hp]) => hca[hp][0]), pairs.map(([, op]) => oca[op][0]));
      return { ...r, rmsd: Math.round(r.rmsd * 100) / 100, n: pairs.length };
    };
    const fits = [];
    const whole = fit(null); if (whole) fits.push({ ...whole, region: 'whole protein' });
    regions.forEach(([label, lo, hi]) => { const f = fit([lo, hi]); if (f) fits.push({ ...f, region: `${label} (${lo}-${hi})` }); });
    if (!fits.length) return null;
    if (fits[0].region === 'whole protein' && fits[0].rmsd <= 3) return fits[0];
    const good = fits.filter(f => f.rmsd <= 3);
    return (good.length ? good : fits).sort((a, b) => b.n - a.n)[0];
  }

  const H = { Accept: 'application/json' };
  const getJSON = async u => { const r = await fetch(u, { headers: H }); if (!r.ok) throw new Error(`${r.status} from ${new URL(u).host}`); return r.json(); };

  async function liveCompute(entry) {
    const [taxid, taxname, common, geneId, symbol, swiss, group] = entry;
    const G = E.genes[gene], hseq = G.seq;
    const say = t => { $('#loading').style.display = 'grid'; $('#loading').textContent = t; };
    say(`Fetching ${symbol} from UniProt…`);
    const res = await getJSON(`https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(`xref:geneid-${geneId}`)}&fields=accession,sequence,organism_name&format=json&size=25`);
    let cands = res.results.map(r => ({ acc: r.primaryAccession, seq: r.sequence.value, latin: r.organism.scientificName }));
    for (const acc of swiss) if (!cands.find(c => c.acc === acc)) {
      try { const r = await getJSON(`https://rest.uniprot.org/uniprotkb/${acc}?fields=accession,sequence,organism_name&format=json`); cands.push({ acc, seq: r.sequence.value, latin: r.organism.scientificName }); } catch (e) { /* skip */ }
    }
    let how = `linked to NCBI Gene ${geneId}`;
    if (!cands.length) {  // many species lack the GeneID cross-reference in UniProt; fall back to gene name within the organism
      const q = `organism_id:${taxid} AND (gene_exact:${symbol} OR gene_exact:${gene})`;
      const r2 = await getJSON(`https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(q)}&fields=accession,sequence,organism_name&format=json&size=25`);
      cands = r2.results.map(r => ({ acc: r.primaryAccession, seq: r.sequence.value, latin: r.organism.scientificName }));
      how = `found by gene name ${symbol} in this organism (UniProt has no link to NCBI Gene ${geneId})`;
    }
    if (!cands.length) throw new Error(`UniProt has no ${symbol} entry for this organism.`);
    cands = cands.sort((a, b) => Math.abs(a.seq.length - hseq.length) - Math.abs(b.seq.length - hseq.length)).slice(0, 8);
    say(`Aligning ${cands.length} candidate sequence${cands.length > 1 ? 's' : ''} to human ${gene}…`);
    await new Promise(r => setTimeout(r, 30));
    const scored = cands.map(c => ({ ...c, ...align(hseq, c.seq) })).sort((a, b) => b.score - a.score);
    // prefer the best-scoring entry that has an AlphaFold model, among entries within 5% of the top score
    say('Fetching the AlphaFold model…');
    let best = scored[0], af = null, pdb = null, afErr = '';
    for (const c of scored.filter(c => c.score >= 0.95 * scored[0].score)) {
      try { af = (await getJSON(`https://alphafold.ebi.ac.uk/api/prediction/${c.acc}`))[0]; pdb = await (await fetch(af.pdbUrl)).text(); best = c; break; }
      catch (e) { afErr = e.message; af = null; }
    }
    const map = best.map, oseq = best.seq;
    const resAt = hp => map[hp] ? oseq[map[hp] - 1] : '-';
    const pairs = Object.keys(map).map(Number), same = pairs.filter(hp => hseq[hp - 1] === resAt(hp));
    const dis = [...new Set(G.disease.map(d => d.pos))];
    const core = gene === 'APP' ? Array.from({ length: 42 }, (_, i) => 672 + i)
      : [...new Set(G.features.filter(f => f.type === 'Transmembrane').flatMap(f => Array.from({ length: f.end - f.start + 1 }, (_, i) => f.start + i)))];
    const active = G.features.filter(f => f.type === 'Active site').map(f => ({ pos: f.start, human: hseq[f.start - 1], org: resAt(f.start) }));
    const drugs = E.drugs.filter(d => d.gene === gene).map(d => ({ pdb: d.pdb, drug: d.drug, total: d.residues.length,
      same: d.residues.filter(p => resAt(p) === hseq[p - 1]).length, diffs: d.residues.filter(p => resAt(p) !== hseq[p - 1]).map(p => `${hseq[p - 1]}${p}${resAt(p)}`) }));
    const rec = { name: common ? common[0].toUpperCase() + common.slice(1) : taxname, group, latin: best.latin, uniprot: best.acc, symbol, gene_id: geneId,
      source: 'NCBI ortholog', live: true,
      note: `Computed live in your browser. ${scored.length} UniProt candidate${scored.length > 1 ? 's' : ''} ${how}; ${best === scored[0] ? 'best alignment to human kept' : 'top match had no AlphaFold model, so the next-best match with one is shown'}. Not pre-checked.`,
      length: oseq.length, identity: +(100 * same.length / pairs.length).toFixed(1), coverage: +(100 * pairs.length / hseq.length).toFixed(1),
      map, seq: oseq, disease_total: dis.length, disease_same: dis.filter(p => resAt(p) === hseq[p - 1]).length,
      core_total: core.length, core_same: core.filter(p => resAt(p) === hseq[p - 1]).length, active, drugs };
    if (af) {
      const humanPdb = await (await fetch(G.human_model)).text();
      const hca = caMap(humanPdb), oca = caMap(pdb);
      const regions = G.features.filter(f => f.type === 'Domain' && f.end - f.start >= 30).map(f => [f.description || 'Domain', f.start, f.end]);
      if (gene === 'APP') regions.push(['amyloid-beta and membrane helix', 672, 723]);
      Object.assign(rec, { model_text: pdb, alphafold: af.entryId, superpose: superpose(hca, oca, map, regions),
        plddt_mean: +(Object.values(oca).reduce((s, v) => s + v[1], 0) / Object.keys(oca).length).toFixed(1) });
      if ((af.uniprotSequence || af.sequence) && (af.uniprotSequence || af.sequence) !== oseq) rec.note += ' AlphaFold model sequence differs from the current UniProt sequence.';
    } else rec.note += ` No AlphaFold model available (${afErr}).`;
    return rec;
  }

  // ---------- search UI ----------
  const types = ['All types', 'Primates', 'Rodents', 'Other mammals', 'Birds', 'Reptiles', 'Amphibians', 'Fish'];
  const sel = document.createElement('select'); sel.id = 'qtype'; sel.setAttribute('aria-label', 'Organism type');
  sel.innerHTML = types.map(t => `<option>${t}</option>`).join('');
  $('#q').after(sel);
  function search() {
    const list = (window.ORTHOLOGS || {})[gene] || [], q = $('#q').value.trim().toLowerCase(), ty = sel.value;
    const hits = list.filter(e => (ty === 'All types' || e[6] === ty) && (!q || e[1].toLowerCase().includes(q) || (e[2] || '').toLowerCase().includes(q)));
    $('#qinfo').textContent = `${hits.length} of ${list.length} organisms with an NCBI ${gene} ortholog`;
    $('#results').innerHTML = (q || ty !== 'All types') ? hits.slice(0, 60).map((e, i) =>
      `<button type="button" class="org" data-i="${list.indexOf(e)}"><i>${esc(e[1])}</i>${e[2] ? ` <span class="muted">${esc(e[2])}</span>` : ''}</button>`).join('') + (hits.length > 60 ? `<span class="muted" style="font-size:.8rem">${hits.length - 60} more; type to narrow</span>` : '') : '';
  }
  $('#q').addEventListener('input', search); sel.addEventListener('change', search);
  $('#genes').addEventListener('click', () => setTimeout(search, 0));
  $('#results').addEventListener('click', async ev => {
    const b = ev.target.closest('.org'); if (!b) return;
    const entry = window.ORTHOLOGS[gene][+b.dataset.i];
    $('#results').querySelectorAll('.org').forEach(x => x.setAttribute('aria-pressed', x === b));
    try { liveOrg = await liveCompute(entry); await render(); }
    catch (e) { $('#loading').style.display = 'grid'; $('#loading').textContent = `Could not build this organism: ${e.message}`; }
  });
  const s = document.createElement('script'); s.src = 'explorer/orthologs.js'; s.onload = search; document.head.append(s);
  window.__live = { align, kabsch, superpose, caMap, liveCompute };
  window.__liveSelfTest = () => {  // human vs human must give 100% and RMSD 0
    const G = E.genes[gene], a = align(G.seq, G.seq);
    return Object.keys(a.map).length === G.seq.length && Object.entries(a.map).every(([h, o]) => +h === o);
  };
})();
