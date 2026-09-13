// fourier-forge smoke test — extracts the inline Fourier engine and verifies DFT invariants.
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const html = fs.readFileSync(path.join(__dirname, 'index.html'), 'utf8');
const m = html.match(/<script id="engine">([\s\S]*?)<\/script>/);
if (!m) { console.error('FAIL: engine script not found'); process.exit(1); }
const ctx = {};
vm.createContext(ctx);
const Fourier = vm.runInContext(m[1] + '\nFourier;', ctx);

let pass = 0, fail = 0;
const T = (name, cond) => {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name); }
};
const near = (a, b, eps) => Math.abs(a - b) < eps;

console.log('FourierForge smoke test');

// heart + star paths (same as UI)
function heartPts(N){
  const pts = [];
  for (let i = 0; i < N; i++){
    const t = i * 2 * Math.PI / N;
    pts.push({ x: 16 * Math.sin(t) ** 3, y: 13 * Math.cos(t) - 5 * Math.cos(2 * t) - 2 * Math.cos(3 * t) - Math.cos(4 * t) });
  }
  return pts;
}
function starPts(N){ // 5-point star, triangle-wave radius (piecewise linear -> rich spectrum)
  const pts = [];
  for (let i = 0; i < N; i++){
    const t = i * 2 * Math.PI / N;
    const ph = (t * 5 / (2 * Math.PI)) % 1;
    const r = 15 - 9 * Math.abs(ph * 2 - 1);
    pts.push({ x: r * Math.sin(t), y: -r * Math.cos(t) });
  }
  return pts;
}

// 1) impulse vector: dft([1,0,0,0]) -> all |X[k]| = 1
{
  const c = Fourier.dft([{ x: 1, y: 0 }, { x: 0, y: 0 }, { x: 0, y: 0 }, { x: 0, y: 0 }]);
  T('impulse: dft([1,0,0,0]) all |X|=1', c.every(x => near(Fourier.magnitude(x), 1, 1e-9)));
}
// 2) constant sequence: N ones -> X[0]=N, rest 0
{
  const c = Fourier.dft(Array(8).fill({ x: 1, y: 0 }));
  T('constant: X[0]=8, rest 0', near(c[0].re, 8, 1e-9) && c.slice(1).every(x => near(Fourier.magnitude(x), 0, 1e-9)));
}
// 3) pure single-frequency circle -> only X[1] = N
{
  const pts = [];
  for (let i = 0; i < 8; i++){ const a = 2 * Math.PI * i / 8; pts.push({ x: Math.cos(a), y: Math.sin(a) }); }
  const c = Fourier.dft(pts);
  T('single-freq circle: only X[1]=8 nonzero', near(c[1].re, 8, 1e-9) && near(c[1].im, 0, 1e-9) && c.filter((x, i) => i !== 1).every(x => near(Fourier.magnitude(x), 0, 1e-9)));
}
// 4) round-trip exact
{
  const p = heartPts(256);
  const err = Fourier.errRms(Fourier.idft(Fourier.dft(p)), p);
  T('round-trip idft(dft(heart)) err < 1e-9 (got ' + err.toExponential(2) + ')', err < 1e-9);
}
// 5) Parseval energy conservation
{
  const p = heartPts(256);
  const r = Fourier.parsevalRatio(p, Fourier.dft(p));
  T('Parseval Σ|X|² = N·Σ|x|² (ratio ' + r.toFixed(12) + ')', near(r, 1, 1e-9));
}
// 6) truncation error monotonically non-increasing in K, exact at K=N (star: rich spectrum)
{
  const p = starPts(256);
  const c = Fourier.dft(p);
  const Ks = [1, 2, 4, 8, 16, 32, 64, 128, 256];
  let ok = true, prev = Infinity;
  const errs = [];
  for (const K of Ks){
    const e = Fourier.errRms(Fourier.partialSum(c, Fourier.topIndices(c, K)), p);
    if (e > prev + 1e-12) ok = false;
    errs.push(K + ':' + e.toFixed(3));
    prev = e;
  }
  T('truncation error non-increasing in K, K=N ≈ 0 (star)', ok && prev < 1e-9);
  console.log('    ' + errs.join('  '));
}
// 7) DC term = centroid
{
  const p = heartPts(256);
  const c = Fourier.dft(p);
  const mx = p.reduce((s, q) => s + q.x, 0) / p.length;
  const my = p.reduce((s, q) => s + q.y, 0) / p.length;
  T('DC term X[0]/N = centroid', near(c[0].re / 256, mx, 1e-9) && near(c[0].im / 256, my, 1e-9));
}
// 8) conjugate symmetry for purely real input: X[N-k] = conj(X[k])
{
  const p = heartPts(64).map(q => ({ x: q.x, y: 0 })); // real-valued sequence
  const c = Fourier.dft(p);
  let ok = true;
  for (let k = 1; k < 32; k++){
    if (!near(c[64 - k].re, c[k].re, 1e-9) || !near(c[64 - k].im, -c[k].im, 1e-9)) ok = false;
  }
  T('conjugate symmetry X[N-k]=conj(X[k]) for real input', ok);
}
// 9) linearity: dft(a+b) = dft(a)+dft(b)
{
  const N = 16;
  const a = [], b = [];
  for (let i = 0; i < N; i++){
    a.push({ x: Math.cos(2 * Math.PI * i / N), y: 0 });
    b.push({ x: 0, y: Math.sin(4 * Math.PI * i / N) });
  }
  const ab = a.map((p, i) => ({ x: p.x, y: b[i].y })); // a real + b imaginary
  const ca = Fourier.dft(a), cb = Fourier.dft(b), cab = Fourier.dft(ab);
  let ok = true;
  for (let k = 0; k < N; k++){
    if (!near(cab[k].re, ca[k].re + cb[k].re, 1e-9) || !near(cab[k].im, ca[k].im + cb[k].im, 1e-9)) ok = false;
  }
  T('linearity dft(a+b)=dft(a)+dft(b)', ok);
}
// 10) determinism
{
  const p = heartPts(128);
  const c1 = Fourier.dft(p), c2 = Fourier.dft(p);
  T('determinism: identical DFT twice', JSON.stringify(c1) === JSON.stringify(c2));
}

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
