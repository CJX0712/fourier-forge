// fourier-forge probe — ASCII plot original vs reconstructed + error table + engine vectors.
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const html = fs.readFileSync(path.join(__dirname, 'index.html'), 'utf8');
const m = html.match(/<script id="engine">([\s\S]*?)<\/script>/);
const ctx = {};
vm.createContext(ctx);
const Fourier = vm.runInContext(m[1] + '\nFourier;', ctx);

function heartPts(N){
  const pts = [];
  for (let i = 0; i < N; i++){
    const t = i * 2 * Math.PI / N;
    pts.push({ x: 16 * Math.sin(t) ** 3, y: 13 * Math.cos(t) - 5 * Math.cos(2 * t) - 2 * Math.cos(3 * t) - Math.cos(4 * t) });
  }
  return pts;
}

function asciiPlot(src, rec, cols, rows){
  const grid = Array.from({ length: rows }, () => Array(cols).fill(' '));
  let minX = 1e9, maxX = -1e9, minY = 1e9, maxY = -1e9;
  for (const p of src){ minX = Math.min(minX, p.x); maxX = Math.max(maxX, p.x); minY = Math.min(minY, p.y); maxY = Math.max(maxY, p.y); }
  const sx = (cols - 2) / (maxX - minX), sy = (rows - 2) / (maxY - minY);
  const put = (p, ch) => {
    const gx = Math.round((p.x - minX) * sx) + 1, gy = Math.round((maxY - p.y) * sy) + 1;
    if (gx >= 0 && gx < cols && gy >= 0 && gy < rows) grid[gy][gx] = ch;
  };
  for (const p of src) put(p, '#');
  if (rec) for (const p of rec) put(p, 'o');
  return grid.map(r => r.join('')).join('\n');
}

const N = 256;
const heart = heartPts(N);
const coefs = Fourier.dft(heart);

console.log('=== probe: original heart vs reconstruction K=3 ===');
{
  const rec = Fourier.partialSum(coefs, Fourier.topIndices(coefs, 3));
  console.log(asciiPlot(heart, rec, 62, 30));
  console.log('errRms K=3: ' + Fourier.errRms(rec, heart).toFixed(3));
}
console.log('\n=== probe: original heart vs reconstruction K=25 ===');
{
  const rec = Fourier.partialSum(coefs, Fourier.topIndices(coefs, 25));
  console.log(asciiPlot(heart, rec, 62, 30));
  console.log('errRms K=25: ' + Fourier.errRms(rec, heart).toFixed(3));
}
console.log('\n=== probe: error table K -> RMS ===');
for (const K of [1, 2, 3, 5, 8, 13, 25, 50, 100, 256]){
  const rec = Fourier.partialSum(coefs, Fourier.topIndices(coefs, K));
  console.log('K=' + String(K).padStart(3) + '  errRms=' + Fourier.errRms(rec, heart).toFixed(4));
}
console.log('\n=== probe: engine vectors ===');
{
  const c = Fourier.dft([{ x: 1, y: 0 }, { x: 0, y: 0 }, { x: 0, y: 0 }, { x: 0, y: 0 }]);
  console.log('impulse dft magnitudes: ' + c.map(x => Fourier.magnitude(x).toFixed(6)).join(' '));
  console.log('parseval ratio: ' + Fourier.parsevalRatio(heart, coefs).toFixed(12));
  console.log('round-trip err: ' + Fourier.errRms(Fourier.idft(coefs), heart).toExponential(3));
}
console.log('\nPROBE OK');
