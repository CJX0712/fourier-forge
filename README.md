# FourierForge · 傅里叶画图机

单文件、零依赖、可离线运行的**傅里叶级数画图机（epicycles）**——把任意闭合路径做 DFT，分解成一串旋转圆链（大圆套小圆），画笔钉在最外层小圆上，转起来就把路径画出来。

- 4 条预设路径：心形 / 星形 / 利萨茹 / 螺旋
- 系数个数 K 可调（1–256）：按幅度取前 K 个频率，实时看误差变化
- 实时显示 RMS 误差与 Parseval 能量比
- 内置自检：冲激/常数/单频向量、往返重构、能量守恒、截断单调性、共轭对称、线性
- 纯 JS 引擎，无外部依赖，打开 `index.html` 即用

## 引擎不变量

| 量 | 公式 | 说明 |
|----|------|------|
| DFT | `X[k] = Σ x[n]·e^(-2πikn/N)` | 复数点序列 |
| IDFT | `x[n] = (1/N) Σ X[k]·e^(2πikn/N)` | 往返误差 < 1e-9 |
| Parseval | `Σ|X|² = N·Σ|x|²` | 能量守恒，比值恒为 1 |
| 截断 | 按幅度取前 K 项 | 误差随 K 单调不增（最小二乘最优） |
| 直流项 | `X[0]/N = 质心` | 路径均点 |
| 共轭对称 | `X[N-k] = conj(X[k])` | 纯实输入 |

## 本地校验

```bash
node _smoke.js   # 10 项不变量测试
node _probe.js   # ASCII 绘图对比 + 误差表 + 引擎向量
```

## 引擎接口

```js
const Fourier = require('./index.html'); // 浏览器内为全局 const Fourier
Fourier.dft(pts);        // [{x,y}] -> [{re,im}]
Fourier.idft(coefs);     // 精确逆变换
Fourier.topIndices(c,K); // 按幅度降序取前 K 个频率下标
Fourier.partialSum(c,idx); // 用选定系数重构路径点
Fourier.errRms(a,b);     // 两组点集的 RMS 距离
Fourier.parsevalRatio(p,c); // 能量比（应恒为 1）
```

## 许可证

MIT — 见 [LICENSE](./LICENSE)。
