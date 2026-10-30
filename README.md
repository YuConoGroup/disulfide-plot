<p align="center">
  <img src="assets/logo.svg" width="120" alt="disulfide-plot logo">
</p>

<p align="center">
  <img src="assets/wordmark.svg" width="720" alt="Disulfide Connectivity Plotter">
</p>

<p align="center">
  <a href="README.md">中文</a> · <a href="README.en.md">English</a>
</p>

# disulfide-plot

把肽序列的二硫键连接画成彩色括号图：序列横排、Cys 标红、每条键一条彩色括号。

## 效果

![两条序列的示例](examples/demo.png)

## 特性

- 纯 Pillow 实现，不需要 matplotlib，也不需要 tkinter / cairo
- 跨度大的键自动放上层
- 多条键分上下轨道，不重叠、不压字
- Cys 红色，不同键用不同颜色
- 输出 PNG，可选 SVG
- 支持 CSV，也支持直接读 BLAST outfmt 6 TSV + 库 FASTA

## 安装

```bash
pip install -r requirements.txt
```

需要 Pillow（`Pillow>=9.0`）。默认字体为 Liberation Serif，换字体用 `--font` / `--bold`。

## 用法

### 1. 内置示例（两条序列）

```bash
python3 disulfide_plot.py --demo --out demo.png
```

### 2. 从 CSV 画

```bash
python3 disulfide_plot.py -i examples/candidates.csv \
  -o candidates.png --title "1-4 / 2-5 / 3-6 connectivity"
```

同时出 SVG：

```bash
python3 disulfide_plot.py -i examples/candidates.csv \
  -o candidates.png --svg candidates.svg
```

### 3. 直接读 BLAST 结果

```bash
python3 disulfide_plot.py \
  --blast /path/to/blast_short.tsv \
  --fasta /path/to/library.fasta \
  --top 20 --cys 6 -o blast_top20.png
```

行标签形如 `u24304 (80.0%/26%)`，括号内是 pident 和 qcovhsp。
只保留 Cys 数等于 `--cys` 的序列，按 `(pident, bitscore)` 降序取前 `--top` 条。

## CSV 格式

```csv
name,sequence,pairs,note
Query,CKGODAPCSKSKDCCGEASMCSKGADGKKTCMIDKLWG*,1-4;2-5;3-6,query
Pu6.9,SCGEEGEGCYTRPCCSLGLKCIGIRSHGGLCREE,1-4;2-5;3-6,
```

- 必须有 `sequence` 列；`name`、`pairs`、`note` 可选
- `pairs` 按 Cys 出现顺序编号，1-based，可写 `1-4;2-5;3-6` 或 `1-4,2-5,3-6`
- 省略 `pairs`：6 个 Cys 默认 `1-4,2-5,3-6`；其他 Cys 数默认首尾配对
- 序列可含 `*`、`O`、`X`，绘图原样显示
- Cys 少于 2 个的行会跳过并打印 `[skip]`

## 参数

| 参数 | 说明 |
|---|---|
| `-i, --input` | 输入 CSV |
| `--blast` | BLAST outfmt 6 TSV（15 列） |
| `--fasta` | 与 `--blast` 配套的库 FASTA |
| `--top` | `--blast` 取前 N 条，默认 20 |
| `--cys` | `--blast` 只保留 N 个 Cys，默认 6 |
| `--pairs` | 覆盖配对，如 `1-4;2-5;3-6` |
| `-o, --out` | 输出 PNG，默认 `disulfide.png` |
| `--svg` | 同时输出 SVG |
| `--title` | 标题 |
| `--font-size` | 字号，默认 34 |
| `--gap` | 轨道间距，默认 15 |
| `--tick` | 竖线额外长度，默认 5 |
| `--row-pad` | 行间留白，默认 14 |
| `--scale` | 超采样倍数，默认 3 |
| `--no-legend` | 不画图例 |
| `--font` / `--bold` | 替换字体 |

## 画法规则

1. 键按端点跨度 `|j-i|` 降序，最大者放 `above` track 0
2. 其余按 `below0 → below1 → above1 → above2 → below2 …` 分轨
3. 竖线只到文字上/下缘，不穿过序列
4. Cys 红色 `#d62728`，其余 `#111111`
5. 键颜色按配对顺序：蓝 / 绿 / 紫 / 橙 / 青 / 洋红

## BLAST 模式

- 读 outfmt 6 第 2 列 `sseqid`、第 3 列 `pident`、第 12 列 `bitscore`、第 15 列 `qcovhsp`
- `pident` 是 HSP 局部一致性，不是全长 identity
- `qcovhsp` 是 query 覆盖度，不是 subject 覆盖度
- 同一 `sseqid` 保留 `(pident, bitscore)` 最大的一条

## 限制

- 序列本身不编码二硫键配对；图按 `pairs` 画，推断配对不代表实验验证
- BLAST 命中常是 8–10 aa 的 HSP，图中画的是整条 subject 的 Cys 配对
- `blastp-short` 的 pident 可能高估约 1.7 个百分点

## 文件结构

```
disulfide-plot/
├── disulfide_plot.py
├── README.md
├── README.en.md
├── requirements.txt
├── assets/
│   ├── logo.svg
│   └── wordmark.svg
└── examples/
    ├── candidates.csv
    ├── candidates.png
    ├── candidates.svg
    ├── blast_top12.png
    └── demo.png
```

## 常见问题

**找不到字体**：默认 Liberation Serif。用 `fc-list | grep -i liberation` 查路径，或用 `--font` 指定。

**中文标题显示方块**：换含中文的字体，如 Noto Sans CJK。

**Cys 数不足**：`pairs` 的编号不能超过序列 Cys 数，否则报 `配对 3-6 超出 Cys 数 4`。

**图太宽**：调小 `--font-size` 或 `--scale`；序列不换行。

## License

MIT
