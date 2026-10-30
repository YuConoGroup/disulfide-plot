#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
disulfide_plot.py — 把肽序列的二硫键连接画成“彩色括号”图。

只依赖 Pillow（PIL）。不需要 matplotlib。
输入 CSV：name,sequence,pairs
  - pairs 可写 "1-4;2-5;3-6" 或 "1-4,2-5,3-6"（1-based，按 Cys 出现顺序编号）
  - 省略 pairs 时默认 1-4,2-5,3-6（6 个 Cys 的 VI/VII 常见配对）
输出：PNG（可选 SVG）

示例：
  python3 disulfide_plot.py --demo --out demo.png
  python3 disulfide_plot.py -i candidates.csv -o fig.png --title "BLAST candidates"
"""
from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont

# ---------------- 默认样式 ----------------
DEFAULT_FONT = "/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf"
DEFAULT_BOLD = "/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf"
CYS_COLOR = "#d62728"
TEXT_COLOR = "#111111"
BG = "#ffffff"
PALETTE = ["#1f4e9c", "#2e8b57", "#7b3fa0", "#d2691e", "#008b8b", "#c71585", "#8b4513", "#4682b4"]


@dataclass
class Record:
    name: str
    sequence: str
    pairs: list[tuple[int, int]]
    note: str = ""


def parse_pairs(s: str | None, n_cys: int) -> list[tuple[int, int]]:
    """把 '1-4;2-5;3-6' 解析为 0-based Cys 索引对。"""
    if not s or not str(s).strip():
        if n_cys == 6:
            return [(0, 3), (1, 4), (2, 5)]
        # 非 6 Cys：默认首尾配
        return [(i, n_cys - 1 - i) for i in range(n_cys // 2)]
    pairs: list[tuple[int, int]] = []
    for tok in re.split(r"[;,|]", str(s)):
        tok = tok.strip()
        if not tok:
            continue
        m = re.fullmatch(r"(\d+)\s*[-–—:]\s*(\d+)", tok)
        if not m:
            raise ValueError(f"无法解析配对: {tok!r}（应形如 1-4）")
        a, b = int(m.group(1)) - 1, int(m.group(2)) - 1
        if not (0 <= a < n_cys and 0 <= b < n_cys):
            raise ValueError(f"配对 {tok} 超出 Cys 数 {n_cys}")
        if a == b:
            raise ValueError(f"配对 {tok} 两端相同")
        pairs.append((a, b))
    return pairs


def read_csv(path: str) -> list[Record]:
    records: list[Record] = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError("CSV 没有表头")
        # 容忍列名大小写/空格
        field = {k.strip().lower(): k for k in reader.fieldnames}
        seq_key = field.get("sequence") or field.get("seq") or field.get("peptide")
        name_key = field.get("name") or field.get("id") or field.get("header")
        pair_key = field.get("pairs") or field.get("disulfides") or field.get("ss")
        note_key = field.get("note") or field.get("source")
        if not seq_key:
            raise ValueError("CSV 必须有 sequence 列")
        for i, row in enumerate(reader, 2):
            seq = (row.get(seq_key) or "").strip().upper()
            if not seq:
                continue
            n_cys = seq.count("C")
            if n_cys < 2:
                print(f"[skip] 第 {i} 行 {seq!r}: Cys < 2")
                continue
            name = (row.get(name_key) or f"seq{i}").strip() if name_key else f"seq{i}"
            pairs = parse_pairs(row.get(pair_key) if pair_key else None, n_cys)
            note = (row.get(note_key) or "").strip() if note_key else ""
            records.append(Record(name, seq, pairs, note))
    return records


def records_from_blast(
    tsv_path: str,
    fasta_path: str,
    *,
    top: int = 20,
    cys: int = 6,
    pairs: str | None = None,
) -> list[Record]:
    """从 BLAST outfmt 6 TSV + 全库 FASTA 自动生成绘图记录。

    TSV 需包含：sseqid, pident, bitscore, qcovhsp（C-004/C-006 的 15 列格式）。
    按 pident 降序、bitscore 降序取前 top 条；只保留 Cys 数 == cys 的序列。
    """
    # 读 FASTA：header 第一个 token 作为 key
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()

    rows: dict[str, tuple[float, float, str]] = {}
    with open(tsv_path, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            sseqid, pident, bits, qcov = parts[1], parts[2], parts[11], parts[14]
            try:
                p, b = float(pident), float(bits)
            except ValueError:
                continue
            # 同名取最佳
            if sseqid not in rows or (p, b) > (rows[sseqid][0], rows[sseqid][1]):
                rows[sseqid] = (p, b, qcov)

    ordered = sorted(rows.items(), key=lambda kv: (kv[1][0], kv[1][1]), reverse=True)
    out: list[Record] = []
    for sseqid, (pident, bits, qcov) in ordered:
        seq = seqs.get(sseqid)
        if not seq:
            # sseqid 可能带描述，尝试用第一段
            seq = seqs.get(sseqid.split("|")[0])
        if not seq:
            continue
        if seq.count("C") != cys:
            continue
        short = sseqid.split("|")[0]
        rec_pairs = parse_pairs(pairs, cys)
        out.append(Record(f"{short} ({pident:.1f}%/{qcov}%)", seq, rec_pairs, sseqid))
        if len(out) >= top:
            break
    return out


def demo_records() -> list[Record]:
    data = [
        ("Query", "CKGODAPCSKSKDCCGEASMCSKGADGKKTCMIDKLWG*", ""),
        ("Pu6.9", "SCGEEGEGCYTRPCCSLGLKCIGIRSHGGLCREE", ""),
    ]
    return [Record(n, s, parse_pairs(p, s.count("C")), "") for n, s, p in data]


def assign_slots(bonds: list[tuple[int, int, int]], n_above: int, n_below: int):
    """
    给每条键分配 (side, track)。规则：
      - 跨度最大者放 above track 0；跨度相同则右侧端点更靠后者优先
        （3 条键时正好是 3-6 在上，1-4/2-5 在下两层，与参考图一致）；
      - 其余按 below0, below1, above1, above2, below2... 展开。
    bonds 元素为 (left_idx, right_idx, original_pair_index)。
    """
    ordered = sorted(bonds, key=lambda b: (abs(b[1] - b[0]), b[1]), reverse=True)
    slots = []
    ai = bi = 0
    for idx, _ in enumerate(ordered):
        if idx == 0:
            side, track = "above", ai
            ai += 1
        else:
            k = (idx - 1) % 4
            if k in (0, 1):
                side, track = "below", bi
                bi += 1
            else:
                side, track = "above", ai
                ai += 1
        slots.append((side, track))
    return list(zip(ordered, slots)), ai, bi


def draw_record(
    draw: ImageDraw.ImageDraw,
    rec: Record,
    x0: float,
    baseline: float,
    font: ImageFont.FreeTypeFont,
    size: int,
    scale: float,
    gap: float,
    tick: float,
    above_tracks: int,
    below_tracks: int,
) -> None:
    """在给定 baseline 上画一行序列 + 该行的二硫键括号。"""
    # 每个字符的 x 起点与宽度
    xs: list[float] = []
    widths: list[float] = []
    x = x0
    for ch in rec.sequence:
        w = font.getlength(ch)
        xs.append(x)
        widths.append(w)
        x += w

    cys_pos = [i for i, ch in enumerate(rec.sequence) if ch == "C"]
    # 端点：Cys 字符中心
    cx = {idx: xs[idx] + widths[idx] / 2 for idx in cys_pos}

    # 文本上下边界
    ascent, descent = font.getmetrics()
    text_top = baseline - ascent
    text_bottom = baseline + descent

    # 先画键（在文字下面层），再画文字，避免线压字
    bonds = []
    for k, (a, b) in enumerate(rec.pairs):
        if a < len(cys_pos) and b < len(cys_pos):
            i, j = cys_pos[a], cys_pos[b]
            bonds.append((min(i, j), max(i, j), k))
    slots, _, _ = assign_slots(bonds, above_tracks, below_tracks)
    for (i, j, k), (side, track) in slots:
        color = PALETTE[k % len(PALETTE)]
        x1, x2 = cx[i], cx[j]
        if x1 > x2:
            x1, x2 = x2, x1
        if side == "above":
            y_line = text_top - (track + 1) * gap - tick
            y_text = text_top - 2 * scale
        else:
            y_line = text_bottom + (track + 1) * gap + tick
            y_text = text_bottom + 2 * scale
        lw = max(2, int(round(1.6 * scale)))
        # 水平线 + 两端短竖线（括号）
        draw.line([(x1, y_line), (x2, y_line)], fill=color, width=lw)
        draw.line([(x1, y_line), (x1, y_text)], fill=color, width=lw)
        draw.line([(x2, y_line), (x2, y_text)], fill=color, width=lw)

    # 文字：Cys 红色
    for i, ch in enumerate(rec.sequence):
        color = CYS_COLOR if ch == "C" else TEXT_COLOR
        draw.text((xs[i], baseline), ch, font=font, fill=color, anchor="ls")


def render(
    records: list[Record],
    out_png: str,
    *,
    font_path: str = DEFAULT_FONT,
    bold_path: str = DEFAULT_BOLD,
    font_size: int = 34,
    title: str | None = None,
    legend: bool = True,
    margin: int = 18,
    gap: int = 15,
    tick: int = 5,
    row_pad: int = 14,
    scale: int = 3,
    dpi_meta: int = 300,
) -> tuple[int, int]:
    """渲染 PNG，返回 (宽, 高)。scale 为超采样倍数。"""
    font = ImageFont.truetype(font_path, font_size * scale)
    bold = ImageFont.truetype(bold_path, font_size * scale)
    small = ImageFont.truetype(font_path, max(12, int(font_size * 0.62)) * scale)
    small_bold = ImageFont.truetype(bold_path, max(12, int(font_size * 0.62)) * scale)

    # 前缀宽度（name:）
    prefix_w = max(font.getlength(r.name + ":") for r in records) if records else 0
    seq_x0 = margin * scale + prefix_w + 10 * scale

    # 每行的上下轨道数
    per_row = []
    for rec in records:
        cys_pos = [i for i, ch in enumerate(rec.sequence) if ch == "C"]
        bonds = []
        for k, (a, b) in enumerate(rec.pairs):
            if a < len(cys_pos) and b < len(cys_pos):
                i, j = cys_pos[a], cys_pos[b]
                bonds.append((min(i, j), max(i, j), k))
        _, na, nb = assign_slots(bonds, 99, 99)
        per_row.append((na, nb))

    ascent, descent = font.getmetrics()
    text_h = (ascent + descent)
    max_seq_w = max((sum(font.getlength(ch) for ch in r.sequence) for r in records), default=0)
    width = int(seq_x0 + max_seq_w + margin * scale)

    title_h = 0
    if title:
        tw = bold.getlength(title) * 1.0
        title_h = int((ascent + descent) * 0.9 + 8 * scale)

    legend_h = 0
    if legend:
        legend_h = int((ascent + descent) * 0.7 + 10 * scale)

    heights = []
    for na, nb in per_row:
        above = na * gap * scale + tick * scale + 8 * scale
        below = nb * gap * scale + tick * scale + 8 * scale
        heights.append(int(above + text_h + below + row_pad * scale))

    height = int(margin * 2 * scale + title_h + sum(heights) + legend_h)

    img = Image.new("RGB", (width, height), BG)
    draw = ImageDraw.Draw(img)

    y = margin * scale
    if title:
        draw.text((margin * scale, y + ascent), title, font=bold, fill=TEXT_COLOR, anchor="ls")
        y += title_h

    for rec, (na, nb), h in zip(records, per_row, heights):
        above = na * gap * scale + tick * scale + 8 * scale
        baseline = y + above + ascent
        # 名称
        draw.text((margin * scale, baseline), rec.name + ":", font=font, fill=TEXT_COLOR, anchor="ls")
        draw_record(draw, rec, seq_x0, baseline, font, font_size, scale, gap * scale, tick * scale, na, nb)
        y += h

    if legend:
        ly = y + (ascent + descent) * 0.4
        lx = margin * scale
        # 用第一条记录的配对做图例
        if records:
            rec = records[0]
            for idx, (a, b) in enumerate(rec.pairs):
                color = PALETTE[idx % len(PALETTE)]
                draw.line([(lx, ly), (lx + 26 * scale, ly)], fill=color, width=max(2, int(1.6 * scale)))
                draw.text((lx + 32 * scale, ly + ascent * 0.35), f"{a+1}-{b+1}", font=small, fill=TEXT_COLOR, anchor="ls")
                lx += 95 * scale

    # 下采样
    img = img.resize((max(1, width // scale), max(1, height // scale)), Image.LANCZOS)
    os.makedirs(os.path.dirname(os.path.abspath(out_png)) or ".", exist_ok=True)
    img.save(out_png, dpi=(dpi_meta, dpi_meta))
    return img.size


def render_svg(records: list[Record], out_svg: str, **kwargs) -> None:
    """可选的简单 SVG 导出（1x 坐标，便于 AI/Inkscape 编辑）。"""
    # 为简洁，SVG 复用 PNG 的布局计算但重画；不追求与 PNG 逐像素一致。
    from xml.sax.saxutils import escape

    font_size = kwargs.get("font_size", 34)
    margin = kwargs.get("margin", 18)
    gap = kwargs.get("gap", 15)
    tick = kwargs.get("tick", 5)
    row_pad = kwargs.get("row_pad", 14)
    title = kwargs.get("title")
    legend = kwargs.get("legend", True)

    # 用 PIL 字体量宽度（与 PNG 同字体）
    font = ImageFont.truetype(kwargs.get("font_path", DEFAULT_FONT), font_size)
    bold = ImageFont.truetype(kwargs.get("bold_path", DEFAULT_BOLD), font_size)
    prefix_w = max(font.getlength(r.name + ":") for r in records)
    seq_x0 = margin + prefix_w + 10
    max_seq_w = max(sum(font.getlength(ch) for ch in r.sequence) for r in records)
    width = int(seq_x0 + max_seq_w + margin)

    ascent, descent = font.getmetrics()
    text_h = ascent + descent
    rows = []
    for rec in records:
        cys_pos = [i for i, ch in enumerate(rec.sequence) if ch == "C"]
        bonds = []
        for k, (a, b) in enumerate(rec.pairs):
            if a < len(cys_pos) and b < len(cys_pos):
                i, j = cys_pos[a], cys_pos[b]
                bonds.append((min(i, j), max(i, j), k))
        _, na, nb = assign_slots(bonds, 99, 99)
        rows.append((rec, na, nb))
    title_h = (ascent + descent) * 0.9 + 8 if title else 0
    legend_h = (ascent + descent) * 0.7 + 10 if legend else 0
    heights = [na * gap + tick + 8 + text_h + nb * gap + tick + 8 + row_pad for _, na, nb in rows]
    height = int(margin * 2 + title_h + sum(heights) + legend_h)

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">']
    parts.append(f'<rect width="{width}" height="{height}" fill="{BG}"/>')
    y = margin
    if title:
        parts.append(f'<text x="{margin}" y="{y+ascent}" font-family="Liberation Serif" font-size="{font_size}" fill="{TEXT_COLOR}">{escape(title)}</text>')
        y += title_h
    for row_idx, (rec, na, nb) in enumerate(rows):
        above = na * gap + tick + 8
        baseline = y + above + ascent
        parts.append(f'<text x="{margin}" y="{baseline}" font-family="Liberation Serif" font-size="{font_size}" fill="{TEXT_COLOR}">{escape(rec.name)}:</text>')
        x = seq_x0
        xs = []
        for ch in rec.sequence:
            w = font.getlength(ch)
            xs.append(x)
            color = CYS_COLOR if ch == "C" else TEXT_COLOR
            parts.append(f'<text x="{x:.2f}" y="{baseline}" font-family="Liberation Serif" font-size="{font_size}" fill="{color}">{escape(ch)}</text>')
            x += w
        cys_pos = [i for i, ch in enumerate(rec.sequence) if ch == "C"]
        cx = {i: xs[i] + font.getlength(rec.sequence[i]) / 2 for i in cys_pos}
        bonds = []
        for k, (a, b) in enumerate(rec.pairs):
            if a < len(cys_pos) and b < len(cys_pos):
                i, j = cys_pos[a], cys_pos[b]
                bonds.append((min(i, j), max(i, j), k))
        slots, _, _ = assign_slots(bonds, 99, 99)
        text_top = baseline - ascent
        text_bottom = baseline + descent
        for (i, j, k), (side, track) in slots:
            color = PALETTE[k % len(PALETTE)]
            x1, x2 = cx[i], cx[j]
            if x1 > x2:
                x1, x2 = x2, x1
            if side == "above":
                y_line = text_top - (track + 1) * gap - tick
                y_text = text_top - 2
            else:
                y_line = text_bottom + (track + 1) * gap + tick
                y_text = text_bottom + 2
            parts.append(f'<path d="M{x1:.2f},{y_text:.2f} L{x1:.2f},{y_line:.2f} L{x2:.2f},{y_line:.2f} L{x2:.2f},{y_text:.2f}" fill="none" stroke="{color}" stroke-width="1.6"/>')
        y += heights[row_idx]
    if legend and records:
        ly = y + (ascent + descent) * 0.4
        lx = margin
        rec = records[0]
        for idx, (a, b) in enumerate(rec.pairs):
            color = PALETTE[idx % len(PALETTE)]
            parts.append(f'<line x1="{lx}" y1="{ly}" x2="{lx+26}" y2="{ly}" stroke="{color}" stroke-width="1.6"/>')
            parts.append(f'<text x="{lx+32}" y="{ly+ascent*0.35}" font-family="Liberation Serif" font-size="{max(12,int(font_size*0.62))}" fill="{TEXT_COLOR}">{a+1}-{b+1}</text>')
            lx += 95
    parts.append("</svg>")
    with open(out_svg, "w", encoding="utf-8") as f:
        f.write("\n".join(parts))


def main(argv: Iterable[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="画肽序列的二硫键连接（彩色括号）")
    ap.add_argument("-i", "--input", help="CSV: name,sequence,pairs[,note]")
    ap.add_argument("--blast", help="BLAST outfmt6 TSV（需 15 列）")
    ap.add_argument("--fasta", help="BLAST 库 FASTA（配合 --blast）")
    ap.add_argument("--top", type=int, default=20, help="--blast 时取前 N 条（默认 20）")
    ap.add_argument("--cys", type=int, default=6, help="--blast 时只保留 N 个 Cys 的序列")
    ap.add_argument("--pairs", default=None, help="配对，如 '1-4;2-5;3-6'；默认 6 Cys 用 1-4/2-5/3-6")
    ap.add_argument("-o", "--out", default="disulfide.png", help="输出 PNG")
    ap.add_argument("--svg", help="同时输出 SVG")
    ap.add_argument("--demo", action="store_true", help="使用内置示例")
    ap.add_argument("--title", default=None)
    ap.add_argument("--font", default=DEFAULT_FONT)
    ap.add_argument("--bold", default=DEFAULT_BOLD)
    ap.add_argument("--font-size", type=int, default=34)
    ap.add_argument("--margin", type=int, default=18)
    ap.add_argument("--gap", type=int, default=15)
    ap.add_argument("--tick", type=int, default=5)
    ap.add_argument("--row-pad", type=int, default=14)
    ap.add_argument("--no-legend", action="store_true")
    ap.add_argument("--scale", type=int, default=3, help="超采样倍数（默认 3）")
    args = ap.parse_args(argv)

    if args.demo:
        records = demo_records()
    elif args.blast:
        if not args.fasta:
            ap.error("--blast 需要同时给 --fasta")
        records = records_from_blast(
            args.blast, args.fasta, top=args.top, cys=args.cys, pairs=args.pairs
        )
    elif args.input:
        records = read_csv(args.input)
    else:
        ap.error("需要 --input 或 --demo")
    if not records:
        ap.error("没有可画的记录")

    size = render(
        records,
        args.out,
        font_path=args.font,
        bold_path=args.bold,
        font_size=args.font_size,
        title=args.title,
        legend=not args.no_legend,
        margin=args.margin,
        gap=args.gap,
        tick=args.tick,
        row_pad=args.row_pad,
        scale=args.scale,
    )
    print(f"[ok] {args.out}  {size[0]}x{size[1]} px")
    if args.svg:
        render_svg(
            records,
            args.svg,
            font_path=args.font,
            bold_path=args.bold,
            font_size=args.font_size,
            title=args.title,
            legend=not args.no_legend,
            margin=args.margin,
            gap=args.gap,
            tick=args.tick,
            row_pad=args.row_pad,
        )
        print(f"[ok] {args.svg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
