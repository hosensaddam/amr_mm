// Render a JSON report specification to .docx (Times New Roman 12 pt, 1-inch margins, A4).
// Usage: node render_docx.js spec.json out.docx
// Block types: h1/h2/h3 {text}, p {text | runs:[{text,bold,italics}]}, bullets {items:[...]},
// table {columns:[...], rows:[[...]], widths:[dxa...], caption, note}, image {path, width_mm, caption}, note {text}, pagebreak
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell, WidthType,
  ShadingType, AlignmentType, ImageRun, PageBreak, LevelFormat, BorderStyle, Footer, PageNumber,
  TableOfContents,
} = require("docx");

const spec = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
const FONT = "Times New Roman";
const PAGE_W = 11906, MARGIN = 1440, TEXT_W = PAGE_W - 2 * MARGIN; // A4, 1 inch margins
const border = { style: BorderStyle.SINGLE, size: 4, color: "808080" };
const borders = { top: border, bottom: border, left: border, right: border };

function runs(b, size) {
  if (b.runs) return b.runs.map(r => new TextRun({ text: String(r.text), bold: !!r.bold, italics: !!r.italics, font: FONT, size: size || 24 }));
  return [new TextRun({ text: String(b.text), font: FONT, size: size || 24, bold: !!b.bold, italics: !!b.italics })];
}
function para(b) {
  return new Paragraph({ children: runs(b), spacing: { after: 120, line: spec.line_spacing || 360 }, alignment: AlignmentType.JUSTIFIED });
}
function caption(text) {
  return new Paragraph({ children: [new TextRun({ text, font: FONT, size: 20, italics: true })], spacing: { before: 60, after: 200 } });
}
function table(b) {
  const n = b.columns.length;
  let widths = b.widths;
  if (!widths) { const w = Math.floor(TEXT_W / n); widths = Array(n).fill(w); widths[n - 1] = TEXT_W - w * (n - 1); }
  const tw = widths.reduce((a, c) => a + c, 0);
  const fs_ = b.font_size || 18;
  const cell = (t, i, hdr) => new TableCell({
    borders, width: { size: widths[i], type: WidthType.DXA },
    shading: hdr ? { fill: "D9E2F3", type: ShadingType.CLEAR, color: "auto" } : undefined,
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    children: String(t ?? "").split("\n").map(line => new Paragraph({ children: [new TextRun({ text: line, font: FONT, size: fs_, bold: hdr })] })),
  });
  const out = [];
  if (b.caption) out.push(new Paragraph({ children: [new TextRun({ text: b.caption, font: FONT, size: 22, bold: true })], spacing: { before: 200, after: 80 }, keepNext: true }));
  out.push(new Table({
    width: { size: tw, type: WidthType.DXA }, columnWidths: widths,
    rows: [new TableRow({ tableHeader: true, children: b.columns.map((c, i) => cell(c, i, true)) }),
      ...b.rows.map(r => new TableRow({ children: r.map((c, i) => cell(c, i, false)) }))],
  }));
  if (b.note) out.push(caption(b.note)); else out.push(new Paragraph({ children: [], spacing: { after: 120 } }));
  return out;
}
function image(b) {
  const buf = fs.readFileSync(b.path);
  // PNG size from IHDR
  const w = buf.readUInt32BE(16), h = buf.readUInt32BE(20);
  const maxpx = Math.min((b.width_mm || 160) / 25.4 * 96, TEXT_W / 1440 * 96);
  const sc = maxpx / w;
  return [new Paragraph({ alignment: AlignmentType.CENTER, keepNext: true, children: [new ImageRun({ type: "png", data: buf, transformation: { width: Math.round(w * sc), height: Math.round(h * sc) }, altText: { title: b.caption || "figure", description: b.caption || "figure", name: path.basename(b.path) } })] }),
    caption(b.caption || "")];
}
// Markdown-style **bold** inside strings
function md(text, size, extra) {
  const parts = String(text).split(/(\*\*[^*]+\*\*)/g).filter(x => x.length);
  return parts.map(t => t.startsWith("**") ? new TextRun({ text: t.slice(2, -2), bold: true, font: FONT, size: size || 24, ...(extra || {}) }) : new TextRun({ text: t, font: FONT, size: size || 24, ...(extra || {}) }));
}
const children = [];
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 }, children: [new TextRun({ text: spec.title, font: FONT, size: 32, bold: true })] }));
if (spec.subtitle) children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 }, children: [new TextRun({ text: spec.subtitle, font: FONT, size: 24, italics: true })] }));
for (const m of spec.meta || []) children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 60 }, children: md(m, 22) }));
children.push(new Paragraph({ children: [], spacing: { after: 200 } }));
let tabN = 0;
for (const b of spec.blocks) {
  if (b.h1) children.push(new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun({ text: b.h1 })] }));
  else if (b.h2) children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, children: [new TextRun({ text: b.h2 })] }));
  else if (b.h3) children.push(new Paragraph({ heading: HeadingLevel.HEADING_3, children: [new TextRun({ text: b.h3 })] }));
  else if (b.p !== undefined) children.push(new Paragraph({ children: md(b.p), spacing: { after: 120, line: spec.line_spacing || 360 }, alignment: AlignmentType.JUSTIFIED }));
  else if (b.note) children.push(new Paragraph({ children: md(b.note, 22), spacing: { before: 120, after: 120 }, shading: { fill: "F2F2F2", type: ShadingType.CLEAR, color: "auto" }, border: { left: { style: BorderStyle.SINGLE, size: 12, color: "4472C4", space: 6 } } }));
  else if (b.bullets) for (const it of b.bullets) children.push(new Paragraph({ numbering: { reference: "bul", level: 0 }, spacing: { after: 60, line: spec.line_spacing || 360 }, children: md(it) }));
  else if (b.table) {
    const t = b.table; const n = t.columns.length;
    const fr = t.widths || Array(n).fill(1 / n); const sum = fr.reduce((a, c) => a + c, 0);
    let widths = fr.map(x => Math.floor(TEXT_W * x / sum)); widths[n - 1] = TEXT_W - widths.slice(0, -1).reduce((a, c) => a + c, 0);
    if (b.title) children.push(new Paragraph({ children: md(b.title, 22), spacing: { before: 200, after: 80 }, keepNext: true }));
    const fs_ = t.size || 18;
    const cell = (v, i, hdr) => new TableCell({ borders, width: { size: widths[i], type: WidthType.DXA },
      shading: hdr ? { fill: "D9E2F3", type: ShadingType.CLEAR, color: "auto" } : undefined, margins: { top: 30, bottom: 30, left: 60, right: 60 },
      children: String(v ?? "").split("\n").map(line => new Paragraph({ children: [new TextRun({ text: line, font: FONT, size: fs_, bold: hdr })] })) });
    children.push(new Table({ width: { size: TEXT_W, type: WidthType.DXA }, columnWidths: widths,
      rows: [new TableRow({ tableHeader: true, children: t.columns.map((c, i) => cell(c, i, true)) }), ...t.rows.map(r => new TableRow({ cantSplit: true, children: r.map((c, i) => cell(c, i, false)) }))] }));
    children.push(new Paragraph({ children: [], spacing: { after: 160 } }));
  }
  else if (b.image) {
    const buf = fs.readFileSync(b.image); const w = buf.readUInt32BE(16), h = buf.readUInt32BE(20);
    const maxpx = TEXT_W / 1440 * 96; const sc = Math.min(maxpx, b.width_px || maxpx) / w;
    children.push(new Paragraph({ alignment: AlignmentType.CENTER, keepNext: true, spacing: { before: 120 }, children: [new ImageRun({ type: "png", data: buf, transformation: { width: Math.round(w * sc), height: Math.round(h * sc) }, altText: { title: "figure", description: String(b.caption).slice(0, 200), name: path.basename(b.image) } })] }));
    children.push(new Paragraph({ children: md(b.caption || "", 20), spacing: { before: 60, after: 240 }, alignment: AlignmentType.JUSTIFIED }));
  }
  else if (b.pagebreak) children.push(new Paragraph({ children: [new PageBreak()] }));
}
const doc = new Document({
  creator: spec.creator || "Updated_Research_Analysis", title: spec.title,
  styles: {
    default: { document: { run: { font: FONT, size: 24 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 28, bold: true, font: FONT }, paragraph: { spacing: { before: 300, after: 120 }, outlineLevel: 0, keepNext: true } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 26, bold: true, font: FONT }, paragraph: { spacing: { before: 220, after: 100 }, outlineLevel: 1, keepNext: true } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 24, bold: true, italics: true, font: FONT }, paragraph: { spacing: { before: 160, after: 80 }, outlineLevel: 2, keepNext: true } },
    ],
  },
  numbering: { config: [{ reference: "bul", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 270 } } } }] }] },
  sections: [{
    properties: { page: { size: { width: PAGE_W, height: 16838 }, margin: { top: MARGIN, right: MARGIN, bottom: MARGIN, left: MARGIN } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ font: FONT, size: 18, children: [spec.footer || "", "  |  Page ", PageNumber.CURRENT] })] })] }) },
    children,
  }],
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(process.argv[3], b); console.log("written", process.argv[3]); });
