const fs = require("fs");
const path = require("path");
const {
  AlignmentType,
  BorderStyle,
  Document,
  Footer,
  Header,
  ImageRun,
  Packer,
  PageBreak,
  PageNumber,
  Paragraph,
  ShadingType,
  Table,
  TableCell,
  TableRow,
  TextRun,
  VerticalAlign,
  WidthType,
} = require("docx");

const ROOT = __dirname;
const OUTPUT = path.join(ROOT, "Manuscript_Q1_rewrite.docx");
const SUPP_OUTPUT = path.join(ROOT, "Supplementary_Information_Q1_rewrite.docx");

const A4_WIDTH = 11906;
const A4_HEIGHT = 16838;
const LEFT_RIGHT = 1008;
const CONTENT_WIDTH = A4_WIDTH - 2 * LEFT_RIGHT;
const border = { style: BorderStyle.SINGLE, size: 1, color: "B8C2CC" };
const borders = { top: border, bottom: border, left: border, right: border, insideHorizontal: border, insideVertical: border };

function clean(text) {
  return text.replace(/\*\*/g, "").replace(/`/g, "").trim();
}

function runs(text, opts = {}) {
  const chunks = [];
  const re = /\*\*(.+?)\*\*/g;
  let cursor = 0;
  let match;
  while ((match = re.exec(text)) !== null) {
    if (match.index > cursor) chunks.push(new TextRun({ text: text.slice(cursor, match.index), ...opts }));
    chunks.push(new TextRun({ text: match[1], bold: true, ...opts }));
    cursor = re.lastIndex;
  }
  if (cursor < text.length) chunks.push(new TextRun({ text: text.slice(cursor), ...opts }));
  return chunks.length ? chunks : [new TextRun({ text, ...opts })];
}

function para(text, options = {}) {
  return new Paragraph({
    spacing: { after: 110, line: 276 },
    ...options,
    children: runs(text, { font: "Arial", size: 20 }),
  });
}

function heading(text, level) {
  const sizes = { 1: 28, 2: 24, 3: 22 };
  return new Paragraph({
    spacing: { before: level === 1 ? 260 : 190, after: 100 },
    pageBreakBefore: text === "Tables" || text === "Figure legends" || text === "References",
    children: [new TextRun({ text, bold: true, font: "Arial", size: sizes[level] || 22 })],
  });
}

function markdownTable(lines) {
  const body = lines
    .filter((line, index) => index !== 1)
    .map((line) => line.split("|").slice(1, -1).map((cell) => clean(cell)));
  const columns = body[0].length;
  const widths = Array(columns).fill(Math.floor(CONTENT_WIDTH / columns));
  widths[columns - 1] += CONTENT_WIDTH - widths.reduce((a, b) => a + b, 0);
  return new Table({
    width: { size: CONTENT_WIDTH, type: WidthType.DXA },
    columnWidths: widths,
    borders,
    rows: body.map((row, rowIndex) => new TableRow({
      children: row.map((cell, index) => new TableCell({
        width: { size: widths[index], type: WidthType.DXA },
        borders,
        verticalAlign: VerticalAlign.TOP,
        margins: { top: 70, bottom: 70, left: 80, right: 80 },
        shading: rowIndex === 0 ? { fill: "E8EEF4", type: ShadingType.CLEAR } : undefined,
        children: [new Paragraph({
          spacing: { after: 0, line: 190 },
          children: [new TextRun({ text: cell, font: "Arial", size: 14, bold: rowIndex === 0 })],
        })],
      })),
    })),
  });
}

function imageParagraph(relativeImage, caption) {
  const imagePath = path.join(ROOT, relativeImage.replaceAll("/", path.sep));
  if (!fs.existsSync(imagePath)) throw new Error(`Missing figure: ${imagePath}`);
  const ext = path.extname(imagePath).slice(1).toLowerCase();
  const buffer = fs.readFileSync(imagePath);
  if (ext !== "png") throw new Error(`Only PNG figures are supported by this builder: ${imagePath}`);
  const imageWidth = buffer.readUInt32BE(16);
  const imageHeight = buffer.readUInt32BE(20);
  const outputWidth = 620;
  const outputHeight = Math.round(outputWidth * imageHeight / imageWidth);
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 80, after: 180 },
    children: [new ImageRun({
      type: ext === "jpg" ? "jpeg" : ext,
      data: buffer,
      transformation: { width: outputWidth, height: outputHeight },
      altText: { title: caption, description: caption, name: path.basename(imagePath) },
    })],
  });
}

function parseMarkdown(markdownPath, isSupplement = false) {
  const lines = fs.readFileSync(markdownPath, "utf8").replace(/\r/g, "").split("\n");
  const children = [];
  let i = 0;
  let firstTitle = true;
  while (i < lines.length) {
    const line = lines[i];
    if (!line.trim()) {
      i += 1;
      continue;
    }
    if (line.startsWith("![")) {
      const match = line.match(/^!\[(.*?)\]\((.*?)\)$/);
      if (match) children.push(imageParagraph(match[2], match[1]));
      i += 1;
      continue;
    }
    if (line.startsWith("|")) {
      const tableLines = [];
      while (i < lines.length && lines[i].startsWith("|")) {
        tableLines.push(lines[i]);
        i += 1;
      }
      children.push(markdownTable(tableLines));
      children.push(new Paragraph({ spacing: { after: 80 }, children: [] }));
      continue;
    }
    const headingMatch = line.match(/^(#{1,3})\s+(.*)$/);
    if (headingMatch) {
      const level = headingMatch[1].length;
      const title = clean(headingMatch[2]);
      if (level === 1 && firstTitle) {
        children.push(new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 240 },
          children: [new TextRun({ text: title, bold: true, font: "Arial", size: 34 })],
        }));
        firstTitle = false;
      } else {
        children.push(heading(title, level));
      }
      i += 1;
      continue;
    }
    if (line.startsWith("**Figure ") || line.startsWith("**Supplementary Figure")) {
      children.push(new Paragraph({ children: [new PageBreak()] }));
      children.push(new Paragraph({
        spacing: { before: 80, after: 100, line: 250 },
        children: runs(line, { font: "Arial", size: 19 }),
      }));
      i += 1;
      continue;
    }
    const paragraphLines = [line.trim()];
    i += 1;
    while (i < lines.length && lines[i].trim() && !lines[i].startsWith("#") && !lines[i].startsWith("|") && !lines[i].startsWith("![")) {
      paragraphLines.push(lines[i].trim());
      i += 1;
    }
    const text = paragraphLines.join(" ");
    const metadata = text.startsWith("**Article type:") || text.startsWith("**Corresponding author:") || text.startsWith("**Running title:") || text.startsWith("**Keywords:");
    children.push(para(text, metadata ? { alignment: AlignmentType.LEFT, spacing: { after: 80, line: 240 } } : {}));
  }
  return children;
}

function makeDocument(markdownName, outputPath, isSupplement = false) {
  const children = parseMarkdown(path.join(ROOT, markdownName), isSupplement);
  const headerText = isSupplement ? "Supplementary Information" : "Chromosome 16 gain in ETV6::RUNX1 ALL";
  const doc = new Document({
    creator: "Shanghai Children's Hospital",
    title: isSupplement ? "Supplementary Information" : "Chromosome 16 gain in ETV6::RUNX1-positive childhood acute lymphoblastic leukemia",
    description: "Journal-neutral Q1-ready manuscript package",
    styles: {
      default: { document: { run: { font: "Arial", size: 20 }, paragraph: { spacing: { after: 110, line: 276 } } } },
    },
    sections: [{
      properties: {
        page: { size: { width: A4_WIDTH, height: A4_HEIGHT }, margin: { top: 1008, right: LEFT_RIGHT, bottom: 1008, left: LEFT_RIGHT } },
      },
      headers: {
        default: new Header({ children: [new Paragraph({
          border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: "7A8793", space: 1 } },
          children: [new TextRun({ text: headerText, font: "Arial", size: 16, color: "4E5963" })],
        })] }),
      },
      footers: {
        default: new Footer({ children: [new Paragraph({
          alignment: AlignmentType.RIGHT,
          children: [new TextRun({ text: "Page ", font: "Arial", size: 16, color: "4E5963" }), new TextRun({ children: [PageNumber.CURRENT], font: "Arial", size: 16, color: "4E5963" })],
        })] }),
      },
      children,
    }],
  });
  return Packer.toBuffer(doc).then((buffer) => fs.writeFileSync(outputPath, buffer));
}

Promise.all([
  makeDocument("Manuscript_Q1_rewrite.md", OUTPUT, false),
  makeDocument("Supplementary_Information_Q1_rewrite.md", SUPP_OUTPUT, true),
]).catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
