// 06_build_report.js — builds outputs/Data_Collection_and_QC_Report.docx from harmonized/report_data.json
const fs = require("fs");
const path = require("path");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, ShadingType, HeadingLevel,
        AlignmentType, LevelFormat, BorderStyle, Footer, PageNumber } = require("docx");

const B = path.resolve(__dirname, "..");
const D = JSON.parse(fs.readFileSync(path.join(B, "harmonized", "report_data.json"), "utf8"));
const S = D.summ;
const FONT = "Times New Roman";
const W = 9026; // A4 text width in DXA with 1-inch margins

const p = (text, opts = {}) => new Paragraph({ spacing: { after: 120, line: 276 }, ...opts,
  children: (Array.isArray(text) ? text : [text]).map(t => typeof t === "string" ? new TextRun({ text: t, font: FONT, size: 22 }) : t) });
const b = (t) => new TextRun({ text: t, bold: true, font: FONT, size: 22 });
const h1 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_1, keepNext: true, spacing: { before: 280, after: 120 }, children: [new TextRun({ text: t, font: FONT })] });
const h2 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_2, keepNext: true, spacing: { before: 200, after: 100 }, children: [new TextRun({ text: t, font: FONT })] });
const bullet = (t) => new Paragraph({ numbering: { reference: "bul", level: 0 }, spacing: { after: 60, line: 264 },
  children: (Array.isArray(t) ? t : [t]).map(x => typeof x === "string" ? new TextRun({ text: x, font: FONT, size: 22 }) : x) });
const note = (t) => new Paragraph({ spacing: { after: 160 }, children: [new TextRun({ text: t, font: FONT, size: 18, italics: true })] });
const caption = (t) => new Paragraph({ spacing: { before: 120, after: 60 }, keepNext: true, children: [new TextRun({ text: t, font: FONT, size: 20, bold: true })] });

function table(rows, widths) {
  const total = widths.reduce((a, c) => a + c, 0);
  const scaled = widths.map(w => Math.floor(w * W / total));
  scaled[scaled.length - 1] += W - scaled.reduce((a, c) => a + c, 0);
  const border = { style: BorderStyle.SINGLE, size: 4, color: "BFBFBF" };
  const borders = { top: border, bottom: border, left: border, right: border };
  return new Table({
    width: { size: W, type: WidthType.DXA }, columnWidths: scaled,
    rows: rows.map((r, i) => new TableRow({ tableHeader: i === 0, cantSplit: true,
      children: r.map((c, j) => new TableCell({ width: { size: scaled[j], type: WidthType.DXA }, borders,
        margins: { top: 40, bottom: 40, left: 80, right: 80 },
        shading: i === 0 ? { fill: "1F4E78", type: ShadingType.CLEAR, color: "auto" } : (i % 2 === 0 ? { fill: "F2F2F2", type: ShadingType.CLEAR, color: "auto" } : undefined),
        children: [new Paragraph({ children: [new TextRun({ text: String(c), font: FONT, size: 18, bold: i === 0, color: i === 0 ? "FFFFFF" : "000000" })] })] })) }))
  });
}

const children = [
  new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 80 }, children: [new TextRun({ text: "Public Metagenomic Data Collection and Quality-Control Report", font: FONT, size: 34, bold: true, color: "1F4E78" })] }),
  p([new TextRun({ text: "Book chapter: Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping (Academic Press, Elsevier)", font: FONT, size: 22, italics: true })]),
  p("Collection date: 31 October 2026. Prepared by an AI assistant (Claude) for author review. Every count below was calculated from the collected records by the scripts in the package; no sequence reads were downloaded and no statistical analysis was run."),

  h1("Summary"),
  p(`Six public repositories were queried with 236 logged queries, yielding ${S["Records retrieved and screened"]} screened records and ${S["Distinct candidate records after deduplication"]} distinct candidates. All ${S["INSDC candidate studies verified against ENA"]} INSDC candidate studies were checked against their official ENA run records. ${S["Studies/projects included"]} studies or projects were included, ${S["Studies/projects as separate contextual datasets"]} kept as separate contextual (food-chain and related) datasets, and ${S["Studies/projects pending verification"]} left pending.`),
  p(`The harmonized master holds ${S["Unique samples in sample master"]} unique repository samples, of which ${S["Included unique samples"]} are included public shotgun metagenomes (${S["Included samples (DS01 counted by biological sample)"]} when the global sewage samples are counted by biological sample). They come from ${S["Countries represented (included samples)"]} countries and ${S["Cities represented (distinct City|Country pairs, source-labelled city, country known)"]} source-labelled cities (counted as city–country pairs), with ${S["Additional locality strings not standardized to cities"]} further locality strings not yet standardized. Gene-level AMR tables exist for ${S["Included samples with gene-level AMR table"]} included samples and a phylum-level taxonomic table for ${S["Included samples with taxonomy table (phylum)"]}; every other sample would need reprocessing from raw reads.`),

  h1("Part A. Verified findings"),
  h2("1. Repositories investigated and queries executed"),
  p("Each database received the same 35 search concepts in six query groups (urban microbiome; urban AMR and resistome; wastewater and sewage; environmental microbiomes; food environments; geographic searches including Dhaka, Delhi and Karachi). Exact strings, URLs, dates, filters and counts are in the Search_Query_Log sheet and Public_Data_Source_Log.csv."),
  caption("Table 1. Records returned, screened and retained by database"),
  table(D.qs, [24, 9, 11, 11, 11, 11, 11]),
  note("Relevant = passed the documented programmatic title/metadata screen. Verified = the official study or project record was retrieved. Retained = included or separate contextual after review. Europe PMC and Zenodo records are leads only and were not verified individually."),
  p("Two repository behaviours shaped the results. The MGnify API v2 'search' parameter matched each query as an exact phrase, so all 35 multi-word queries returned zero; 26 adapted single-keyword and biome-lineage queries were added and logged separately. Zenodo combines words with OR by default, so its queries were rerun with all words required."),

  h2("2. Studies and accessions identified"),
  caption("Table 2. Deduplicated candidate records by type and decision"),
  table(D.invs, [30, 14, 14, 18, 24]),
  note("ERP015409 was merged with PRJEB13831 (the same study under its secondary accession). 105 INSDC candidates had no public runs and 223 had no shotgun metagenomic runs."),

  h2("3. Datasets verified"),
  caption("Table 3. Core datasets named in the brief"),
  table(D.core, [9, 27, 20, 9, 11, 7, 17]),
  p("A further 106 INSDC studies and three MG-RAST series (GLUSEEN urban soils, Swiss wastewater and river bacteriomes, Polish WWTP resistome series) were included from the discovery searches. Their classification rests on title review and ENA records, and their sample metadata have not been checked against publications."),

  h2("4. Data actually downloaded"),
  bullet("Zenodo 10.5281/zenodo.17590646 (CC BY 4.0): Indian wastewater processed archive, 318 MB; MD5 matches the Zenodo record. Contains Bracken phylum/genus/species tables, the RGI (CARD) gene table and sample metadata."),
  bullet("Zenodo 10.5281/zenodo.6919377 (CC BY 4.0; Martiny et al.): metadata, ResFinder annotation, diversity and ARG tables (ARG table 445 MB). Rows were extracted for the 1,004 runs in the selected studies that the resource covers (825 of them in included samples)."),
  bullet("Munk et al. 2022 supplementary package and full text (Europe PMC); Singh et al. 2026 full text; MetaSUB metadata (GitHub commit 1d738748); FARME DB files (≈18 MB); CARD 4.0.2 data archive (4.4 MB)."),
  bullet("Metadata only: ENA run reports for 593 studies, extended run/sample metadata for 132 studies, NCBI BioSample attributes for 2,139 samples and MG-RAST records for 585 metagenomes."),
  bullet([b("Not downloaded: "), new TextRun({ text: "raw FASTQ for any dataset (about 24.1 TB compressed across included samples), the 6.4 GB rRNA count table, MGnify analyses, and publisher source-data files (HTTP 403 from the workspace).", font: FONT, size: 22 })]),

  h2("5. Exclusions"),
  caption("Table 4. Samples excluded within included studies"),
  table(D.ex, [70, 30]),
  p("At study level, 437 INSDC and 29 MG-RAST candidates were excluded: no public shotgun runs, host-associated, agricultural, marine or remote settings, process-engineering or method studies, plastisphere or enrichment libraries. Each has its reason in the Study_Inventory and Exclusion_Log sheets."),

  h2("6. Cities, countries and South Asian coverage"),
  p(`The proposal's 35 cities are a target, not a dataset property. Verified coverage is broader but uneven: the global sewage dataset alone holds 242 distinct city labels (243 city–country pairs) in 101 countries, consistent with the 243 cities and 101 countries reported by Munk et al. The three named South Asian cities are present: ${S["Proposal target cities present (Dhaka / Delhi / Karachi)"]} included samples. Karachi appears only in the global sewage dataset, and MetaSUB contains no South Asian city.`),
  caption("Table 5. Included samples from South Asian countries"),
  table(D.sa, [40, 30, 30]),

  h2("7. Environmental sample types"),
  caption("Table 6. Included samples by standardized sample type"),
  table(D.st, [26, 13, 13, 14, 17, 17]),
  note("Food-related studies (840 samples) are held separately and are not counted above."),

  h2("8. Taxonomic and AMR information"),
  bullet("Taxonomy: 447 Indian samples have Kraken 2/Bracken profiles (phylum in the workbook; genus and species in the source files). 625 samples have genus- and phylum-level diversity indices from the Zenodo 6919377 resource. 156 MG-RAST metagenomes have online annotations that were not downloaded."),
  bullet("AMR: 447 samples with RGI read mapping against CARD and 625 with KMA read mapping against ResFinder, together 316,983 gene-level rows. The two pipelines differ in database, naming and count definition and are kept separate. Counts are raw (not normalized for gene length or depth)."),
  bullet("Functional profiles: none downloaded; the worksheet is an empty, labelled template."),
  bullet("Reference resources: FARME DB is online but unchanged since March 2019 (20,724 DNA and 48,178 protein sequences from 29 functional-metagenomic projects; non-commercial reuse). CARD 4.0.2 (August 2026) is current; commercial use is restricted."),

  h2("9. Incomplete metadata"),
  caption("Table 7. Missing values among included samples"),
  table(D.mm, [40, 30, 30]),
  note("City is blank for samples whose location is a locality string rather than a source-labelled city; Locality_As_Reported holds that text. Environmental_Matrix (env_medium) is rarely submitted."),

  h2("10. Quality-control results"),
  caption("Table 8. Quality-control checks on the sample master"),
  table(D.qc, [14, 48, 14, 24]),
  p("The single failed check is a source error: 86 MetaSUB records place Kyiv in Poland, and 10 more Kyiv records have no country that could be taken from the same city without repeating the error. All 96 have country set to 'Not yet verified', with original values kept. Separately, 81 included samples also carry amplicon runs; these are listed apart and excluded from analysis fields. Other source problems handled the same way: 60 MetaSUB dates of 1970-01-01 (an epoch placeholder), one ENA date with year 0216, and the Kochi/Cochin spelling variants in the sewage dataset."),

  h2("11. Source-verification checks"),
  caption("Table 9. Claims checked against official records"),
  table(D.ver, [13, 24, 49, 14]),

  h2("12. Suitability for subsequent analysis"),
  caption("Table 10. Readiness of the proposed analyses"),
  table(D.rd, [22, 16, 62]),

  h1("Part B. Research objectives the data cannot currently support"),
  bullet("Cross-city or cross-study statistical comparison of diversity or resistomes: profiles come from different pipelines and most datasets have raw reads only."),
  bullet("Environmental association analyses: no pollution, climate, sanitation or antibiotic-use covariates have been joined to samples."),
  bullet("Predictive modelling of disease, outbreaks or transmission: no outcome data are matched to samples by place and time."),
  bullet("Food-system and food-business implications: the food datasets are heterogeneous, mostly not urban-located and carry no food-safety outcomes."),
  bullet("Spatial hotspot mapping: coordinates are often approximate or city-level; only observed-point maps are defensible now."),

  h1("Part C. Recommendations"),
  bullet("Agree a manageable raw-read subset (for example, DS01 plus DS02 plus the South Asian studies, stratified by matrix) and reprocess it uniformly with one taxonomic classifier and CARD 4.0.2, so that cross-study comparison becomes possible."),
  bullet("Download the 6.4 GB rRNA table from Zenodo 6919377 if genus-level taxonomy is needed for the 625 covered samples before reprocessing."),
  bullet("Retrieve manually the Munk et al. Source Data and the Singh et al. Supplementary Data S7 (site coordinates) and source data."),
  bullet("Build a documented lookup table standardizing the 625 locality strings to cities, reviewed by the authors."),
  bullet("Screen the 242 accession-linked publications and 178 Zenodo records, and resolve the 25 pending studies."),
  bullet("Assemble city-level covariates (population density, sanitation coverage, climate, antibiotic use) from authoritative sources before any association analysis."),

  h1("Part D. Unresolved issues requiring author decisions"),
  bullet("Confirmation of all study-level eligibility decisions (recorded in scripts/curation_decisions.py)."),
  bullet("Whether FARME DB is the proposal's 'FARME data management' resource."),
  bullet("The CARD release used by Singh et al. 2026 (not stated in any source accessed)."),
  bullet("The publication linked to the Dhaka project PRJNA1399457 (none found in repository records; the DOI in the brief was not verified)."),
  bullet("Approval of any large download or computing work."),
  bullet("Whether DS01 analyses should use biological samples (798 in Munk et al.) rather than repository accessions."),
  p([new TextRun({ text: "Responsible AI use: ", font: FONT, size: 22, bold: true }), new TextRun({ text: "this data package and report were prepared with AI assistance. The authors should verify the records and decisions before relying on them and disclose AI use as required by Elsevier policy.", font: FONT, size: 22 })]),
];

const doc = new Document({
  creator: "AI assistant (Claude) for the chapter authors",
  title: "Public Metagenomic Data Collection and QC Report",
  styles: { default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 28, bold: true, color: "1F4E78", font: FONT }, paragraph: { outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 24, bold: true, color: "2E75B6", font: FONT }, paragraph: { outlineLevel: 1 } }] },
  numbering: { config: [{ reference: "bul", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 270 } } } }] }] },
  sections: [{ properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: ["Page ", PageNumber.CURRENT], font: FONT, size: 18 })] })] }) },
    children }],
});
Packer.toBuffer(doc).then(buf => { fs.writeFileSync(path.join(B, "outputs", "Data_Collection_and_QC_Report.docx"), buf); console.log("written"); });
