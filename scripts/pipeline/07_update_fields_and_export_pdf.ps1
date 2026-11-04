# Step 07 - Let Microsoft Word update every field (SEQ captions, REF cross-references, TOC,
# List of Figures, List of Tables, PAGE/NUMPAGES), save the DOCX and export a PDF that matches it.
# Usage (from this folder):  powershell -ExecutionPolicy Bypass -File 07_update_fields_and_export_pdf.ps1
$ErrorActionPreference = "Stop"
$docs = Join-Path (Split-Path (Split-Path $PSScriptRoot)) "06_Research_Documents"
$docx = Join-Path $docs "Supplementary_Material_Final.docx"
$pdf  = Join-Path $docs "Supplementary_Material_Final.pdf"
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
try {
    $doc = $word.Documents.Open($docx, $false, $false, $false)
    Write-Output ("{0:T}  opened" -f (Get-Date))
    $doc.Repaginate()
    for ($pass = 1; $pass -le 2; $pass++) {           # two passes: numbers first, then page numbers in lists
        [void]$doc.Fields.Update()
        foreach ($t in $doc.TablesOfContents) { [void]$t.Update() }
        foreach ($t in $doc.TablesOfFigures)  { [void]$t.Update() }
        foreach ($s in $doc.Sections) {
            foreach ($h in $s.Headers) { [void]$h.Range.Fields.Update() }
            foreach ($f in $s.Footers) { [void]$f.Range.Fields.Update() }
        }
        $doc.Repaginate()
    }
    $text = $doc.Content.Text
    $errors = ([regex]::Matches($text, "Error! (Reference source not found|Bookmark not defined|No table of)")).Count
    Write-Output ("Fields: {0}; TOCs: {1}; TablesOfFigures: {2}; Pages: {3}; field errors: {4}" -f `
        $doc.Fields.Count, $doc.TablesOfContents.Count, $doc.TablesOfFigures.Count, `
        $doc.ComputeStatistics(2), $errors)
    $doc.Save()
    Write-Output ("{0:T}  DOCX saved; exporting PDF" -f (Get-Date))
    # ExportAsFixedFormat(OutputFileName, ExportFormat=17 PDF, OpenAfterExport, OptimizeFor=0 print,
    #   Range=0 whole document, From, To, Item=0 content, IncludeDocProps, KeepIRM,
    #   CreateBookmarks=1 headings, DocStructureTags, BitmapMissingFonts, UseISO19005_1)
    # DocStructureTags is $false: tagging 100+ large tables makes Word's export extremely slow.
    $doc.ExportAsFixedFormat($pdf, 17, $false, 0, 0, 1, 1, 0, $true, $true, 1, $false, $true, $false)
    Write-Output ("{0:T}  PDF exported" -f (Get-Date))
    $doc.Close($false)
} finally {
    $word.Quit()
    [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($word)
}
Write-Output "Saved $docx and $pdf"
