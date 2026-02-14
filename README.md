# VCF Filter

A high-performance Rust library for filtering VCF (Variant Call Format) files using an expressive filter DSL. Designed for bioinformatics workflows that need to filter variants based on quality scores, annotations, and clinical significance.

## Features

- **Expressive filter syntax** — Boolean logic, comparisons, string matching
- **Structured annotation support** — Access SnpEff ANN, LOF, NMD subfields by name
- **Wildcard matching** — `ANN[*].Gene_Name == "BRCA1"` matches any annotation
- **Auto-detection** — Parses INFO field metadata from VCF headers
- **Zero-copy parsing** — Efficient processing of large VCF files

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
vcf-filter = "0.1.1"
```

## Quick Start

```rust
use vcf_filter::FilterEngine;

// Parse the VCF header to extract INFO field definitions
let header = concat!(
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n",
    "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID'\">"
);

let engine = FilterEngine::new(header).unwrap();

// Filter a VCF data row
let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30;ANN=G|missense|HIGH|BRCA1|ENSG123";

// Evaluate filters
assert!(engine.evaluate("QUAL > 30", row).unwrap());
assert!(engine.evaluate(r#"FILTER == "PASS""#, row).unwrap());
assert!(engine.evaluate("DP >= 30", row).unwrap());
assert!(engine.evaluate(r#"ANN[0].Gene_Name == "BRCA1""#, row).unwrap());
```

## Filter Syntax

### Comparison Operators

| Operator | Example | Description |
|----------|---------|-------------|
| `==` | `FILTER == "PASS"` | Equal |
| `!=` | `CLNSIG != "Benign"` | Not equal |
| `>` | `QUAL > 30` | Greater than |
| `<` | `DP < 100` | Less than |
| `>=` | `QUAL >= 30` | Greater than or equal |
| `<=` | `DP <= 50` | Less than or equal |
| `contains` | `CLNDN contains "cancer"` | Substring match |

### Logical Operators

| Operator | Example | Description |
|----------|---------|-------------|
| `&&` | `QUAL > 30 && DP >= 10` | Logical AND |
| `\|\|` | `CLNSIG == "Benign" \|\| CLNSIG == "Likely_benign"` | Logical OR |
| `!` | `!exists(LOF)` | Logical NOT |
| `()` | `(A \|\| B) && C` | Grouping |

### Field Access

```rust
// Built-in VCF columns
"CHROM == \"chr1\""      // Chromosome
"POS > 100000"           // Position
"ID == \"rs123\""        // Variant ID
"REF == \"A\""           // Reference allele
"ALT == \"G\""           // Alternate allele
"QUAL >= 30"             // Quality score
"FILTER == \"PASS\""     // Filter status

// FORMAT/sample fields (genotype data)
"GT == \"0/0\""          // Homozygous reference
"GT == \"0/1\""          // Heterozygous
"GT == \"1/1\""          // Homozygous alternate

// INFO fields
"DP >= 30"               // Read depth
"AF > 0.01"              // Allele frequency

// Structured annotations (index access)
"ANN[0].Gene_Name"           // First annotation's gene
"ANN[0].Annotation_Impact"   // First annotation's impact

// Wildcard access (any match)
"ANN[*].Annotation_Impact == \"HIGH\""  // Any annotation has HIGH impact
"ANN[*].Gene_Name == \"BRCA1\""         // Any annotation for BRCA1
```

### Functions

| Function | Example | Description |
|----------|---------|-------------|
| `exists()` | `exists(CLNSIG)` | True if field is present and not missing |

## Examples

### Filter by Quality and Depth

```rust
let filter = "QUAL >= 30 && DP >= 10";
if engine.evaluate(filter, row)? {
    println!("Variant passes quality filters");
}
```

### Filter by Clinical Significance

```rust
// Keep pathogenic or likely pathogenic variants
let filter = r#"CLNSIG == "Pathogenic" || CLNSIG == "Likely_pathogenic""#;
```

### Filter by Annotation Impact

```rust
// Keep variants with HIGH or MODERATE impact
let filter = r#"ANN[*].Annotation_Impact == "HIGH" || ANN[*].Annotation_Impact == "MODERATE""#;
```

### Complex Filter with Grouping

```rust
// High quality, passing variants in BRCA genes
let filter = r#"QUAL > 30 && FILTER == "PASS" && (ANN[*].Gene_Name == "BRCA1" || ANN[*].Gene_Name == "BRCA2")"#;
```

### Filter by Disease Association

```rust
// Variants associated with cancer
let filter = r#"CLNDN contains "cancer" || CLNDN contains "carcinoma""#;
```

### Efficient Batch Processing

```rust
// Pre-parse filter and rows for better performance
let expr = engine.parse_filter("QUAL > 30 && DP >= 10")?;

for line in vcf_rows {
    let row = engine.parse_row(line)?;
    if engine.evaluate_parsed(&expr, &row)? {
        // Process matching variant
    }
}
```

## Supported Annotation Formats

The library automatically detects structured annotations from VCF header descriptions:

### SnpEff ANN Format
```
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
```

Access subfields by name (normalized: spaces→`_`, dots→`_`, slashes→`_`):
- `ANN[0].Allele`
- `ANN[0].Annotation_Impact`
- `ANN[0].Gene_Name`
- `ANN[0].HGVS_c`

### SnpEff LOF/NMD Format
```
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
```

## API Reference

### FilterEngine

```rust
impl FilterEngine {
    /// Create a new engine from VCF header
    pub fn new(header: &str) -> Result<Self>;
    
    /// Evaluate a filter against a row
    pub fn evaluate(&self, filter: &str, row: &str) -> Result<bool>;
    
    /// Parse a row for reuse
    pub fn parse_row(&self, row: &str) -> Result<VcfRow>;
    
    /// Parse a filter for reuse
    pub fn parse_filter(&self, filter: &str) -> Result<Expr>;
    
    /// Evaluate pre-parsed filter against pre-parsed row
    pub fn evaluate_parsed(&self, expr: &Expr, row: &VcfRow) -> Result<bool>;
    
    /// Get INFO field metadata
    pub fn info_map(&self) -> &InfoMap;
}
```

### Error Handling

```rust
use vcf_filter::{FilterEngine, VcfFilterError};

match engine.evaluate(filter, row) {
    Ok(true) => println!("Match"),
    Ok(false) => println!("No match"),
    Err(VcfFilterError::FilterParseError(msg)) => eprintln!("Invalid filter: {}", msg),
    Err(VcfFilterError::RowParseError(msg)) => eprintln!("Invalid row: {}", msg),
    Err(e) => eprintln!("Error: {}", e),
}
```

## License

MIT

## Author

This library was vibe coded by Michael Simmons using GitHub Copilot with Claude Opus 4.5.

## Contributing

Contributions welcome! Please see the [Copilot Instructions](.github/copilot-instructions.md) for architecture details and development guidelines.
