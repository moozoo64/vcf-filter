//! # VCF Filter Library
//!
//! A high-performance library for filtering VCF (Variant Call Format) files
//! using flexible filter expressions.
//!
//! ## Features
//!
//! - Parse VCF headers to extract INFO field metadata
//! - Parse structured annotations (ANN, LOF, NMD) with automatic subfield detection
//! - Filter expressions with comparison, logical, and containment operators
//! - Array access with indexing (`[0]`) and wildcard (`[*]`) support
//!
//! ## Example
//!
//! ```rust
//! use vcf_filter::FilterEngine;
//!
//! // VCF header with INFO field definitions (SnpEff ANN format)
//! let header = concat!(
//!     "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n",
//!     "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID'\">"
//! );
//!
//! // Initialize the filter engine with the header
//! let engine = FilterEngine::new(header).unwrap();
//!
//! // A VCF data row with ANN annotation
//! let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30;ANN=G|missense|HIGH|BRCA1|ENSG123";
//!
//! // Evaluate filters
//! assert!(engine.evaluate("QUAL > 30", row).unwrap());
//! assert!(engine.evaluate(r#"FILTER == "PASS""#, row).unwrap());
//! assert!(engine.evaluate("DP >= 30", row).unwrap());
//! assert!(engine.evaluate(r#"ANN[0].Gene_Name == "BRCA1""#, row).unwrap());
//! ```
//!
//! ## Filter Expression Syntax
//!
//! ### Comparison Operators
//! - `==` Equal
//! - `!=` Not equal
//! - `>` Greater than
//! - `<` Less than
//! - `>=` Greater than or equal
//! - `<=` Less than or equal
//! - `contains` String containment
//!
//! ### Logical Operators
//! - `&&` Logical AND
//! - `||` Logical OR
//! - `!` Logical NOT
//!
//! ### Field Access
//! - `QUAL` - Built-in VCF column
//! - `FILTER` - Filter status
//! - `DP` - INFO field
//! - `ANN[0].Gene_Name` - First annotation's gene name
//! - `ANN[*].Annotation_Impact` - Any annotation's impact (wildcard)
//!
//! ### Functions
//! - `exists(field)` - Check if a field exists

/// Embedded README.md documentation
const README: &str = include_str!("../README.md");

/// Returns the embedded README.md documentation.
///
/// # Example
///
/// ```rust
/// use vcf_filter::docs;
///
/// let documentation = docs();
/// println!("{}", documentation);
/// ```
pub fn docs() -> &'static str {
    README
}

pub mod error;
pub mod eval;
pub mod filter;
pub mod header;
pub mod row;
pub mod value;

pub use error::{Result, VcfFilterError};
pub use filter::{AccessPart, BinaryOp, Expr, UnaryOp};
pub use header::{InfoField, InfoMap, InfoNumber, InfoType};
pub use row::VcfRow;
pub use value::Value;

use crate::eval::evaluate;
use crate::filter::parse_filter;
use crate::header::parse_header;
use crate::row::parse_row;

/// The main filter engine for evaluating VCF filters.
///
/// Create an instance with `FilterEngine::new(header)` and then use
/// `evaluate(filter, row)` to test rows against filter expressions.
#[derive(Debug, Clone)]
pub struct FilterEngine {
    /// Parsed INFO field metadata from the header.
    info_map: InfoMap,
}

impl FilterEngine {
    /// Create a new FilterEngine from a VCF header.
    ///
    /// # Arguments
    ///
    /// * `header` - The VCF header string containing ##INFO lines
    ///
    /// # Returns
    ///
    /// A new `FilterEngine` instance ready to evaluate filters.
    ///
    /// # Example
    ///
    /// ```rust
    /// use vcf_filter::FilterEngine;
    ///
    /// let header = r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">"#;
    /// let engine = FilterEngine::new(header).unwrap();
    /// ```
    pub fn new(header: &str) -> Result<Self> {
        let info_map = parse_header(header)?;
        Ok(Self { info_map })
    }

    /// Evaluate a filter expression against a VCF row.
    ///
    /// # Arguments
    ///
    /// * `filter` - The filter expression string (e.g., `"QUAL > 30"`)
    /// * `row` - A single VCF data row (tab-separated)
    ///
    /// # Returns
    ///
    /// `true` if the row matches the filter, `false` otherwise.
    ///
    /// # Example
    ///
    /// ```rust
    /// use vcf_filter::FilterEngine;
    ///
    /// let header = r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">"#;
    /// let engine = FilterEngine::new(header).unwrap();
    ///
    /// let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30";
    /// assert!(engine.evaluate("QUAL > 30", row).unwrap());
    /// ```
    pub fn evaluate(&self, filter: &str, row: &str) -> Result<bool> {
        let parsed_row = parse_row(row, &self.info_map)?;
        let expr = parse_filter(filter).map_err(|errs| {
            VcfFilterError::FilterParseError(
                errs.into_iter()
                    .map(|e| e.to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
            )
        })?;
        let result = evaluate(&expr, &parsed_row, &self.info_map)?;
        Ok(result.as_bool().unwrap_or(false))
    }

    /// Parse a VCF row without evaluating a filter.
    ///
    /// Useful for inspecting row data or performing multiple evaluations
    /// against the same row.
    ///
    /// # Arguments
    ///
    /// * `row` - A single VCF data row (tab-separated)
    ///
    /// # Returns
    ///
    /// A parsed `VcfRow` structure.
    pub fn parse_row(&self, row: &str) -> Result<VcfRow> {
        parse_row(row, &self.info_map)
    }

    /// Parse a filter expression without evaluating it.
    ///
    /// Useful for validating filter syntax or caching parsed expressions.
    ///
    /// # Arguments
    ///
    /// * `filter` - The filter expression string
    ///
    /// # Returns
    ///
    /// A parsed `Expr` AST.
    pub fn parse_filter(&self, filter: &str) -> Result<Expr> {
        parse_filter(filter).map_err(|errs| {
            VcfFilterError::FilterParseError(
                errs.into_iter()
                    .map(|e| e.to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
            )
        })
    }

    /// Evaluate a pre-parsed filter expression against a pre-parsed row.
    ///
    /// This is more efficient when evaluating the same filter against
    /// many rows, or the same row against many filters.
    ///
    /// # Arguments
    ///
    /// * `expr` - The parsed filter expression
    /// * `row` - The parsed VCF row
    ///
    /// # Returns
    ///
    /// `true` if the row matches the filter, `false` otherwise.
    pub fn evaluate_parsed(&self, expr: &Expr, row: &VcfRow) -> Result<bool> {
        let result = evaluate(expr, row, &self.info_map)?;
        Ok(result.as_bool().unwrap_or(false))
    }

    /// Get the INFO field metadata map.
    ///
    /// Useful for inspecting what fields are available and their types.
    pub fn info_map(&self) -> &InfoMap {
        &self.info_map
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FULL_HEADER: &str = r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Aggregate germline classification for this single variant; multiple values are separated by a vertical bar">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">"#;

    const REAL_ROW: &str = "chr1\t186308857\trs3737940\tA\tG\t37.5\tPASS\tANN=G|synonymous_variant|LOW|PRG4|ENSG00000116690|transcript|ENST00000445192.6|protein_coding|7/13|c.3138A>G|p.Pro1046Pro|3183/5044|3138/4215|1046/1404||,G|synonymous_variant|LOW|PRG4|ENSG00000116690|transcript|ENST00000367483.8|protein_coding|6/12|c.3015A>G|p.Pro1005Pro|3060/4921|3015/4092|1005/1363||,G|synonymous_variant|LOW|PRG4|ENSG00000116690|transcript|ENST00000367485.4|protein_coding|5/11|c.2859A>G|p.Pro953Pro|2904/4765|2859/3936|953/1311||,G|synonymous_variant|LOW|PRG4|ENSG00000116690|transcript|ENST00000635041.1|protein_coding|6/12|c.3009A>G|p.Pro1003Pro|3054/4915|3009/4086|1003/1361||,G|downstream_gene_variant|MODIFIER|TPR|ENSG00000047410|transcript|ENST00000367478.9|protein_coding||c.*5114T>C|||||2795|,G|downstream_gene_variant|MODIFIER|PRG4|ENSG00000116690|transcript|ENST00000533951.5|protein_coding||c.*2379A>G|||||2379|WARNING_TRANSCRIPT_NO_STOP_CODON,G|downstream_gene_variant|MODIFIER|PRG4|ENSG00000116690|transcript|ENST00000367482.8|protein_coding||c.*175A>G|||||175|WARNING_TRANSCRIPT_INCOMPLETE,G|downstream_gene_variant|MODIFIER|RNU6-1240P|ENSG00000202025|transcript|ENST00000365155.1|snRNA||n.*2968T>C|||||2968|;CLNSIG=Benign;CLNDN=not_provided|Camptodactyly-arthropathy-coxa_vara-pericarditis_syndrome|not_specified";

    #[test]
    fn test_engine_creation() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        assert!(engine.info_map().contains_key("ANN"));
        assert!(engine.info_map().contains_key("CLNSIG"));
    }

    #[test]
    fn test_real_row_qual() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        assert!(engine.evaluate("QUAL > 30", REAL_ROW).unwrap());
        assert!(engine.evaluate("QUAL < 40", REAL_ROW).unwrap());
        assert!(engine.evaluate("QUAL == 37.5", REAL_ROW).unwrap());
    }

    #[test]
    fn test_real_row_filter() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        assert!(engine.evaluate(r#"FILTER == "PASS""#, REAL_ROW).unwrap());
    }

    #[test]
    fn test_real_row_clnsig() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        assert!(engine.evaluate(r#"CLNSIG == "Benign""#, REAL_ROW).unwrap());
    }

    #[test]
    fn test_real_row_ann_first() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        assert!(
            engine
                .evaluate(r#"ANN[0].Annotation == "synonymous_variant""#, REAL_ROW)
                .unwrap()
        );
        assert!(
            engine
                .evaluate(r#"ANN[0].Annotation_Impact == "LOW""#, REAL_ROW)
                .unwrap()
        );
        assert!(
            engine
                .evaluate(r#"ANN[0].Gene_Name == "PRG4""#, REAL_ROW)
                .unwrap()
        );
    }

    #[test]
    fn test_real_row_ann_wildcard() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        // Any annotation has MODIFIER impact
        assert!(
            engine
                .evaluate(r#"ANN[*].Annotation_Impact == "MODIFIER""#, REAL_ROW)
                .unwrap()
        );
        // Any annotation is for TPR gene
        assert!(
            engine
                .evaluate(r#"ANN[*].Gene_Name == "TPR""#, REAL_ROW)
                .unwrap()
        );
        // No annotation has HIGH impact
        assert!(
            !engine
                .evaluate(r#"ANN[*].Annotation_Impact == "HIGH""#, REAL_ROW)
                .unwrap()
        );
    }

    #[test]
    fn test_real_row_chrom_pos() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        assert!(engine.evaluate(r#"CHROM == "chr1""#, REAL_ROW).unwrap());
        assert!(engine.evaluate("POS > 186000000", REAL_ROW).unwrap());
    }

    #[test]
    fn test_real_row_complex_filter() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        // Combined filter: high quality, passed, and benign
        assert!(
            engine
                .evaluate(
                    r#"QUAL > 30 && FILTER == "PASS" && CLNSIG == "Benign""#,
                    REAL_ROW
                )
                .unwrap()
        );
    }

    #[test]
    fn test_real_row_contains() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        assert!(
            engine
                .evaluate(r#"CLNDN contains "Camptodactyly""#, REAL_ROW)
                .unwrap()
        );
    }

    #[test]
    fn test_parse_and_reuse() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();

        // Parse filter once
        let expr = engine.parse_filter("QUAL > 30").unwrap();

        // Parse row once
        let row = engine.parse_row(REAL_ROW).unwrap();

        // Evaluate multiple times efficiently
        assert!(engine.evaluate_parsed(&expr, &row).unwrap());
    }

    #[test]
    fn test_missing_field() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        // LOF doesn't exist in this row
        assert!(!engine.evaluate("exists(LOF)", REAL_ROW).unwrap());
        assert!(engine.evaluate("!exists(LOF)", REAL_ROW).unwrap());
    }

    #[test]
    fn test_existing_field() {
        let engine = FilterEngine::new(FULL_HEADER).unwrap();
        assert!(engine.evaluate("exists(CLNSIG)", REAL_ROW).unwrap());
        assert!(engine.evaluate("exists(ANN)", REAL_ROW).unwrap());
    }

    #[test]
    fn test_doctest_scenario() {
        // This matches the doctest example exactly
        let header = r#"
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID'">
"#;
        let engine = FilterEngine::new(header).unwrap();

        // Check that ANN subfields were parsed
        let ann_field = engine.info_map().get("ANN").unwrap();
        let subfields = ann_field.subfields.as_ref().unwrap();
        assert_eq!(subfields[0], "Allele");
        assert_eq!(subfields[1], "Annotation");
        assert_eq!(subfields[2], "Annotation_Impact");
        assert_eq!(subfields[3], "Gene_Name");
        assert_eq!(subfields[4], "Gene_ID");

        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30;ANN=G|missense|HIGH|BRCA1|ENSG123";

        assert!(engine.evaluate("QUAL > 30", row).unwrap());
        assert!(engine.evaluate(r#"FILTER == "PASS""#, row).unwrap());
        assert!(engine.evaluate("DP >= 30", row).unwrap());
        assert!(
            engine
                .evaluate(r#"ANN[0].Gene_Name == "BRCA1""#, row)
                .unwrap()
        );
    }
}
