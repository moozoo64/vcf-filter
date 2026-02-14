//! VCF row parser using nom.
//!
//! Parses individual VCF data rows into structured `VcfRow` objects,
//! including parsing of INFO fields and structured annotations like ANN.

use std::collections::HashMap;

use crate::error::{Result, VcfFilterError};
use crate::header::{InfoField, InfoMap, InfoNumber, InfoType};
use crate::value::Value;

/// A parsed VCF data row.
#[derive(Debug, Clone)]
pub struct VcfRow {
    /// Chromosome (CHROM column).
    pub chrom: String,
    /// Position (POS column).
    pub pos: u64,
    /// Variant ID (ID column).
    pub id: Option<String>,
    /// Reference allele (REF column).
    pub ref_allele: String,
    /// Alternate allele(s) (ALT column).
    pub alt_alleles: Vec<String>,
    /// Quality score (QUAL column).
    pub qual: Option<f64>,
    /// Filter status (FILTER column).
    pub filter: Vec<String>,
    /// INFO fields parsed into values.
    pub info: HashMap<String, Value>,
    /// FORMAT fields (sample genotype data like GT, DP, GQ).
    pub format: HashMap<String, Value>,
}

/// A single annotation from a structured field like ANN.
pub type Annotation = HashMap<String, String>;

impl VcfRow {
    /// Get a value from the row by field name.
    ///
    /// Supports built-in fields (CHROM, POS, REF, ALT, QUAL, FILTER, ID)
    /// and INFO fields.
    pub fn get(&self, field: &str) -> Value {
        match field {
            "CHROM" => Value::String(self.chrom.clone()),
            "POS" => Value::Number(self.pos as f64),
            "ID" => self
                .id
                .as_ref()
                .map(|s| Value::String(s.clone()))
                .unwrap_or(Value::Missing),
            "REF" => Value::String(self.ref_allele.clone()),
            "ALT" => {
                if self.alt_alleles.len() == 1 {
                    Value::String(self.alt_alleles[0].clone())
                } else {
                    Value::Array(
                        self.alt_alleles
                            .iter()
                            .map(|s| Value::String(s.clone()))
                            .collect(),
                    )
                }
            }
            "QUAL" => self.qual.map(Value::Number).unwrap_or(Value::Missing),
            "FILTER" => {
                if self.filter.len() == 1 {
                    Value::String(self.filter[0].clone())
                } else {
                    Value::Array(
                        self.filter
                            .iter()
                            .map(|s| Value::String(s.clone()))
                            .collect(),
                    )
                }
            }
            _ => {
                // Check INFO fields first, then FORMAT fields.
                // INFO is the primary namespace for filter fields and should
                // not be shadowed by FORMAT fields with the same name (e.g., DP).
                if let Some(val) = self.info.get(field) {
                    val.clone()
                } else {
                    self.format.get(field).cloned().unwrap_or(Value::Missing)
                }
            }
        }
    }
}

/// Parse INFO field values based on their type.
fn parse_info_value(raw: &str, field: &InfoField) -> Value {
    // Handle structured fields with subfields (like ANN)
    if let Some(ref subfield_names) = field.subfields {
        // Split by comma for multiple annotations
        let annotations: Vec<Value> = raw
            .split(',')
            .map(|ann| {
                let parts: Vec<&str> = ann.split('|').collect();
                let mut map = HashMap::new();
                for (i, name) in subfield_names.iter().enumerate() {
                    if let Some(val) = parts.get(i) {
                        map.insert(name.clone(), val.to_string());
                    }
                }
                // Convert to a nested Value structure
                Value::Array(
                    subfield_names
                        .iter()
                        .map(|name| {
                            map.get(name)
                                .map(|v| Value::String(v.clone()))
                                .unwrap_or(Value::Missing)
                        })
                        .collect(),
                )
            })
            .collect();

        // Store as array of arrays, but also keep the raw parsed data accessible
        return Value::Array(annotations);
    }

    // Handle based on type and number
    match (&field.number, &field.field_type) {
        (InfoNumber::Count(1), InfoType::Integer) => raw
            .parse::<i64>()
            .map(|n| Value::Number(n as f64))
            .unwrap_or(Value::String(raw.to_string())),
        (InfoNumber::Count(1), InfoType::Float) => raw
            .parse::<f64>()
            .map(Value::Number)
            .unwrap_or(Value::String(raw.to_string())),
        (InfoNumber::Flag, _) => Value::Bool(true),
        (_, InfoType::Integer) => {
            // Multiple integers
            let values: Vec<Value> = raw
                .split(',')
                .map(|s| {
                    s.parse::<i64>()
                        .map(|n| Value::Number(n as f64))
                        .unwrap_or(Value::String(s.to_string()))
                })
                .collect();
            if values.len() == 1 {
                values.into_iter().next().unwrap()
            } else {
                Value::Array(values)
            }
        }
        (_, InfoType::Float) => {
            // Multiple floats
            let values: Vec<Value> = raw
                .split(',')
                .map(|s| {
                    s.parse::<f64>()
                        .map(Value::Number)
                        .unwrap_or(Value::String(s.to_string()))
                })
                .collect();
            if values.len() == 1 {
                values.into_iter().next().unwrap()
            } else {
                Value::Array(values)
            }
        }
        _ => {
            // String or unknown - check for multiple values
            if raw.contains(',') && !raw.contains('|') {
                Value::Array(
                    raw.split(',')
                        .map(|s| Value::String(s.to_string()))
                        .collect(),
                )
            } else {
                Value::String(raw.to_string())
            }
        }
    }
}

/// Parse INFO field when no metadata is available.
fn parse_info_value_unknown(raw: &str) -> Value {
    // Try to parse as number
    if let Ok(n) = raw.parse::<f64>() {
        return Value::Number(n);
    }

    // Check for multiple values
    if raw.contains(',') && !raw.contains('|') {
        return Value::Array(
            raw.split(',')
                .map(|s| Value::String(s.to_string()))
                .collect(),
        );
    }

    Value::String(raw.to_string())
}

/// Parse the INFO column into a map of field names to values.
fn parse_info_column(info_str: &str, info_map: &InfoMap) -> HashMap<String, Value> {
    let mut result = HashMap::new();

    if info_str == "." {
        return result;
    }

    for field in info_str.split(';') {
        if field.is_empty() {
            continue;
        }

        if let Some((key, value)) = field.split_once('=') {
            let parsed_value = if let Some(field_meta) = info_map.get(key) {
                parse_info_value(value, field_meta)
            } else {
                parse_info_value_unknown(value)
            };
            result.insert(key.to_string(), parsed_value);
        } else {
            // Flag field (no value)
            result.insert(field.to_string(), Value::Bool(true));
        }
    }

    result
}

/// Parse a single VCF data row.
///
/// # Arguments
///
/// * `row` - A single line from the VCF file (tab-separated)
/// * `info_map` - The INFO field metadata from the header
///
/// # Returns
///
/// A parsed `VcfRow` structure.
pub fn parse_row(row: &str, info_map: &InfoMap) -> Result<VcfRow> {
    let fields: Vec<&str> = row.split('\t').collect();

    if fields.len() < 8 {
        return Err(VcfFilterError::RowParseError(format!(
            "Expected at least 8 columns, got {}",
            fields.len()
        )));
    }

    let chrom = fields[0].to_string();

    let pos = fields[1]
        .parse::<u64>()
        .map_err(|e| VcfFilterError::RowParseError(format!("Invalid POS: {}", e)))?;

    let id = if fields[2] == "." {
        None
    } else {
        Some(fields[2].to_string())
    };

    let ref_allele = fields[3].to_string();

    let alt_alleles: Vec<String> = if fields[4] == "." {
        vec![]
    } else {
        fields[4].split(',').map(|s| s.to_string()).collect()
    };

    let qual = if fields[5] == "." {
        None
    } else {
        fields[5].parse::<f64>().ok()
    };

    let filter: Vec<String> = if fields[6] == "." {
        vec![]
    } else {
        fields[6].split(';').map(|s| s.to_string()).collect()
    };

    let info = parse_info_column(fields[7], info_map);

    // Parse FORMAT and sample columns if present (columns 9 and 10+)
    let format = if fields.len() >= 10 {
        parse_format_columns(fields[8], fields[9])
    } else {
        HashMap::new()
    };

    Ok(VcfRow {
        chrom,
        pos,
        id,
        ref_allele,
        alt_alleles,
        qual,
        filter,
        info,
        format,
    })
}

/// Parse FORMAT and sample columns into a HashMap.
///
/// FORMAT column contains colon-separated field names (e.g., "GT:DP:GQ"),
/// and sample column contains corresponding colon-separated values (e.g., "0/1:30:99").
fn parse_format_columns(format_str: &str, sample_str: &str) -> HashMap<String, Value> {
    let mut result = HashMap::new();

    let format_keys: Vec<&str> = format_str.split(':').collect();
    let sample_values: Vec<&str> = sample_str.split(':').collect();

    for (i, key) in format_keys.iter().enumerate() {
        if let Some(value) = sample_values.get(i) {
            let val = if *value == "." {
                Value::Missing
            } else {
                Value::String(value.to_string())
            };
            result.insert(key.to_string(), val);
        }
    }

    result
}

/// Helper to access a subfield from a structured annotation.
///
/// # Arguments
///
/// * `row` - The parsed VCF row
/// * `field` - The INFO field name (e.g., "ANN")
/// * `index` - The annotation index (0-based)
/// * `subfield` - The subfield name (e.g., "Gene_Name")
/// * `info_map` - The header metadata
///
/// # Returns
///
/// The value at the specified subfield, or Missing if not found.
pub fn get_annotation_subfield(
    row: &VcfRow,
    field: &str,
    index: usize,
    subfield: &str,
    info_map: &InfoMap,
) -> Value {
    let field_meta = match info_map.get(field) {
        Some(f) => f,
        None => return Value::Missing,
    };

    let subfield_names = match &field_meta.subfields {
        Some(s) => s,
        None => return Value::Missing,
    };

    let subfield_index = match subfield_names.iter().position(|s| s == subfield) {
        Some(i) => i,
        None => return Value::Missing,
    };

    let annotations = match row.info.get(field) {
        Some(Value::Array(arr)) => arr,
        _ => return Value::Missing,
    };

    let annotation = match annotations.get(index) {
        Some(Value::Array(arr)) => arr,
        _ => return Value::Missing,
    };

    annotation
        .get(subfield_index)
        .cloned()
        .unwrap_or(Value::Missing)
}

/// Get all values for a subfield across all annotations (for wildcard access).
pub fn get_all_annotation_subfields(
    row: &VcfRow,
    field: &str,
    subfield: &str,
    info_map: &InfoMap,
) -> Vec<Value> {
    let field_meta = match info_map.get(field) {
        Some(f) => f,
        None => return vec![],
    };

    let subfield_names = match &field_meta.subfields {
        Some(s) => s,
        None => return vec![],
    };

    let subfield_index = match subfield_names.iter().position(|s| s == subfield) {
        Some(i) => i,
        None => return vec![],
    };

    let annotations = match row.info.get(field) {
        Some(Value::Array(arr)) => arr,
        _ => return vec![],
    };

    annotations
        .iter()
        .filter_map(|ann| {
            if let Value::Array(arr) = ann {
                arr.get(subfield_index).cloned()
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::parse_header;

    const HEADER: &str = r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">"#;

    #[test]
    fn test_parse_simple_row() {
        let info_map = parse_header(HEADER).unwrap();
        let row = "chr1\t12345\trs123\tA\tG\t30.5\tPASS\tEND=12400";
        let parsed = parse_row(row, &info_map).unwrap();

        assert_eq!(parsed.chrom, "chr1");
        assert_eq!(parsed.pos, 12345);
        assert_eq!(parsed.id, Some("rs123".to_string()));
        assert_eq!(parsed.ref_allele, "A");
        assert_eq!(parsed.alt_alleles, vec!["G"]);
        assert_eq!(parsed.qual, Some(30.5));
        assert_eq!(parsed.filter, vec!["PASS"]);
        assert_eq!(parsed.info.get("END"), Some(&Value::Number(12400.0)));
    }

    #[test]
    fn test_parse_row_with_ann() {
        let info_map = parse_header(HEADER).unwrap();
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tANN=G|missense_variant|HIGH|BRCA1|ENSG123|transcript|ENST456|protein_coding|1/10|c.100A>G|p.Thr34Ala|100/500|100/400|34/133||";
        let parsed = parse_row(row, &info_map).unwrap();

        // Check annotation access
        let gene = get_annotation_subfield(&parsed, "ANN", 0, "Gene_Name", &info_map);
        assert_eq!(gene, Value::String("BRCA1".to_string()));

        let impact = get_annotation_subfield(&parsed, "ANN", 0, "Annotation_Impact", &info_map);
        assert_eq!(impact, Value::String("HIGH".to_string()));
    }

    #[test]
    fn test_parse_row_with_multiple_annotations() {
        let info_map = parse_header(HEADER).unwrap();
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tANN=G|missense|HIGH|BRCA1|E1|t|T1|pc|1|c.1|p.1|1|1|1||,G|synonymous|LOW|BRCA2|E2|t|T2|pc|2|c.2|p.2|2|2|2||";
        let parsed = parse_row(row, &info_map).unwrap();

        let gene1 = get_annotation_subfield(&parsed, "ANN", 0, "Gene_Name", &info_map);
        assert_eq!(gene1, Value::String("BRCA1".to_string()));

        let gene2 = get_annotation_subfield(&parsed, "ANN", 1, "Gene_Name", &info_map);
        assert_eq!(gene2, Value::String("BRCA2".to_string()));

        // Test wildcard access
        let all_genes = get_all_annotation_subfields(&parsed, "ANN", "Gene_Name", &info_map);
        assert_eq!(all_genes.len(), 2);
    }

    #[test]
    fn test_get_builtin_fields() {
        let info_map = parse_header(HEADER).unwrap();
        let row = "chr1\t12345\trs123\tA\tG\t30.5\tPASS\t.";
        let parsed = parse_row(row, &info_map).unwrap();

        assert_eq!(parsed.get("CHROM"), Value::String("chr1".to_string()));
        assert_eq!(parsed.get("POS"), Value::Number(12345.0));
        assert_eq!(parsed.get("QUAL"), Value::Number(30.5));
        assert_eq!(parsed.get("FILTER"), Value::String("PASS".to_string()));
    }

    #[test]
    fn test_info_field_takes_precedence_over_format_field() {
        let header = r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">"#;
        let info_map = parse_header(header).unwrap();
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP\t0/1:8";
        let parsed = parse_row(row, &info_map).unwrap();

        assert_eq!(parsed.get("DP"), Value::Number(30.0));
    }

    #[test]
    fn test_format_field_used_when_info_missing() {
        let info_map = parse_header(HEADER).unwrap();
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15";
        let parsed = parse_row(row, &info_map).unwrap();

        assert_eq!(parsed.get("DP"), Value::String("15".to_string()));
    }
}
