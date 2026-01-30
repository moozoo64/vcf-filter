//! VCF header parser.
//!
//! Parses ##INFO lines to extract field metadata, including subfield names
//! for structured annotations like ANN, LOF, and NMD.

use std::collections::HashMap;

use crate::error::Result;

/// The number of values an INFO field can have.
#[derive(Debug, Clone, PartialEq)]
pub enum InfoNumber {
    /// A fixed number of values.
    Count(usize),
    /// One value per alternate allele (A).
    PerAltAllele,
    /// One value per possible genotype (G).
    PerGenotype,
    /// One value per allele including reference (R).
    PerAllele,
    /// Variable number of values (.).
    Variable,
    /// Flag type (0 values, presence indicates true).
    Flag,
}

/// The data type of an INFO field.
#[derive(Debug, Clone, PartialEq)]
pub enum InfoType {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

/// Metadata for a single INFO field parsed from the header.
#[derive(Debug, Clone)]
pub struct InfoField {
    /// The field identifier (e.g., "ANN", "DP").
    pub id: String,
    /// The number of values this field can have.
    pub number: InfoNumber,
    /// The data type of the field.
    pub field_type: InfoType,
    /// The description from the header.
    pub description: String,
    /// Subfield names for structured fields (e.g., ANN).
    /// Extracted from the description if it contains a format specification.
    pub subfields: Option<Vec<String>>,
}

/// Map of INFO field ID to its metadata.
pub type InfoMap = HashMap<String, InfoField>;

/// Parse the Number attribute from an INFO line.
fn parse_number(s: &str) -> InfoNumber {
    match s {
        "A" => InfoNumber::PerAltAllele,
        "G" => InfoNumber::PerGenotype,
        "R" => InfoNumber::PerAllele,
        "." => InfoNumber::Variable,
        "0" => InfoNumber::Flag,
        _ => InfoNumber::Count(s.parse().unwrap_or(1)),
    }
}

/// Parse the Type attribute from an INFO line.
fn parse_type(s: &str) -> InfoType {
    match s {
        "Integer" => InfoType::Integer,
        "Float" => InfoType::Float,
        "Flag" => InfoType::Flag,
        "Character" => InfoType::Character,
        _ => InfoType::String,
    }
}

/// Parse the key=value pairs from inside the INFO angle brackets.
fn parse_info_attrs(content: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    let mut remaining = content;

    while !remaining.is_empty() {
        // Find the key
        let eq_pos = match remaining.find('=') {
            Some(p) => p,
            None => break,
        };
        let key = remaining[..eq_pos].trim();
        remaining = &remaining[eq_pos + 1..];

        // Parse the value
        let value = if remaining.starts_with('"') {
            // Quoted value - find the closing quote
            remaining = &remaining[1..];
            let end_quote = remaining.find('"').unwrap_or(remaining.len());
            let val = &remaining[..end_quote];
            remaining = &remaining[(end_quote + 1).min(remaining.len())..];
            // Skip comma if present
            if remaining.starts_with(',') {
                remaining = &remaining[1..];
            }
            val.to_string()
        } else {
            // Unquoted value - find comma or end
            let comma_pos = remaining.find(',').unwrap_or(remaining.len());
            let val = &remaining[..comma_pos];
            remaining = if comma_pos < remaining.len() {
                &remaining[comma_pos + 1..]
            } else {
                ""
            };
            val.to_string()
        };

        attrs.insert(key.to_string(), value);
    }

    attrs
}

/// Extract subfield names from a description string.
///
/// Looks for patterns like:
/// - "Format: 'Gene_Name | Gene_ID | ...'"
/// - "'Allele | Annotation | Annotation_Impact | ...'"
fn extract_subfields(description: &str) -> Option<Vec<String>> {
    // Look for content between single quotes that contains pipe separators
    let start = description.find('\'')?;
    let end = description[start + 1..].find('\'')? + start + 1;
    let format_str = &description[start + 1..end];

    // Check if it looks like a pipe-separated format
    if !format_str.contains('|') {
        return None;
    }

    let subfields: Vec<String> = format_str
        .split('|')
        .map(|s| {
            s.trim()
                .replace(' ', "_")
                .replace('.', "_")
                .replace('/', "_")
        })
        .filter(|s| !s.is_empty())
        .collect();

    if subfields.len() >= 2 {
        Some(subfields)
    } else {
        None
    }
}

/// Parse a single ##INFO line.
fn parse_info_line(line: &str) -> Option<InfoField> {
    let line = line.strip_prefix("##INFO=<")?;
    let line = line.strip_suffix('>')?;

    let attrs = parse_info_attrs(line);

    let id = attrs.get("ID")?.clone();
    let number_str = attrs.get("Number")?;
    let type_str = attrs.get("Type")?;
    let description = attrs.get("Description").cloned().unwrap_or_default();

    let number = parse_number(number_str);
    let field_type = parse_type(type_str);

    let subfields = extract_subfields(&description);

    Some(InfoField {
        id,
        number,
        field_type,
        description,
        subfields,
    })
}

/// Parse all ##INFO lines from a VCF header string.
///
/// # Arguments
///
/// * `header` - The full VCF header as a string (all lines starting with ##)
///
/// # Returns
///
/// A map of INFO field IDs to their metadata.
pub fn parse_header(header: &str) -> Result<InfoMap> {
    let mut info_map = HashMap::new();

    for line in header.lines() {
        let line = line.trim();
        if line.starts_with("##INFO=<") {
            if let Some(field) = parse_info_line(line) {
                info_map.insert(field.id.clone(), field);
            }
        }
    }

    Ok(info_map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_info() {
        let header = r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position">"#;
        let map = parse_header(header).unwrap();

        let field = map.get("END").unwrap();
        assert_eq!(field.id, "END");
        assert_eq!(field.number, InfoNumber::Count(1));
        assert_eq!(field.field_type, InfoType::Integer);
        assert!(field.subfields.is_none());
    }

    #[test]
    fn test_parse_ann_with_subfields() {
        let header = r#"##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">"#;
        let map = parse_header(header).unwrap();

        let field = map.get("ANN").unwrap();
        assert_eq!(field.id, "ANN");
        assert_eq!(field.number, InfoNumber::Variable);
        assert_eq!(field.field_type, InfoType::String);

        let subfields = field.subfields.as_ref().unwrap();
        assert_eq!(subfields[0], "Allele");
        assert_eq!(subfields[1], "Annotation");
        assert_eq!(subfields[2], "Annotation_Impact");
        assert_eq!(subfields[3], "Gene_Name");
    }

    #[test]
    fn test_parse_lof_with_subfields() {
        let header = r#"##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">"#;
        let map = parse_header(header).unwrap();

        let field = map.get("LOF").unwrap();
        let subfields = field.subfields.as_ref().unwrap();
        assert_eq!(subfields.len(), 4);
        assert_eq!(subfields[0], "Gene_Name");
        assert_eq!(subfields[1], "Gene_ID");
    }

    #[test]
    fn test_parse_multiple_info_lines() {
        let header = r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">"#;

        let map = parse_header(header).unwrap();
        assert_eq!(map.len(), 3);
        assert!(map.contains_key("DP"));
        assert!(map.contains_key("AF"));
        assert!(map.contains_key("CLNSIG"));

        assert_eq!(map.get("AF").unwrap().number, InfoNumber::PerAltAllele);
    }
}
