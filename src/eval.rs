//! Expression evaluator for VCF filter expressions.
//!
//! Evaluates parsed filter ASTs against VCF row data.

use crate::error::{Result, VcfFilterError};
use crate::filter::{AccessPart, BinaryOp, Expr, UnaryOp};
use crate::header::InfoMap;
use crate::row::{VcfRow, get_all_annotation_subfields, get_annotation_subfield};
use crate::value::Value;

/// Evaluate a filter expression against a VCF row.
///
/// # Arguments
///
/// * `expr` - The parsed filter expression AST
/// * `row` - The parsed VCF row
/// * `info_map` - The header metadata for resolving field types
///
/// # Returns
///
/// The result of evaluating the expression as a `Value`.
pub fn evaluate(expr: &Expr, row: &VcfRow, info_map: &InfoMap) -> Result<Value> {
    match expr {
        Expr::Number(n) => Ok(Value::Number(*n)),
        Expr::String(s) => Ok(Value::String(s.clone())),
        Expr::Bool(b) => Ok(Value::Bool(*b)),
        Expr::Var(parts) => resolve_variable(parts, row, info_map),
        Expr::Binary(left, op, right) => evaluate_binary(left, op, right, row, info_map),
        Expr::Unary(op, inner) => evaluate_unary(op, inner, row, info_map),
        Expr::Exists(parts) => {
            let value = resolve_variable(parts, row, info_map)?;
            Ok(Value::Bool(!value.is_missing()))
        }
    }
}

/// Resolve a variable access path to a value.
fn resolve_variable(parts: &[AccessPart], row: &VcfRow, info_map: &InfoMap) -> Result<Value> {
    if parts.is_empty() {
        return Ok(Value::Missing);
    }

    // Get the base field name
    let field_name = match &parts[0] {
        AccessPart::Field(name) => name,
        _ => {
            return Err(VcfFilterError::EvaluationError(
                "Variable must start with a field name".to_string(),
            ));
        }
    };

    // Namespace-qualified access: INFO.<field> or FORMAT.<field>
    if field_name == "INFO" || field_name == "FORMAT" {
        if parts.len() < 2 {
            return Ok(Value::Missing);
        }

        let namespaced_field = match &parts[1] {
            AccessPart::Field(name) => name,
            _ => return Ok(Value::Missing),
        };

        return Ok(resolve_with_base(
            Some(field_name),
            namespaced_field,
            &parts[2..],
            row,
            info_map,
        ));
    }

    Ok(resolve_with_base(
        None,
        field_name,
        &parts[1..],
        row,
        info_map,
    ))
}

/// Resolve field access against a specific namespace.
///
/// `namespace` values:
/// - `None`: unqualified lookup (`row.get` semantics)
/// - `Some("INFO")`: strict INFO lookup
/// - `Some("FORMAT")`: strict FORMAT lookup
fn resolve_with_base(
    namespace: Option<&str>,
    field_name: &str,
    access_parts: &[AccessPart],
    row: &VcfRow,
    info_map: &InfoMap,
) -> Value {
    let base_value = match namespace {
        Some("INFO") => row.info.get(field_name).cloned().unwrap_or(Value::Missing),
        Some("FORMAT") => row
            .format
            .get(field_name)
            .cloned()
            .unwrap_or(Value::Missing),
        _ => row.get(field_name),
    };

    if access_parts.is_empty() {
        return base_value;
    }

    // Handle structured field access (e.g., ANN[0].Gene_Name)
    let mut current_index: Option<usize> = None;
    let mut is_wildcard = false;
    let mut subfield_name: Option<String> = None;

    for part in access_parts {
        match part {
            AccessPart::Index(i) => {
                current_index = Some(*i);
            }
            AccessPart::Wildcard => {
                is_wildcard = true;
            }
            AccessPart::Field(name) => {
                subfield_name = Some(name.clone());
            }
        }
    }

    // If we have a subfield access
    if let Some(ref subfield) = subfield_name {
        // FORMAT namespace does not support annotation-style subfield access.
        if namespace == Some("FORMAT") {
            return Value::Missing;
        }

        if is_wildcard {
            // Return array of all matching subfield values
            let values = get_all_annotation_subfields(row, field_name, subfield, info_map);
            return Value::Array(values);
        } else if let Some(idx) = current_index {
            // Return specific index's subfield
            return get_annotation_subfield(row, field_name, idx, subfield, info_map);
        }
    }

    // Array access without subfield
    if let Some(idx) = current_index {
        if let Value::Array(arr) = base_value {
            return arr.get(idx).cloned().unwrap_or(Value::Missing);
        }
    }

    Value::Missing
}

/// Evaluate a binary operation.
fn evaluate_binary(
    left: &Expr,
    op: &BinaryOp,
    right: &Expr,
    row: &VcfRow,
    info_map: &InfoMap,
) -> Result<Value> {
    let left_val = evaluate(left, row, info_map)?;
    let right_val = evaluate(right, row, info_map)?;

    // Handle wildcard comparisons (array on left side)
    if let Value::Array(ref arr) = left_val {
        let result = match op {
            BinaryOp::Eq => arr.iter().any(|v| values_equal(v, &right_val)),
            BinaryOp::NotEq => arr.iter().all(|v| !values_equal(v, &right_val)),
            BinaryOp::Contains => arr.iter().any(|v| value_contains(v, &right_val)),
            _ => {
                // For numeric comparisons, check if any match
                arr.iter()
                    .any(|v| compare_values(v, op, &right_val).unwrap_or(false))
            }
        };
        return Ok(Value::Bool(result));
    }

    match op {
        BinaryOp::Eq => Ok(Value::Bool(values_equal(&left_val, &right_val))),
        BinaryOp::NotEq => Ok(Value::Bool(!values_equal(&left_val, &right_val))),
        BinaryOp::Lt => Ok(Value::Bool(compare_values(&left_val, op, &right_val)?)),
        BinaryOp::Gt => Ok(Value::Bool(compare_values(&left_val, op, &right_val)?)),
        BinaryOp::LtEq => Ok(Value::Bool(compare_values(&left_val, op, &right_val)?)),
        BinaryOp::GtEq => Ok(Value::Bool(compare_values(&left_val, op, &right_val)?)),
        BinaryOp::Contains => Ok(Value::Bool(value_contains(&left_val, &right_val))),
        BinaryOp::And => {
            let left_bool = value_to_bool(&left_val)?;
            if !left_bool {
                return Ok(Value::Bool(false));
            }
            let right_bool = value_to_bool(&right_val)?;
            Ok(Value::Bool(right_bool))
        }
        BinaryOp::Or => {
            let left_bool = value_to_bool(&left_val)?;
            if left_bool {
                return Ok(Value::Bool(true));
            }
            let right_bool = value_to_bool(&right_val)?;
            Ok(Value::Bool(right_bool))
        }
    }
}

/// Evaluate a unary operation.
fn evaluate_unary(op: &UnaryOp, inner: &Expr, row: &VcfRow, info_map: &InfoMap) -> Result<Value> {
    let val = evaluate(inner, row, info_map)?;

    match op {
        UnaryOp::Not => {
            let bool_val = value_to_bool(&val)?;
            Ok(Value::Bool(!bool_val))
        }
    }
}

/// Check if two values are equal.
fn values_equal(left: &Value, right: &Value) -> bool {
    match (left, right) {
        (Value::String(l), Value::String(r)) => l == r,
        (Value::Number(l), Value::Number(r)) => (l - r).abs() < f64::EPSILON,
        (Value::Bool(l), Value::Bool(r)) => l == r,
        (Value::Missing, Value::Missing) => true,
        // Try numeric comparison if one is a string that looks like a number
        (Value::String(s), Value::Number(n)) | (Value::Number(n), Value::String(s)) => s
            .parse::<f64>()
            .map(|sn| (sn - n).abs() < f64::EPSILON)
            .unwrap_or(false),
        _ => false,
    }
}

/// Check if left contains right (string containment).
fn value_contains(left: &Value, right: &Value) -> bool {
    match (left, right) {
        (Value::String(l), Value::String(r)) => l.contains(r.as_str()),
        _ => false,
    }
}

/// Compare two values with a comparison operator.
fn compare_values(left: &Value, op: &BinaryOp, right: &Value) -> Result<bool> {
    if matches!(left, Value::Missing) || matches!(right, Value::Missing) {
        return Ok(false);
    }

    let left_num = left.as_number();
    let right_num = right.as_number();

    match (left_num, right_num) {
        (Some(l), Some(r)) => {
            let result = match op {
                BinaryOp::Lt => l < r,
                BinaryOp::Gt => l > r,
                BinaryOp::LtEq => l <= r,
                BinaryOp::GtEq => l >= r,
                _ => {
                    return Err(VcfFilterError::EvaluationError(format!(
                        "Unexpected operator in compare_values: {:?}",
                        op
                    )));
                }
            };
            Ok(result)
        }
        _ => {
            // String comparison for non-numeric values
            match (left, right) {
                (Value::String(l), Value::String(r)) => {
                    let result = match op {
                        BinaryOp::Lt => l < r,
                        BinaryOp::Gt => l > r,
                        BinaryOp::LtEq => l <= r,
                        BinaryOp::GtEq => l >= r,
                        _ => {
                            return Err(VcfFilterError::EvaluationError(format!(
                                "Unexpected operator in compare_values: {:?}",
                                op
                            )));
                        }
                    };
                    Ok(result)
                }
                _ => Err(VcfFilterError::TypeMismatch {
                    left: left.type_name().to_string(),
                    right: right.type_name().to_string(),
                }),
            }
        }
    }
}

/// Convert a value to a boolean.
fn value_to_bool(val: &Value) -> Result<bool> {
    match val {
        Value::Bool(b) => Ok(*b),
        Value::Missing => Ok(false),
        Value::String(s) => Ok(!s.is_empty()),
        Value::Number(n) => Ok(*n != 0.0),
        Value::Array(arr) => Ok(!arr.is_empty()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::parse_filter;
    use crate::header::parse_header;
    use crate::row::parse_row;

    const HEADER: &str = r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="Disease name">"#;

    fn eval_filter(filter: &str, row_str: &str, header: &str) -> bool {
        let info_map = parse_header(header).unwrap();
        let row = parse_row(row_str, &info_map).unwrap();
        let expr = parse_filter(filter).unwrap();
        let result = evaluate(&expr, &row, &info_map).unwrap();
        result.as_bool().unwrap_or(false)
    }

    #[test]
    fn test_qual_comparison() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30";
        assert!(eval_filter("QUAL > 30", row, HEADER));
        assert!(eval_filter("QUAL >= 50", row, HEADER));
        assert!(!eval_filter("QUAL > 50", row, HEADER));
    }

    #[test]
    fn test_filter_comparison() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30";
        assert!(eval_filter(r#"FILTER == "PASS""#, row, HEADER));
        assert!(!eval_filter(r#"FILTER == "FAIL""#, row, HEADER));
    }

    #[test]
    fn test_info_field() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30";
        assert!(eval_filter("DP >= 30", row, HEADER));
        assert!(eval_filter("DP == 30", row, HEADER));
        assert!(!eval_filter("DP > 30", row, HEADER));
    }

    #[test]
    fn test_missing_numeric_field_comparison_is_false() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\t.";
        assert!(!eval_filter("DP > 12", row, HEADER));
        assert!(!eval_filter("DP <= 12", row, HEADER));
    }

    #[test]
    fn test_ann_subfield_access() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tANN=G|missense_variant|HIGH|BRCA1|ENSG123|transcript|ENST456|protein_coding|1/10|c.100A>G|p.Thr34Ala|100/500|100/400|34/133||";
        assert!(eval_filter(r#"ANN[0].Gene_Name == "BRCA1""#, row, HEADER));
        assert!(eval_filter(
            r#"ANN[0].Annotation_Impact == "HIGH""#,
            row,
            HEADER
        ));
        assert!(!eval_filter(r#"ANN[0].Gene_Name == "BRCA2""#, row, HEADER));
    }

    #[test]
    fn test_ann_wildcard_access() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tANN=G|missense|LOW|BRCA1|E1|t|T1|pc|1|c.1|p.1|1|1|1||,G|synonymous|HIGH|BRCA2|E2|t|T2|pc|2|c.2|p.2|2|2|2||";
        // Any annotation has HIGH impact
        assert!(eval_filter(
            r#"ANN[*].Annotation_Impact == "HIGH""#,
            row,
            HEADER
        ));
        // Any annotation is for BRCA1
        assert!(eval_filter(r#"ANN[*].Gene_Name == "BRCA1""#, row, HEADER));
        // No annotation has MODERATE impact
        assert!(!eval_filter(
            r#"ANN[*].Annotation_Impact == "MODERATE""#,
            row,
            HEADER
        ));
    }

    #[test]
    fn test_logical_and() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30";
        assert!(eval_filter(r#"QUAL > 30 && DP >= 30"#, row, HEADER));
        assert!(!eval_filter(r#"QUAL > 60 && DP >= 30"#, row, HEADER));
    }

    #[test]
    fn test_logical_or() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30";
        assert!(eval_filter(r#"QUAL > 60 || DP >= 30"#, row, HEADER));
        assert!(!eval_filter(r#"QUAL > 60 || DP > 30"#, row, HEADER));
    }

    #[test]
    fn test_not_operator() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30";
        assert!(eval_filter(r#"!exists(CLNSIG)"#, row, HEADER));
        assert!(!eval_filter(r#"!exists(DP)"#, row, HEADER));
    }

    #[test]
    fn test_exists_function() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30;CLNSIG=Benign";
        assert!(eval_filter("exists(DP)", row, HEADER));
        assert!(eval_filter("exists(CLNSIG)", row, HEADER));
        assert!(!eval_filter("exists(AF)", row, HEADER));
    }

    #[test]
    fn test_dp_numeric_comparison_from_format_when_info_missing() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15";
        assert!(eval_filter("DP > 12", row, HEADER));
        assert!(!eval_filter("DP > 20", row, HEADER));
    }

    #[test]
    fn test_namespaced_info_and_format_resolution() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP\t0/1:15";

        // Unqualified DP keeps default (INFO-first) behavior
        assert!(eval_filter("DP == 30", row, HEADER));

        // Namespaced fields resolve strictly
        assert!(eval_filter("INFO.DP == 30", row, HEADER));
        assert!(eval_filter("FORMAT.DP == 15", row, HEADER));
        assert!(!eval_filter("FORMAT.DP == 30", row, HEADER));
    }

    #[test]
    fn test_namespaced_exists_behavior() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15";

        assert!(!eval_filter("exists(INFO.DP)", row, HEADER));
        assert!(eval_filter("exists(FORMAT.DP)", row, HEADER));
    }

    #[test]
    fn test_namespaced_info_annotation_access() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tANN=G|missense_variant|HIGH|BRCA1|ENSG123|transcript|ENST456|protein_coding|1/10|c.100A>G|p.Thr34Ala|100/500|100/400|34/133||,G|synonymous_variant|LOW|BRCA2|ENSG456|transcript|ENST789|protein_coding|2/10|c.101A>G|p.Thr35Ala|101/500|101/400|35/133||";

        assert!(eval_filter(
            r#"INFO.ANN[0].Gene_Name == "BRCA1""#,
            row,
            HEADER
        ));
        assert!(eval_filter(
            r#"INFO.ANN[*].Annotation_Impact == "HIGH""#,
            row,
            HEADER
        ));
    }

    #[test]
    fn test_format_subfield_chain_is_missing() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15";
        assert!(!eval_filter("FORMAT.DP.X == 15", row, HEADER));
    }

    #[test]
    fn test_contains_operator() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tCLNDN=Breast_cancer_familial";
        assert!(eval_filter(r#"CLNDN contains "cancer""#, row, HEADER));
        assert!(!eval_filter(r#"CLNDN contains "diabetes""#, row, HEADER));
    }

    #[test]
    fn test_complex_expression() {
        let row = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30;ANN=G|missense|HIGH|BRCA1|E1|t|T1|pc|1|c.1|p.1|1|1|1||";
        assert!(eval_filter(
            r#"(QUAL > 30 && FILTER == "PASS") && ANN[0].Annotation_Impact == "HIGH""#,
            row,
            HEADER
        ));
    }
}
