//! Filter expression parser using chumsky.
//!
//! Parses filter expressions like:
//! - `QUAL > 30`
//! - `FILTER == "PASS"`
//! - `ANN[0].Gene_Name == "BRCA1"`
//! - `ANN[*].Annotation_Impact == "HIGH"` (any match)
//! - `DP > 10 && QUAL >= 30`
//! - `CLNSIG == "Benign" || CLNSIG == "Likely_benign"`

use chumsky::prelude::*;

/// Binary operators for comparisons and logic.
#[derive(Debug, Clone, PartialEq)]
pub enum BinaryOp {
    // Comparison
    Eq,       // ==
    NotEq,    // !=
    Lt,       // <
    Gt,       // >
    LtEq,     // <=
    GtEq,     // >=
    Contains, // contains (string contains)

    // Logical
    And, // &&
    Or,  // ||
}

/// Unary operators.
#[derive(Debug, Clone, PartialEq)]
pub enum UnaryOp {
    Not, // !
}

/// Part of a variable access path.
#[derive(Debug, Clone, PartialEq)]
pub enum AccessPart {
    /// A field name (e.g., "ANN", "Gene_Name").
    Field(String),
    /// An array index (e.g., [0], [1]).
    Index(usize),
    /// Wildcard array access (e.g., [*] - matches any).
    Wildcard,
}

/// A filter expression AST node.
#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    /// A numeric literal.
    Number(f64),
    /// A string literal.
    String(String),
    /// A boolean literal.
    Bool(bool),
    /// A variable reference with access path.
    /// e.g., `ANN[0].Gene_Name` becomes `Var([Field("ANN"), Index(0), Field("Gene_Name")])`
    Var(Vec<AccessPart>),
    /// A binary operation.
    Binary(Box<Expr>, BinaryOp, Box<Expr>),
    /// A unary operation.
    Unary(UnaryOp, Box<Expr>),
    /// Check if a field exists (is not missing).
    Exists(Vec<AccessPart>),
}

impl Expr {
    /// Create a simple variable expression.
    pub fn var(name: &str) -> Self {
        Expr::Var(vec![AccessPart::Field(name.to_string())])
    }
}

/// Create the filter expression parser.
pub fn parser() -> impl Parser<char, Expr, Error = Simple<char>> {
    recursive(|full_expr| {
        // Number literal
        let number = text::int(10)
            .chain::<char, _, _>(just('.').chain(text::digits(10)).or_not().flatten())
            .collect::<String>()
            .map(|s| Expr::Number(s.parse().unwrap()))
            .padded();

        // String literal (double-quoted)
        let string = just('"')
            .ignore_then(filter(|c| *c != '"').repeated())
            .then_ignore(just('"'))
            .collect::<String>()
            .map(Expr::String)
            .padded();

        // Boolean literal
        let boolean = choice((
            text::keyword("true").to(Expr::Bool(true)),
            text::keyword("false").to(Expr::Bool(false)),
        ))
        .padded();

        // Identifier (field name)
        let ident = text::ident().padded();

        // Array index: [0], [1], etc.
        let array_index = just('[')
            .ignore_then(
                just('*')
                    .to(AccessPart::Wildcard)
                    .or(text::int(10).map(|s: String| AccessPart::Index(s.parse().unwrap()))),
            )
            .then_ignore(just(']'));

        // Field access: .FieldName
        let field_access = just('.')
            .ignore_then(text::ident().padded())
            .map(|s: String| AccessPart::Field(s));

        // Variable with optional access chain: ANN[0].Gene_Name
        let variable = ident
            .map(|s: String| AccessPart::Field(s))
            .then(choice((array_index.clone(), field_access.clone())).repeated())
            .map(|(first, rest)| {
                let mut parts = vec![first];
                parts.extend(rest);
                Expr::Var(parts)
            })
            .padded();

        // exists(field) function
        let exists_fn = text::keyword("exists")
            .padded()
            .ignore_then(
                just('(')
                    .padded()
                    .ignore_then(
                        text::ident()
                            .padded()
                            .map(|s: String| AccessPart::Field(s))
                            .then(choice((array_index, field_access)).repeated())
                            .map(|(first, rest)| {
                                let mut parts = vec![first];
                                parts.extend(rest);
                                parts
                            }),
                    )
                    .then_ignore(just(')').padded()),
            )
            .map(Expr::Exists);

        // Parenthesized expression (uses full_expr recursively)
        let paren_expr = just('(')
            .padded()
            .ignore_then(full_expr)
            .then_ignore(just(')').padded());

        // Atoms: literals, variables, or parenthesized expressions
        let atom = choice((exists_fn, boolean, number, string, paren_expr, variable));

        // Unary operators (!)
        let unary = just('!')
            .padded()
            .repeated()
            .then(atom)
            .foldr(|_op, expr| Expr::Unary(UnaryOp::Not, Box::new(expr)));

        // Comparison operators
        let cmp_op = choice((
            just("==").to(BinaryOp::Eq),
            just("!=").to(BinaryOp::NotEq),
            just("<=").to(BinaryOp::LtEq),
            just(">=").to(BinaryOp::GtEq),
            just("<").to(BinaryOp::Lt),
            just(">").to(BinaryOp::Gt),
            text::keyword("contains").to(BinaryOp::Contains),
        ))
        .padded();

        // Comparison expressions
        let comparison = unary
            .clone()
            .then(cmp_op.then(unary).repeated())
            .foldl(|left, (op, right)| Expr::Binary(Box::new(left), op, Box::new(right)));

        // Logical AND (&&)
        let and_op = just("&&").padded().to(BinaryOp::And);
        let and_expr = comparison
            .clone()
            .then(and_op.then(comparison).repeated())
            .foldl(|left, (op, right)| Expr::Binary(Box::new(left), op, Box::new(right)));

        // Logical OR (||)
        let or_op = just("||").padded().to(BinaryOp::Or);
        and_expr
            .clone()
            .then(or_op.then(and_expr).repeated())
            .foldl(|left, (op, right)| Expr::Binary(Box::new(left), op, Box::new(right)))
    })
    .then_ignore(end())
}

/// Parse a filter expression string into an AST.
pub fn parse_filter(filter: &str) -> Result<Expr, Vec<Simple<char>>> {
    parser().parse(filter)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_comparison() {
        let expr = parse_filter("QUAL > 30").unwrap();
        assert!(matches!(expr, Expr::Binary(_, BinaryOp::Gt, _)));
    }

    #[test]
    fn test_parse_string_comparison() {
        let expr = parse_filter(r#"FILTER == "PASS""#).unwrap();
        assert!(matches!(expr, Expr::Binary(_, BinaryOp::Eq, _)));
    }

    #[test]
    fn test_parse_array_access() {
        let expr = parse_filter(r#"ANN[0].Gene_Name == "BRCA1""#).unwrap();
        if let Expr::Binary(left, BinaryOp::Eq, _) = expr {
            if let Expr::Var(parts) = *left {
                assert_eq!(parts.len(), 3);
                assert_eq!(parts[0], AccessPart::Field("ANN".to_string()));
                assert_eq!(parts[1], AccessPart::Index(0));
                assert_eq!(parts[2], AccessPart::Field("Gene_Name".to_string()));
            } else {
                panic!("Expected Var");
            }
        } else {
            panic!("Expected Binary");
        }
    }

    #[test]
    fn test_parse_wildcard_access() {
        let expr = parse_filter(r#"ANN[*].Annotation_Impact == "HIGH""#).unwrap();
        if let Expr::Binary(left, BinaryOp::Eq, _) = expr {
            if let Expr::Var(parts) = *left {
                assert_eq!(parts[1], AccessPart::Wildcard);
            } else {
                panic!("Expected Var");
            }
        } else {
            panic!("Expected Binary");
        }
    }

    #[test]
    fn test_parse_logical_and() {
        let expr = parse_filter("QUAL > 30 && DP >= 10").unwrap();
        assert!(matches!(expr, Expr::Binary(_, BinaryOp::And, _)));
    }

    #[test]
    fn test_parse_logical_or() {
        let expr = parse_filter(r#"CLNSIG == "Benign" || CLNSIG == "Likely_benign""#).unwrap();
        assert!(matches!(expr, Expr::Binary(_, BinaryOp::Or, _)));
    }

    #[test]
    fn test_parse_complex_expression() {
        let expr = parse_filter(
            r#"(QUAL > 30 && FILTER == "PASS") || ANN[*].Annotation_Impact == "HIGH""#,
        )
        .unwrap();
        assert!(matches!(expr, Expr::Binary(_, BinaryOp::Or, _)));
    }

    #[test]
    fn test_parse_not_operator() {
        let expr = parse_filter(r#"!exists(DP)"#).unwrap();
        assert!(matches!(expr, Expr::Unary(UnaryOp::Not, _)));
    }

    #[test]
    fn test_parse_exists_function() {
        let expr = parse_filter("exists(CLNSIG)").unwrap();
        assert!(matches!(expr, Expr::Exists(_)));
    }

    #[test]
    fn test_parse_contains() {
        let expr = parse_filter(r#"CLNDN contains "BRCA""#).unwrap();
        assert!(matches!(expr, Expr::Binary(_, BinaryOp::Contains, _)));
    }

    #[test]
    fn test_parse_boolean_literal() {
        let expr = parse_filter("true").unwrap();
        assert_eq!(expr, Expr::Bool(true));
    }

    #[test]
    fn test_parse_number_literal() {
        let expr = parse_filter("42.5").unwrap();
        assert_eq!(expr, Expr::Number(42.5));
    }
}
