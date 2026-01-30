# VCF Filter - Copilot Instructions

## Project Overview

A Rust library for filtering VCF (Variant Call Format) bioinformatics files using expressive filter DSL. Parses VCF headers for metadata, then evaluates filter expressions against data rows.

## Architecture

```
FilterEngine (lib.rs)           <- Public API: new(header) → evaluate(filter, row)
    ├── header.rs               <- Parses ##INFO lines, extracts subfield names from descriptions
    ├── filter.rs               <- Chumsky 0.9 parser: filter DSL → Expr AST
    ├── row.rs                  <- Parses tab-separated VCF rows into VcfRow
    ├── eval.rs                 <- Evaluates Expr against VcfRow using InfoMap
    ├── value.rs                <- Value enum: String|Number|Bool|Array|Missing
    └── error.rs                <- thiserror-based VcfFilterError variants
```

**Data flow**: Header string → `InfoMap` (field metadata) → Row string + InfoMap → `VcfRow` → Filter string → `Expr` AST → `evaluate()` → `Value::Bool`

## Dependencies

- **chumsky 0.9**: Filter expression parser (NOT nom - header/row parsing uses pure Rust string operations)
- **thiserror 2**: Error type derivation
- **pretty_assertions** (dev): Better test diffs

### Why chumsky 0.9 (not 0.12+)

Chumsky 0.9 is intentionally pinned despite newer versions existing:

1. **API stability** — Chumsky 1.0+ underwent major breaking changes (0.10 → 0.11 → 0.12) with different `Parser` trait signatures, combinator names, and error handling
2. **Documentation** — Most tutorials and examples online use 0.9 patterns
3. **Simpler recursion** — The 0.9 `recursive()` pattern is more straightforward for expression grammars
4. **Battle-tested** — 0.9 is widely used and stable

Upgrading to 0.12+ would require rewriting the entire parser with new combinator syntax.

Note: We intentionally avoided nom for header/row parsing. VCF format is simple enough that manual string parsing is cleaner and avoids nom 8's breaking API changes.

## Key Patterns

### Parser Construction (filter.rs)

Uses Chumsky 0.9 recursive parser combinator pattern. The entire parser is wrapped in `recursive()` to support parenthesized expressions containing full filter expressions:

```rust
pub fn parser() -> impl Parser<char, Expr, Error = Simple<char>> {
    recursive(|full_expr| {
        // Define atoms, operators, precedence levels inside the recursive block
        // paren_expr uses full_expr to allow nested expressions
        let paren_expr = just('(').padded().ignore_then(full_expr).then_ignore(just(')').padded());
        // ... build up to or_expr
    })
    .then_ignore(end())
}
```

When adding new syntax:

```rust
// Pattern: define atom parsers, combine with choice(), apply precedence via recursive()
let new_feature = text::keyword("keyword").ignore_then(...).map(Expr::NewVariant);
let atom = choice((number, string, boolean, exists_fn, new_feature, variable, paren_expr));
```

### Header Parsing (header.rs)

Uses pure Rust string parsing (no nom). Key function `parse_info_attrs()` handles quoted values in `##INFO=<...>` lines.

Subfield extraction from descriptions looks for pipe-separated formats in single quotes:

- `"Format: 'Gene_Name | Gene_ID | ...'"` → extracts subfield names
- Subfield names are normalized: spaces→`_`, dots→`_`, slashes→`_`

The `InfoNumber` enum includes a `Flag` variant for `Number=0` fields.

### Structured Annotation Access

VCF annotations (ANN, LOF, NMD) use pipe-separated subfields extracted from header descriptions. Access pattern: `ANN[index].SubfieldName` or `ANN[*].SubfieldName` for wildcard.

Internally, annotations are stored as `Value::Array(Vec<Value::Array>)` - outer array is annotations, inner array is subfields indexed by position.

### Value Comparison (eval.rs)

- Wildcard `[*]` returns `Value::Array`, comparisons use `any()` semantics for `==`, `all()` for `!=`
- `Missing` values propagate through comparisons (similar to SQL NULL)
- Type coercion: strings parse to numbers when compared with numeric literals
- `contains` operator only works on String values

### Error Handling

All public APIs return `Result<T, VcfFilterError>`. Use specific error variants:

```rust
VcfFilterError::HeaderParseError(_)   // Invalid ##INFO syntax
VcfFilterError::FilterParseError(_)   // Invalid filter DSL
VcfFilterError::EvaluationError(_)    // Runtime evaluation failure
VcfFilterError::UnknownField(_)       // Field not in header
VcfFilterError::InvalidIndex { .. }   // Array index out of bounds
VcfFilterError::TypeMismatch { .. }   // Incompatible types in comparison
```

## Build & Test

```bash
cargo build                    # Build library and example binary
cargo test                     # Run all tests (47 total: 44 unit + 3 doctests)
cargo run                      # Run main.rs demo showing filter evaluation
```

Tests are co-located in each module plus integration tests in `lib.rs` using real VCF data constants `FULL_HEADER` and `REAL_ROW`.

## Filter DSL Reference

| Syntax     | Example                    | Notes                                 |
| ---------- | -------------------------- | ------------------------------------- |
| Comparison | `QUAL > 30`, `DP >= 10`    | Works on Number values                |
| Equality   | `FILTER == "PASS"`         | String comparison                     |
| Logical    | `A && B`, `A \|\| B`, `!A` | Short-circuit evaluation              |
| Contains   | `CLNDN contains "cancer"`  | Substring match                       |
| Exists     | `exists(LOF)`              | True if field present and non-missing |
| Index      | `ANN[0].Gene_Name`         | First annotation's subfield           |
| Wildcard   | `ANN[*].Impact == "HIGH"`  | Any annotation matches                |
| Grouping   | `(A \|\| B) && C`          | Parentheses for precedence            |

### Built-in Fields

Access without INFO lookup: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`

## Adding New Features

1. **New operator**: Add variant to `BinaryOp`/`UnaryOp` in filter.rs, parser rule in `cmp_op`, eval case in `evaluate_binary()`
2. **New function**: Follow `exists()` pattern - add `Expr::FunctionName(args)` variant, parser rule, eval case
3. **New field type**: Extend `Value` enum, add `as_X()` converter, update `type_name()` for errors

## Common Pitfalls

1. **Chumsky recursive parser**: The entire expression parser must be inside `recursive()` for parentheses to work with full expressions, not just atoms
2. **Doctests with multiline strings**: Use `concat!()` macro instead of raw strings in doc comments to avoid whitespace issues
3. **ANN subfield count**: Header description must have ≥2 pipe-separated subfields to be recognized as structured
4. **Filter keyword clashes**: `contains` is a keyword - identifiers like `contains_count` work, but `contains` alone triggers the operator
