//! Command-line VCF filter tool.
//!
//! Usage: vcf-filter -filter <expression>
//!
//! Example:
//!   zcat test.vcf.gz | vcf-filter -filter "QUAL > 30 && exists(CLNSIG)" | bgzip -c > out.vcf.gz

use std::io::{self, BufRead, Write};
use vcf_filter::FilterEngine;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    // Parse arguments
    let filter_expr = match parse_args(&args) {
        Ok(expr) => expr,
        Err(msg) => {
            eprintln!("{}", msg);
            std::process::exit(1);
        }
    };

    if let Err(e) = run_filter(&filter_expr) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn parse_args(args: &[String]) -> Result<String, String> {
    if args.len() < 3 {
        return Err(format!(
            "Usage: {} -filter <expression>\n\n\
             Example:\n  \
             zcat test.vcf.gz | {} -filter \"QUAL > 30 && exists(CLNSIG)\" | bgzip -c > out.vcf.gz",
            args[0], args[0]
        ));
    }

    if args[1] != "-filter" && args[1] != "--filter" {
        return Err(format!(
            "Unknown option: {}. Use -filter <expression>",
            args[1]
        ));
    }

    Ok(args[2].clone())
}

fn run_filter(filter_expr: &str) -> Result<(), Box<dyn std::error::Error>> {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut stdout_lock = stdout.lock();

    let mut header_lines = Vec::new();
    let mut engine: Option<FilterEngine> = None;
    let mut passed = 0u64;
    let mut total = 0u64;

    for line_result in stdin.lock().lines() {
        let line = line_result?;

        if line.starts_with('#') {
            // Accumulate header lines
            header_lines.push(line.clone());
            writeln!(stdout_lock, "{}", line)?;

            // When we hit #CHROM, we have the full header
            if line.starts_with("#CHROM") {
                let header_str = header_lines.join("\n");
                engine = Some(FilterEngine::new(&header_str)?);
            }
        } else {
            // Data row
            let eng = engine
                .as_ref()
                .ok_or("No VCF header found before data rows")?;

            total += 1;
            if eng.evaluate(filter_expr, &line)? {
                passed += 1;
                writeln!(stdout_lock, "{}", line)?;
            }
        }
    }

    eprintln!("vcf-filter: {}/{} variants passed filter", passed, total);
    Ok(())
}
