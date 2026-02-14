use std::fs;
use vcf_filter::FilterEngine;

fn main() {
    let content = fs::read_to_string("samples/sample.vcf").expect("Failed to read sample.vcf");
    let lines: Vec<&str> = content.lines().collect();

    // Find header lines and data lines
    let header_end = lines.iter().position(|l| l.starts_with("#CHROM")).unwrap();
    let header = lines[..=header_end].join("\n");
    let data_lines = &lines[header_end + 1..];

    // Create filter engine
    let engine = FilterEngine::new(&header).expect("Failed to parse header");

    // Test GT == "0/0" filter
    println!("Testing filter: GT == \"0/0\"");
    println!("----------------------------");

    for line in data_lines {
        if line.is_empty() {
            continue;
        }
        let result = engine.evaluate("GT == \"0/0\"", line);
        match result {
            Ok(val) => {
                // Get position info for context
                let fields: Vec<&str> = line.split('\t').collect();
                let chrom = fields[0];
                let pos = fields[1];
                let id = fields[2];
                let gt = fields.get(9).unwrap_or(&"N/A");
                println!("{}:{} ({}) GT={} -> {}", chrom, pos, id, gt, val);
            }
            Err(e) => println!("Error: {:?}", e),
        }
    }
}
