//! Example usage of the vcf-filter library.

use vcf_filter::{FilterEngine, docs};

fn main() {
    // Example VCF header with INFO field definitions
    let header = r#"
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="Disease name">
"#;

    // Create the filter engine
    let engine = FilterEngine::new(header).expect("Failed to parse header");

    // Example VCF row
    let row = "chr1\t186308857\trs3737940\tA\tG\t37.5\tPASS\tDP=30;ANN=G|synonymous_variant|LOW|PRG4|ENSG00000116690|transcript|ENST00000445192.6|protein_coding|7/13|c.3138A>G|p.Pro1046Pro|3183/5044|3138/4215|1046/1404||;CLNSIG=Benign";

    // Test various filters
    let filters = [
        "QUAL > 30",
        r#"FILTER == "PASS""#,
        r#"CLNSIG == "Benign""#,
        r#"ANN[0].Gene_Name == "PRG4""#,
        r#"ANN[0].Annotation_Impact == "LOW""#,
        r#"QUAL > 30 && FILTER == "PASS""#,
        "DP >= 30",
        "exists(CLNSIG)",
        "!exists(LOF)",
    ];

    println!("VCF Filter Engine Demo");
    println!("======================\n");
    println!("Row: {}...\n", &row[..80]);

    for filter in &filters {
        match engine.evaluate(filter, row) {
            Ok(result) => println!("  {} => {}", filter, result),
            Err(e) => println!("  {} => ERROR: {}", filter, e),
        }
    }

    // Demonstrate the docs() function
    println!("\n\nEmbedded Documentation Preview");
    println!("===============================");
    let documentation = docs();
    let preview: String = documentation
        .lines()
        .take(10)
        .collect::<Vec<_>>()
        .join("\n");
    println!("{}", preview);
    println!(
        "...\n(Total documentation length: {} bytes)",
        documentation.len()
    );
}
