use needletail::parse_fastx_file;
use flate2::write::GzEncoder;
use flate2::Compression;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::prelude::*;

pub fn tab_to_map(
    filename: String,
    query: &str,
    _accept: &str,
) -> std::collections::HashMap<std::string::String, String> {
    let mut map = HashMap::new();
    let f = File::open(filename).expect("classification file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        let h: Vec<&str> = v[0].split(' ').collect();
        if v[1].contains(query)
        /* && (v[4] == accept)*/
        {
            map.insert(String::from(h[0]), String::from(v[1]));
        }
    }
    map
}

#[allow(unused_assignments)]
pub fn read_filter_pe(
    class_map: std::collections::HashMap<std::string::String, String>,
    filenames: Vec<&str>,
    query: &str,
    prefix: &str,
    exclude: bool,
) {
    let query_cleaned = query.replace(" ", "_");

    // Counters
    let mut excluded: usize = 0;
    let mut included: usize = 0;

    // Output gzipped FASTQs
    let fq1 = File::create(format!("{}_{}_R1.fq.gz", prefix, query_cleaned))
        .expect("could not create R1!");
    let mut gz1 = GzEncoder::new(fq1, Compression::default());

    let fq2 = File::create(format!("{}_{}_R2.fq.gz", prefix, query_cleaned))
        .expect("could not create R2!");
    let mut gz2 = GzEncoder::new(fq2, Compression::default());

    // needletail readers for PE files (R1, R2)
    let mut reader1 = parse_fastx_file(filenames[0])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[0]));
    let mut reader2 = parse_fastx_file(filenames[1])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[1]));

    loop {
        let rec1_opt = reader1.next();
        let rec2_opt = reader2.next();

        let (rec1, rec2) = match (rec1_opt, rec2_opt) {
            (Some(Ok(r1)), Some(Ok(r2))) => (r1, r2),
            _ => break, // EOF or error in either
        };

        // Require qualities (FASTQ)
        let (q1, q2) = match (rec1.qual(), rec2.qual()) {
            (Some(q1), Some(q2)) => (q1, q2),
            _ => continue,
        };

        // Header lines: just "@ID"
        let id1 = String::from_utf8(rec1.id().to_vec()).expect("Non UTF-8 id in R1");
        let header1 = format!("@{}", id1);

        let id2 = String::from_utf8(rec2.id().to_vec()).expect("Non UTF-8 id in R2");
        let header2 = format!("@{}", id2);

        // Sequences
        let seq1 = String::from_utf8(rec1.seq().to_vec()).expect("Non UTF-8 sequence in R1");
        let seq2 = String::from_utf8(rec2.seq().to_vec()).expect("Non UTF-8 sequence in R2");

        // Qualities
        let qual1 = String::from_utf8(q1.to_vec()).expect("Non UTF-8 quality in R1");
        let qual2 = String::from_utf8(q2.to_vec()).expect("Non UTF-8 quality in R2");

        // Key is first token of header1 (like before)
        let h_tokens: Vec<&str> = header1.split(' ').collect();
        let key = h_tokens[0]; // e.g. "@READ_ID"

        if exclude {
            if !class_map.contains_key(key) {
                gz1.write_all(format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes())
                    .expect("could not write R1!");
                gz2.write_all(format!("{}\n{}\n+\n{}\n", header2, seq2, qual2).as_bytes())
                    .expect("could not write R2!");
                excluded += 1;
            }
        } else {
            if class_map.contains_key(key) {
                gz1.write_all(format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes())
                    .expect("could not write R1!");
                gz2.write_all(format!("{}\n{}\n+\n{}\n", header2, seq2, qual2).as_bytes())
                    .expect("could not write R2!");
                included += 1;
            }
        }
    }

    gz1.finish().expect("Could not close new R1 file");
    gz2.finish().expect("Could not close new R2 file");

    if exclude {
        eprintln!(
            "Excluded {} read pairs with classification containing '{}' from output files",
            excluded, query
        );
    } else {
        eprintln!(
            "Wrote {} read-pairs with classification containing '{}' to output files",
            included, query
        );
    }
}

//quicker to adjust the function to se
#[allow(unused_assignments)]
pub fn read_filter_se(
    class_map: std::collections::HashMap<std::string::String, String>,
    filenames: Vec<&str>,
    query: &str,
    prefix: &str,
    exclude: bool,
) {
    let query_cleaned = query.replace(" ", "_");

    let mut excluded: usize = 0;
    let mut included: usize = 0;

    let fq1 =
        File::create(format!("{}_{}.fq.gz", prefix, query_cleaned)).expect("could not create R1!");
    let mut gz1 = GzEncoder::new(fq1, Compression::default());

    // Single-end FASTQ (or FASTA, but we require qualities)
    let mut reader1 = parse_fastx_file(filenames[0])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[0]));

    while let Some(record_res) = reader1.next() {
        let record = match record_res {
            Ok(r) => r,
            Err(_) => continue,
        };

        let qual = match record.qual() {
            Some(q) => q,
            None => continue, // this function is explicitly quality-based FASTQ
        };

        // Header: just "@ID"
        let id1 = String::from_utf8(record.id().to_vec()).expect("Non UTF-8 id");
        let header1 = format!("@{}", id1);

        let seq1 = String::from_utf8(record.seq().to_vec()).expect("Non UTF-8 sequence");
        let qual1 = String::from_utf8(qual.to_vec()).expect("Non UTF-8 quality");

        let h_tokens: Vec<&str> = header1.split(' ').collect();
        let key = h_tokens[0]; // e.g. "@READ_ID"

        if exclude {
            if !class_map.contains_key(key) {
                gz1.write_all(format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes())
                    .expect("Could not write forward read(-s) to file");
                excluded += 1;
            }
        } else {
            if class_map.contains_key(key) {
                gz1.write_all(format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes())
                    .expect("Could not write reverse read(-s) to file");
                included += 1;
            }
        }
    }

    gz1.finish().expect("Could not close new read file");

    if exclude {
        eprintln!(
            "Excluded {} read pairs  with classification containing '{}' from output files",
            excluded, query
        );
    } else {
        eprintln!(
            "Wrote {} read-pairs with classification containing '{}' to output files",
            included, query
        );
    }
}