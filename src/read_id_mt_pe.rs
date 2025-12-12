use bit_vec_serde::BitVec;
use xxh3;
use fnv;
use crate::kmer;
use probability::prelude::*;
use rayon::prelude::*;
//use rayon::ThreadPoolBuilder;
use crate::seq;
use needletail::parse_fastx_file;
use std;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::sync::Arc;
use std::time::SystemTime;

pub fn false_prob_map(
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
) -> fnv::FnvHashMap<usize, f64> {
    //first change a colors to accessions map to a accessions to colors map
    let mut accession_color = HashMap::new();
    for (color, accession) in colors_accession {
        accession_color.insert(accession, color);
    }
    let mut false_positive_p = fnv::FnvHashMap::default();
    for (key, value) in ref_kmers_in {
        let color = accession_color.get(key).unwrap();
        false_positive_p.insert(
            *color.to_owned(),
            false_prob(bloom_size as f64, num_hash as f64, *value as f64),
        );
    }
    false_positive_p
}

#[inline]
pub fn bitwise_and(vector_of_bitvectors: &Vec<&BitVec>) -> BitVec {
    let mut first = vector_of_bitvectors[0].to_owned();
    if vector_of_bitvectors.len() == 1 {
        first
    } else {
        for i in 1..vector_of_bitvectors.len() {
            let j = i as usize;
            first.intersect(&vector_of_bitvectors[j]);
        }
        first
    }
}

//change a vector of strings to a comma delimited string
#[inline]
pub fn vec_strings_to_string(vector_in: &Vec<String>) -> String {
    let mut comma_separated = String::new();
    for s in vector_in {
        comma_separated.push_str(&s.to_string());
        comma_separated.push_str(",");
    }
    comma_separated.pop();
    comma_separated
}

pub fn search_index_classic(
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>,
    map: &fnv::FnvHashSet<std::string::String>,
    bloom_size: usize,
    num_hash: usize,
    no_hits_num: usize,
) -> fnv::FnvHashMap<usize, usize> {
    let mut report = fnv::FnvHashMap::default();
    let empty_bitvec = BitVec::from_elem(no_hits_num, false);
    for k in map {
        let mut kmer_slices = Vec::new();
        for i in 0..num_hash {
            let bit_index =
                xxh3::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
            let bi = bit_index as usize;
            if !bigsi_map.contains_key(&bi) || bigsi_map[&bi] == empty_bitvec {
                break;
            } else {
                kmer_slices.push(bigsi_map.get(&bi).unwrap());
            }
        }
        if kmer_slices.len() < num_hash {
            *report.entry(no_hits_num).or_insert(0) += 1;
            break;
        } else {
            let first = bitwise_and(&kmer_slices);
            let mut color = 0;
            for item in first {
                if item {
                    *report.entry(color).or_insert(0) += 1;
                }
                color += 1;
            }
        }
    }
    report
}

pub fn search_index(
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>,
    map: &fnv::FnvHashSet<std::string::String>,
    bloom_size: usize,
    num_hash: usize,
    no_hits_num: usize,
    start_sample: usize,
) -> fnv::FnvHashMap<usize, usize> {
    let mut report = fnv::FnvHashSet::default();
    let mut final_report = fnv::FnvHashMap::default();
    let mut counter = 0;
    for k in map {
        if counter < start_sample {
            let mut kmer_slices = Vec::new();
            for i in 0..num_hash {
                let bit_index =
                    xxh3::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
                match bigsi_map.get(&(bit_index as usize)) {
                    Some(bit_array) => kmer_slices.push(bit_array),
                    None => break,
                };
            }
            if kmer_slices.len() < num_hash {
                *final_report.entry(no_hits_num).or_insert(0) += 1;
                break;
            } else {
                let first = bitwise_and(&kmer_slices);
                let mut color: usize = 0;
                for item in &first {
                    if item {
                        report.insert(color);
                        *final_report.entry(color).or_insert(0) += 1;
                    }
                    color += 1;
                }
            }
        } else {
            let mut kmer_slices = Vec::new();
            for i in 0..num_hash {
                let bit_index =
                    xxh3::hash64_with_seed(&k.as_bytes(), i as u64) % bloom_size as u64;
                match bigsi_map.get(&(bit_index as usize)) {
                    Some(bit_array) => kmer_slices.push(bit_array),
                    None => break,
                };
            }
            if kmer_slices.len() < num_hash {
                *final_report.entry(no_hits_num).or_insert(0) += 1;
                break;
            } else {
                let first = bitwise_and(&kmer_slices);
                for item in &report {
                    if first[*item] {
                        *final_report.entry(*item).or_insert(0) += 1;
                    }
                }
            }
        }
        counter += 1;
    }
    final_report
}

#[inline]
pub fn not_fp_signicant(
    observations: usize,
    child_fp: &fnv::FnvHashMap<usize, f64>,
    fp_correct: f64,
    taxon: usize,
    taxon_hits: usize,
) -> bool {
    let p_false = child_fp.get(&taxon).unwrap();
    let distribution = Binomial::new(observations, *p_false);
    let critical_value = observations as f64 * p_false;
    let mpf = distribution.mass(taxon_hits);
    ((taxon_hits as f64) < critical_value)
        || (((taxon_hits as f64) > critical_value) && (mpf >= fp_correct))
}

//kmer poll classification plus raw count data as output
// 1. filter hits below false prob threshold
// 2. find tophits
// 3 if only one tophit, report with unique flag
pub fn kmer_poll_plus<'a>(
    report: &fnv::FnvHashMap<usize, usize>,
    kmer_length: usize,
    colors_accession: &'a fnv::FnvHashMap<usize, String>,
    child_fp: &fnv::FnvHashMap<usize, f64>,
    no_hits_num: usize,
    fp_correct: f64,
) -> (String, usize, usize, &'a str, usize) {
    let mut count_vec: Vec<_> = report.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    if (count_vec[0].0 == &no_hits_num) && (count_vec.len() == 1) {
        return (
            "no_hits".to_string(),
            0 as usize,
            kmer_length,
            "accept",
            0 as usize,
        );
    };
    let mut significant_hits = Vec::new();
    for t in &count_vec {
        if t.0 == &no_hits_num {
            continue;
        } else if not_fp_signicant(kmer_length, child_fp, fp_correct, *t.0, *t.1) {
            continue;
        } else {
            significant_hits.push(t.to_owned());
        }
    }
    if significant_hits.len() == 0 {
        (
            "no_significant_hits".to_string(),
            0 as usize,
            kmer_length,
            "reject",
            0 as usize,
        )
    } else {
        //significant_hits.sort_by(|a, b| b.1.cmp(a.1));
        let first_tophit = colors_accession[&significant_hits[0].0].to_owned();
        let mut top_hits = Vec::new();
        for h in &significant_hits {
            if h.1 == significant_hits[0].1 {
                top_hits.push(colors_accession[&h.0].to_owned())
            }
        }
        if top_hits.len() == 1 {
            (
                first_tophit.to_owned(),
                significant_hits[0].1.to_owned(),
                kmer_length,
                "accept",
                top_hits.len(),
            )
        } else {
            (
                vec_strings_to_string(&top_hits),
                significant_hits[0].1.to_owned(),
                kmer_length,
                "reject",
                top_hits.len(),
            )
        }
    }
}

pub struct SeqRead {
    pub id: String,       //id including >
    pub seq: Vec<String>, // sequence
}

impl SeqRead {
    pub fn new() -> SeqRead {
        SeqRead {
            id: String::new(),
            seq: Vec::new(),
        }
    }
}

pub struct SeqReadstr<'a> {
    pub id: &'a str,       //id including >
    pub seq: Vec<&'a str>, // sequence
}

impl<'a> SeqReadstr<'a> {
    pub fn new() -> SeqReadstr<'a> {
        SeqReadstr {
            id: "",
            seq: Vec::new(),
        }
    }
}

#[allow(unused_assignments)]
pub fn parallel_vec(
    vec: &Vec<(String, Vec<String>)>,
    bigsi: &fnv::FnvHashMap<usize, BitVec>,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    no_hits_num: usize,
    k: usize,
    m: usize,
    d: usize,
    fp_correct: f64,
    start_sample: usize,
) -> std::vec::Vec<(std::string::String, String, usize, usize, String, usize)> {
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);
    let my_bigsi: Arc<&fnv::FnvHashMap<usize, BitVec>> = Arc::new(bigsi);
    let false_positive_p_arc: Arc<fnv::FnvHashMap<usize, f64>> = Arc::new(false_positive_p);
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            if r.1[0].len() < k {
                (
                    r.0.to_owned(),
                    "too_short".to_string(),
                    0 as usize,
                    0 as usize,
                    "accept".to_string(),
                    0 as usize,
                )
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector_skip_n_set(&r.1, k, d)
                } else {
                    kmer::minimerize_vector_skip_n_set(&r.1, k, m, d)
                };
                let report = if start_sample == 0 {
                    search_index_classic(&child_bigsi, &map, bloom_size, num_hash, no_hits_num)
                } else {
                    search_index(
                        &child_bigsi,
                        &map,
                        bloom_size,
                        num_hash,
                        no_hits_num,
                        start_sample,
                    )
                };
                if report.is_empty() {
                    (
                        r.0.to_owned(),
                        "no_hits".to_string(),
                        0 as usize,
                        map.len(),
                        "accept".to_string(),
                        0 as usize,
                    )
                } else {
                    let classification = kmer_poll_plus(
                        &report,
                        map.len(),
                        &colors_accession,
                        &child_fp,
                        no_hits_num,
                        fp_correct,
                    );
                    (
                        r.0.to_owned(),
                        classification.0,
                        classification.1.to_owned(),
                        classification.2.to_owned(),
                        classification.3.to_owned(),
                        classification.4,
                    )
                }
            }
        })
        .collect();
    out_vec
}

#[allow(unused_assignments)]
pub fn parallel_vec_str(
    vec: &Vec<(&str, Vec<&str>)>,
    bigsi: &fnv::FnvHashMap<usize, BitVec>,
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    no_hits_num: usize,
    k: usize,
    m: usize,
    d: usize,
    fp_correct: f64,
    start_sample: usize,
) -> std::vec::Vec<(std::string::String, String, usize, usize, String, usize)> {
    let false_positive_p = false_prob_map(colors_accession, ref_kmers_in, bloom_size, num_hash);
    let my_bigsi: Arc<&fnv::FnvHashMap<usize, BitVec>> = Arc::new(bigsi);
    let false_positive_p_arc: Arc<fnv::FnvHashMap<usize, f64>> = Arc::new(false_positive_p);
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let child_bigsi = my_bigsi.clone();
            let child_fp = false_positive_p_arc.clone();
            if r.1[0].len() < k {
                (
                    r.0.to_owned(),
                    "too_short".to_string(),
                    0 as usize,
                    0 as usize,
                    "accept".to_string(),
                    0 as usize,
                )
            } else {
                let map = if m == 0 {
                    kmer::kmerize_vector_skip_n_set_str(&r.1, k, d)
                } else {
                    kmer::minimerize_vector_skip_n_set_str(&r.1, k, m, d)
                };
                let report = if start_sample == 0 {
                    search_index_classic(&child_bigsi, &map, bloom_size, num_hash, no_hits_num)
                } else {
                    search_index(
                        &child_bigsi,
                        &map,
                        bloom_size,
                        num_hash,
                        no_hits_num,
                        start_sample,
                    )
                };
                if report.is_empty() {
                    (
                        r.0.to_owned(),
                        "no_hits".to_string(),
                        0 as usize,
                        map.len(),
                        "accept".to_string(),
                        0 as usize,
                    )
                } else {
                    let classification = kmer_poll_plus(
                        &report,
                        map.len(),
                        &colors_accession,
                        &child_fp,
                        no_hits_num,
                        fp_correct,
                    );
                    (
                        r.0.to_owned(),
                        classification.0,
                        classification.1.to_owned(),
                        classification.2.to_owned(),
                        classification.3.to_owned(),
                        classification.4,
                    )
                }
            }
        })
        .collect();
    out_vec
}

#[allow(unused_assignments)]
pub fn stream_fasta(
    filenames: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>, //has to be an Arc ?
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize,
    _t: usize, //threads (unused)
    d: usize,  //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    start_sample: usize,
) {
    let no_hits_num: usize = colors_accession.len();
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");

    let search_time = SystemTime::now();
    let mut read_count: usize = 0;

    let mut vec: Vec<(String, Vec<String>)> = Vec::with_capacity(b);

    // needletail handles gzip/plain FASTA/FASTQ
    let mut reader = parse_fastx_file(filenames[0])
        .unwrap_or_else(|_| panic!("Error reading FASTA/FASTQ file {}", filenames[0]));

    while let Some(record_res) = reader.next() {
        let record = record_res.expect("Invalid FASTA/FASTQ record");

        // mimic "id including >"
        let id = format!(
            ">{}",
            String::from_utf8(record.id().to_vec()).expect("Non UTF-8 id")
        );

        let seq_bytes = record.seq();
        if seq_bytes.is_empty() {
            continue;
        }
        let seq_str = String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence");

        vec.push((id, vec![seq_str]));

        if vec.len() >= b {
            let c = parallel_vec(
                &vec,
                bigsi_map,
                colors_accession,
                ref_kmers_in,
                bloom_size,
                num_hash,
                no_hits_num,
                k,
                m,
                d,
                fp_correct,
                start_sample,
            );
            read_count += c.len();
            eprint!(" {} reads classified\r", read_count);
            for id in c {
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        id.0, id.1, id.2, id.3, id.4, id.5
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
            vec.clear();
        }
    }

    // process remaining records
    if !vec.is_empty() {
        let c = parallel_vec(
            &vec,
            bigsi_map,
            colors_accession,
            ref_kmers_in,
            bloom_size,
            num_hash,
            no_hits_num,
            k,
            m,
            d,
            fp_correct,
            start_sample,
        );
        read_count += c.len();
        for id in c {
            file.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    id.0, id.1, id.2, id.3, id.4, id.5
                )
                .as_bytes(),
            )
            .expect("could not write results!");
        }
    }

    eprint!(" {} reads classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            eprintln!("Error: {:?}", e);
        }
    }
}
/*
#[allow(unused_assignments)]
pub fn stream_fasta_str(
    filenames: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>, //has to be an Arc ?
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize,
    t: usize, //threads
    d: usize, //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    start_sample: usize,
)
{
    let f = File::open(filenames[0]).expect("file not found");
    let iter1 = io::BufReader::new(f).lines();
    let mut vec = Vec::with_capacity(b);
    let no_hits_num: usize = colors_accession.len();
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let mut sub_string = String::new();
    ThreadPoolBuilder::new()
        .num_threads(t)
        .build_global()
        .expect("Can't initialize ThreadPoolBuilder");
    let search_time = SystemTime::now();
    let mut count = 0;
    let mut read_count = 0;
    let mut fasta = SeqReadstr::new();
    for line in iter1 {
        let l = line.unwrap();
        if count == 0 {
            fasta.id = &l;
        } else {
            if l.contains('>') {
                if sub_string.len() > 0 {
                    fasta.seq = vec![&sub_string];
                    vec.push((fasta.id, fasta.seq));
                    fasta.id = &l;
                    sub_string.clear();
                }
            } else {
                sub_string.push_str(&l);
            }
        }
        count += 1;
        if vec.len() % b == 0 {
            let c = parallel_vec_str(
                &vec,
                bigsi_map,
                colors_accession,
                ref_kmers_in,
                bloom_size,
                num_hash,
                no_hits_num,
                k,
                m,
                d,
                fp_correct,
                start_sample,
            );
            read_count += c.len();
            eprint!(" {} reads classified\r", read_count);
            for id in c {
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        id.0, id.1, id.2, id.3, id.4, id.5
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
            vec.clear();
        }
    }
    fasta.seq = vec![&sub_string];
    vec.push((fasta.id, fasta.seq));
    let c = parallel_vec_str(
        &vec,
        bigsi_map,
        colors_accession,
        ref_kmers_in,
        bloom_size,
        num_hash,
        no_hits_num,
        k,
        m,
        d,
        fp_correct,
        start_sample,
    );
    read_count += c.len();
    for id in c {
        file.write_all(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                id.0, id.1, id.2, id.3, id.4, id.5
            )
            .as_bytes(),
        )
        .expect("could not write results!");
    }
    eprint!(" {} reads classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
}
*/

pub fn false_prob(m: f64, k: f64, n: f64) -> f64 {
    let e = std::f64::consts::E;
    (1.0 - e.powf(-((k * (n + 0.5)) / (m - 1.0)))).powf(k)
}

#[allow(unused_assignments)]
pub fn per_read_stream_pe(
    filenames: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>, //has to be an Arc
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize,
    _t: usize, //threads (unused)
    d: usize,  //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    qual_offset: u8,
    start_sample: usize,
) {
    let search_time = SystemTime::now();
    let mut vec: Vec<(String, Vec<String>)> = Vec::with_capacity(b);
    let no_hits_num: usize = colors_accession.len();
    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");
    let mut read_count: usize = 0;

    let mut reader1 = parse_fastx_file(filenames[0])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[0]));
    let mut reader2 = parse_fastx_file(filenames[1])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[1]));

    loop {
        let rec1_opt = reader1.next();
        let rec2_opt = reader2.next();

        let (rec1, rec2) = match (rec1_opt, rec2_opt) {
            (Some(Ok(r1)), Some(Ok(r2))) => (r1, r2),
            _ => break, // stop when either file ends or an error occurs
        };

        let (q1, q2) = match (rec1.qual(), rec2.qual()) {
            (Some(q1), Some(q2)) => (q1, q2),
            _ => continue, // skip pairs without qualities
        };

        // ID from R1; prefix with '@' to match FASTQ style
        let id = format!(
            "@{}",
            String::from_utf8(rec1.id().to_vec()).expect("Non UTF-8 id")
        );

        let seq1_bytes = rec1.seq();
        let seq2_bytes = rec2.seq();
        if seq1_bytes.is_empty() && seq2_bytes.is_empty() {
            continue;
        }

        let seq1_str =
            String::from_utf8(seq1_bytes.to_vec()).expect("Non UTF-8 sequence in R1");
        let seq2_str =
            String::from_utf8(seq2_bytes.to_vec()).expect("Non UTF-8 sequence in R2");
        let qual1_str =
            String::from_utf8(q1.to_vec()).expect("Non UTF-8 quality in R1");
        let qual2_str =
            String::from_utf8(q2.to_vec()).expect("Non UTF-8 quality in R2");

        let masked1 = seq::qual_mask(&seq1_str, &qual1_str, qual_offset);
        let masked2 = seq::qual_mask(&seq2_str, &qual2_str, qual_offset);

        vec.push((id, vec![masked1, masked2]));

        if vec.len() >= b {
            let c = parallel_vec(
                &vec,
                bigsi_map,
                colors_accession,
                ref_kmers_in,
                bloom_size,
                num_hash,
                no_hits_num,
                k,
                m,
                d,
                fp_correct,
                start_sample,
            );
            read_count += c.len();
            eprint!("{} read pairs classified\r", read_count);
            for id in c {
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        id.0, id.1, id.2, id.3, id.4, id.5
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
            vec.clear();
        }
    }

    // process remaining read pairs
    if !vec.is_empty() {
        let c = parallel_vec(
            &vec,
            bigsi_map,
            colors_accession,
            ref_kmers_in,
            bloom_size,
            num_hash,
            no_hits_num,
            k,
            m,
            d,
            fp_correct,
            start_sample,
        );
        read_count += c.len();
        eprint!("{} read pairs classified\r", read_count);
        for id in c {
            file.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    id.0, id.1, id.2, id.3, id.4, id.5
                )
                .as_bytes(),
            )
            .expect("could not write results!");
        }
    }

    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} read pairs in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            eprintln!("Error: {:?}", e);
        }
    }
}

#[allow(unused_assignments)]
pub fn per_read_stream_se(
    filenames: Vec<&str>,
    bigsi_map: &fnv::FnvHashMap<usize, BitVec>, //has to be an Arc ?
    colors_accession: &fnv::FnvHashMap<usize, String>,
    ref_kmers_in: &fnv::FnvHashMap<String, usize>,
    bloom_size: usize,
    num_hash: usize,
    k: usize,
    m: usize,  //0 == no m, otherwise minimizer
    _t: usize, //threads
    d: usize,  //downsample factor
    fp_correct: f64,
    b: usize,
    prefix: &str,
    qual_offset: u8,
    start_sample: usize,
) {
    let search_time = SystemTime::now();
    let mut vec: Vec<(String, Vec<String>)> = Vec::with_capacity(b);
    let no_hits_num: usize = colors_accession.len();
    let mut read_count: usize = 0;

    let mut file =
        File::create(format!("{}_reads.txt", prefix)).expect("could not create outfile!");

    // needletail reader
    let mut reader = parse_fastx_file(filenames[0])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[0]));

    while let Some(record_res) = reader.next() {
        let record = record_res.expect("Invalid FASTQ record");

        let qual = match record.qual() {
            Some(q) => q,
            None => continue, // this function is qual-based
        };

        let id = format!(
            "@{}",
            String::from_utf8(record.id().to_vec()).expect("Non UTF-8 id")
        );

        let seq_bytes = record.seq();
        if seq_bytes.is_empty() {
            continue;
        }

        let seq_str = String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence");
        let qual_str = String::from_utf8(qual.to_vec()).expect("Non UTF-8 quality");

        let masked = seq::qual_mask(&seq_str, &qual_str, qual_offset);

        vec.push((id, vec![masked]));

        if vec.len() >= b {
            let c = parallel_vec(
                &vec,
                bigsi_map,
                colors_accession,
                ref_kmers_in,
                bloom_size,
                num_hash,
                no_hits_num,
                k,
                m,
                d,
                fp_correct,
                start_sample,
            );
            read_count += c.len();
            eprint!("{} reads classified\r", read_count);
            for id in c {
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        id.0, id.1, id.2, id.3, id.4, id.5
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
            vec.clear();
        }
    }

    // last block
    if !vec.is_empty() {
        let c = parallel_vec(
            &vec,
            bigsi_map,
            colors_accession,
            ref_kmers_in,
            bloom_size,
            num_hash,
            no_hits_num,
            k,
            m,
            d,
            fp_correct,
            start_sample,
        );
        read_count += c.len();
        eprint!("{} reads classified\r", read_count);
        for id in c {
            file.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    id.0, id.1, id.2, id.3, id.4, id.5
                )
                .as_bytes(),
            )
            .expect("could not write results!");
        }
    }

    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            eprintln!("Error: {:?}", e);
        }
    }
}
