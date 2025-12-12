use fnv;
use crate::seq;
use std;
use std::cmp;
use needletail::parse_fastx_file;
use packed_seq::PackedSeqVec;
use packed_seq::SeqVec;
use simd_minimizers::canonical_minimizers;
use seq_hash::NtHasher;

pub fn read_fasta(filename: String) -> Vec<String> {
    // Use needletail so that both plain and compressed FASTA/FASTQ are supported.
    let mut reader = parse_fastx_file(&filename)
        .unwrap_or_else(|_| panic!("Error reading FASTA/FASTQ file {}", filename));
    let mut vec: Vec<String> = Vec::new();

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid FASTA/FASTQ record");
        let seq = record.seq();
        if seq.len() > 1 {
            vec.push(String::from_utf8(seq.to_vec()).expect("Non UTF-8 sequence"));
        }
    }

    vec
}

pub fn read_fasta_mf(filename: String) -> (Vec<String>, Vec<String>) {
    // Same as read_fasta, but also returns the header (label) for each sequence.
    // Needletail automatically handles gzip-compressed input.
    let mut reader = parse_fastx_file(&filename)
        .unwrap_or_else(|_| panic!("Error reading FASTA/FASTQ file {}", filename));

    let mut labels: Vec<String> = Vec::new();
    let mut seqs: Vec<String> = Vec::new();

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid FASTA/FASTQ record");
        // record.id() does NOT include the leading '>'.
        labels.push(String::from_utf8(record.id().to_vec()).expect("Non UTF-8 id"));
        let seq = record.seq();
        if seq.len() > 1 {
            seqs.push(String::from_utf8(seq.to_vec()).expect("Non UTF-8 sequence"));
        }
    }

    (labels, seqs)
}

#[inline]
pub fn kmerize_vector(
    v: &Vec<String>,
    k: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    for l in v {
        if l.len() < k {
            continue;
        } else {
            let length_l = l.len();
            if length_l < k {
                continue;
            } else {
                let l_r = revcomp(&l);
                for i in (0..l.len() - k + 1).step_by(d) {
                    if seq::has_no_n(l[i..i + k].as_bytes()) {
                        if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                            let count = map
                                .entry(l[i..i + k].to_string().to_uppercase())
                                .or_insert(0);
                            *count += 1;
                        } else {
                            let count = map
                                .entry(
                                    l_r[length_l - (i + k)..length_l - i]
                                        .to_string()
                                        .to_uppercase(),
                                )
                                .or_insert(0);
                            *count += 1;
                        }
                    }
                }
            }
        }
    }
    map
}

#[inline]
pub fn kmerize_vector_set(
    v: &Vec<String>,
    k: usize,
    d: usize,
) -> fnv::FnvHashSet<std::string::String> {
    let mut set = fnv::FnvHashSet::default();
    for l in v {
        if l.len() < k {
            continue;
        } else {
            let length_l = l.len();
            if length_l < k {
                continue;
            } else {
                let l_r = revcomp(&l);
                for i in (0..l.len() - k + 1).step_by(d) {
                    if seq::has_no_n(l[i..i + k].as_bytes()) {
                        if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                            set.insert(l[i..i + k].to_string().to_uppercase());
                        } else {
                            set.insert(
                                l_r[length_l - (i + k)..length_l - i]
                                    .to_string()
                                    .to_uppercase(),
                            );
                        }
                    }
                }
            }
        }
    }
    set
}

#[inline]
pub fn kmerize_vector_uppercase(
    v: Vec<String>,
    k: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    for l in v {
        let length_l = l.len();
        let l_r = revcomp(&l);
        for i in (0..l.len() - k + 1).step_by(d) {
            //if i % d == 0 {
            if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                *count += 1;
            } else {
                let count = map
                    .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                    .or_insert(0);
                *count += 1;
            }
            //}
        }
    }
    map
}

#[inline]
pub fn kmerize_vector_skip_n(
    v: &Vec<String>,
    k: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();
    for l in v {
        let length_l = l.len();
        let l_r = revcomp(&l);
        for i in (0..l.len() - k + 1).step_by(d) {
            //if i % d == 0 {
            if seq::has_no_n(l[i..i + k].as_bytes()) {
                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                    let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                    *count += 1;
                } else {
                    let count = map
                        .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                        .or_insert(0);
                    *count += 1;
                }
            } else {
                continue;
            }
            //}
        }
    }
    map
}

#[inline]
pub fn kmerize_vector_skip_n_set(
    v: &Vec<String>,
    k: usize,
    d: usize,
) -> fnv::FnvHashSet<std::string::String> {
    let mut map = fnv::FnvHashSet::default();
    for l in v {
        let length_l = l.len();
        let l_r = revcomp(&l);
        for i in (0..l.len() - k + 1).step_by(d) {
            if seq::has_no_n(l[i..i + k].as_bytes()) {
                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                    map.insert(l[i..i + k].to_string());
                } else {
                    map.insert(l_r[length_l - (i + k)..length_l - i].to_string());
                }
            } else {
                continue;
            }
        }
    }
    map
}

#[inline]
pub fn kmerize_vector_skip_n_set_str(
    v: &Vec<&str>,
    k: usize,
    d: usize,
) -> fnv::FnvHashSet<std::string::String> {
    let mut map = fnv::FnvHashSet::default();
    for l in v {
        let length_l = l.len();
        let l_r = revcomp(&l);
        for i in (0..l.len() - k + 1).step_by(d) {
            if seq::has_no_n(l[i..i + k].as_bytes()) {
                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                    map.insert(l[i..i + k].to_string());
                } else {
                    map.insert(l_r[length_l - (i + k)..length_l - i].to_string());
                }
            } else {
                continue;
            }
        }
    }
    map
}

#[inline]
pub fn kmerize_string(l: String, k: usize) -> Option<fnv::FnvHashMap<std::string::String, usize>> {
    let mut map = fnv::FnvHashMap::default();
    let length_l = l.len();
    let l_r = revcomp(&l);
    if length_l < k {
        None
    } else {
        Some({
            for i in 0..l.len() - k + 1 {
                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                    let count = map
                        .entry(l[i..i + k].to_string().to_uppercase())
                        .or_insert(0);
                    *count += 1;
                } else {
                    let count = map
                        .entry(
                            l_r[length_l - (i + k)..length_l - i]
                                .to_string()
                                .to_uppercase(),
                        )
                        .or_insert(0);
                    *count += 1;
                }
            }
            map
        })
    }
}

pub fn minimerize_vector(
    v: Vec<String>,
    k: usize,
    m: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();

    for l in v {
        let len = l.len();
        if len < k || len < m {
            continue;
        }

        let bytes = l.as_bytes();
        let packed = PackedSeqVec::from_ascii(bytes);

        // Number of m-mers
        let n_mers = bytes.len().saturating_sub(m).saturating_add(1);
        if n_mers == 0 {
            continue;
        }

        // Window size in m-mer space corresponding to a k-long window
        let w = if k > m { k - m + 1 } else { 1 };

        // Canonical nucleotide hasher
        let hasher: NtHasher<true, 7> = NtHasher::new(m);

        let mut minimizer_positions: Vec<u32> = Vec::new();

        let _vals: Vec<u64> = canonical_minimizers(m, w)
            .hasher(&hasher)
            .run(packed.as_slice(), &mut minimizer_positions)
            .values_u64()
            .collect();

        for (idx, &pos_u32) in minimizer_positions.iter().enumerate() {
            // respect downsampling stride d in window index space
            if d > 1 && (idx % d != 0) {
                continue;
            }

            let pos = pos_u32 as usize;
            if pos + m > len {
                continue;
            }

            let minmer = &l[pos..pos + m];

            // This version did NOT skip N originally, so we keep all minimizers.
            *map.entry(minmer.to_string()).or_insert(0) += 1;
        }
    }

    map
}

pub fn minimerize_vector_skip_n(
    v: &Vec<String>,
    k: usize,
    m: usize,
    d: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();

    for l in v {
        let len = l.len();
        if len < k || len < m {
            continue;
        }

        let bytes = l.as_bytes();
        let packed = PackedSeqVec::from_ascii(bytes);

        let n_mers = bytes.len().saturating_sub(m).saturating_add(1);
        if n_mers == 0 {
            continue;
        }

        let w = if k > m { k - m + 1 } else { 1 };

        let hasher: NtHasher<true, 7> = NtHasher::new(m);

        let mut minimizer_positions: Vec<u32> = Vec::new();

        let _vals: Vec<u64> = canonical_minimizers(m, w)
            .hasher(&hasher)
            .run(packed.as_slice(), &mut minimizer_positions)
            .values_u64()
            .collect();

        for (idx, &pos_u32) in minimizer_positions.iter().enumerate() {
            if d > 1 && (idx % d != 0) {
                continue;
            }

            let pos = pos_u32 as usize;
            if pos + m > len {
                continue;
            }

            let minmer = &l[pos..pos + m];

            // Skip minimizers containing N/n
            if minmer
                .as_bytes()
                .iter()
                .any(|&b| b == b'N' || b == b'n')
            {
                continue;
            }

            *map.entry(minmer.to_ascii_uppercase()).or_insert(0) += 1;
        }
    }

    map
}

pub fn minimerize_vector_skip_n_set(
    v: &Vec<String>,
    k: usize,
    m: usize,
    d: usize,
) -> fnv::FnvHashSet<std::string::String> {
    let mut set = fnv::FnvHashSet::default();

    for l in v {
        let len = l.len();
        if len < k || len < m {
            continue;
        }

        let bytes = l.as_bytes();
        let packed = PackedSeqVec::from_ascii(bytes);

        let n_mers = bytes.len().saturating_sub(m).saturating_add(1);
        if n_mers == 0 {
            continue;
        }

        let w = if k > m { k - m + 1 } else { 1 };

        let hasher: NtHasher<true, 7> = NtHasher::new(m);

        let mut minimizer_positions: Vec<u32> = Vec::new();

        let _vals: Vec<u64> = canonical_minimizers(m, w)
            .hasher(&hasher)
            .run(packed.as_slice(), &mut minimizer_positions)
            .values_u64()
            .collect();

        for (idx, &pos_u32) in minimizer_positions.iter().enumerate() {
            if d > 1 && (idx % d != 0) {
                continue;
            }

            let pos = pos_u32 as usize;
            if pos + m > len {
                continue;
            }

            let minmer = &l[pos..pos + m];

            if minmer
                .as_bytes()
                .iter()
                .any(|&b| b == b'N' || b == b'n')
            {
                continue;
            }

            set.insert(minmer.to_ascii_uppercase());
        }
    }

    set
}

pub fn minimerize_vector_skip_n_set_str(
    v: &Vec<&str>,
    k: usize,
    m: usize,
    d: usize,
) -> fnv::FnvHashSet<std::string::String> {
    let mut set = fnv::FnvHashSet::default();

    for l in v {
        // l: &&str; deref to &str
        let l: &str = *l;
        let len = l.len();
        if len < k || len < m {
            continue;
        }

        let bytes = l.as_bytes();
        let packed = PackedSeqVec::from_ascii(bytes);

        let n_mers = bytes.len().saturating_sub(m).saturating_add(1);
        if n_mers == 0 {
            continue;
        }

        let w = if k > m { k - m + 1 } else { 1 };

        let hasher: NtHasher<true, 7> = NtHasher::new(m);

        let mut minimizer_positions: Vec<u32> = Vec::new();

        let _vals: Vec<u64> = canonical_minimizers(m, w)
            .hasher(&hasher)
            .run(packed.as_slice(), &mut minimizer_positions)
            .values_u64()
            .collect();

        for (idx, &pos_u32) in minimizer_positions.iter().enumerate() {
            if d > 1 && (idx % d != 0) {
                continue;
            }

            let pos = pos_u32 as usize;
            if pos + m > len {
                continue;
            }

            let minmer = &l[pos..pos + m];

            if minmer
                .as_bytes()
                .iter()
                .any(|&b| b == b'N' || b == b'n')
            {
                continue;
            }

            set.insert(minmer.to_ascii_uppercase());
        }
    }

    set
}

pub fn kmers_from_fq(filename: String, k: usize) -> fnv::FnvHashMap<String, usize> {
    // Use needletail to read (possibly gzipped) FASTQ/FASTA and extract canonical k-mers.
    let mut map = fnv::FnvHashMap::default();
    let mut reader = parse_fastx_file(&filename)
        .unwrap_or_else(|_| panic!("Error reading FASTQ/FASTA file {}", filename));

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid FASTQ/FASTA record");
        let seq = record.seq();
        if seq.len() < k {
            continue;
        }

        // The original implementation worked on ASCII strings; we keep that behaviour.
        let l = String::from_utf8(seq.to_vec()).expect("Non UTF-8 sequence");
        let length_l = l.len();
        if length_l < k {
            continue;
        }
        let k1 = k / 2;

        for i in 0..=(length_l - k) {
            let end = i + k;
            let substring = &l[i..end];
            if substring.contains('N') {
                continue;
            }

            let substring_rev = revcomp(substring);

            // Choose canonical orientation (lexicographically smallest half, then tie-break).
            let canonical = match substring[..k1].cmp(&substring_rev[..k1]) {
                cmp::Ordering::Less => substring.to_string(),
                cmp::Ordering::Greater => substring_rev,
                cmp::Ordering::Equal => {
                    if substring[k1..] <= substring_rev[k1..] {
                        substring.to_string()
                    } else {
                        substring_rev
                    }
                }
            };

            let counter = map.entry(canonical).or_insert(0);
            *counter += 1;
        }
    }

    map
}

pub fn kmers_from_fq_qual(
    filename: String,
    k: usize,
    _d: usize,
    qual_offset: u8,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();

    // needletail handles plain / gz FASTQ automatically
    let mut reader = parse_fastx_file(&filename)
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filename));

    let d = 1usize; // original code hard-coded d = 1

    while let Some(record_res) = reader.next() {
        let record = record_res.expect("Invalid FASTQ record");

        // This function is quality-based; skip records without qualities (e.g., FASTA).
        let qual = match record.qual() {
            Some(q) => q,
            None => continue,
        };

        let seq_bytes = record.seq();
        if seq_bytes.is_empty() {
            continue;
        }

        let seq_str = String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence");
        let qual_str = String::from_utf8(qual.to_vec()).expect("Non UTF-8 quality string");

        // Apply your existing quality mask
        let masked = seq::qual_mask(&seq_str, &qual_str, qual_offset);
        let l = masked;
        let length_l = l.len();
        if length_l < k {
            continue;
        }

        let l_r = revcomp(&l);

        for i in 0..=length_l - k {
            if i % d != 0 {
                continue;
            }
            if !seq::has_no_n(l[i..i + k].as_bytes()) {
                continue;
            }

            if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                *count += 1;
            } else {
                let count = map
                    .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                    .or_insert(0);
                *count += 1;
            }
        }
    }

    map
}

pub fn kmers_from_fq_minimizer(
    filename: String,
    k: usize,
    m: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();

    // needletail handles gz / plain FASTQ/FASTA
    let mut reader = parse_fastx_file(&filename)
        .unwrap_or_else(|_| panic!("Error reading FASTQ/FASTA file {}", filename));

    while let Some(record_res) = reader.next() {
        let record = record_res.expect("Invalid FASTA/FASTQ record");
        let seq_bytes = record.seq();
        if seq_bytes.len() < m || seq_bytes.len() < k {
            continue;
        }

        // Keep using String so the rest of your code style is consistent
        let seq_str = String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence");
        let l = seq_str;
        if l.len() < m || l.len() < k {
            continue;
        }

        let bytes = l.as_bytes();
        let packed = PackedSeqVec::from_ascii(bytes);

        // number of m-mers
        let n_mers = bytes.len() - m + 1;
        if n_mers == 0 {
            continue;
        }

        // approximate "per k-mer minimizer" by using window size in m-mer space
        let w = if k > m { k - m + 1 } else { 1 };

        // canonical nucleotide hasher
        let hasher: NtHasher<true, 7> = NtHasher::new(m);

        let mut minimizer_positions: Vec<u32> = Vec::new();

        let _vals: Vec<u64> = canonical_minimizers(m, w)
            .hasher(&hasher)
            .run(packed.as_slice(), &mut minimizer_positions)
            .values_u64()
            .collect();

        for &pos_u32 in &minimizer_positions {
            let pos = pos_u32 as usize;
            if pos + m > l.len() {
                continue;
            }

            let minmer = &l[pos..pos + m];

            // optional but usually desirable: skip minimizers with N/n
            if minmer
                .as_bytes()
                .iter()
                .any(|&b| b == b'N' || b == b'n')
            {
                continue;
            }

            *map.entry(minmer.to_string()).or_insert(0) += 1;
        }
    }

    map
}

pub fn kmers_fq_pe(
    filenames: Vec<&str>,
    k: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();

    for filename in filenames {
        let mut reader = parse_fastx_file(filename)
            .unwrap_or_else(|_| panic!("Error reading FASTQ/FASTA file {}", filename));

        while let Some(record_res) = reader.next() {
            let record = record_res.expect("Invalid FASTA/FASTQ record");

            let seq_bytes = record.seq();
            if seq_bytes.len() < k {
                continue;
            }

            let l = String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence");
            let length_l = l.len();
            if length_l < k {
                continue;
            }

            let l_r = revcomp(&l);

            for i in 0..=length_l - k {
                if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                    let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                    *count += 1;
                } else {
                    let count = map
                        .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                        .or_insert(0);
                    *count += 1;
                }
            }
        }
    }

    map
}

pub fn kmers_fq_pe_qual(
    filenames: Vec<&str>,
    k: usize,
    d: usize,
    qual_offset: u8,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();

    // filenames[0] = R1, filenames[1] = R2
    let mut reader1 = parse_fastx_file(filenames[0])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[0]));
    let mut reader2 = parse_fastx_file(filenames[1])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[1]));

    // helper closure to process *one* read (either R1 or R2)
    let mut process_read = |seq_bytes: &[u8], qual_bytes: &[u8]| {
        let seq_str =
            String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence in PE read");
        let qual_str =
            String::from_utf8(qual_bytes.to_vec()).expect("Non UTF-8 quality string in PE read");

        let masked = seq::qual_mask(&seq_str, &qual_str, qual_offset);
        let l = masked;
        let length_l = l.len();
        if length_l < k {
            return;
        }

        let l_r = revcomp(&l);

        for i in 0..=length_l - k {
            if i % d != 0 {
                continue;
            }
            if !seq::has_no_n(l[i..i + k].as_bytes()) {
                continue;
            }

            if l[i..i + k] < l_r[length_l - (i + k)..length_l - i] {
                let count = map.entry(l[i..i + k].to_string()).or_insert(0);
                *count += 1;
            } else {
                let count = map
                    .entry(l_r[length_l - (i + k)..length_l - i].to_string())
                    .or_insert(0);
                *count += 1;
            }
        }
    };

    loop {
        let rec1_opt = reader1.next();
        let rec2_opt = reader2.next();

        let (rec1, rec2) = match (rec1_opt, rec2_opt) {
            (Some(Ok(r1)), Some(Ok(r2))) => (r1, r2),
            _ => break, // stop when either file ends or error
        };

        let (q1, q2) = match (rec1.qual(), rec2.qual()) {
            (Some(q1), Some(q2)) => (q1, q2),
            _ => continue, // skip pairs without qualities
        };

        process_read(&rec1.seq(), q1);
        process_read(&rec2.seq(), q2);
    }

    map
}

pub fn kmers_fq_pe_minimizer(
    filenames: Vec<&str>,
    k: usize,
    m: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map = fnv::FnvHashMap::default();

    for filename in filenames {
        let mut reader = parse_fastx_file(filename)
            .unwrap_or_else(|_| panic!("Error reading FASTQ/FASTA/PE file {}", filename));

        while let Some(record_res) = reader.next() {
            let record = record_res.expect("Invalid FASTA/FASTQ record");
            let seq_bytes = record.seq();
            if seq_bytes.len() < m || seq_bytes.len() < k {
                continue;
            }

            let seq_str = String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence");
            let l = seq_str;
            if l.len() < m || l.len() < k {
                continue;
            }

            let bytes = l.as_bytes();
            let packed = PackedSeqVec::from_ascii(bytes);

            let n_mers = bytes.len() - m + 1;
            if n_mers == 0 {
                continue;
            }

            let w = if k > m { k - m + 1 } else { 1 };

            let hasher: NtHasher<true, 7> = NtHasher::new(m);

            let mut minimizer_positions: Vec<u32> = Vec::new();

            let _vals: Vec<u64> = canonical_minimizers(m, w)
                .hasher(&hasher)
                .run(packed.as_slice(), &mut minimizer_positions)
                .values_u64()
                .collect();

            for &pos_u32 in &minimizer_positions {
                let pos = pos_u32 as usize;
                if pos + m > l.len() {
                    continue;
                }

                let minmer = &l[pos..pos + m];

                if minmer
                    .as_bytes()
                    .iter()
                    .any(|&b| b == b'N' || b == b'n')
                {
                    continue;
                }

                *map.entry(minmer.to_string()).or_insert(0) += 1;
            }
        }
    }

    map
}

pub fn kmers_fq_pe_minimizer_qual(
    filenames: Vec<&str>,
    k: usize,
    m: usize,
    d: usize,
    qual_offset: u8,
) -> fnv::FnvHashMap<std::string::String, usize> {
    use needletail::parse_fastx_file;

    let mut map = fnv::FnvHashMap::default();

    // filenames[0] = R1, filenames[1] = R2
    let mut reader1 = parse_fastx_file(filenames[0])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[0]));
    let mut reader2 = parse_fastx_file(filenames[1])
        .unwrap_or_else(|_| panic!("Error reading FASTQ file {}", filenames[1]));

    // shared closure: process one read (seq + qual)
    let mut process_read = |seq_bytes: &[u8], qual_bytes: &[u8]| {
        if seq_bytes.len() < m || seq_bytes.len() < k {
            return;
        }

        let seq_str =
            String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence (PE)");
        let qual_str =
            String::from_utf8(qual_bytes.to_vec()).expect("Non UTF-8 quality string (PE)");

        let masked = crate::seq::qual_mask(&seq_str, &qual_str, qual_offset);
        if masked.len() < m || masked.len() < k {
            return;
        }

        let l = masked;
        let bytes = l.as_bytes();
        let packed = PackedSeqVec::from_ascii(bytes);

        let n_mers = bytes.len() - m + 1;
        if n_mers == 0 {
            return;
        }

        let w = if k > m { k - m + 1 } else { 1 };
        let hasher: NtHasher<true, 7> = NtHasher::new(m);

        let mut minimizer_positions: Vec<u32> = Vec::new();

        let _vals: Vec<u64> = canonical_minimizers(m, w)
            .hasher(&hasher)
            .run(packed.as_slice(), &mut minimizer_positions)
            .values_u64()
            .collect();

        for (_idx, &pos_u32) in minimizer_positions.iter().enumerate() {
            let pos = pos_u32 as usize;
            if pos + m > l.len() {
                continue;
            }

            // downsample in *position* space or in window index space; here: position
            if d > 1 && (pos % d != 0) {
                continue;
            }

            let minmer = &l[pos..pos + m];

            if minmer
                .as_bytes()
                .iter()
                .any(|&b| b == b'N' || b == b'n')
            {
                continue;
            }

            *map.entry(minmer.to_string()).or_insert(0) += 1;
        }
    };

    loop {
        let rec1_opt = reader1.next();
        let rec2_opt = reader2.next();

        let (rec1, rec2) = match (rec1_opt, rec2_opt) {
            (Some(Ok(r1)), Some(Ok(r2))) => (r1, r2),
            _ => break, // stop when either file ends or an error occurs
        };

        let (q1, q2) = match (rec1.qual(), rec2.qual()) {
            (Some(q1), Some(q2)) => (q1, q2),
            _ => continue,
        };

        // NOTE: seq() / qual() return Cow<[u8]>, so use .as_ref() to get &[u8]
        let s1 = rec1.seq();
        let s2 = rec2.seq();

        process_read(s1.as_ref(), q1.as_ref());
        process_read(s2.as_ref(), q2.as_ref());
    }

    map
}

pub fn kmers_from_fq_minimizer_qual(
    filename: String,
    k: usize,
    m: usize,
    _d: usize,
    qual_offset: u8,
) -> fnv::FnvHashMap<std::string::String, usize> {
    // Map: minimizer (String) -> count
    let mut map = fnv::FnvHashMap::default();

    // needletail handles plain / gzipped FASTQ/FASTA
    let mut reader = parse_fastx_file(&filename)
        .unwrap_or_else(|_| panic!("Error reading FASTQ/FASTA file {}", filename));

    while let Some(record_res) = reader.next() {
        let record = record_res.expect("Invalid FASTA/FASTQ record");

        // We only do quality-based masking when qualities exist (FASTQ).
        // For FASTA (no qualities), we skip – same semantics as "this function is for *qual*".
        let qual = match record.qual() {
            Some(q) => q,
            None => continue,
        };

        let seq_bytes = record.seq();
        if seq_bytes.is_empty() {
            continue;
        }

        // Convert to String so we can reuse your existing qual_mask logic.
        let seq_str = String::from_utf8(seq_bytes.to_vec()).expect("Non UTF-8 sequence");
        let qual_str = String::from_utf8(qual.to_vec()).expect("Non UTF-8 quality string");

        // Apply your quality masking (this returns a String).
        let masked = seq::qual_mask(&seq_str, &qual_str, qual_offset);
        let l = masked;

        // Need at least k for k-mers and at least m for m-mer minimizers.
        if l.len() < k || l.len() < m {
            continue;
        }

        let bytes = l.as_bytes();

        // Pack sequence for SIMD minimizer computation.
        let packed_seq = PackedSeqVec::from_ascii(bytes);

        // Number of m-mers in this masked read.
        let n_mers = bytes.len() - m + 1;
        if n_mers == 0 {
            continue;
        }

        // Window size in m-mers.
        // To approximate the original behavior (one minimizer per k-mer),
        // use w = k - m + 1, but guard against k < m.
        let w = if k > m { k - m + 1 } else { 1 };

        // Canonical NtHasher for nucleotides.
        let hasher: NtHasher<true, 7> = NtHasher::new(m);

        // Positions (in m-mer coordinates) of minimizers.
        let mut minimizer_positions: Vec<u32> = Vec::new();

        // Drive the SIMD minimizer pipeline; we don't need hash values themselves,
        // but we must collect them to consume the iterator.
        let _vals: Vec<u64> = canonical_minimizers(m, w)
            .hasher(&hasher)
            .run(packed_seq.as_slice(), &mut minimizer_positions)
            .values_u64()
            .collect();

        // For each minimizer position, extract the m-mer string and count it.
        for &pos_u32 in &minimizer_positions {
            let pos = pos_u32 as usize;
            if pos + m > l.len() {
                continue;
            }

            let minmer = &l[pos..pos + m];

            // Preserve the original "ignore N" semantics.
            if minmer
                .as_bytes()
                .iter()
                .any(|&b| b == b'N' || b == b'n')
            {
                continue;
            }

            let entry = map.entry(minmer.to_string()).or_insert(0);
            *entry += 1;
        }
    }

    map
}

pub fn clean_map(
    map: fnv::FnvHashMap<std::string::String, usize>,
    t: usize,
) -> fnv::FnvHashMap<std::string::String, usize> {
    let mut map_clean = fnv::FnvHashMap::default();
    for (key, value) in map {
        if value > t {
            map_clean.insert(key, value);
        }
    }
    map_clean
}

pub fn revcomp(dna: &str) -> String {
    let mut rc_dna = String::with_capacity(dna.len());
    for c in dna.chars().rev() {
        rc_dna.push(switch_base(&c))
    }
    rc_dna
}

fn switch_base(c: &char) -> char {
    match c {
        'a' => 't',
        'c' => 'g',
        't' => 'a',
        'g' => 'c',
        'u' => 'a',
        'n' => 'n',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        'N' => 'N',
        _ => 'N',
    }
}

//auto cutoff inference from Zam Iqbal's Cortex
pub fn auto_cutoff(map: fnv::FnvHashMap<std::string::String, usize>) -> usize {
    let mut histo_map = fnv::FnvHashMap::default();
    let mut max_cov = 0;
    for (_key, value) in &map {
        if value > &max_cov {
            max_cov = *value;
        }
        *histo_map.entry(value).or_insert(0) += 1;
    }
    //eprintln!("length histo {}", histo_map.len());
    //not in original code of Zam: do not filter when mean coverage is close to 1
    let mut sum = 0;
    for (i, p) in &histo_map {
        sum += *i * p;
    }
    let num_kmers_total = map.len();
    let total_mean: f64 = sum as f64 / num_kmers_total as f64;
    if total_mean < 1.5 {
        0
    } else {
        let mut count_vec = Vec::with_capacity(max_cov);
        for c in 1..max_cov {
            if !histo_map.contains_key(&c) {
                count_vec.push((c, 0));
            } else {
                count_vec.push((c, *histo_map.get(&c).unwrap()));
            }
        }
        /*
        let mut count_vec: Vec<_> = histo_map.iter().collect();
        count_vec.sort_by(|&(a, _), &(b, _)| a.cmp(&b));
        */
        let mut coverages: Vec<usize> = Vec::with_capacity(count_vec.len());
        for v in count_vec {
            coverages.push(v.1);
        }
        //first pseudo-derivative
        let mut d1 = Vec::new();
        for i in 1..coverages.len() - 1 {
            d1.push(coverages[i] as f64 / coverages[i + 1] as f64);
        }
        //second pseudo-derivative
        let mut d2 = Vec::new();
        for i in 0..d1.len() - 1 {
            d2.push(d1[i] / d1[i + 1]);
        }
        let mut first_pos_d1 = 0;
        let mut first_pos_d2 = 0;
        let threshold: f64 = 1.0;
        for (i, p) in d1.iter().enumerate() {
            if p < &threshold {
                first_pos_d1 = i + 1;
                break;
            }
        }
        for (i, p) in d2.iter().enumerate() {
            if p < &threshold {
                first_pos_d2 = i + 1;
                break;
            }
        }
        //estimate coverage (mean), exclude singleton k-mers
        let mut bigsum = 0;
        for (i, p) in coverages[1..].iter().enumerate() {
            bigsum += i * p;
        }
        let num_kmers: usize = coverages[1..].iter().sum();
        let mean: f64 = bigsum as f64 / num_kmers as f64;
        if (first_pos_d1 > 0) && ((first_pos_d1 as f64) < (mean * 0.75)) {
            first_pos_d1
        } else if first_pos_d2 > 0 {
            first_pos_d2
        } else {
            cmp::max(1, (mean / 2.0).ceil() as usize)
        }
    }
}

/* pseudominimizer based on murmurhash
pub fn find_minimizer(kmer: &str, m: usize) -> String {
    let mut minimizer = kmer[0..m].to_string();
    let mut value = murmur_hash64a(minimizer.as_bytes(), 0);
    for i in 1..kmer.len() - m + 1 {
        let alt_value = murmur_hash64a(kmer[i..i + m].as_bytes(), 0);
        if alt_value < value {
            minimizer = kmer[i..i + m].to_string();
            value = alt_value;
        }
    }
    minimizer
}*/
/*
pub fn find_minimizer(kmer: &str, m: usize) -> String {
    let kmer_r = revcomp(&kmer);
    let length = kmer_r.len();
    let mut minimizers = Vec::new();
    for i in 1..kmer.len() - m + 1 {
        minimizers.push(kmer[i..i + m].to_string());
        minimizers.push(kmer_r[length - (i + m)..length - i].to_string());
    }
    minimizers.sort();
    minimizers[0].to_string()
}
*/
pub fn find_minimizer(seq: &str, m: usize) -> String {
    // Edge cases: empty or shorter than m
    if m == 0 || seq.len() < m {
        return String::new();
    }

    // Pack the sequence (ACGT) for SIMD minimizer
    let bytes = seq.as_bytes();
    let packed_seq = PackedSeqVec::from_ascii(bytes);

    // Number of k-mers in the sequence
    let n_kmers = bytes.len() - m + 1;

    // Single window that spans the whole sequence in terms of k-mers.
    // This gives you one *global* canonical minimizer for the sequence.
    let w = n_kmers;

    // Nucleotide hasher for canonical minimizers
    let hasher = <seq_hash::NtHasher>::new(m);

    // Positions of minimizers will be pushed into this vector
    let mut minimizer_positions: Vec<u32> = Vec::new();

    // Run the SIMD minimizer scan; we don't actually need the hash values,
    // only the positions, but we must drive the iterator.
    let _vals: Vec<u64> = canonical_minimizers(m, w)
        .hasher(&hasher)
        .run(packed_seq.as_slice(), &mut minimizer_positions)
        .values_u64()
        .collect();

    if minimizer_positions.is_empty() {
        return String::new();
    }

    // Take the first (and only, given our window choice) minimizer position
    let pos = minimizer_positions[0] as usize;
    seq[pos..pos + m].to_string()
}
