// use ahash::{AHashMap, AHashSet, HashMap, HashMapExt};
use ahash::{AHashMap, AHashSet, HashMapExt};
use bio::io::{
    fasta,
    fastq::{self, Record, Records},
};
use std::{
    collections::HashMap,
    fs,
    hash::BuildHasherDefault,
    hash::Hash,
    io::BufReader,
    path::{Path, PathBuf},
};
//use bio::io::fastq;
use brotli::{CompressorReader, CompressorWriter};
use clap::{builder::OsStr, Parser, Subcommand};
use colored::*;
use flate2::read::GzDecoder;
use flate2::write::ZlibEncoder;
use flate2::Compression;
use itertools::Itertools;
use nohash_hasher::NoHashHasher;
use nohash_hasher::{BuildNoHashHasher, IntMap};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fs::File;
use std::io;

const OUTPUT_EXT: &str = "prmap";

/*
   TODO: Implement smarter species array, "zeroes" should just be a u64, with
   TODO: the first bit either set or unset. Species u64 has the opposite and
   TODO: only keeps 63 bits for species.
   TODO: Change Vec with ThinVec
*/

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    // #[arg(
    //     short,
    //     long,
    //     value_name = "INT",
    //     help = "A positive integer value less than 33."
    // )]
    // kmersize: usize,
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Map read against an indexed database
    Map {
        #[arg(
            value_name = "FASTQ",
            help = "Currently only supports one pair of fastq files."
        )]
        input: Vec<PathBuf>,

        #[arg(
            short, long, value_name = "PREFIX",
            help = format!("Output prefix. If a directory is given, the file\n\
                            name will be the same as the input file name with an\n\
                            added extension. If a file name is given it will \n\
                            also get added the extension.\n\
                            The extension is: <pair no>_{}.",OUTPUT_EXT)
        )]
        output: Option<PathBuf>,

        #[arg(
            short,
            long,
            value_name = "PATH",
            help = "Path to the template file created with the index subcommand"
        )]
        template: PathBuf,

        #[arg(
            short,
            long,
            value_name = "INT",
            help = "A positive integer value less than 33."
        )]
        kmersize: usize,
    },

    /// Create and indexed database
    Index {
        #[arg(
            short,
            long,
            value_name = "INT",
            help = "A positive integer value less than 33."
        )]
        kmersize: usize,

        #[arg(
            short,
            long,
            value_name = "TSV",
            help = "A GTDB formatted metadata file containing metada on each\n\
                    of the reference genomes provided."
        )]
        metadata: PathBuf,

        #[arg(
            short,
            long,
            value_name = "DIR",
            help = "Path to GTDB files, stored in the GTDB directory structure"
        )]
        gtdb: PathBuf,

        #[arg(
            short,
            long,
            value_name = "PATH",
            help = "Path to where index files should be stored. The main \n\
                    index will be stored as <PATH>, a taxonomy file will be \n\
                    stored as <PATH>_tax. Both files will be stored as \n\
                    compressed binary files."
        )]
        output: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();

    // let max_kmer_bit = max_bits(&cli.kmersize * 2); // TODO twice in code
    // let _kmer_overflow_bits = usize::MAX - max_kmer_bit; // TODO twice in code

    match &cli.command {
        Some(Commands::Map {
            input,
            output,
            template,
            kmersize,
        }) => {
            let mut file = File::open(&template).unwrap();
            let mut file_decomp_input = brotli::Decompressor::new(&mut file, 4096);
            let db: AHashMap<usize, Vec<u64>> =
                bincode::deserialize_from(&mut file_decomp_input).unwrap();

            let mut tax_filename = template.file_name().unwrap().to_owned();
            tax_filename.push("_tax");
            let mut tax_group_file_path: PathBuf = template.clone();
            tax_group_file_path.set_file_name(tax_filename);

            let mut tax_group_file = File::open(&tax_group_file_path).unwrap();
            let tax_group_file_decomp = brotli::Decompressor::new(&mut tax_group_file, 4096);
            let tax_groups: Vec<String> = bincode::deserialize_from(tax_group_file_decomp).unwrap();

            // TODO DEBUG print
            println!("tax_groups: {:?}", tax_groups);

            calc_prob_from_input(&kmersize, &db, &input, &tax_groups);
        }
        Some(Commands::Index {
            kmersize,
            metadata,
            gtdb,
            output,
        }) => {
            println!("Get metadata");
            let metadata_map = match load_metadata(metadata, gtdb) {
                Ok(metadata_loaded) => metadata_loaded,
                Err(_) => unreachable!(),
            };

            // TODO Debug print
            // println!("{:?}", metadata_map);

            // TODO: Mem debug
            let tax_groups: Vec<String> = get_unique_tax_groups(&metadata_map);
            let mut tax_group_file = output.clone();
            let derived_tax_filename = format!(
                "{}{}",
                tax_group_file.file_name().unwrap().to_str().unwrap(),
                "_tax"
            );
            tax_group_file.set_file_name(derived_tax_filename);
            serialize_compress_write(&tax_group_file, &tax_groups);
            // println!("Tax groups: {:?}", tax_groups);

            create_db(&metadata_map, &kmersize, &tax_groups, &output);
        }
        None => {}
    }
}

fn get_unique_tax_groups(metadata_map: &AHashMap<String, MetaEntry>) -> Vec<String> {
    metadata_map
        .values()
        .map(|entry| entry.gtdb_tax.to_string())
        .unique()
        .sorted_unstable()
        .collect()
}

fn calc_prob_from_input(
    k: &usize,
    db: &AHashMap<usize, Vec<u64>>,
    input: &Vec<PathBuf>,
    tax_groups: &Vec<String>,
) -> Result<(), std::io::Error> {
    let max_kmer_bit = max_bits(k * 2);
    let kmer_overflow_bits = usize::MAX - max_kmer_bit;

    // prior is 1 divided by number of species in database + 1 for unkown.
    let prior: f64 = 1.0 / (tax_groups.len() + 1) as f64;

    // TODO debug print
    // println!("Prior: {}", prior);

    let mut record_counter: usize = 0;
    // +1 for unkown species
    let mut species_sum_probs: Vec<f64> = vec![prior; tax_groups.len() + 1];

    // TODO debug print
    // println!("species sum prob 0: {:?}", species_sum_probs);

    // let mut fq_records: Vec<Records<BufReader<GzDecoder<&mut File>>>> = vec![];

    // for fq_file_path in input {
    //     let mut fq_file = fs::File::open(fq_file_path).unwrap();
    //     let decoder = GzDecoder::new(&mut fq_file);
    //     fq_records.push(fastq::Reader::new(decoder).records());
    // }

    let mut fq_records: Vec<Records<BufReader<GzDecoder<File>>>> = input
        .iter()
        .map(|path| fs::File::open(path).unwrap())
        .map(|file| GzDecoder::new(file))
        .map(|decoder| fastq::Reader::new(decoder).records())
        .collect();

    // TODO DEBUG
    // println!("Here comes:3");
    while let Some(Ok(record)) = fq_records[0].next() {
        record_counter += 1;
        let mut kmers: Vec<usize> = vec![]; // !Should this be a tuple?

        kmers.extend(get_kmers(
            record.seq(),
            k,
            &max_kmer_bit,
            &kmer_overflow_bits,
        ));
        // TODO DEBUG
        // println!("Here comes:2");
        match fq_records.get_mut(1) {
            Some(fq_record_rev) => {
                if let Some(Ok(record_rev)) = fq_record_rev.next() {
                    kmers.extend(get_kmers(
                        record_rev.seq(),
                        k,
                        &max_kmer_bit,
                        &kmer_overflow_bits,
                    ));
                } else {
                    eprintln!(
                        "{} FASTQ files do not seem to contain an equal amount of records.",
                        "Error:".bold().red()
                    );
                    std::process::exit(1);
                }
            }
            None => (),
        };

        let read_species_distr: Vec<f64> = annotate_kmers(&kmers, db, prior, tax_groups.len());
        // TODO DEBUG print
        // println!("Here comes:");
        // println!("read_species_distr: {:?}", read_species_distr);
        for (i, prob) in read_species_distr.iter().enumerate() {
            species_sum_probs[i] += *prob;
        }
    }

    for sum_prob in species_sum_probs.iter_mut() {
        *sum_prob /= record_counter as f64;
    }

    // Todo remove debug print
    println!("species sum prob {:?}", species_sum_probs);
    println!("Norm value: {}", record_counter);

    Ok(())
}

fn annotate_kmers(
    kmers: &Vec<usize>,
    db: &AHashMap<usize, Vec<u64>>,
    prior: f64,
    tax_groups_length: usize,
) -> Vec<f64> {
    // index 0 in sample counts are unkown kmers.
    // Indexes >1 are species from species vec.
    let mut sample_counts: Vec<usize> = vec![0; tax_groups_length + 1];

    // TODO debug var
    let mut debug = false;
    // println!("debug set to true");

    for kmer in kmers {
        match db.get(kmer) {
            Some(species_list) => {
                // TODO debug
                if debug == true {
                    println!("Species list: {:?}", species_list);
                }

                for (int_index, species_int) in species_list.iter().enumerate() {
                    let mut bit_index = 0;
                    let mut cache: u64 = species_int.clone(); // Guess ints are copied so no need to clone
                    while cache != 0 {
                        // Uneven cache == 1. Uneven cache -> rightmost bit is 1
                        if (cache & 1) == 1 {
                            sample_counts[int_index + bit_index + 1] += 1;
                        }
                        cache >>= 1;
                        bit_index += 1;
                    }
                }
            }
            None => {
                sample_counts[0] += 1;
            }
        }
        // TODO debug var
        debug = false;
    }

    // for kmer in kmers {
    //     match db.get(kmer) {
    //         Some(species_list) => {
    //             for (int_index, species_int) in species_list.iter().enumerate() {
    //                 let mut bit_index = 0;
    //                 let mut cache: u64 = species_int.clone();
    //                 while cache != 0 {
    //                     cache >>= 1;
    //                     if cache != *species_int {
    //                         sample_counts[int_index + bit_index] += 1;
    //                     }
    //                     bit_index += 1;
    //                 }
    //             }
    //         }
    //         None => {
    //             sample_counts[0] += 1;
    //         }
    //     }
    // }

    // println!("Sample counts: {:?}", &sample_counts);

    calc_probs(&sample_counts, prior)
}

fn calc_probs(sample_counts: &Vec<usize>, prior: f64) -> Vec<f64> {
    let mut sample_distr: Vec<f64> = Vec::with_capacity(sample_counts.len());
    let count_sum: f64 = sample_counts.iter().map(|v| *v as f64 + prior).sum();

    for count in sample_counts {
        let prob: f64 = (*count as f64 + prior) / count_sum as f64;
        sample_distr.push(prob);
    }
    sample_distr
}

fn create_db(
    metadata_map: &AHashMap<String, MetaEntry>,
    k: &usize,
    tax_groups: &Vec<String>,
    output: &PathBuf,
) {
    let max_kmer_bit = max_bits(k * 2);
    let kmer_overflow_bits = usize::MAX - max_kmer_bit;

    //let max_tax_bit_index = (tax_groups.len() as f32 / 64.0).ceil() as usize;
    let max_tax_bit_index = (tax_groups.len() as f32 / 8.0).ceil() as usize; //  ! REMEMBER TO SET BIT SIZE

    // let mut db: AHashMap<usize, AHashSet<u32>> = AHashMap::with_capacity(1000000);
    // let mut db: AHashMap<usize, Vec<u8>> = AHashMap::with_capacity(100000000);

    // let mut db: AHashMap<usize, Vec<u64>> = AHashMap::new();
    let mut db: AHashMap<usize, Vec<u8>> = AHashMap::with_capacity(100_000_000);
    // let mut db: AHashMap<usize, Vec<u64>, BuildHasherDefault<NoHashHasher<u8>>> =
    //     AHashMap::with_hasher(BuildHasherDefault::default());

    // let mut db: HashMap<usize, Vec<u64>, BuildHasherDefault<NoHashHasher<u8>>> =
    //     HashMap::with_capacity_and_hasher(100_000_000, BuildHasherDefault::default());
    // let mut db: HashMap<usize, Vec<u64>, BuildHasherDefault<NoHashHasher<u8>>> =
    //     HashMap::with_hasher(BuildHasherDefault::default());

    // let mut db: IntMap<usize, Vec<u64>> = IntMap::default();
    // let mut db: IntMap<usize, Vec<u64>> = IntMap::default();

    // let mut db: FxHashMap<usize, Vec<u64>> = FxHashMap::with_capacity(100_000_000);

    // let mut db: AHashMap<usize, Vec<u8>> = AHashMap::new(); //  ! REMEMBER TO SET BIT SIZE

    // TODO remove debug print
    println!("DB size start: {:?}", std::mem::size_of_val(&db));
    let mut shared_kmers = 0;

    for metadata in metadata_map.values() {
        // TODO: Not necessarry to find index for each file, only each species
        let tax_8bit_index: (usize, usize) = get_tax_8bit_index(metadata, tax_groups);
        // let tax_8bit_index: (usize, usize) = get_tax_8bit_index(metadata, tax_groups);  //  ! REMEMBER TO SET BIT SIZE
        let mut file = fs::File::open(&metadata.path).unwrap(); // TODO std::io::Error

        let decoder = GzDecoder::new(&mut file);
        let records = fasta::Reader::new(decoder).records();

        // TODO Debug code
        // let mut db_kmer_count = 0;
        // let mut db_cap_count = 0;
        // let mut db_read_count = 0;

        for record in records {
            let kmers: AHashSet<usize> =
                get_unique_kmers(record.unwrap().seq(), k, &max_kmer_bit, &kmer_overflow_bits);

            // TODO rremove
            // db_read_count += 1;
            // db_kmer_count += kmers.len();
            // db_cap_count += kmers.capacity();

            for kmer in kmers {
                // db.entry(kmer).or_insert(vec![0u64; max_tax_bit_index])[tax_8bit_index.0] |=
                //     1 << tax_8bit_index.1;
                db.entry(kmer).or_insert(vec![0u8; max_tax_bit_index])[tax_8bit_index.0] |=  //  ! REMEMBER TO SET BIT SIZE
                    1 << tax_8bit_index.1;

                //if db.contains_key(&kmer) {
                //    shared_kmers += 1;
                //}
            }
        }

        // ! Memory calc - Do this again with "queue" or "smallvec"
        // let test_vec = vec![0u64; max_tax_bit_index];
        // println!("Vec: {:?}", test_vec);
        // println!("Size: {}", std::mem::size_of_val(&test_vec));
        // for i in test_vec {
        //     println!("\tentry size: {:?}", std::mem::size_of_val(&i));
        // }

        // println!(
        //     "kmers pr read: {}",
        //     db_kmer_count as f64 / db_read_count as f64
        // );

        // println!(
        //     "capacity pr read: {}",
        //     db_cap_count as f64 / db_read_count as f64
        // );

        // TODO: Delete old code
        // let dummy: () = records.map(|record| {
        //     get_unique_kmers(record.unwrap().seq(), k, &max_kmer_bit, &kmer_overflow_bits)
        //         .iter()
        //         .map(|kmer| {
        //             db.entry(*kmer)
        //                 .and_modify(|taxids| taxids.extend([metadata.ncbi_taxid]))
        //                 .or_insert(AHashSet::from([metadata.ncbi_taxid]));
        //         });
        // });

        // for kmer_set in kmer_sets {
        //     for kmer in kmer_set {
        //         db.entry(kmer)
        //             .and_modify(|taxids| taxids.extend([metadata.ncbi_taxid]))
        //             .or_insert(AHashSet::from([metadata.ncbi_taxid]));
        //     }
        // }
    }

    // TODO MEM debug

    // TODO remove debug print
    // println!("DB size end:\t{:?}", std::mem::size_of_val(&db));
    // for (kmer, tax_vec) in &db {
    //     println!("{:?}\t{:?}", kmer, tax_vec);
    //     println!("Size of vec:\t{:?}", std::mem::size_of_val(tax_vec));
    //     println!("Size of key:\t{:?}", std::mem::size_of_val(kmer));
    //     let one: Vec<u8> = vec![0u8; 1];
    //     println!("Size vec size 1:\t{:?}", std::mem::size_of_val(&one));
    //     let empty: Vec<u8> = vec![0u8; 0];
    //     println!("Size empty vec:\t{:?}", std::mem::size_of_val(&empty));

    //     let test_ar: [u64; 4] = [0u64; 4];
    //     let mut test_box: Box<[u64; 4]> = Box::new([0u64; 4]);
    //     test_box[0] = 12u64;
    //     test_box[1] = 122u64;
    //     test_box[2] = 12332u64;
    //     test_box[3] = 1234322u64;
    //     println!("Size of ar:\t{:?}", std::mem::size_of_val(&test_ar));
    //     println!("Size of box:\t{:?}", std::mem::size_of_val(&test_box));
    //     break;
    // }
    // println!("Shared kmers: {}", &shared_kmers);

    // TODO: Uncomment index write
    serialize_compress_write(output, &db);

    println!("Number of kmers: {}", db.len());
}

fn serialize_compress_write<T>(output: &PathBuf, object: &T)
where
    T: Serialize,
{
    let mut out_file = File::create(output).unwrap();
    let mut comp = brotli::CompressorWriter::new(
        &mut out_file,
        4096,      // buffer size
        6 as u32,  // quality 0-11
        20 as u32, // lg_window_size
    );
    bincode::serialize_into(&mut comp, object);
}

fn get_tax_8bit_index(meta_entry: &MetaEntry, tax_groups: &Vec<String>) -> (usize, usize) {
    let tax_index = tax_groups
        .iter()
        .find_position(|tax| *tax == &meta_entry.gtdb_tax);

    match tax_index {
        Some((index, _)) => {
            // let sub_index = index % 64;
            let sub_index = index % 8; //  ! REMEMBER TO SET BIT SIZE

            // let bit_index = (index - sub_index) / 64;
            let bit_index = (index - sub_index) / 8; //  ! REMEMBER TO SET BIT SIZE

            return (bit_index, sub_index);
        }
        None => unreachable!(),
    }
}

fn nuc2int(nuc: &u8) -> Result<usize, usize> {
    match nuc {
        b'A' => Ok(0),
        b'a' => Ok(0),
        b'T' => Ok(3),
        b't' => Ok(3),
        b'G' => Ok(2),
        b'g' => Ok(2),
        b'C' => Ok(1),
        b'c' => Ok(1),
        _ => Err(99),
    }
}

/*
    Get the maximum unsigned integer available given b number of bits.
*/
fn max_bits(b: usize) -> usize {
    if (usize::BITS as usize) == b {
        usize::MAX
    } else {
        (1 << b) - 1
    }
}

fn rev_kmer(kmer: usize, k: usize) -> usize {
    let mut v = kmer;
    let mut r = kmer;

    let mut s = 8 * core::mem::size_of::<usize>() - 2;
    v >>= 2;

    while v != 0 {
        // Quit early if there are no more 1s to shift in
        r <<= 2; // Make room for the next significant pair of bits
        r |= v & 2; // Add the bit to the reverse variable at second pos
        r |= v & 1; // Add the bit to the reverse variable at first pos
        v >>= 2; // Go to the next significant pair of bits
        s -= 2; // Decrement the leftover bit count
    }

    // Shift the reversal to the correct position and return the reversal with
    // respect to kmer size.

    r <<= s;

    r >> ((usize::BITS as usize) - k * 2)
}

fn rev_comp(kmer: usize, overflow: usize, k: usize) -> usize {
    let out_kmer = !kmer - overflow;
    rev_kmer(out_kmer, k)
}

fn get_kmers(
    nuc_string: &[u8],
    k: &usize,
    max_kmer_bit: &usize,
    kmer_overflow_bits: &usize,
) -> Vec<usize> {
    let mut kmer_vec: Vec<usize> = Vec::with_capacity(500);
    let mut kmer: usize = 0;
    let mut ignore_count = 0;

    for (nuc_string_index, nuc) in nuc_string.iter().enumerate() {
        kmer <<= 2;

        match nuc2int(nuc) {
            Ok(nuc_int) => kmer += nuc_int,
            Err(_) => {
                ignore_count += k // Ignore nucs as long as kmer contains a non-nuc
            }
        }

        while kmer > *max_kmer_bit {
            kmer -= max_kmer_bit + 1
        }

        if ignore_count > 0 {
            ignore_count -= 1;
            continue;
        }

        if nuc_string_index >= *k {
            let kmer_rev_comp = rev_comp(kmer, *kmer_overflow_bits, *k);
            if kmer <= kmer_rev_comp {
                kmer_vec.push(kmer);
            } else {
                kmer_vec.push(kmer_rev_comp);
            }
        }
    }

    kmer_vec
}

fn get_unique_kmers(
    nuc_string: &[u8],
    k: &usize,
    max_kmer_bit: &usize,
    kmer_overflow_bits: &usize,
) -> AHashSet<usize> {
    let mut kmer_set: AHashSet<usize> = AHashSet::with_capacity(500000);
    let mut kmer: usize = 0;
    let mut ignore_count = 0;

    // TODO remove debug print
    // println!("Seq size: {:?}", nuc_string.len());

    for (nuc_string_index, nuc) in nuc_string.iter().enumerate() {
        kmer <<= 2;

        match nuc2int(nuc) {
            Ok(nuc_int) => kmer += nuc_int,
            Err(_) => {
                ignore_count += k // Ignore nucs as long as kmer contains a non-nuc
            }
        }

        while kmer > *max_kmer_bit {
            kmer -= max_kmer_bit + 1
        }

        if ignore_count > 0 {
            ignore_count -= 1;
            continue;
        }

        if nuc_string_index >= *k {
            let kmer_rev_comp = rev_comp(kmer, *kmer_overflow_bits, *k);
            if kmer <= kmer_rev_comp {
                kmer_set.insert(kmer);
            } else {
                kmer_set.insert(kmer_rev_comp);
            }
        }
    }

    kmer_set
}

#[derive(Debug)]
struct MetaEntry {
    filename: String,
    path: PathBuf,
    ncbi_taxid: u32,
    gtdb_tax: String,
}

fn load_metadata(
    meta_file: &PathBuf,
    gtdb_dir: &PathBuf,
) -> Result<AHashMap<String, MetaEntry>, Box<dyn Error>> {
    // Build the CSV reader
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(meta_file)?;

    type Record = std::collections::HashMap<String, String>;

    // Iterate over each record.
    let mut metadata = AHashMap::new();
    let mut record_count = 0;
    let mut record_stored_count = 0;
    for result in rdr.deserialize() {
        // The iterator yields Result<StringRecord, Error>, so we check the
        // error here.
        let record: Record = result?;
        record_count += 1;
        match create_metadata_record(&record, &gtdb_dir) {
            Some(valid_entry) => {
                let filename: String = valid_entry.filename.clone();
                metadata.insert(filename, valid_entry);
                record_stored_count += 1;
            }
            None => (),
        }
    }
    println!("Records in metadata file: {}", record_count);
    println!("Records with existing files: {}", record_stored_count);
    Ok(metadata)
}

fn create_metadata_record(
    csv_record: &std::collections::HashMap<String, String>,
    gtdb_dir: &PathBuf,
) -> Option<MetaEntry> {
    let accession: Vec<&str> = csv_record["accession"].trim().split(&['_', '.']).collect();
    
    if accession.len() < 4 {
        println!("Weird accession no.: {:?}", accession);
        return None;
    }

    let filename = format!(
        "{}_{}.{}_genomic.fna.gz",
        accession[1], accession[2], accession[3]
    );

    // GTDB paths looks like:
    // <gtdb_dir>GCA/123/456/789/GCA_123456789.1_genomic.fna.gz
    let path: PathBuf = gtdb_dir.join(
        [
            accession[1],
            accession[2].get(..3).unwrap(),
            accession[2].get(3..6).unwrap(),
            accession[2].get(6..).unwrap(),
            &filename,
        ]
        .iter()
        .collect::<PathBuf>(),
    );

    if path.is_file() {
        let ncbi_taxid: u32 = csv_record["ncbi_species_taxid"].trim().parse().unwrap();

        let gtdb_tax_full = csv_record["gtdb_taxonomy"].trim().to_string();
        let gtdb_tax = get_gtdb_taxonomy(&gtdb_tax_full, 's');

        let metadata_record = MetaEntry {
            filename,
            path,
            ncbi_taxid,
            gtdb_tax,
        };

        return Some(metadata_record);
    }

    return None;
}

fn get_gtdb_taxonomy(tax_record: &String, tax_level: char) -> String {
    let tax_pattern = &format!("{}__", tax_level);
    for gtdb_tax in tax_record.split(';') {
        if gtdb_tax.starts_with(tax_pattern) {
            return gtdb_tax[3..].to_string();
        }
    }
    unreachable!();
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn test_nuc2int() {
        let result_a = nuc2int(&b'a');
        let result_ua = nuc2int(&b'A');
        let result_c = nuc2int(&b'c');
        let result_uc = nuc2int(&b'C');
        let result_g = nuc2int(&b'g');
        let result_ug = nuc2int(&b'G');
        let result_t = nuc2int(&b't');
        let result_ut = nuc2int(&b'T');
        let result_err = nuc2int(&b'N');

        assert_eq!(result_a, Ok(0));
        assert_eq!(result_ua, Ok(0));
        assert_eq!(result_c, Ok(1));
        assert_eq!(result_uc, Ok(1));
        assert_eq!(result_g, Ok(2));
        assert_eq!(result_ug, Ok(2));
        assert_eq!(result_t, Ok(3));
        assert_eq!(result_ut, Ok(3));
        assert_eq!(result_err, Err(99));
    }

    #[test]
    fn test_max_bits() {
        let result = max_bits(2 * 16);
        assert_eq!(result, 4_294_967_295);

        let result_max64bit = max_bits(64);
        assert_eq!(result_max64bit, 18_446_744_073_709_551_615);
    }

    #[test]
    fn test_rev_kmer() {
        let k = 5;
        // kmer:     ACTTG -> 00 01 11 11 10 -> 126
        // rev_kmer: GTTCA -> 10 11 11 01 00 -> 756
        let result = rev_kmer(126, k);
        assert_eq!(result, 756);
    }

    #[test]
    fn test_rev_comp() {
        let k = 5;
        // kmer:          ACTTG -> 00 01 11 11 10 -> 126
        // rev_comp_kmer: CAAGT -> 01 00 00 10 11 -> 267
        let max_kmer_bit = max_bits(k * 2);
        let overflow_bits = usize::MAX - max_kmer_bit;
        let result = rev_comp(126, overflow_bits, k);
        assert_eq!(result, 267);
    }

    #[test]
    fn test_get_unique_kmers() {
        let k = 5;
        /*
         * string: ACTTGTC
         *      bit repr: 00 01 11 11 10 11 01
         *
         * kmers: ACTTG CTTGT TTGTC
         *      ACTTG bit: 00 01 11 11 10 -->  126
         *      CTTGT bit: 01 11 11 10 11 -->  507
         *      TTGTC bit: 11 11 10 11 01 --> 1005
         * rev_comp_kmers: CAAGT ACAAG GACAA
         *      CAAGT bit: 01 00 00 10 11 -->  267
         *      ACAAG bit: 00 01 00 00 10 -->   66
         *      GACAA bit: 10 00 01 00 00 -->  528
         *
         * Expect kmers extracted: 126, 66, 528
         */
        let nuc_string: &[u8] = &[b'A', b'C', b'T', b'T', b'G', b'T', b'C'];

        let max_kmer_bit = max_bits(k * 2);
        let kmer_overflow_bits = usize::MAX - max_kmer_bit;

        //let mut kmer_set: AHashSet<usize> = AHashSet::new();

        let kmer_set = get_unique_kmers(nuc_string, &k, &max_kmer_bit, &kmer_overflow_bits);

        assert_eq!(kmer_set.contains(&126), true);
        assert_eq!(kmer_set.contains(&66), true);
        assert_eq!(kmer_set.contains(&528), true);
        assert_eq!(kmer_set.contains(&507), false);
        assert_eq!(kmer_set.contains(&1005), false);
        assert_eq!(kmer_set.contains(&999999), false);
        assert_eq!(kmer_set.len(), 3);
    }
}
