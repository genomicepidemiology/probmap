// TODO: Not necessarry to find index for each file, only each species

use ahash::{AHashMap, AHashSet, HashMapExt};
use bio::io::{
    fasta,
    fastq::{self, Record, Records},
};
use brotli::{CompressorReader, CompressorWriter};
use clap::{builder::OsStr, Parser, Subcommand};
use colored::*;
use fast_log::{plugin::file, Config};
use flate2::read::GzDecoder;
use flate2::write::ZlibEncoder;
use flate2::Compression;
use itertools::{enumerate, Itertools};
use log::{error, info, warn};
use nohash_hasher::NoHashHasher;
use nohash_hasher::{BuildNoHashHasher, IntMap};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use smallvec::{smallvec, SmallVec};
use std::{
    collections::HashMap,
    error::Error,
    fs,
    fs::File,
    hash::BuildHasherDefault,
    hash::Hash,
    io::BufReader,
    path::{Path, PathBuf},
    sync::mpsc,
    thread::{self, JoinHandle},
};
use std::{io::BufWriter, usize};

const IO_BLOCKSIZE: usize = 16_777_216;
const OUTPUT_EXT: &str = "probmap";
static DB_PREFIX: [usize; 16] = [0, 15, 4, 11, 8, 7, 9, 3, 5, 14, 1, 10, 13, 2, 6, 12];

// Optimal distro for pairs
// static DB_PREFIX: [usize; 16] = [0, 15, 4, 11, 8, 7, 12, 3, 1, 14, 5, 10, 9, 6, 12, 2];

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

        #[arg(long, value_name = "INT", help = "Number of thread to spawn (max).")]
        threads: usize,
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

        #[arg(
            short,
            long,
            value_name = "INT",
            help = "Number of thread to spawn (max)."
        )]
        threads: usize,
    },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Map {
            input,
            output,
            template,
            threads,
        }) => {
            println!("Opening file: {:?}", &template);
            let mut info_file = File::open(&template).unwrap();
            let index_info: IndexInfo = bincode::deserialize_from(info_file).unwrap();
            // let mut file_decomp_input = brotli::Decompressor::new(&mut file, 4096);
            println!("{:?}", index_info);

            // TODO handle no output given or output is a directory.
            let log_file = create_file_path(&output.clone().unwrap(), "", "log");
            fast_log::init(
                Config::new()
                    .file(log_file.to_str().unwrap())
                    .chan_len(Some(10000)),
            )
            .unwrap();

            // let db: AHashMap<usize, Vec<u8>> =
            //     bincode::deserialize_from(&mut file_decomp_input).unwrap();

            // let mut tax_filename = template.file_name().unwrap().to_owned();
            // tax_filename.push("_tax");
            // let mut tax_group_file_path: PathBuf = template.clone();
            // tax_group_file_path.set_file_name(tax_filename);

            // let mut tax_group_file = File::open(&tax_group_file_path).unwrap();
            // let tax_group_file_decomp = brotli::Decompressor::new(&mut tax_group_file, 4096);
            // let tax_groups: Vec<String> = bincode::deserialize_from(tax_group_file_decomp).unwrap();

            let mut tax_file_path: PathBuf = template.parent().unwrap().to_owned();
            tax_file_path.push(index_info.tax_file);
            let mut tax_file = File::open(tax_file_path).unwrap();
            let mut tax_groups: Vec<String> = bincode::deserialize_from(tax_file).unwrap();
            tax_groups.insert(0, "Unknown".to_string());

            // TODO DEBUG print
            // println!("tax_groups: {:?}", tax_groups);

            let species_probs = match calc_prob_from_input(
                &index_info.k,
                &index_info.index_files,
                &input,
                &tax_groups,
                &threads,
            ) {
                Ok(species_sum_probs) => species_sum_probs,
                Err(_) => {
                    eprintln!("{} Uknown.", "Error:".bold().red());
                    std::process::exit(1);
                }
            };

            let most_probable_indeces = get_largest_values(&species_probs, 20);
            println!("most_probable_indeces: {:?}", &most_probable_indeces);
            //let mut most_probable_species: Vec<(&String, f64)> = vec![];
            let most_probable_species: Vec<(&String, f64)> = most_probable_indeces
                .iter()
                .map(|(prob, index)| (&tax_groups[*index], *prob))
                .collect();

            for (species, prob) in most_probable_species {
                log::info!("{}\t{}", prob, species);
            }
        }
        Some(Commands::Index {
            kmersize,
            metadata,
            gtdb,
            output,
            threads,
        }) => {
            let tax_group_file = create_file_path(&output, "tax", "probmap");
            let log_file = create_file_path(&output, "", "log");
            let info_file = create_file_path(&output, "", "probmap");

            fast_log::init(
                Config::new()
                    .file(log_file.to_str().unwrap())
                    .chan_len(Some(10000)),
            )
            .unwrap();

            log::info!("Loading metadata.");

            let metadata_map = match load_metadata(metadata, gtdb) {
                Ok(metadata_loaded) => metadata_loaded,
                Err(_) => unreachable!(),
            };

            log::info!("Metadata Loaded.");

            let tax_groups: Vec<String> = get_unique_tax_groups(&metadata_map);
            serialize_compress_write(&tax_group_file, &tax_groups);

            let index_files: Vec<PathBuf> =
                create_db(&metadata_map, *kmersize, &tax_groups, &output, threads);

            let index_info = IndexInfo {
                k: *kmersize,
                tax_file: tax_group_file
                    .file_name()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .to_owned(),
                index_files,
            };
            serialize_compress_write(&info_file, &index_info);

            println!("{:?}", index_info);
        }
        None => {}
    }
}

fn create_file_path(base_file: &PathBuf, append: &str, ext: &str) -> PathBuf {
    let mut file_path = base_file.clone();

    if append != "" {
        let derived_filename = format!(
            "{}_{}",
            file_path.file_name().unwrap().to_str().unwrap(),
            append
        );
        file_path.set_file_name(derived_filename);
    }

    file_path.set_extension(ext);
    file_path
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
    // db: &AHashMap<usize, Vec<u8>>,
    index_files: &Vec<PathBuf>,
    input: &Vec<PathBuf>,
    tax_groups: &Vec<String>,
    nthreads: &usize,
) -> Result<Vec<f64>, std::io::Error> {
    log::info!("Start mapping.");
    // NOTES:
    // let mut db: AHashMap<usize, Vec<u16>>

    let max_kmer_bit = max_bits(k * 2);
    let kmer_overflow_bits = usize::MAX - max_kmer_bit;

    // prior is 1 divided by number of species in database + 1 for unkown.
    let prior: f64 = 1.0 / (tax_groups.len()) as f64;

    let mut record_counter: usize = 0;
    // +1 for unkown species
    // let mut species_sum_probs: Vec<f64> = vec![prior; tax_groups.len() + 1];
    let mut species_sum_probs: Vec<f64> = vec![0.0; tax_groups.len()];

    // let mut kmers: AHashSet<usize> = AHashSet::new();
    // TODO: Remove this temp prefix vec. Should be found in info file!
    let temp_db_prefix: Vec<SmallVec<[usize; 16]>> = vec![
        smallvec![0],
        smallvec![15],
        smallvec![4],
        smallvec![11],
        smallvec![8],
        smallvec![7],
        smallvec![9],
        smallvec![3],
        smallvec![5],
        smallvec![14],
        smallvec![1],
        smallvec![10],
        smallvec![13],
        smallvec![2],
        smallvec![6],
        smallvec![12],
    ];

    for (i, file_path) in enumerate(index_files) {
        log::info!("Mapping to file {:?}", &file_path);
        let db_prefix_chunk = &temp_db_prefix[i];

        let index_file = File::open(file_path).unwrap();
        let kmer_db: AHashMap<usize, Vec<u16>> = bincode::deserialize_from(index_file).unwrap();

        let mut fq_records: Vec<Records<BufReader<GzDecoder<File>>>> = input
            .iter()
            .map(|path| fs::File::open(path).unwrap())
            .map(|file| GzDecoder::new(file))
            .map(|decoder| fastq::Reader::new(decoder).records())
            .collect();

        let mut kmers: Vec<usize> = vec![];

        while let Some(fastq_chunk) = load_chunk_of_fastq(&mut fq_records, 1_000_000) {
            if i == 0 {
                record_counter += fastq_chunk.len();
                // log::info!("Mapping: {}", &record_counter);
            }

            kmers.extend(
                fastq_chunk
                    .par_iter()
                    .map(|record| get_kmers(&record, k, &max_kmer_bit, &kmer_overflow_bits))
                    .flatten()
                    .filter(|kmer| !ignore_kmer(*kmer, db_prefix_chunk))
                    .collect::<Vec<usize>>(),
            );
        }

        let total_kmers = kmers.len();
        log::info!("Kmers total {}", total_kmers);

        let species_distr: Vec<f64> =
            annotate_kmers(kmers, total_kmers, &kmer_db, prior, tax_groups.len());

        for (i, prob) in species_distr.iter().enumerate() {
            species_sum_probs[i] += *prob;
        }
    }

    for sum_prob in species_sum_probs.iter_mut() {
        *sum_prob /= index_files.len() as f64;
    }

    // Todo remove debug print
    log::info!("End mapping.");
    // println!("species sum prob {:?}", species_sum_probs);
    println!("Reads value: {}", record_counter);

    Ok(species_sum_probs)
}

fn get_largest_values(data: &Vec<f64>, count: usize) -> Vec<(f64, usize)> {
    // Create a vector of tuples with values and their original indices
    let mut sorted_data_with_indices: Vec<(f64, usize)> =
        data.iter().enumerate().map(|(i, &val)| (val, i)).collect();

    // Sort the vector of tuples in descending order based on values
    // let mut sorted_data_with_indices = data_with_indices.clone();
    sorted_data_with_indices.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    // Take the first 'count' elements
    let largest_values_with_indices: Vec<(f64, usize)> =
        sorted_data_with_indices.into_iter().take(count).collect();

    largest_values_with_indices
}

fn load_chunk_of_fastq(
    fq_records: &mut Vec<Records<BufReader<GzDecoder<File>>>>,
    n: usize,
) -> Option<Vec<Vec<u8>>> {
    let mut fastq_chunk: Vec<Vec<u8>> = vec![];
    let mut record_counter = 0;

    while let Some(Ok(record)) = fq_records[0].next() {
        let seq = record.seq().into_iter().map(|i| *i).collect_vec();
        fastq_chunk.push(seq);

        match fq_records.get_mut(1) {
            Some(fq_record_rev) => {
                if let Some(Ok(record_rev)) = fq_record_rev.next() {
                    let seq_rev = record_rev.seq().into_iter().map(|i| *i).collect_vec();
                    fastq_chunk.push(seq_rev);
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

        record_counter += 1;
        if record_counter == n {
            break;
        }
    }

    if fastq_chunk.len() == 0 {
        log::info!("No more records");
        return None;
    }

    Some(fastq_chunk)
}

fn check_elements_exist<T>(input: T)
where
    T: IntoParallelIterator<Item = usize>,
{
    input.into_par_iter();
}

fn annotate_kmers<T>(
    // kmers: &AHashSet<usize>,
    kmers: T,
    total_kmer_count: usize,
    //db: &AHashMap<usize, Vec<u8>>,
    db: &AHashMap<usize, Vec<u16>>,
    prior: f64,
    tax_groups_length: usize,
) -> Vec<f64>
where
    T: IntoParallelIterator<Item = usize>,
{
    // index 0 in sample counts are unkown kmers.
    // Indexes >=1 are species from species vec.
    let mut sample_counts: Vec<usize> = vec![0; tax_groups_length];

    let tax_lists: Vec<_> = kmers
        .into_par_iter()
        .filter_map(|kmer| match &db.get(&kmer) {
            Some(tax_list) => Some(tax_list.clone()),
            None => None,
        })
        .collect();

    let known_kmers = tax_lists.len();

    for tax_index in tax_lists.into_iter().flatten() {
        sample_counts[*tax_index as usize + 1] += 1;
    }

    let unknown_kmers = total_kmer_count - known_kmers;
    sample_counts[0] += unknown_kmers;

    // for kmer in kmers {
    //     match db.get(&kmer) {
    //         // ! ERROR Ignores unknown kmers
    //         None => (),
    //         Some(tax_list) => {
    //             for tax_index in tax_list {
    //                 sample_counts[*tax_index as usize + 1] += 1;
    //             }
    //         }
    //     }
    // }

    // 8bit sparse but complete tax vector implementation
    // for kmer in kmers {
    //     match db.get(kmer) {
    //         Some(tax_bit_encoded_list) => {
    //             let tax_list = extract_tax_from_list(tax_bit_encoded_list);
    //             for tax_index in &tax_list {
    //                 sample_counts[tax_index + 1] += 1; // + 1 as 0 index is unknown.
    //             }
    //         }
    //         None => {
    //             sample_counts[0] += 1;
    //         }
    //     }
    // }

    calc_probs(&sample_counts, prior)
}

fn extract_tax_from_list(tax_list: &Vec<u8>) -> Vec<usize> {
    const BYTE_SIZE: usize = 8;
    let mut tax_index_count = 0;
    let mut zero_count_switch = false;
    let mut tax_output: Vec<usize> = vec![];

    for tax_entry in tax_list {
        match tax_entry {
            0 => {
                // Entry is 0 byte (not a zero counter)
                if zero_count_switch == false {
                    tax_index_count += BYTE_SIZE;
                    zero_count_switch = true;
                // Entry is a zero bit counter but it is zero
                } else {
                    zero_count_switch = false;
                }
            }
            not_zero => {
                // Entry is a zero bit counter and not zero
                if zero_count_switch == true {
                    tax_index_count += *not_zero as usize;
                    zero_count_switch = false;
                // Entry is a non-zero byte
                } else {
                    push_tax_from_bit_vec_to_tax_vec(*tax_entry, tax_index_count, &mut tax_output);
                    tax_index_count += BYTE_SIZE;
                }
            }
        }
    }

    tax_output
}

fn push_tax_from_bit_vec_to_tax_vec(
    tax_entry: u8,
    tax_index_count: usize,
    tax_output: &mut Vec<usize>,
) {
    /*
     * Outputs a list of zero-indexed tax.
     */
    let mut cache: u8 = tax_entry;
    let mut bit_index = tax_index_count;
    while cache != 0 {
        // Uneven cache == 1. Uneven cache -> rightmost bit is 1
        if (cache & 1) == 1 {
            tax_output.push(bit_index);
        }
        cache >>= 1;
        bit_index += 1;
    }
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
    k: usize,
    tax_groups: &Vec<String>,
    output: &PathBuf,
    nthreads: &usize,
) -> Vec<PathBuf> {
    log::info!("Start create_db with {} taxes.", tax_groups.len());
    const SPLITS: usize = 16;

    let db_prefix_chunks = split_db_prefix(SPLITS);

    let mut index_files: Vec<PathBuf> = vec![];

    let max_kmer_bit = max_bits(k * 2);
    let kmer_overflow_bits = usize::MAX - max_kmer_bit;

    log::info!("Will split index in {} chunks", db_prefix_chunks.len());
    for (db_index, db_prefixes) in enumerate(db_prefix_chunks) {
        log::info!("Creating index for index split: {:?}", &db_prefixes);
        let mut db: AHashMap<usize, Vec<u16>> = AHashMap::with_capacity(100_000_000);

        let (sender, receiver) = mpsc::sync_channel::<(usize, AHashSet<usize>)>(64);
        let mut thread_handles: Vec<JoinHandle<()>> = vec![];

        let worklists = split_tax_groups_into_chunks(nthreads, metadata_map);

        log::info!("Spawning {} threads", worklists.len());
        for (i, worklist) in enumerate(&worklists) {
            log::info!("\tthread {}: {:?}", i, worklist);
        }

        for worklist in worklists {
            let sender = sender.clone();
            let tax_groups = tax_groups.clone();
            let kmer_filter = db_prefixes.clone();

            thread_handles.push(thread::spawn(move || {
                log::info!(
                    "Start worklist thread with the kmer filter {:?}",
                    kmer_filter
                );
                for metadata in worklist {
                    let tax_index: usize = get_tax_index(&metadata, &tax_groups);
                    let file_path = metadata.path.clone();

                    let mut file = fs::File::open(&file_path).unwrap(); // TODO std::io::Error

                    let decoder = GzDecoder::new(&mut file);
                    let records = fasta::Reader::new(decoder).records();

                    // let mut kmers: AHashSet<usize> = AHashSet::new();

                    for record in records {
                        // kmers.extend(get_unique_kmers(
                        //     record.unwrap().seq(),
                        //     &k,
                        //     &max_kmer_bit,
                        //     &kmer_overflow_bits,
                        // ));
                        sender.send((
                            tax_index,
                            get_unique_kmers(
                                record.unwrap().seq(),
                                &k,
                                &max_kmer_bit,
                                &kmer_overflow_bits,
                                &kmer_filter,
                            ),
                        ));
                    }
                    log::info!(
                        "worklist thread finished: {:?}",
                        &file_path.file_name().unwrap()
                    );
                    // sender.send((tax_index, kmers));
                }
            }));
        }
        drop(sender);

        let mut length_log = 1000_000;

        for (tax_index, kmers) in receiver {
            for kmer in kmers {
                annotate_kmer_with_tax(&mut db, kmer, tax_index)
            }
            if length_log < db.len() {
                length_log = 1000_000 + db.len();
                log::info!("Reciever got {} kmers.", db.len());
            }
        }

        for handle in thread_handles {
            handle.join().unwrap();
        }

        log::info!(
            "Finished creating split with prefixes {:?} with {} kmers.",
            &db_prefixes,
            db.len()
        );

        // let mut output_split_filename = output.file_name().unwrap().to_owned();
        // output_split_filename.push(db_index.to_string());
        // let mut output_split = output.clone();
        // output_split.set_file_name(output_split_filename);

        let output_split = create_file_path(&output, &db_index.to_string(), "probmap");
        serialize_compress_write(&output_split, &db);

        // TODO DEBUG
        // for key in db.keys() {
        //     log::info!("{}", kmer_int_2_kmer_str(*key, k));
        // }

        log::info!(
            "Wrote database with prefixes {:?} to: {:?}",
            &db_prefixes,
            &output_split.as_os_str()
        );

        index_files.push(output_split);
    }

    index_files
}

fn annotate_kmer_with_tax(db: &mut AHashMap<usize, Vec<u16>>, kmer: usize, tax_index: usize) {
    match db.get_mut(&kmer) {
        Some(tax_vec) => {
            // let new_tax_vec = insert_index_into_tax_vec(&tax_vec, tax_index);
            // db.insert(kmer, new_tax_vec);
            tax_vec.push(tax_index as u16);
        }
        None => {
            //db.insert(kmer, new_tax_vec(tax_index));
            db.insert(kmer, vec![tax_index as u16]);
        }
    }
}

fn split_db_prefix(splits: usize) -> Vec<SmallVec<[usize; 16]>> {
    let mut splits_vec: Vec<SmallVec<[usize; 16]>> = vec![];
    let prefix_size = DB_PREFIX.len();
    let mut chunk_size = prefix_size / splits;
    if prefix_size % splits != 0 {
        chunk_size += 1;
    }

    for i in (0..prefix_size).step_by(chunk_size) {
        let mut prefix_chunk: SmallVec<[usize; 16]> = smallvec![0; chunk_size];
        // let mut prefix_chunk = [0usize; chunk_size];
        for j in 0..chunk_size {
            prefix_chunk[j] = DB_PREFIX[i + j];
        }
        splits_vec.push(prefix_chunk);
    }

    splits_vec
}

fn split_kmer_db_into_chunks(splits: usize, db: &AHashMap<usize, Vec<u16>>) -> Vec<Vec<usize>> {
    let mut worklists: Vec<Vec<usize>> = vec![];

    let mut chunk_size = db.len() / splits;
    if db.len() % splits != 0 {
        chunk_size += 1;
    }

    let mut counter = 0;
    for chunk in &db.keys().chunks(chunk_size) {
        worklists.push(chunk.cloned().collect_vec());
    }

    worklists
}

fn split_tax_groups_into_chunks<'a>(
    splits: &usize,
    metadata_map: &AHashMap<String, MetaEntry>,
) -> Vec<Vec<MetaEntry>> {
    let mut splits = *splits;
    let mut worklists: Vec<Vec<MetaEntry>> = vec![];
    let metadata_size = metadata_map.len();
    if metadata_size < splits {
        splits = metadata_size;
    }

    let mut chunk_size = metadata_size / splits;
    if metadata_size % splits != 0 {
        chunk_size += 1;
    }

    let mut worklist: Vec<MetaEntry> = vec![];
    let mut current_chunk = 0;
    for metadata in metadata_map.values() {
        if current_chunk < chunk_size {
            current_chunk += 1;
            worklist.push(metadata.clone());
        } else if current_chunk == chunk_size {
            worklists.push(worklist);
            worklist = vec![];
            worklist.push(metadata.clone());
            current_chunk = 1;
        }
    }
    worklists.push(worklist);

    worklists
}

fn insert_index_tax_bit_in_byte(index_max: usize, tax_index: usize, tax_byte: u8) -> u8 {
    let index_diff = index_max - tax_index;
    let zeroes = 8 - index_diff;
    let new_byte = 1 << zeroes;
    return tax_byte | new_byte;
}

fn copy_values_between_vectors(
    src_index: usize,
    src_vec: &Vec<u8>,
    dest_index: usize,
    dest_vec: &mut Vec<u8>,
) {
    /*
     * Copies values from the given index of a source vector to destination
     * vector with the destination index offset given.
     * If destination vector is too short the last values will be pushed.
     */
    for i in 0..(src_vec.len() - src_index) {
        if (dest_index + i) < dest_vec.len() {
            dest_vec[dest_index + i] = src_vec[src_index + i];
        } else {
            dest_vec.push(src_vec[src_index + i])
        }
    }
}

fn insert_or_push_to_vec(vec: &mut Vec<u8>, index: usize, value: u8) {
    if vec.len() > index {
        vec[index] = value;
    } else {
        vec.push(value);
    }
}

fn insert_index_into_tax_vec(tax_vec: &Vec<u8>, tax_index: usize) -> Vec<u8> {
    const BYTE_SIZE: usize = 8;
    let mut new_tax_vec: Vec<u8> = vec![0u8; tax_vec.len()];
    let mut new_tax_vec_i: usize = 0;
    let mut tax_vec_offset: usize = 0;
    let mut index_counter: usize = 0;
    let mut zero_count_switch = false;

    for i in 0..tax_vec.len() {
        let tax_entry = tax_vec[i + tax_vec_offset];
        new_tax_vec[new_tax_vec_i] = tax_entry;

        // Entry is 0 byte (not a zero counter)
        if tax_entry == 0 && zero_count_switch == false {
            index_counter += BYTE_SIZE;
            zero_count_switch = true;

            if index_counter > tax_index {
                new_tax_vec[new_tax_vec_i] =
                    insert_index_tax_bit_in_byte(index_counter, tax_index, 0);
                new_tax_vec_i += 1;
                // Zero byte followed by non-zero counter.
                if tax_vec[i + tax_vec_offset + 1] != 0 {
                    let new_count = tax_vec[i + tax_vec_offset + 1] - BYTE_SIZE as u8;
                    insert_or_push_to_vec(&mut new_tax_vec, new_tax_vec_i, 0);
                    new_tax_vec_i += 1;
                    insert_or_push_to_vec(&mut new_tax_vec, new_tax_vec_i, new_count);

                    tax_vec_offset += 2;
                    new_tax_vec_i += 1;

                    index_counter += BYTE_SIZE + new_count as usize;

                // Zero byte followed by counter==0. As zero byte is no
                // longer zero, the counter should be removed. Done by
                // skipping the index.
                } else {
                    tax_vec_offset += 2;
                    new_tax_vec.pop();
                }
                copy_values_between_vectors(
                    i + tax_vec_offset,
                    &tax_vec,
                    new_tax_vec_i,
                    &mut new_tax_vec,
                );

                break;
            }

        // Entry is a zero bit counter but it is zero
        } else if tax_entry == 0 && zero_count_switch == true {
            zero_count_switch = false;

        // Entry is a zero bit counter and not zero
        } else if zero_count_switch == true {
            index_counter += tax_entry as usize;
            zero_count_switch = false;

            if index_counter > tax_index {
                let post_zeroes = index_counter - tax_index - BYTE_SIZE + 1;
                let pre_zeroes = (tax_index - (index_counter - tax_entry as usize)) + BYTE_SIZE;

                // Prev index to counter is always a zero byte.
                new_tax_vec_i -= 1;

                // Insert zeroes before index byte and index byte
                insert_zeroes_in_tax_vec(pre_zeroes, &mut new_tax_vec, &mut new_tax_vec_i, true);
                new_tax_vec_i += 1;
                tax_vec_offset += 1;
                // Insert zeroes after index byte
                insert_zeroes_in_tax_vec(post_zeroes, &mut new_tax_vec, &mut new_tax_vec_i, false);

                copy_values_between_vectors(
                    i + tax_vec_offset,
                    &tax_vec,
                    new_tax_vec_i,
                    &mut new_tax_vec,
                );

                break;
            }

        // Entry is a non-zero byte
        } else {
            index_counter += BYTE_SIZE;
            if index_counter > tax_index {
                new_tax_vec[new_tax_vec_i] =
                    insert_index_tax_bit_in_byte(index_counter, tax_index, tax_entry);
                new_tax_vec_i += 1;
                tax_vec_offset += 1;
                copy_values_between_vectors(
                    i + tax_vec_offset,
                    &tax_vec,
                    new_tax_vec_i,
                    &mut new_tax_vec,
                );
                break;
            }
        }

        new_tax_vec_i += 1;
    }

    // Index to insert is larger than the largest existing index
    if tax_index > index_counter {
        let mut rest = tax_index - index_counter;

        // Inserts 0, 255 pairs
        let zero_chunk_count = rest / 263;
        for _ in 0..zero_chunk_count {
            new_tax_vec.push(0);
            new_tax_vec.push(255);
        }

        // Insert zero count < 255
        rest -= zero_chunk_count * 263;
        let mut tax_vec_i = new_tax_vec.len();
        insert_zeroes_in_tax_vec(rest, &mut new_tax_vec, &mut tax_vec_i, true);
    }

    new_tax_vec
}

fn insert_zeroes_in_tax_vec(
    zeroes: usize,
    new_tax_vec: &mut Vec<u8>,
    new_tax_vec_i: &mut usize,
    include_index_byte: bool,
) {
    /*
     * zeroes must be <= 263
     * if include_index_byte is false then zeroes % 8 == 0 must be true.
     */

    const BYTE_SIZE: usize = 8;
    let mut remaining_zeroes = zeroes;

    if remaining_zeroes >= BYTE_SIZE {
        let zero_count = ((remaining_zeroes - BYTE_SIZE) / 8) * 8;
        insert_or_push_to_vec(new_tax_vec, *new_tax_vec_i, 0);
        *new_tax_vec_i += 1;
        insert_or_push_to_vec(new_tax_vec, *new_tax_vec_i, zero_count as u8);
        *new_tax_vec_i += 1;
        remaining_zeroes = remaining_zeroes - zero_count - BYTE_SIZE;
    }

    if include_index_byte {
        // Create and insert/push the final byte
        let entry_byte: u8 = 1 << remaining_zeroes;
        insert_or_push_to_vec(new_tax_vec, *new_tax_vec_i, entry_byte);
    }
}

fn new_tax_vec(tax_index: usize) -> Vec<u8> {
    const BYTE_SIZE: usize = 8;
    // Create a tax vector with no annotated tax up till the provided
    // tax index.
    // This means first entry is zero, given the provided tax is >7.
    // The following entries then describes non-annotated indices in
    // chunks of maximum 255 + 8 = 263.

    // After a zero will ALWAYS follow the number of next zeroes until
    // next annotated index. This number can be zero itself.

    // 0000 0000, 0000 0000, 1000 0000 = 0, 0, 128       => 15
    // 0000 0000, 0000 0000, 1000 0010 = 0, 0, 130       => 9, 15
    // 0000 0010, 1000 0010 = 2, 130 => 1, 9, 15

    // ex1: 1, ex2: 1, ex3: 0, ex4: 0
    let zero_vec_count = (tax_index as f32 / 263_f32).floor() as usize;
    // ex1: 135, ex2: 7, ex3: 4, ex4: 15
    let zeroes_left = tax_index - zero_vec_count * 263;
    // ex1: 128, ex2: 0, ex3: 0, ex4: 8
    let zero_count = (zeroes_left / BYTE_SIZE) * BYTE_SIZE;
    let zero_count_vec_size: usize;
    if zero_count > 0 {
        zero_count_vec_size = 2;
    } else {
        zero_count_vec_size = 0;
    }
    // ex1: 7, ex2: 7, ex3: 4, ex4: 7
    let zero_bits = zeroes_left - zero_count;
    // ex1: 5, ex2: 3, ex3: 1, ex4: 3
    let mut tax_vec = vec![0u8; zero_vec_count * 2 + zero_count_vec_size + 1];

    // ex1: [0, 255, 0, 0, 0]
    // ex2: [0, 255, 0]
    // ex3: [0]
    // ex4: [0, 0, 0]
    let zero_count_index = zero_vec_count * 2;
    for i in (0..zero_count_index).step_by(2) {
        tax_vec[i] = 0;
        tax_vec[i + 1] = 255;
    }

    // ex1: [0, 255, 0, 120, 0]
    // ex2: [0, 255, 0]
    // ex3: [0]
    // ex4: [0, 0, 0]
    let last_byte_index;
    if zero_count > 8 {
        tax_vec[zero_count_index] = 0;
        tax_vec[zero_count_index + 1] = (zero_count - 8) as u8;
        last_byte_index = zero_count_index + 2;
    } else if zero_count == 8 {
        tax_vec[zero_count_index] = 0;
        tax_vec[zero_count_index + 1] = 0;
        last_byte_index = zero_count_index + 2;
    } else {
        last_byte_index = zero_count_index;
    }

    let last_byte: u8 = 1 << zero_bits;
    tax_vec[last_byte_index] = last_byte;

    tax_vec
}

fn serialize_compress_write<T>(output: &PathBuf, object: &T)
where
    T: Serialize,
{
    let out_file = File::create(output).unwrap();
    let mut buffer = BufWriter::with_capacity(IO_BLOCKSIZE, out_file);
    // let mut comp = brotli::CompressorWriter::new(
    //     &mut buffer,
    //     4096,      // buffer size
    //     6 as u32,  // quality 0-11
    //     20 as u32, // lg_window_size
    // );
    bincode::serialize_into(&mut buffer, object);
}

fn get_tax_index(meta_entry: &MetaEntry, tax_groups: &Vec<String>) -> usize {
    let tax_index = tax_groups
        .iter()
        .find_position(|tax| *tax == &meta_entry.gtdb_tax);

    match tax_index {
        Some((index, _)) => index,
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

fn int2nuc(n: usize) -> Result<&'static str, ()> {
    match n {
        0 => Ok("A"),
        1 => Ok("C"),
        2 => Ok("G"),
        3 => Ok("T"),
        _ => unreachable!(),
    }
}

fn kmer_int_2_kmer_str(kmer_int: usize, k: usize) -> String {
    let bit_mask: usize = 3;
    let mut kmer_int_cache = kmer_int.clone();
    let mut kmer_list: Vec<&str> = vec![];
    for _ in 0..k {
        let nuc_int = kmer_int_cache & bit_mask;
        kmer_list.push(int2nuc(nuc_int).unwrap());
        kmer_int_cache >>= 2;
    }
    kmer_list.reverse();
    String::from_iter(kmer_list)
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
    filter: &SmallVec<[usize; 16]>,
) -> AHashSet<usize> {
    let mut kmer_set: AHashSet<usize> = AHashSet::with_capacity(500_000);
    let mut kmer: usize = 0;
    let mut ignore_count = 0;

    for (nuc_string_index, nuc) in nuc_string.iter().enumerate() {
        kmer <<= 2;

        match nuc2int(nuc) {
            Ok(nuc_int) => {
                kmer += nuc_int;
            }
            Err(_) => {
                // Ignore nucs as long as kmer contains a non-nuc
                ignore_count += k
            }
        }

        while kmer > *max_kmer_bit {
            kmer -= max_kmer_bit + 1
        }

        if ignore_count > 0 {
            ignore_count -= 1;
            continue;
        }

        if nuc_string_index >= *k - 1 {
            let kmer_rev_comp = rev_comp(kmer, *kmer_overflow_bits, *k);
            if kmer <= kmer_rev_comp {
                if ignore_kmer(kmer, &filter) {
                    continue;
                };
                kmer_set.insert(kmer);
            } else {
                if ignore_kmer(kmer_rev_comp, &filter) {
                    continue;
                };
                kmer_set.insert(kmer_rev_comp);
            }
        }
    }

    kmer_set
}

fn ignore_kmer(kmer: usize, filter: &SmallVec<[usize; 16]>) -> bool {
    // Check if combination of the first two nucleotides are in
    // filter
    let comparison_bit = kmer & 15;
    for filter_bit in filter {
        if comparison_bit == *filter_bit {
            return false;
        }
    }
    true
}

#[derive(Debug, Clone)]
struct MetaEntry {
    filename: String,
    path: PathBuf,
    ncbi_taxid: u32,
    gtdb_tax: String,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct IndexInfo {
    k: usize,
    tax_file: String,
    index_files: Vec<PathBuf>,
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
    use itertools::enumerate;

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
        let filter: SmallVec<[usize; 16]> = SmallVec::from_iter(DB_PREFIX);

        let kmer_set =
            get_unique_kmers(nuc_string, &k, &max_kmer_bit, &kmer_overflow_bits, &filter);

        assert_eq!(kmer_set.len(), 3);
        assert_eq!(kmer_set.contains(&126), true);
        assert_eq!(kmer_set.contains(&66), true);
        assert_eq!(kmer_set.contains(&528), true);
        assert_eq!(kmer_set.contains(&507), false);
        assert_eq!(kmer_set.contains(&1005), false);
        assert_eq!(kmer_set.contains(&999999), false);

        let filter2: SmallVec<[usize; 16]> = smallvec![14, 11];
        let kmer_set2 =
            get_unique_kmers(nuc_string, &k, &max_kmer_bit, &kmer_overflow_bits, &filter2);
        assert_eq!(kmer_set2.len(), 1);
        assert_eq!(kmer_set2.contains(&126), true);
        // 507 is in filter, but is > 66 which is not in filter, hence both
        // should be ignored.
        assert_eq!(kmer_set2.contains(&507), false);
        assert_eq!(kmer_set2.contains(&66), false);
    }

    #[test]
    fn test_new_tax_vec() {
        /*
         * Example 1: [0, 255, 0, 120, 64]
         * Would annotate tax a index 263 + 8 + 120 + 7 - 1 = 397
         * Above 64 = 0100 0000
         *
         * Example 2: [0, 255, 128] = 263 + 8 - 1 = 270
         * 128 = 1000 0000
         *
         * Example 3: [16] = b00010000 => 4
         *
         * Example 4: [0, 0, 128] = 8 + 0 + 8 - 1 = 15
         *
         * Example 5: [0, 255, 0, 255, 32] = 263 + 263 + 5 - 1 = 530
         * 16 = 0001 0000
         */

        let tax_index_ex1 = 397;
        let tax_index_ex2 = 270;
        let tax_index_ex3 = 4;
        let tax_index_ex4 = 15;
        let tax_index_ex5 = 530;
        let tax_index_ex6 = 13;

        let tax_vec_ex1 = new_tax_vec(tax_index_ex1);
        assert_eq!(tax_vec_ex1, [0, 255, 0, 120, 64]);

        let tax_vec_ex2 = new_tax_vec(tax_index_ex2);
        assert_eq!(tax_vec_ex2, [0, 255, 128]);

        let tax_vec_ex3 = new_tax_vec(tax_index_ex3);
        assert_eq!(tax_vec_ex3, [16]);

        let tax_vec_ex4 = new_tax_vec(tax_index_ex4);
        assert_eq!(tax_vec_ex4, [0, 0, 128]);

        let tax_vec_ex5 = new_tax_vec(tax_index_ex5);
        assert_eq!(tax_vec_ex5, [0, 255, 0, 255, 16]);

        let tax_vec_ex6 = new_tax_vec(tax_index_ex6);
        assert_eq!(tax_vec_ex6, [0, 0, 32]);
    }

    #[test]
    fn test_insert_or_push_to_vec() {
        let mut test_vec = vec![0u8; 3];

        insert_or_push_to_vec(&mut test_vec, 1, 112);
        assert_eq!(test_vec, [0, 112, 0]);

        insert_or_push_to_vec(&mut test_vec, 3, 27);
        assert_eq!(test_vec, [0, 112, 0, 27]);

        // TODO This should be made to fail somehow.
        insert_or_push_to_vec(&mut test_vec, 20, 20);
        assert_eq!(test_vec, [0, 112, 0, 27, 20]);
    }

    #[test]
    fn test_copy_values_between_vectors() {
        let src_vec: Vec<u8> = vec![0, 0, 10, 4, 5, 6, 7, 8];
        let mut dest_vec: Vec<u8> = vec![3, 4, 5, 6];
        let src_index = 3;
        let dest_index = 1;
        copy_values_between_vectors(src_index, &src_vec, dest_index, &mut dest_vec);
        assert_eq!(dest_vec, [3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn test_insert_zeroes_in_tax_vec() {
        let mut tax_vec_1 = vec![1, 2];
        let mut tax_vec_i_1 = 2;
        insert_zeroes_in_tax_vec(170, &mut tax_vec_1, &mut tax_vec_i_1, true);
        assert_eq!(tax_vec_1, [1, 2, 0, 160, 4]);

        let mut tax_vec_2 = vec![1, 2];
        let mut tax_vec_i_2 = 2;
        insert_zeroes_in_tax_vec(6, &mut tax_vec_2, &mut tax_vec_i_2, true);
        assert_eq!(tax_vec_2, [1, 2, 64]);

        let mut tax_vec_3 = vec![0, 255, 0, 120, 0];
        let mut tax_vec_i_3 = 2;
        insert_zeroes_in_tax_vec(25, &mut tax_vec_3, &mut tax_vec_i_3, true);
        assert_eq!(tax_vec_3, [0, 255, 0, 16, 2]);

        let mut dest_vec_1 = vec![2, 3, 0, 0];
        let mut dest_vec_i_1 = 2;
        insert_zeroes_in_tax_vec(40, &mut dest_vec_1, &mut dest_vec_i_1, false);
        assert_eq!(dest_vec_1, [2, 3, 0, 32]);

        let mut dest_vec_3 = vec![2, 3];
        let mut dest_vec_i_3 = 1;
        insert_zeroes_in_tax_vec(0, &mut dest_vec_3, &mut dest_vec_i_3, false);
        assert_eq!(dest_vec_3, [2, 3]);

        let mut dest_vec_4 = vec![0, 255, 0, 16, 2];
        let mut dest_vec_i_4 = 5;
        insert_zeroes_in_tax_vec(88, &mut dest_vec_4, &mut dest_vec_i_4, false);
        assert_eq!(dest_vec_4, [0, 255, 0, 16, 2, 0, 80]);
    }

    #[test]
    fn test_insert_index_into_tax_vec() {
        /*
         * Insertion cases:
         * 1. Insertion into existing byte.
         * 2. Insertion into zero byte followed by non-zero counter.
         * 3. Insertion into counter more than 8 indeces away from zero byte.
         * 4. Insertion into counter less than 8 indeces away from zero byte.
         * 5. Insertion into zero byte followed by counter==0.
         * 6. Insertion after greatest index.
         *
         * 1. Insert index 392 into tax_vec with index 397 annotated.
         * [0, 255, 0, 120, 64] -> 263 + 8 + 120 + 7 - 1 -> 397
         * 64 = 0100 0000 => 0100 0010 = 66
         * [0, 255, 0, 120, 66] -> 392, 397
         *
         * 2. Insert index 265 into tax_vec with index 397 annotated.
         * [0, 255, 0, 120, 64] -> 263 + 8 + 120 + 7 - 1 -> 397
         * [0, 255, 4, 0, 112, 64] -> 265, 397
         * 4 = 0000 0100
         *
         * 3. Insert index 288 into tax_vec with index 397 annotated.
         * [0, 255, 0, 120, 64]
         * [0, 255, 0, 16, 2, 0, 88, 64] -> 288, 397
         *  263 + 128 = 391
         *  391 - 288 = 103 = post_zeroes
         *  288 - (255+8+8)
         *
         * 4. Insert index 272 into tax_vec with index 397 annotated.
         * [0, 255, 0, 120, 64]
         * [0, 255, 0, 0, 2, 0, 104, 64] -> 272, 397
         *
         * 5. Insert index 270 into tax_vec with index 272 annotated.
         * [0, 255, 0, 0, 2] -> 272
         * [0, 255, 128, 2] > 270, 272
         *
         * 6. Insert index 420 into tax_vec with index 397 annotated.
         * [0, 255, 0, 120, 64] -> 397
         * [0, 255, 0, 120, 64, 0, 8, 32] -> 397, 420
         *
         * 7. Insert index 720 into tax_vec with index 397 annotated.
         * [0, 255, 0, 120, 64] -> 397
         * [0, 255, 0, 120, 64, 0, 255, 0, 48, 4] -> 397, 720
         *
         * 8. Insert index 5 into tax_vec with index 5, 11, 22, 30 annotated.
         * [32, 8, 128, 64] -> 5, 11, 22, 30
         *
         */

        let tax_vec_1: Vec<u8> = vec![0, 255, 0, 120, 64];
        let tax_vec_new_1 = insert_index_into_tax_vec(&tax_vec_1, 392);
        assert_eq!(tax_vec_new_1, [0, 255, 0, 120, 66]);

        let tax_vec_2: Vec<u8> = vec![0, 255, 0, 120, 64];
        let tax_vec_new_2 = insert_index_into_tax_vec(&tax_vec_2, 265);
        assert_eq!(tax_vec_new_2, [0, 255, 4, 0, 112, 64]);

        let tax_vec_3: Vec<u8> = vec![0, 255, 0, 120, 64];
        let tax_vec_new_3 = insert_index_into_tax_vec(&tax_vec_3, 288);
        assert_eq!(tax_vec_new_3, [0, 255, 0, 16, 2, 0, 88, 64]);

        let tax_vec_4: Vec<u8> = vec![0, 255, 0, 120, 64];
        let tax_vec_new_4 = insert_index_into_tax_vec(&tax_vec_4, 272);
        assert_eq!(tax_vec_new_4, [0, 255, 0, 0, 2, 0, 104, 64]);

        let tax_vec_5: Vec<u8> = vec![0, 255, 0, 0, 2];
        let tax_vec_new_5 = insert_index_into_tax_vec(&tax_vec_5, 270);
        assert_eq!(tax_vec_new_5, [0, 255, 128, 2]);

        let tax_vec_6: Vec<u8> = vec![0, 255, 0, 120, 64];
        let tax_vec_new_6 = insert_index_into_tax_vec(&tax_vec_6, 420);
        assert_eq!(tax_vec_new_6, [0, 255, 0, 120, 64, 0, 8, 32]);

        let tax_vec_7: Vec<u8> = vec![0, 255, 0, 120, 64];
        let tax_vec_new_7 = insert_index_into_tax_vec(&tax_vec_7, 720);
        assert_eq!(tax_vec_new_7, [0, 255, 0, 120, 64, 0, 255, 0, 48, 4]);

        let tax_vec_8: Vec<u8> = vec![32, 8, 128, 64];
        let tax_vec_new_8 = insert_index_into_tax_vec(&tax_vec_8, 5);
        assert_eq!(tax_vec_new_8, [32, 8, 128, 64]);
    }

    #[test]
    fn test_push_tax_from_bit_vec_to_tax_vec() {
        let tax_index_count = 100;
        let mut tax_output_1: Vec<usize> = vec![];
        let mut tax_output_2: Vec<usize> = vec![];
        let mut tax_output_3: Vec<usize> = vec![];
        let entry_1 = 0b_0000_0000; // 0
        let entry_2 = 0b_0000_0001; // 1
        let entry_3 = 0b_0100_1001; // 1, 4, 7

        push_tax_from_bit_vec_to_tax_vec(entry_1, tax_index_count, &mut tax_output_1);
        assert_eq!(tax_output_1, []);

        push_tax_from_bit_vec_to_tax_vec(entry_2, tax_index_count, &mut tax_output_2);
        assert_eq!(tax_output_2, [100]);

        push_tax_from_bit_vec_to_tax_vec(entry_3, tax_index_count, &mut tax_output_3);
        assert_eq!(tax_output_3, [100, 103, 106]);
    }

    #[test]
    fn test_extract_tax_from_list() {
        let tax_vec_1: Vec<u8> = vec![0, 255, 0, 120, 66];
        let tax_list_1 = extract_tax_from_list(&tax_vec_1);
        assert_eq!(tax_list_1, [392, 397]);

        let tax_vec_2: Vec<u8> = vec![0, 255, 4, 0, 112, 64];
        let tax_list_2 = extract_tax_from_list(&tax_vec_2);
        assert_eq!(tax_list_2, [265, 397]);

        let tax_vec_3: Vec<u8> = vec![0, 255, 0, 16, 2, 0, 88, 64];
        let tax_list_3 = extract_tax_from_list(&tax_vec_3);
        assert_eq!(tax_list_3, [288, 397]);

        let tax_vec_4: Vec<u8> = vec![0, 255, 0, 0, 2, 0, 104, 64];
        let tax_list_4 = extract_tax_from_list(&tax_vec_4);
        assert_eq!(tax_list_4, [272, 397]);

        let tax_vec_5: Vec<u8> = vec![0, 255, 128, 2];
        let tax_list_5 = extract_tax_from_list(&tax_vec_5);
        assert_eq!(tax_list_5, [270, 272]);

        let tax_vec_6: Vec<u8> = vec![0, 255, 0, 120, 64, 0, 8, 32];
        let tax_list_6 = extract_tax_from_list(&tax_vec_6);
        assert_eq!(tax_list_6, [397, 420]);

        let tax_vec_7: Vec<u8> = vec![0, 255, 0, 120, 64, 0, 255, 0, 48, 4];
        let tax_list_7 = extract_tax_from_list(&tax_vec_7);
        assert_eq!(tax_list_7, [397, 720]);
    }

    #[test]
    fn test_calc_probs() {
        let prior: f64 = 0.5;
        let sample_counts: Vec<usize> = vec![0, 3, 0, 0, 4, 0];
        let probs: Vec<f64> = calc_probs(&sample_counts, prior);
        assert_eq!(probs, [0.05, 0.35, 0.05, 0.05, 0.45, 0.05]);
    }

    // #[test]
    // fn test_annotate_kmers() {
    //     let database: AHashMap<usize, Vec<u8>> = AHashMap::from([
    //         (1000usize, vec![0, 255, 0, 120, 66]),          // 392, 397
    //         (1001usize, vec![0, 255, 4, 0, 112, 64]),       // 265, 397
    //         (1002usize, vec![0, 255, 0, 120, 66]),          // 392, 397
    //         (1003usize, vec![0, 255, 0, 16, 2, 0, 88, 64]), // 288, 397
    //     ]);
    //     let kmers: Vec<usize> = vec![1000, 1003, 1002, 1000, 666];
    //     let tax_groups_length = 399;
    //     let prior: f64 = 0.1;

    //     let annotated_kmers = annotate_kmers(&kmers, &database, prior, tax_groups_length);

    //     let mut sample_counts: Vec<usize> = vec![0; 399 + 1];
    //     // + 1 due to unkown category at index 0
    //     sample_counts[392 + 1] = 3;
    //     sample_counts[397 + 1] = 4;
    //     sample_counts[288 + 1] = 1;
    //     sample_counts[0] = 1;
    //     let expected_result: Vec<f64> = calc_probs(&sample_counts, prior);

    //     for (i, entry) in enumerate(annotated_kmers) {
    //         let entry_mult = (entry * 1000f64).round();
    //         let expect_mult = (expected_result[i] * 1000f64).round();
    //         assert_eq!(entry_mult, expect_mult);
    //     }
    // }

    #[test]
    fn test_split_db_prefix() {
        let prefixes = split_db_prefix(10);
        assert_eq!(prefixes[0][0], 0);
        assert_eq!(prefixes[0][1], 15);
        assert_eq!(prefixes[1][0], 4);
        assert_eq!(prefixes[1][1], 11);

        let prefixes = split_db_prefix(1);
        assert_eq!(prefixes[0][0], 0);
        assert_eq!(prefixes[0][3], 11);
        assert_eq!(prefixes[0][15], 12);

        let prefixes = split_db_prefix(16);
        assert_eq!(prefixes[0][0], 0);
        assert_eq!(prefixes[3][0], 11);
        assert_eq!(prefixes[15][0], 12);
    }

    #[test]
    fn test_split_kmer_db_into_chunks() {
        let kmer_db: AHashMap<usize, Vec<u16>> = AHashMap::from([
            (0, vec![1, 2, 3]),
            (1, vec![4, 5, 6]),
            (10, vec![1]),
            (20, vec![2, 3]),
            (350, vec![3]),
        ]);
        let worklists_3_pairs = split_kmer_db_into_chunks(3, &kmer_db);
        let worklists_2_triples = split_kmer_db_into_chunks(2, &kmer_db);
        let worklists_5_singletons = split_kmer_db_into_chunks(5, &kmer_db);

        assert_eq!(worklists_3_pairs.len(), 3);
        assert_eq!(worklists_3_pairs[0].len(), 2);
        assert_eq!(worklists_3_pairs[2].len(), 1);

        assert_eq!(worklists_2_triples.len(), 2);
        assert_eq!(worklists_2_triples[0].len(), 3);
        assert_eq!(worklists_2_triples[1].len(), 2);

        assert_eq!(worklists_5_singletons.len(), 5);
        assert_eq!(worklists_5_singletons[0].len(), 1);
    }

    #[test]
    fn test_kmer_int_2_kmer_str() {
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
        assert_eq!(kmer_int_2_kmer_str(126, k), "ACTTG");
        assert_eq!(kmer_int_2_kmer_str(507, k), "CTTGT");
        assert_eq!(kmer_int_2_kmer_str(1005, k), "TTGTC");
        assert_eq!(kmer_int_2_kmer_str(267, k), "CAAGT");
        assert_eq!(kmer_int_2_kmer_str(66, k), "ACAAG");
        assert_eq!(kmer_int_2_kmer_str(528, k), "GACAA");
    }
}
