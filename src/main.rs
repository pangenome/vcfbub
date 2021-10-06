use clap::{App, Arg};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufReader, BufWriter};
//use std::String;
//use fnv::{FnvHashMap, FnvHashSet};
use std::collections::{HashMap, HashSet};
use vcf::*;

fn get_header(infile: &str) -> vcf::VCFHeader {
    VCFReader::new(BufReader::new(MultiGzDecoder::new(
        File::open(infile).unwrap(),
    )))
    .unwrap()
    .header()
    .clone()
}

fn for_each_line_in_vcf(infile: &str, mut callback: impl FnMut(&mut VCFRecord, usize)) {
    let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(
        File::open(infile).unwrap(),
    )))
    .unwrap();

    // verify that we have a file we can operate on
    if !reader.header().info_list().any(|x| x == b"LV") {
        panic!("[vcfbub] error: Input VCF must have LV (snarl level) annotations.")
    }
    if !reader.header().info_list().any(|x| x == b"PS") {
        panic!("[vcfbub] error: Input VCF must have PS (parent snarl) annotations.")
    }

    // prepare VCFRecord object
    let mut vcf_record = VCFRecord::new(reader.header().clone());

    let mut idx = 0;

    while reader.next_record(&mut vcf_record).unwrap() {
        callback(&mut vcf_record, idx);
        idx += 1;
    }
}

fn get_level(vcf_record: &VCFRecord) -> i32 {
    let x = vcf_record.info(b"LV").unwrap().clone();
    let y = String::from_utf8(x[0].clone()).unwrap();
    y.parse::<i32>().unwrap()
}

fn get_snarl_id(vcf_record: &VCFRecord) -> String {
    String::from_utf8(vcf_record.id[0].clone()).unwrap()
}

fn get_parent_snarl(vcf_record: &VCFRecord) -> String {
    match vcf_record.info(b"PS") {
        Some(y) => String::from_utf8(y[0].clone()).unwrap(),
        None => "".to_string(),
    }
}

fn get_bubble_length(vcf_record: &VCFRecord) -> usize {
    let mut alt_lengths = Vec::new();
    alt_lengths.push(vcf_record.reference.len());
    for alt in vcf_record.alternative.iter() {
        alt_lengths.push(alt.len());
    }
    let biggest = alt_lengths.iter().max().unwrap();
    let smallest = alt_lengths.iter().min().unwrap();
    biggest - smallest
}

fn get_max_allele_length(vcf_record: &VCFRecord) -> usize {
    let mut alt_lengths = Vec::new();
    alt_lengths.push(vcf_record.reference.len());
    for alt in vcf_record.alternative.iter() {
        alt_lengths.push(alt.len());
    }
    alt_lengths.iter().max().unwrap().clone()
}

fn main() {
    let matches = App::new("vcfbub")
        .version("0.1")
        .author("Erik Garrison <erik.garrison@gmail.com>")
        .about(concat!(
            "\nFilter vg deconstruct output using variant nesting information.\n",
            "This uses the snarl tree decomposition to describe the nesting of variant bubbles.\n",
            "Nesting information must be given in LV (level) and PS (parent snarl) tags."
        ))
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("FILE")
                .help("Filter this input VCF file.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("max-level")
                .short("l")
                .long("max-level")
                .value_name("LEVEL")
                .help("Filter sites with LV > LEVEL.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("max-length")
                .short("b")
                .long("max-length")
                .value_name("LENGTH")
                .help("Filter sites whose longest allele is greater than LENGTH.")
                .takes_value(true),
        )
        .get_matches();

    let infile = matches.value_of("input").unwrap();

    let max_level = match matches.value_of("max-level") {
        Some(l) => l.parse::<i32>().unwrap(),
        None => i32::MAX,
    };

    let max_length = match matches.value_of("max-length") {
        Some(l) => l.parse::<usize>().unwrap(),
        None => usize::MAX,
    };

    let vcf_header = get_header(infile);

    // make a pass to check the levels
    // and decide what to filter
    // recording things by index in the file
    let mut parent_bubble: HashMap<String, String> = HashMap::new();
    let mut popped_bubbles: HashSet<String> = HashSet::new();
    for_each_line_in_vcf(infile, |vcf_record, idx| {
        let level = get_level(vcf_record);
        let max_allele_length = get_max_allele_length(vcf_record);
        if level > max_level || max_allele_length > max_length {
            popped_bubbles.insert(get_snarl_id(vcf_record));
        }
        let snarl_id = get_snarl_id(vcf_record);
        let parent_snarl = get_parent_snarl(vcf_record);
        parent_bubble.insert(snarl_id, parent_snarl);
    });

    // setup writer
    let handle = std::io::stdout();
    let buf = BufWriter::new(handle);
    let mut writer = VCFWriter::new(buf, &vcf_header).unwrap();
    for_each_line_in_vcf(infile, |vcf_record, idx| {
        let snarl_id = get_snarl_id(vcf_record);
        if !popped_bubbles.contains(&snarl_id) {
            writer.write_record(&vcf_record).unwrap()
        }
    });

    // read one record
    //reader.next_record(&mut vcf_record).unwrap();
}
