use clap::{App, Arg};
//use std::String;
//use fnv::{FnvHashMap, FnvHashSet};
use std::collections::HashSet;
use rust_htslib::bcf::{Reader, Writer, Read, Record, Header, Format};


fn for_each_line_in_vcf(infile: &str, mut callback: impl FnMut(&mut Record, usize)) {
    let mut reader = Reader::from_path(infile).expect("Error opening VCF");

    // verify that we have a file we can operate on
    //if !reader.header().info_list().any(|x| x == b"LV") {
    //    panic!("[vcfbub] error: Input VCF must have LV (snarl level) annotations.")
    //}
    //if !reader.header().info_list().any(|x| x == b"PS") {
    //    panic!("[vcfbub] error: Input VCF must have PS (parent snarl) annotations.")
    // }

    for (idx, vcf_record) in reader.records().enumerate() {
        callback(&mut vcf_record.unwrap(), idx);
    }
}

fn get_level(vcf_record: &Record) -> i32 {
    let lv_info = vcf_record.info(b"LV");
    let lv_opt = lv_info.integer().expect("Could not parse LV");
    let lv_int = match lv_opt {
        Some(lv_opt) => {
            let lv_array = *lv_opt;
            lv_array[0]
        },
        None => 0
    };
    lv_int as i32
}

fn get_snarl_id(vcf_record: &Record) -> String {
    String::from_utf8(vcf_record.id().clone()).unwrap()
}

fn get_parent_snarl(vcf_record: &Record) -> String {
    match vcf_record.info(b"PS").string() {
        Ok(ps_opt) => {
            match ps_opt {
                Some(ps_array) => Some(String::from_utf8_lossy(&*ps_array[0]).to_string()).unwrap(),
                None => "".to_string(),
            }
        },
        Err(_) => "".to_string(),
    }    
}

/*
fn get_bubble_length(vcf_record: &Record) -> usize {
    let mut alt_lengths = Vec::new();
    alt_lengths.push(vcf_record.reference.len());
    for alt in vcf_record.alternative.iter() {
        alt_lengths.push(alt.len());
    }
    let biggest = alt_lengths.iter().max().unwrap();
    let smallest = alt_lengths.iter().min().unwrap();
    biggest - smallest
}
*/

fn get_max_allele_length(vcf_record: &Record) -> usize {
    let mut allele_lengths = Vec::new();
    for allele in vcf_record.alleles().iter() {
        allele_lengths.push(allele.len());
    }
    allele_lengths.iter().max().unwrap().clone()
}

fn get_ref_allele_length(vcf_record: &Record) -> usize {
    vcf_record.alleles()[0].len()
}

fn main() {
    let matches = App::new("vcfbub")
        .version("0.1.0")
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
                .required(true)
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
            Arg::with_name("max-allele-length")
                .short("a")
                .long("max-allele-length")
                .value_name("LENGTH")
                .help("Filter sites whose max allele length is greater than LENGTH.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("min-allele-length")
                .short("A")
                .long("min-allele-length")
                .value_name("LENGTH")
                .help("Filter sites whose max allele length is less than LENGTH.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("max-ref-length")
                .short("r")
                .long("max-ref-length")
                .value_name("LENGTH")
                .help("Filter sites whose reference allele is longer than LENGTH.")
                .takes_value(true),
                
        )
        .arg(
            Arg::with_name("min-ref-length")
                .short("R")
                .long("min-ref-length")
                .value_name("LENGTH")
                .help("Filter sites whose reference allele is less than LENGTH.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("debug")
                .short("d")
                .long("debug")
                .help("Print debugging information about which sites are removed.")
                .takes_value(false),
        )
        .get_matches();

    let infile = matches.value_of("input").unwrap();

    let max_level = match matches.value_of("max-level") {
        Some(l) => l.parse::<i32>().unwrap(),
        None => i32::MAX,
    };

    let max_allele_length = match matches.value_of("max-allele-length") {
        Some(l) => l.parse::<usize>().unwrap(),
        None => usize::MAX,
    };
    let min_allele_length = match matches.value_of("min-allele-length") {
        Some(l) => l.parse::<usize>().unwrap(),
        None => usize::MAX,
    };

    let max_ref_length = match matches.value_of("max-ref-length") {
        Some(l) => l.parse::<usize>().unwrap(),
        None => usize::MAX,
    };
    let min_ref_length = match matches.value_of("min-ref-length") {
        Some(l) => l.parse::<usize>().unwrap(),
        None => usize::MAX,
    };

    let print_debug = matches.is_present("debug");

    // make a pass to check the levels
    // and decide what to filter
    // recording things by index in the file
    let mut keep_record: Vec<bool> = Vec::new();
    let mut popped_bubbles: HashSet<String> = HashSet::new();
    for_each_line_in_vcf(infile, |vcf_record, _| {
        let level = get_level(vcf_record);
        let max_length = get_max_allele_length(vcf_record);
        let ref_length = get_ref_allele_length(vcf_record);
        let mut keep = true;
        if level > max_level {
            keep = false;
        }
        // if max_length > max_allele_length || ref_length > max_ref_length {
        if max_length > max_allele_length || max_length < min_allele_length || ref_length > max_ref_length || ref_length < min_ref_length {
            if print_debug {
                eprintln!(
                    "popped {} {}",
                    get_snarl_id(vcf_record),
                    vcf_record.pos()
                );
            }
            popped_bubbles.insert(get_snarl_id(vcf_record));
            keep = false;
        }
        keep_record.push(keep);
    });

    // setup writer
    let in_vcf = Reader::from_path(infile).expect("Error opening VCF");
    let in_header = in_vcf.header();
    let in_samples = in_vcf.header().samples();    
    let mut out_header = Header::from_template_subset(in_header, &[]).unwrap();
    for sample in in_samples {
        out_header.push_sample(sample);
    }
    let mut writer = Writer::from_stdout(&out_header, true, Format::Vcf).unwrap();
    for_each_line_in_vcf(infile, |vcf_record, idx| {
        let snarl_id = get_snarl_id(vcf_record);
        let parent_snarl = get_parent_snarl(vcf_record);
        if keep_record[idx] // we keep records when they're marked to keep
            // or if they're marked don't-keep but their parent was popped and they weren't popped
            || (!keep_record[idx]
                && popped_bubbles.contains(&parent_snarl)
                && !popped_bubbles.contains(&snarl_id))
        {
	    writer.write(&vcf_record).unwrap()
        }
    });

    // read one record
    //reader.next_record(&mut vcf_record).unwrap();
}
