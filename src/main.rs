use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use vcf::*;
use clap::{App, Arg};

fn main() {
    let matches = App::new("vcflevel")
        .version("0.1")
        .author("Erik Garrison <erik.garrison@gmail.com>")
        .about("Filter vg deconstruct output using snarl tree annotations.")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("FILE")
                .help("Filter this input VCF file.")
                .takes_value(true),
        )
        .get_matches();

    let infile = matches.value_of("input").unwrap();

    let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(
        File::open(infile).unwrap(),
    )))
    .unwrap();

    // prepare VCFRecord object
    let mut vcf_record = VCFRecord::new(reader.header().clone());

    // read one record
    reader.next_record(&mut vcf_record).unwrap();

    let handle = std::io::stdout();
    let buf = BufWriter::new(handle);
    let mut writer = VCFWriter::new(buf, reader.header()).unwrap();

    while reader.next_record(&mut vcf_record).unwrap() {
        writer.write_record(&vcf_record).unwrap();
    }
}
