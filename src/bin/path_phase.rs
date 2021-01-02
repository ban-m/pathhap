use clap::{App, Arg};
use std::io::{BufRead, BufReader};
#[macro_use]
extern crate log;
const SEED: u64 = 24;
fn main() -> std::io::Result<()> {
    let matches = App::new("path_phasing")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Path phasing toolkit")
        .setting(clap::AppSettings::ArgRequiredElseHelp)
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("path")
                .long("path")
                .short("p")
                .value_name("PATHS")
                .takes_value(true)
                .required(true)
                .help("Input path files<TSV>."),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .takes_value(true)
                .default_value(&"1")
                .help("number of threads<INT>"),
        )
        .arg(
            Arg::with_name("max_occ")
                .long("max_occ")
                .takes_value(true)
                .default_value(&"15")
                .help("maximum occurence of nodes<INT>"),
        )
        .get_matches();
    let paths = get_input(matches.value_of("path").unwrap())?;
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let max_occ: usize = matches
        .value_of("max_occ")
        .and_then(|num| num.parse().ok())
        .unwrap();
    if let Some(sub_m) = matches.subcommand().1 {
        let level = match sub_m.occurrences_of("verbose") {
            0 => "warn",
            1 => "info",
            2 => "debug",
            _ => "trace",
        };
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    debug!("Begin path phasing");
    let result = path_phasing::phase(&paths, max_occ, SEED);
    for (id, hap) in result {
        println!("{}\t{}", id, hap);
    }
    Ok(())
}

fn get_input<P: AsRef<std::path::Path>>(
    path: P,
) -> std::io::Result<Vec<(String, Vec<(u64, u64)>)>> {
    std::fs::File::open(path).map(BufReader::new).map(|rdr| {
        rdr.lines()
            .filter_map(|line| line.ok())
            .map(|line| {
                let mut line = line.split('\t');
                let id = line.next().unwrap().to_string();
                let path: Vec<_> = line
                    .map(|elm| {
                        let mut elm = elm.split(':');
                        let node: u64 = elm.next().unwrap().parse().unwrap();
                        let cluster: u64 = elm.next().unwrap().parse().unwrap();
                        (node, cluster)
                    })
                    .collect();
                (id, path)
            })
            .collect()
    })
}
