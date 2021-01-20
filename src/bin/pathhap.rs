use clap::{App, Arg};
use std::io::{BufRead, BufReader};

fn main() -> std::io::Result<()> {
    let matches = App::new("pathhap")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Path phasing CLI")
        .setting(clap::AppSettings::ArgRequiredElseHelp)
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("input")
                .long("input")
                .short("i")
                .value_name("PATHS")
                .takes_value(true)
                .required(true)
                .help("Input FASTA file."),
        )
        .arg(
            Arg::with_name("mode")
                .long("mode")
                .takes_value(true)
                .default_value(&"Phase")
                .possible_values(&["Phase", "Hap"])
                .help("Phasing(path-based) or Haplotyping(node-based)"),
        )
        .arg(
            Arg::with_name("param")
                .long("param")
                .short("p")
                .takes_value(true)
                .default_value(&"15")
                .help(
                    "Parameter. Max coverage for Phase module, max length of path for Hap module.",
                ),
        )
        .get_matches();
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    let paths = match matches.value_of("input") {
        Some(res) => get_input_file(res)?,
        None => panic!(),
    };
    let param: usize = match matches.value_of("param").unwrap().parse::<usize>() {
        Ok(res) => res,
        Err(w) => panic!(
            "{} is not valid. {:?}",
            matches.value_of("param").unwrap(),
            w
        ),
    };
    let reuslt = path_phasing::phase(&paths, param);
    let mut result: Vec<_> = reuslt.into_iter().collect();
    result.sort_by_key(|x| x.1);
    println!("ID\tHap");
    for (id, hap) in result {
        println!("{}\t{}", id, hap);
    }
    Ok(())
}

fn get_input_file(path: &str) -> std::io::Result<Vec<(String, Vec<(u64, u64)>)>> {
    std::fs::File::open(path).map(BufReader::new).map(|rdr| {
        rdr.lines()
            .filter_map(|x| x.ok())
            .filter_map(|line| {
                let mut line = line.split('\t');
                let id = line.next()?.to_string();
                let path: Vec<_> = line
                    .map(|unit| {
                        let mut nc = unit.split(':');
                        let node: u64 = match nc.next() {
                            Some(res) => match res.parse::<u64>() {
                                Ok(res) => res,
                                Err(w) => panic!("Error. Could not parse:{}, {:?}", unit, w),
                            },
                            None => panic!("Error. {} is too short.", unit),
                        };
                        let cluster: u64 = match nc.next() {
                            Some(res) => match res.parse::<u64>() {
                                Ok(res) => res,
                                Err(w) => panic!("Error. Could not parse:{},{:?}", unit, w),
                            },
                            None => panic!("Error. {} is too short.", unit),
                        };
                        if let Some(_) = nc.next() {
                            panic!("{} too long.", unit);
                        }
                        (node, cluster)
                    })
                    .collect();
                Some((id, path))
            })
            .collect()
    })
}
