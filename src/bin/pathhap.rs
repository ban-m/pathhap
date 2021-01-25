use clap::{App, Arg, SubCommand};
use std::io::{BufRead, BufReader};

fn subcommand_phasing() -> App<'static, 'static> {
    SubCommand::with_name("phasing")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Phasing by exact phasing algorithm on graphs.")
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
                .value_name("INPUT")
                .takes_value(true)
                .required(true)
                .help("Input file(TSV). See README.md for the format."),
        )
        .arg(
            Arg::with_name("threads")
                .long("threads")
                .takes_value(true)
                .default_value("1")
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("max_coverage")
                .long("max_cov")
                .short("m")
                .takes_value(true)
                .default_value("20")
                .help("Maximum number of paths on a boundary nodes"),
        )
}

fn subcommand_haplotyping() -> App<'static, 'static> {
    SubCommand::with_name("haplotyping")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Phasing by exact haplotyping algorithm on graphs.")
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
                .value_name("INPUT")
                .takes_value(true)
                .required(true)
                .help("Input file(TSV). See README.md for the format."),
        )
        .arg(
            Arg::with_name("threads")
                .long("threads")
                .takes_value(true)
                .default_value("1")
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("max_prefetch")
                .long("max_pref")
                .short("m")
                .takes_value(true)
                .default_value("20")
                .help("Maximum number of prefetch nodes on a boundary nodes. Roughtly eqaul to the maximum path length."),
        )
        .arg(
            Arg::with_name("error_rate")
                .long("error_rate")
                .short("e")
                .takes_value(true)
                .default_value("0.05")
                .help("Error rate of the path")
        )
}

fn main() -> std::io::Result<()> {
    let matches = App::new("pathhap")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Path phasing CLI")
        .setting(clap::AppSettings::ArgRequiredElseHelp)
        .subcommand(subcommand_haplotyping())
        .subcommand(subcommand_phasing())
        .get_matches();
    let paths = {
        let sub_m = matches.subcommand().1.unwrap();
        let level = match sub_m.occurrences_of("verbose") {
            0 => "warn",
            1 => "info",
            2 => "debug",
            _ => "trace",
        };
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
        if let Some(threads) = sub_m.value_of("threads") {
            match threads.parse::<usize>() {
                Ok(res) => rayon::ThreadPoolBuilder::new()
                    .num_threads(res)
                    .build_global()
                    .unwrap(),
                Err(why) => panic!("{:?}--mulformed threads number", why),
            }
        }
        match sub_m.value_of("input") {
            Some(res) => get_input_file(res)?,
            None => panic!(),
        }
    };
    let reuslt = match matches.subcommand() {
        ("phasing", Some(sub_m)) => {
            let max_cov = sub_m.value_of("max_coverage").unwrap();
            let max_cov: usize = match max_cov.parse::<usize>() {
                Ok(res) => res,
                Err(why) => panic!(
                    "Couldn't parse {}. Please input a valid value(Integer):{:?}",
                    max_cov, why
                ),
            };
            path_phasing::phase(&paths, max_cov)
        }
        ("haplotyping", Some(sub_m)) => {
            let max_len = sub_m.value_of("max_prefetch").unwrap();
            let max_len: usize = match max_len.parse() {
                Ok(res) => res,
                Err(why) => panic!(
                    "Couldn't parse {}. Please input a valid value(Integer):{:?}",
                    max_len, why
                ),
            };
            let error_rate = sub_m.value_of("error_rate").unwrap();
            let error_rate: f64 = match error_rate.parse() {
                Ok(res) => res,
                Err(why) => panic!(
                    "Couldn't parse {}. Please input a valid value(float):{:?}",
                    error_rate, why
                ),
            };
            path_phasing::haplotyping(&paths, max_len, error_rate)
        }
        _ => unreachable!(),
    };
    let mut result: Vec<_> = reuslt.into_iter().collect();
    result.sort_by_key(|x| x.1);
    println!("ID\tHap");
    for (id, hap) in result {
        println!("{}\t{}", id, hap);
    }
    Ok(())
}

type PathWithID = (String, Vec<(u64, u64)>);
fn get_input_file(path: &str) -> std::io::Result<Vec<PathWithID>> {
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
                        if nc.next().is_some() {
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
