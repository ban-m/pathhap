use clap::{App, Arg, SubCommand};
use log::{debug, error};
use std::collections::HashMap;
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
                .value_name("THREADS")
                .takes_value(true)
                .default_value("1")
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("max_coverage")
                .long("max_cov")
                .short("m")
                .value_name("COVERAGE")
                .takes_value(true)
                .default_value("20")
                .help("Maximum number of paths on each boundary"),
        )
        .arg(
            Arg::with_name("min_coverage")
                .long("min_cov")
                .value_name("COVERAGE")
                .takes_value(true)
                .default_value("8")
                .help("Minimum number of paths on each node."),
        )
        .arg(
            Arg::with_name("retry")
                .long("retry")
                .value_name("RETRY")
                .takes_value(true)
                .default_value("7")
                .help("Execute algorithm in [RETRY] times, then pick the best one."),
        )
        .arg(
            Arg::with_name("subsample")
                .long("subsample")
                .value_name("SUBSAMPLE")
                .takes_value(true)
                .default_value("100000")
                .help("Takes top [SUBSAMPLE] bipartition on from each boundary."),
        )
        .arg(
            Arg::with_name("seed")
                .long("seed")
                .value_name("SEED")
                .takes_value(true)
                .default_value("24")
                .help("Seed for a pseudo random number generator."),
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
                .help("Input file(CSV). See README.md for the format."),
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
            let max_cov: usize = sub_m.value_of("max_coverage").map(parse_to_usize).unwrap();
            let min_cov = sub_m.value_of("min_coverage").map(parse_to_usize).unwrap();
            let retry = sub_m.value_of("retry").map(parse_to_usize).unwrap() as u64;
            let subsample = sub_m.value_of("subsample").map(parse_to_usize);
            let seed = sub_m.value_of("seed").map(parse_to_usize).unwrap();
            let seed = seed as u64;
            use rand::{seq::SliceRandom, Rng, SeedableRng};
            let mut rng: rand_xoshiro::Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
            use rayon::prelude::*;
            let seed = rng.gen::<u64>();
            (0..retry)
                .into_par_iter()
                .map(|i| {
                    let mut paths = paths.clone();
                    let mut rng: rand_xoshiro::Xoroshiro128PlusPlus =
                        SeedableRng::seed_from_u64(seed + i);
                    paths.shuffle(&mut rng);
                    let (lk, asn) =
                        path_phasing::phase_with_lk(&paths, max_cov, min_cov, subsample);
                    let asn: HashMap<String, _> =
                        asn.into_iter().map(|(k, v)| (k.to_string(), v)).collect();
                    debug!("Likelihood\t{}\t{:?}", i, lk);
                    (lk, asn)
                })
                .max_by(|a, b| (a.0).partial_cmp(&(b.0)).unwrap())
                .unwrap()
                .1
        }
        ("haplotyping", Some(sub_m)) => {
            error!("This functionality is under development.");
            error!("Please use `phasing` module instead.");
            let max_len = sub_m.value_of("max_prefetch").map(parse_to_usize).unwrap();
            let error_rate = sub_m.value_of("error_rate").unwrap();
            let error_rate: f64 = match error_rate.parse() {
                Ok(res) => res,
                Err(why) => panic!(
                    "Couldn't parse {}. Please input a valid value(float):{:?}",
                    error_rate, why
                ),
            };
            path_phasing::haplotyping(&paths, max_len, error_rate)
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect()
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
                let mut line = line.split(',');
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

fn parse_to_usize(input: &str) -> usize {
    match input.parse::<usize>() {
        Ok(res) => res,
        Err(why) => panic!(
            "Couldn't parse {}. Please input a valid value(Integer):{:?}",
            input, why
        ),
    }
}
