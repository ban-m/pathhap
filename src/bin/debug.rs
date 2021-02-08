use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
const MAX_LENGTH: usize = 7;
const MIN_LENGTH: usize = 3;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    // template
    let template1: Vec<u64> = vec![
        (0..=90).collect::<Vec<_>>(),
        (80..=90).collect::<Vec<_>>(),
        (80..=90).collect::<Vec<_>>(),
        (91..100).collect::<Vec<_>>(),
        vec![60],
        (100..120).collect::<Vec<_>>(),
    ]
    .concat();
    let template2: Vec<u64> = {
        let mut copied = template1.clone();
        let mut index = 0;
        copied.retain(|_| {
            index += 1;
            !(20..50).contains(&(index - 1))
        });
        copied
    };
    eprintln!("{}\t{}", template1.len(), template2.len());
    let coverage = 10;
    // println!("Read Num:{}", read_num);
    println!("Seed\tLK\tElapsed\tER\tEN\tSub");
    let seed = 4357869001101629047;
    test(&template1, &template2, coverage, 0f64, seed, Some(8192));
}

fn test(
    template1: &[u64],
    template2: &[u64],
    coverage: usize,
    error_rate: f64,
    seed: u64,
    sub: Option<usize>,
) {
    let read_num = coverage * (template1.len() + template2.len()) / 4;
    let mut counter: HashMap<u64, u64> = HashMap::new();
    // (UnitID, OccurenceesSoFar).
    let hap1: Vec<_> = template1
        .iter()
        .map(|&unit| {
            let c = counter.entry(unit).or_default();
            *c += 1;
            (unit, *c - 1)
        })
        .collect();
    let hap2: Vec<_> = template2
        .iter()
        .map(|&unit| {
            let c = counter.entry(unit).or_default();
            *c += 1;
            (unit, *c - 1)
        })
        .collect();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
    let reads: Vec<_> = (0..read_num)
        .map(|x| {
            let id = format!("{}", x);
            let hap = if 2 * x / read_num == 0 { &hap1 } else { &hap2 };
            let length = rng.gen_range(MIN_LENGTH..=MAX_LENGTH);
            let path = sim_read(&hap, &mut rng, length, error_rate, &counter);
            (id, path)
        })
        .collect();
    let start = std::time::Instant::now();
    let rng_seed = rng.gen::<u64>();
    let (lk, result) = (0..1u64)
        .into_iter()
        .map(|i| {
            let mut reads = reads.clone();
            let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(rng_seed + i);
            reads.shuffle(&mut rng);
            let (lk, asn) = path_phasing::phase_with_lk(&reads, 20, 6, sub);
            let asn: HashMap<String, _> =
                asn.into_iter().map(|(x, y)| (x.to_string(), y)).collect();
            eprintln!("{},{:.3}", i, lk);
            (lk, asn)
        })
        .max_by(|x, y| (x.0).partial_cmp(&(y.0)).unwrap())
        .unwrap();
    let elapsed = (std::time::Instant::now() - start).as_millis();
    let errors = reads
        .iter()
        .filter(|(id, _path)| {
            assert!(result.contains_key(id.as_str()), "{}", id);
            let pred = *result.get(id.as_str()).unwrap();
            let i = id.parse::<usize>().unwrap();
            let answer = (2 * i / read_num) as u8;
            if answer != pred {
                true
            } else {
                false
            }
        })
        .count();
    let accuracy = (read_num - errors) as f64 / (reads.len() as f64);
    let (accuracy, _errors) = if 0.5 < accuracy {
        (accuracy, errors)
    } else {
        (1. - accuracy, read_num - errors)
    };
    let sub = sub.unwrap();
    println!(
        "{}\t{:.3}\t{}\t{:.3}\t{}\t{:.3}\t{}",
        seed, lk, coverage, error_rate, accuracy, elapsed, sub
    );
    assert!(error_rate < 0.10);
}

fn sim_read<R: Rng>(
    hap: &[(u64, u64)],
    rng: &mut R,
    length: usize,
    error_rate: f64,
    unit_counts: &HashMap<u64, u64>,
) -> Vec<(u64, u64)> {
    let start_position = rng.gen_range(0..hap.len() - length);
    let mut path: Vec<_> = hap
        .iter()
        .skip(start_position)
        .take(length)
        .map(|&(n, c)| {
            if rng.gen_bool(error_rate) {
                (n, rng.gen::<u64>() % (unit_counts[&n] + 1))
            } else {
                (n, c)
            }
        })
        .collect();
    if rng.gen_bool(0.5) {
        path.reverse();
    }
    path
}
