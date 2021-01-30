use log::debug;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    // Parameters.
    let template_length = 500;
    let duplication_rate = 0.01;
    let hap_dup_rate = 0.001;
    let hap_ins_rate = 0.001;
    let hap_del_rate = 0.001;
    let read_num = 3000;
    let min_length = 3;
    let max_length = 6;
    let error_rate = 0.05;
    let iteration: u64 = 10;
    let sub = None;
    // let sub = Some(10_000);
    let elapsed: Vec<_> = (0..iteration)
        .map(|seed| {
            let start = std::time::Instant::now();
            test(
                template_length,
                (duplication_rate, hap_dup_rate, hap_ins_rate, hap_del_rate),
                (read_num, min_length, max_length, error_rate),
                seed,
                sub,
            );
            let ti = (std::time::Instant::now() - start).as_millis();
            println!("{}", ti);
            ti
        })
        .collect();
    let mean = elapsed.iter().sum::<u128>() as f64 / iteration as f64;
    let sd: f64 =
        elapsed.iter().map(|x| x * x).sum::<u128>() as f64 / iteration as f64 - mean * mean;
    println!("{:.4}(SD:{:.4})", mean, sd);
}

fn test(
    template_length: usize,
    (duplication_rate, hap_dup_rate, hap_ins_rate, hap_del_rate): (f64, f64, f64, f64),
    (read_num, min_length, max_length, error_rate): (usize, usize, usize, f64),
    seed: u64,
    sub: Option<usize>,
) {
    debug!("Start Test module");
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
    let mut counter: HashMap<_, u64> = HashMap::new();
    // (UnitID, OccurenceesSoFar).
    let template: Vec<(_, u64)> = (0..template_length)
        .map(|_| {
            if rng.gen_bool(duplication_rate) {
                // Reuse unit.
                let prev_unit: u64 = (rng.gen::<usize>() % counter.len()) as u64;
                let unit = (prev_unit, counter[&prev_unit]);
                *counter.get_mut(&prev_unit).unwrap() += 1;
                unit
            } else {
                let unit = counter.len() as u64;
                counter.insert(unit, 0);
                (unit, 0)
            }
        })
        .collect();
    // eprintln!("{:?}", template);
    // Let's make diploids.
    let hap1: Vec<_> = template
        .iter()
        .filter_map(|(x, y)| {
            if rng.gen_bool(hap_dup_rate) {
                // Pick some nodes and duplicate it.
                let &(unit, occs) = template.choose(&mut rng).unwrap();
                Some((unit, 2 * (occs + 1)))
            } else if rng.gen_bool(hap_ins_rate) {
                let new_unit = counter.len() as u64;
                counter.insert(new_unit, 0);
                Some((new_unit, 0))
            } else if rng.gen_bool(hap_del_rate) {
                None
            } else {
                Some((*x, 2 * y))
            }
        })
        .collect();
    // eprintln!("hap1:{:?}", hap1);
    let hap2: Vec<_> = template
        .iter()
        .filter_map(|(x, y)| {
            if rng.gen_bool(hap_dup_rate) {
                // Pick some nodes and duplicate it.
                let &(unit, occs) = template.choose(&mut rng).unwrap();
                Some((unit, 2 * (occs + 1) + 1))
            } else if rng.gen_bool(hap_ins_rate) {
                let new_unit = counter.len() as u64;
                counter.insert(new_unit, 0);
                Some((new_unit, 0))
            } else if rng.gen_bool(hap_del_rate) {
                None
            } else {
                Some((*x, 2 * y + 1))
            }
        })
        .collect();
    // eprintln!("hap2:{:?}", hap2);
    debug!("Sim reads.");
    let mut reads: Vec<_> = (0..read_num)
        .map(|x| {
            let id = format!("{}", x);
            let hap = if 2 * x / read_num == 0 { &hap1 } else { &hap2 };
            let path = sim_read(&hap, &mut rng, min_length, max_length, error_rate, &counter);
            (id, path)
        })
        .collect();
    reads.shuffle(&mut rng);
    let (_, result) = path_phasing::phase_with_lk(&reads, 22, sub);
    let hap1 = result.get(&"0");
    let hap2_id = format!("{}", read_num - 1);
    let hap2 = result.get(hap2_id.as_str());
    let errors = reads
        .iter()
        .filter(|(id, _)| {
            assert!(result.contains_key(id.as_str()), "{}", id);
            let pred = result.get(id.as_str());
            let i = id.parse::<usize>().unwrap();
            let answer = if 2 * i / read_num == 0 { hap1 } else { hap2 };
            match answer != pred {
                true => true,
                x => x,
            }
        })
        .count();
    let error_rate = errors as f64 / (reads.len() as f64);
    eprintln!("Error Rate:{:?}", error_rate);
    assert!(error_rate < 0.1);
}

fn sim_read<R: Rng>(
    hap: &[(u64, u64)],
    rng: &mut R,
    min_length: usize,
    max_length: usize,
    error_rate: f64,
    unit_counts: &HashMap<u64, u64>,
) -> Vec<(u64, u64)> {
    let length = rng.gen_range(min_length..=max_length);
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
