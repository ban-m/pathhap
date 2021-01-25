// use path_phasing::haplotype_cc;
use log::debug;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::HashMap;
fn sim_path_error<R: Rng>(
    template: &[(u64, u64)],
    cluster_num: &HashMap<u64, u64>,
    rng: &mut R,
    min_len: usize,
    max_len: usize,
    error_rate: f64,
) -> Vec<(u64, u64)> {
    let length = rng.gen::<usize>() % (max_len + 1 - min_len) + min_len;
    let start_position = rng.gen::<usize>() % (template.len() - length);
    let end_position = start_position + length;
    let paths: Vec<_> = if rng.gen_bool(0.5) {
        template[start_position..end_position].to_vec()
    } else {
        template
            .iter()
            .skip(start_position)
            .take(length)
            .rev()
            .cloned()
            .collect()
    };
    paths
        .into_iter()
        .map(|(node, cluster)| {
            if rng.gen_bool(error_rate) {
                (node, rng.gen_range(0..cluster_num[&node]))
            } else {
                (node, cluster)
            }
        })
        .collect()
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let template1: Vec<_> = (0..10).map(|x| (x, 0)).collect();
    let template2: Vec<_> = (0..10).map(|x| (x, 1)).collect();
    let cluster_num: HashMap<u64, u64> = (0..10).map(|x| (x, 2)).collect();
    let min_len = 3;
    let max_len = 6;
    let err = 0.1;
    let path_num = 10;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
    let mut paths: Vec<_> = (0..path_num)
        .map(|i| {
            let path = if i < path_num / 2 {
                sim_path_error(&template1, &cluster_num, &mut rng, min_len, max_len, err)
            } else {
                sim_path_error(&template2, &cluster_num, &mut rng, min_len, max_len, err)
            };
            (format!("{}", i), path)
        })
        .collect();
    for (i, p) in paths.iter() {
        debug!("{}\t{:?}", i, p);
    }
    paths.shuffle(&mut rng);
    let result = path_phasing::haplotyping(&paths, 20, 0.05);
    let cluster1 = *result.get("0").unwrap();
    let cluster2: String = format!("{}", path_num - 1);
    let cluster2 = *result.get(cluster2.as_str()).unwrap();
    assert_ne!(cluster1, cluster2);
    for i in 0..path_num {
        let id: String = format!("{}", i);
        if i < path_num / 2 {
            assert_eq!(cluster1, result[id.as_str()]);
        } else {
            assert_eq!(cluster2, result[id.as_str()]);
        }
    }
}
