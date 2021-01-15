use path_phasing::phase;
use rand::{seq::SliceRandom, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
fn sim_path<R: Rng>(
    template: &[(u64, u64)],
    rng: &mut R,
    min_len: usize,
    max_len: usize,
) -> Vec<(u64, u64)> {
    let length = rng.gen_range(min_len..max_len + 1);
    let start_position = rng.gen_range(0..template.len() - length);
    let end_position = start_position + length;
    if rng.gen_bool(0.5) {
        template[start_position..end_position].to_vec()
    } else {
        template
            .iter()
            .skip(start_position)
            .take(length)
            .rev()
            .cloned()
            .collect()
    }
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let template1: Vec<(u64, u64)> = vec![0, 1, 2, 3, 4, 5, 8, 9, 10]
        .into_iter()
        .map(|n| (n, 0))
        .collect();
    // eprintln!("{:?}", template1);
    let template2: Vec<(u64, u64)> = vec![0, 1, 2, 3, 6, 7, 8, 9, 10]
        .into_iter()
        .map(|n| (n, 1))
        .collect();
    // eprintln!("{:?}", template2);
    let min_len = 3;
    let max_len = 6;
    let path_num = 10;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
    let mut paths: Vec<_> = (0..path_num)
        .map(|i| {
            let path = if i < path_num / 2 {
                sim_path(&template1, &mut rng, min_len, max_len)
            } else {
                sim_path(&template2, &mut rng, min_len, max_len)
            };
            (format!("{}", i), path)
        })
        .collect();
    paths.shuffle(&mut rng);
    let result = phase(&paths, 20);
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
