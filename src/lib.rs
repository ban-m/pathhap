//! A library to provide path clustering algorithm.
//! It has several low-level functions and one top-level wrapping API.
//! Usually, users only need to check the [phase](phase) function.
//! # Algorithm
//! Under the hood, [phase] function consists of distinct steps:
//! 1. It construct a graph from the given input.
//! 2. If it has several connected components, it solves path clustering problem on each connected component.
//! 3. If a connected component has a node with occurence more than `max_occ`, then it down-samples the path until all nodes has occurence less than or equals to `max_occ`
//! 4. It solves path clustering problem incrementaly.
//! 5. It predicts the clustering of discarded reads.
use log::*;
use std::collections::{HashMap, HashSet};
pub mod em_progressive;
mod exact_phase;
mod find_union;
mod model;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro128StarStar;
/// Phasing API.
/// Phase given paths into two haplotypes.
/// # Example
/// let paths = vec![
///     ("ID1".to_string(), vec![(0, 0), (1, 0), (2, 0), (3, 0)]),
///     ("ID2".to_string(), vec![(0, 0), (1, 0), (2, 0), (3, 0)]),
///     ("ID3".to_string(), vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
///     ("ID4".to_string(), vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
/// ];
/// let reuslt = phase(&paths, 15, 24);
/// assert_eq!(
///     result,
///     vec![
///         ("ID1".to_string(), 0),
///         ("ID2".to_string(), 0),
///         ("ID3".to_string(), 1),
///         ("ID4".to_string(), 1)
///     ]
/// );
/// // More complicated example.
/// let path1 = (
///     "ID1".to_string(),
///     vec![(0, 0), (1, 0), (10, 0), (3, 0), (10, 4)],
/// );
/// let path2 = (
///     "ID2".to_string(),
///     vec![(1, 0), (10, 0), (3, 0), (10, 4), (110, 1)],
/// );
/// let path3 = (
///     "ID3".to_string(),
///     vec![(0, 1), (1, 1), (10, 1), (3, 1), (10, 1), (110, 0)],
/// );
/// let path3 = ("ID4".to_string(), vec![(10, 1), (3, 1), (10, 1)]);
/// let paths = vec![path1, path2, path3, path4];
/// let reuslt = phase(&paths, 15 ,24);
/// assert_eq!(
///     result,
///     vec![
///         ("ID1".to_string(), 0),
///         ("ID2".to_string(), 0),
///         ("ID3".to_string(), 1),
///         ("ID4".to_string(), 1)
///     ]
/// );
pub fn phase<'a>(
    paths: &'a [(String, Vec<(u64, u64)>)],
    max_occ: usize,
    seed: u64,
) -> HashMap<&'a str, u8> {
    const REPEAT_NUM: usize = 1;
    // First, decompose into each connected component.
    let mut result = HashMap::new();
    let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(seed);
    use rand::seq::SliceRandom;
    let mut components = split_paths(paths);
    debug!("NumOfCC\t{}", components.len());
    let mut total_lk = 0.;
    for paths in components.iter_mut() {
        paths.shuffle(&mut rng);
        let (phased_paths, lk) = (0..REPEAT_NUM)
            .map(|_| phase_cc(&paths, max_occ))
            .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
            .unwrap();
        debug!("Maximum likelihood:{:.3}", lk);
        total_lk += lk;
        result.extend(phased_paths);
    }
    debug!("Total log likelihood:{:.3}", total_lk);
    result
}

fn split_paths(paths: &[(String, Vec<(u64, u64)>)]) -> Vec<Vec<&(String, Vec<(u64, u64)>)>> {
    let mut fu = find_union::FindUnion::new(paths.len());
    let path_sketch: Vec<HashSet<_>> = paths
        .iter()
        .map(|p| p.1.iter().map(|x| x.0).collect())
        .collect();
    for (idx, p1) in path_sketch.iter().enumerate() {
        for (jdx, p2) in path_sketch.iter().enumerate().skip(idx + 1) {
            if !p1.is_disjoint(p2) {
                fu.unite(idx, jdx);
            }
        }
    }
    let mut cluster = vec![];
    for idx in 0..paths.len() {
        if fu.find(idx).unwrap() == idx {
            let mut cl = vec![];
            for (jdx, p) in paths.iter().enumerate() {
                if fu.find(jdx).unwrap() == idx {
                    cl.push(p);
                }
            }
            cluster.push(cl);
        }
    }
    assert_eq!(cluster.iter().map(|x| x.len()).sum::<usize>(), paths.len());
    cluster
}

// Phaing connected components.
fn phase_cc<'a>(
    paths: &[&'a (String, Vec<(u64, u64)>)],
    max_occ: usize,
) -> (HashMap<&'a str, u8>, f64) {
    // Re-numbering nodes.
    debug!("Start");
    let idx2id: HashMap<usize, _> = paths.iter().enumerate().map(|(i, p)| (i, &p.0)).collect();
    let nodes: HashMap<u64, usize> = {
        let nodes: HashSet<_> = paths.iter().flat_map(|x| x.1.iter().map(|x| x.0)).collect();
        let mut nodes: Vec<_> = nodes.into_iter().collect();
        nodes.sort();
        nodes
            .into_iter()
            .enumerate()
            .map(|(idx, n)| (n, idx))
            .collect()
    };
    let mut paths: Vec<(usize, Vec<(usize, usize)>)> = paths
        .iter()
        .enumerate()
        .map(|(id, (_, path))| {
            let path: Vec<_> = path.iter().map(|&(n, c)| (nodes[&n], c as usize)).collect();
            (id, path)
        })
        .collect();
    if false {
        let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(24);
        let (result, lk) = em_progressive::em_progressive_clustering(&paths, &mut rng);
        let result: HashMap<_, _> = result
            .iter()
            .map(|&(idx, hap)| (idx2id[&idx].as_str(), hap))
            .collect();

        return (result, lk);
    }
    debug!(
        "Renaming {} paths. number of nodes:{}",
        paths.len(),
        nodes.len()
    );
    // Construct a graph.
    // This is the coverage, or occurence, of a node.
    let mut nodes: Vec<usize> = {
        let mut nodes = vec![0; nodes.len()];
        for (_, path) in paths.iter() {
            for &(n, _) in path.iter() {
                nodes[n] += 1;
            }
        }
        nodes
    };
    let mut removed_paths = vec![];
    paths.sort_by_key(|x| x.1.len());
    loop {
        let (max_node, &max) = nodes.iter().enumerate().max_by_key(|x| x.1).unwrap();
        if max <= max_occ {
            break;
        }
        // Remove shortest path containing max_node.
        let idx: usize = paths
            .iter()
            .position(|x| x.1.iter().any(|x| x.0 == max_node))
            .unwrap();
        let path = paths.remove(idx);
        assert!(path.1.iter().any(|x| x.0 == max_node));
        for &(n, _) in path.1.iter() {
            nodes[n] -= 1;
        }
        removed_paths.push(path);
    }
    debug!("Removed:Remaines={}:{}", removed_paths.len(), paths.len());
    let mut result = exact_phase::exact_phase(&paths);
    // PSEUOD count 1.
    let model = model::Model::new(&paths, &result, 1);
    debug!("Finished");
    for (idx, path) in removed_paths {
        result.push((idx, model.predict_path(&path)));
        paths.push((idx, path));
    }
    assert_eq!(paths.len(), idx2id.len());
    let model = model::Model::new(&paths, &result, 0);
    let likelihoods = paths.iter().map(|(_, x)| model.likelihood(x)).sum::<f64>();
    debug!("LK:{:.4}", likelihoods);
    // result.sort_by_key(|x| x.0);
    let result: HashMap<_, _> = result
        .iter()
        .map(|&(idx, hap)| (idx2id[&idx].as_str(), hap))
        .collect();
    (result, likelihoods)
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro256PlusPlus;
    #[test]
    fn it_works() {
        assert!(true);
    }
    #[test]
    fn phase_test_1() {
        let paths = vec![
            ("ID1".to_string(), vec![(0, 0), (1, 0), (2, 0), (3, 0)]),
            ("ID2".to_string(), vec![(0, 0), (1, 0), (2, 0), (3, 0)]),
            ("ID3".to_string(), vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
            ("ID4".to_string(), vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
        ];
        let result = phase(&paths, 14, 24);
        assert_eq!(result["ID1"], result["ID2"]);
        assert_eq!(result["ID3"], result["ID4"]);
        assert_ne!(result["ID1"], result["ID4"]);
    }
    #[test]
    fn phase_test_2() {
        // More complicated example.
        let path1 = (
            "ID1".to_string(),
            vec![(0, 0), (1, 0), (10, 0), (3, 0), (10, 4)],
        );
        let path2 = (
            "ID2".to_string(),
            vec![(1, 0), (10, 0), (3, 0), (10, 4), (110, 1)],
        );
        let path3 = (
            "ID3".to_string(),
            vec![(0, 1), (1, 1), (10, 1), (3, 1), (10, 1), (110, 0)],
        );
        let path4 = ("ID4".to_string(), vec![(10, 1), (3, 1), (10, 1)]);
        let paths = vec![path1, path2, path3, path4];
        let result = phase(&paths, 14, 24);
        assert_eq!(result["ID1"], result["ID2"]);
        assert_eq!(result["ID3"], result["ID4"]);
        assert_ne!(result["ID1"], result["ID4"]);
    }
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
    #[test]
    fn phase_test_random_linear() {
        let template1: Vec<_> = (0..10).map(|x| (x, 0)).collect();
        let template2: Vec<_> = (0..10).map(|x| (x, 1)).collect();
        let min_len = 3;
        let max_len = 6;
        let path_num = 10;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                let path = if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                };
                (format!("{}", i), path)
            })
            .collect();
        let result = phase(&paths, 20, 24);
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
    #[test]
    fn phase_test_random_loop() {
        let template1: Vec<(u64, u64)> = {
            let mut count: HashMap<_, u64> = HashMap::new();
            let nodes = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 3, 9, 10];
            nodes
                .iter()
                .map(|&node| {
                    let count = count.entry(node).or_default();
                    *count += 1;
                    (node, *count - 1)
                })
                .collect()
        };
        // eprintln!("{:?}", template1);
        let template2: Vec<(u64, u64)> = {
            let mut count: HashMap<_, u64> = HashMap::new();
            let nodes = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 3, 9, 10];
            nodes
                .iter()
                .map(|&node| {
                    let count = count.entry(node).or_insert(2);
                    *count += 1;
                    (node, *count - 1)
                })
                .collect()
        };
        // eprintln!("{:?}", template2);
        let min_len = 3;
        let max_len = 6;
        let path_num = 10;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                let path = if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                };
                (format!("{}", i), path)
            })
            .collect();
        let result = phase(&paths, 20, 24);
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
    #[test]
    fn phase_test_random_branch() {
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
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                let path = if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                };
                (format!("{}", i), path)
            })
            .collect();
        let result = phase(&paths, 20, 24);
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
    #[test]
    fn phase_test_random_loop_branch() {
        let template1: Vec<(u64, u64)> = {
            let mut count: HashMap<_, u64> = HashMap::new();
            let nodes = vec![0, 1, 2, 5, 4, 3, 2, 6, 7, 8, 11, 12, 14];
            nodes
                .iter()
                .map(|&node| {
                    let count = count.entry(node).or_default();
                    *count += 1;
                    (node, *count - 1)
                })
                .collect()
        };
        // eprintln!("{:?}", template1);
        let template2: Vec<(u64, u64)> = {
            let mut count: HashMap<_, u64> = HashMap::new();
            let nodes = vec![0, 1, 2, 5, 4, 3, 2, 6, 9, 10, 11, 12, 13, 11, 12, 14];
            nodes
                .iter()
                .map(|&node| {
                    let count = count.entry(node).or_insert(2);
                    *count += 1;
                    (node, *count - 1)
                })
                .collect()
        };
        // eprintln!("{:?}", template2);
        let min_len = 3;
        let max_len = 6;
        let path_num = 10;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                let path = if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                };
                (format!("{}", i), path)
            })
            .collect();
        let result = phase(&paths, 20, 24);
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
    #[test]
    fn phase_test_random_hard() {
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let mut nodes = vec![];
        for _ in 0..20 {
            nodes.push(rng.gen::<u64>() % 100);
        }
        let template1: Vec<(u64, u64)> = {
            let mut count: HashMap<_, u64> = HashMap::new();
            let mut nodes: Vec<_> = nodes.iter().skip(3).copied().collect();
            nodes.sort();
            nodes
                .iter()
                .map(|&node| {
                    let count = count.entry(node).or_default();
                    *count += 1;
                    (node, *count - 1)
                })
                .collect()
        };
        // eprintln!("{:?}", template1);
        let template2: Vec<(u64, u64)> = {
            let mut count: HashMap<_, u64> = HashMap::new();
            let mut nodes: Vec<_> = nodes;
            nodes.sort();
            nodes
                .iter()
                .map(|&node| {
                    let count = count.entry(node).or_insert(2);
                    *count += 1;
                    (node, *count - 1)
                })
                .collect()
        };
        // eprintln!("{:?}", template2);
        let min_len = 3;
        let max_len = 6;
        let path_num = 10;
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                let path = if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                };
                (format!("{}", i), path)
            })
            .collect();
        let result = phase(&paths, 20, 24);
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
    #[test]
    fn phase_test_random_linear_error() {
        let template1: Vec<_> = (0..10).map(|x| (x, 0)).collect();
        let template2: Vec<_> = (0..10).map(|x| (x, 1)).collect();
        let cluster_num: HashMap<u64, u64> = (0..10).map(|x| (x, 2)).collect();
        let min_len = 3;
        let max_len = 6;
        let err = 0.1;
        let path_num = 10;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                let path = if i < path_num / 2 {
                    sim_path_error(&template1, &cluster_num, &mut rng, min_len, max_len, err)
                } else {
                    sim_path_error(&template2, &cluster_num, &mut rng, min_len, max_len, err)
                };
                (format!("{}", i), path)
            })
            .collect();
        let result = phase(&paths, 20, 24);
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
}
