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
use std::collections::HashMap;
mod clustering;
pub mod em_progressive;
mod find_union;
pub use clustering::haplotype_cc;
pub use clustering::phase_cc;
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
pub fn phase(paths: &[(String, Vec<(u64, u64)>)], max_occ: usize) -> HashMap<&str, u8> {
    phase_with_lk(paths, max_occ).1
}

pub fn phase_with_lk(
    paths: &[(String, Vec<(u64, u64)>)],
    max_occ: usize,
) -> (f64, HashMap<&str, u8>) {
    // First, decompose into each connected component.
    // let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(seed);
    // use rand::seq::SliceRandom;
    debug!("Begin");
    let components = split_paths(paths);
    debug!("NumOfCC\t{}", components.len());
    let mut total_lk = 0.;
    let phasing: HashMap<_, u8> = components
        .iter()
        .map(|paths| {
            let normed_paths = normalize_path(paths);
            let phased_paths = phase_cc(&normed_paths, max_occ);
            // Polish clustering by EM algorthm
            let (phased_paths, lk) = em_progressive::em_clustering(&normed_paths, &phased_paths);
            debug!("Maximum likelihood:{:.3}", lk);
            total_lk += lk;
            phased_paths
                .iter()
                .zip(paths)
                .map(|(&phase, (id, _))| (id.as_str(), phase))
                .collect::<Vec<(&str, u8)>>()
        })
        .fold(HashMap::new(), |mut acc, y| {
            acc.extend(y);
            acc
        });
    debug!("Total log likelihood:{:.3}", total_lk);
    (total_lk, phasing)
}

pub fn haplotyping(
    paths: &[(String, Vec<(u64, u64)>)],
    max_length: usize,
    er: f64,
) -> HashMap<&str, u8> {
    let components = split_paths(paths);
    debug!("NumOfCC\t{}", components.len());
    let mut total_lk = 0.;
    let phasing: HashMap<_, u8> = components
        .iter()
        .map(|paths| {
            let normed_paths = normalize_path(paths);
            // Check error_rate is not too large.
            {
                let mut cluster_num: HashMap<_, usize> = HashMap::new();
                for p in normed_paths.iter() {
                    for &(n, c) in p.iter() {
                        let cl = cluster_num.entry(n).or_default();
                        *cl = (*cl).max(c + 1);
                    }
                }
                for (n, c) in cluster_num {
                    if 1f64 < c as f64 * er {
                        error!("Error rate is too large.");
                        error!("Some node(node {}) has {} clusters.", n, c);
                        error!("Node {} is the {}-th smallest nodes in the dataset.", n, n);
                        panic!("Panic from above error.");
                    }
                }
            }
            let (phased_paths, lk) = haplotype_cc(&normed_paths, max_length, er);
            debug!("Maximum likelihood:{:.3}", lk);
            total_lk += lk;
            phased_paths
                .iter()
                .zip(paths)
                .map(|(&phase, (id, _))| (id.as_str(), phase))
                .collect::<Vec<(&str, u8)>>()
        })
        .fold(HashMap::new(), |mut acc, y| {
            acc.extend(y);
            acc
        });
    debug!("Total log likelihood:{:.3}", total_lk);
    phasing
}

type PathWithID = (String, Vec<(u64, u64)>);
fn split_paths(paths: &[PathWithID]) -> Vec<Vec<&PathWithID>> {
    let mut fu = find_union::FindUnion::new(paths.len());
    let mut path_on_nodes: HashMap<_, Vec<usize>> = HashMap::new();
    for (idx, (_, path)) in paths.iter().enumerate() {
        for &(n, _) in path.iter() {
            path_on_nodes.entry(n).or_default().push(idx);
        }
    }
    for indices in path_on_nodes.values() {
        for pair in indices.windows(2) {
            fu.unite(pair[0], pair[1]).unwrap();
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

// Normalize path. In other words, each path would be converted into a
// format so that it has continuous node id, continuoud cluster number.
// The order of the input path should be the same as the output order.
// To see the example, see test.
fn normalize_path(paths: &[&(String, Vec<(u64, u64)>)]) -> Vec<Vec<(usize, usize)>> {
    let nodes: HashMap<u64, usize> = {
        let mut nodes: Vec<_> = paths
            .iter()
            .flat_map(|(_, p)| p.iter())
            .map(|x| x.0)
            .collect();
        nodes.sort_unstable();
        nodes.dedup();
        nodes.into_iter().enumerate().map(|(i, x)| (x, i)).collect()
    };
    let clusters: Vec<HashMap<u64, usize>> = {
        let mut cluster: Vec<Vec<_>> = vec![vec![]; nodes.len()];
        for (_, path) in paths.iter() {
            for (n, c) in path.iter() {
                cluster[nodes[n]].push(*c);
            }
        }
        cluster.iter_mut().for_each(|xs| {
            xs.sort_unstable();
            xs.dedup();
        });
        cluster
            .into_iter()
            .map(|xs| xs.into_iter().enumerate().map(|(x, y)| (y, x)).collect())
            .collect()
    };
    paths
        .iter()
        .map(|(_, path)| {
            let path: Vec<_> = path
                .iter()
                .map(|(n, c)| {
                    let node = nodes[n];
                    let cluster = clusters[node][c];
                    (node, cluster)
                })
                .collect();
            path
        })
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::seq::SliceRandom;
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro256PlusPlus;
    #[test]
    fn it_works() {
        assert!(true);
    }
    #[test]
    fn normalize_path_test() {
        let paths: Vec<_> = vec![
            ("0".to_string(), vec![(0, 0), (2, 3), (3, 0), (6, 0)]),
            ("1".to_string(), vec![(0, 1), (2, 4), (3, 1), (6, 1)]),
            ("2".to_string(), vec![(6, 3), (3, 5)]),
        ];
        let paths: Vec<_> = paths.iter().collect();
        let normed_paths = normalize_path(&paths);
        assert_eq!(normed_paths[0], vec![(0, 0), (1, 0), (2, 0), (3, 0)]);
        assert_eq!(normed_paths[1], vec![(0, 1), (1, 1), (2, 1), (3, 1)]);
        assert_eq!(normed_paths[2], vec![(3, 2), (2, 2)]);
    }
    #[test]
    fn phase_test_1() {
        let paths = vec![
            ("ID1".to_string(), vec![(0, 0), (1, 0), (2, 0), (3, 0)]),
            ("ID2".to_string(), vec![(0, 0), (1, 0), (2, 0), (3, 0)]),
            ("ID3".to_string(), vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
            ("ID4".to_string(), vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
        ];
        let result = phase(&paths, 14);
        eprintln!("{:?}", result);
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
        let result = phase(&paths, 14);
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
    fn phase_test_random_linear_noerror() {
        let template1: Vec<_> = (0..10).map(|x| (x, 0)).collect();
        let template2: Vec<_> = (0..10).map(|x| (x, 1)).collect();
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
        {
            let mut result: Vec<_> = result.iter().collect();
            result.sort();
            let mut paths = paths.clone();
            paths.sort_by(|x, y| (x.0).cmp(&y.0));
            for ((id, val), (_, p)) in result.iter().zip(paths.iter()) {
                eprintln!("{}\t{}\t{:?}", id, val, p);
            }
        }
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
    #[test]
    fn hap_test_random_hard() {
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
        let result = haplotyping(&paths, 20, 0.05);
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
    fn hap_test_random_linear_error() {
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
        paths.shuffle(&mut rng);
        let result = haplotyping(&paths, 20, 0.05);
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
