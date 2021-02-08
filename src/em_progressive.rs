use log::debug;
use rand::Rng;
use std::collections::{HashMap, HashSet};
const PSEUDO_COUNT: f64 = 1f64;
pub fn em_clustering(paths: &[Vec<(usize, usize)>], init_cluster: &[u8]) -> (Vec<u8>, f64) {
    assert_eq!(paths.len(), init_cluster.len());
    let mut weights: Vec<_> = init_cluster
        .iter()
        .map(|&x| if x == 0 { 1f64 } else { 0f64 })
        .collect();
    let mut model = Model::create_dip(&paths, &weights);
    debug!("Start EM algorithm.");
    let mut lk = paths.iter().map(|path| model.lk(path)).sum::<f64>();
    for i in 0.. {
        debug!("LK\t{}\t{:?}", i, lk);
        for (w, path) in weights.iter_mut().zip(paths.iter()) {
            *w = model.weights_dip(&path);
        }
        model.update_dip(&paths, &weights);
        let next_lk = paths.iter().map(|path| model.lk(path)).sum::<f64>();
        if next_lk - lk < 0.000001 {
            break;
        } else {
            lk = next_lk;
        }
    }
    let assignments: Vec<_> = weights
        .iter()
        .map(|&h1| if 0.5 < h1 { 0 } else { 1 })
        .collect();
    // debug!("ASN\tBefore\tAfter\tPath");
    // for ((b, a), p) in assignments.iter().zip(init_cluster).zip(paths.iter()) {
    //     debug!("ASN\t{}\t{}\t{:?}", b, a, p);
    // }
    (assignments, lk)
}
pub fn em_progressive_clustering<R: Rng>(
    paths: &[(usize, Vec<(usize, usize)>)],
    rng: &mut R,
) -> (Vec<(usize, u8)>, f64) {
    let (nodes, edges) = {
        let nodes = paths
            .iter()
            .filter_map(|x| x.1.iter().map(|x| x.0).max())
            .max()
            .unwrap()
            + 1;
        let mut edges = vec![HashSet::new(); nodes];
        for (_, path) in paths.iter() {
            for w in path.windows(2) {
                let from = w[0].0;
                let to = w[1].0;
                edges[from].insert(to);
                edges[to].insert(from);
            }
        }
        let edges: Vec<_> = edges
            .into_iter()
            .map(|x| {
                let mut x: Vec<_> = x.into_iter().collect();
                x.sort_unstable();
                x
            })
            .collect();
        (nodes, edges)
    };
    let genotypes = {
        let mut genotypes: HashMap<usize, usize> = HashMap::new();
        for (_, path) in paths.iter() {
            for &(n, c) in path.iter() {
                genotypes
                    .entry(n)
                    .and_modify(|x| *x = (*x).max(c + 1))
                    .or_insert(c + 1);
            }
        }
        genotypes
    };
    let clustering_order = sort_graph(nodes, &edges);
    debug!("Sorted.");
    debug!("{:?}", clustering_order);
    // These weights would be erased if it is not the weights of reads in R(v_0);
    // TODO: Make some random initialization.
    let mut weights_on_reads: Vec<Vec<_>> = (0..paths.len())
        .map(|_| {
            let w = rng.gen_range(0f64..=1f64);
            vec![w, 1f64 - w]
        })
        .collect();
    let mut clustered_units = HashSet::new();
    let selected_paths: HashSet<_> = (0..paths.len()).collect();
    for (i, node) in clustering_order.into_iter().enumerate() {
        debug!("Clustering {}-th node({})", i, node);
        clustered_units.insert(node);
        if clustered_units.len() == 1 {
            for &index in selected_paths.iter() {
                if let Some(&(_, c)) = paths[index].1.iter().find(|&&(n, _)| n == node) {
                    weights_on_reads[index].iter_mut().for_each(|x| *x = 0f64);
                    weights_on_reads[index][c] = 1f64;
                }
            }
            continue;
        } else {
            clustering_on(
                &mut clustered_units,
                &genotypes,
                &mut weights_on_reads,
                paths,
                &selected_paths,
            );
        }
    }
    let model = Model::create_on(&clustered_units, &genotypes, &weights_on_reads, paths);
    let likelihoods = paths.iter().map(|(_, path)| model.lk(path)).sum::<f64>();
    debug!("Likelihoods:{}", likelihoods);
    let assignments: Vec<_> = weights_on_reads
        .iter()
        .zip(paths.iter())
        .map(|(ws, &(id, _))| if ws[1] < ws[0] { (id, 0) } else { (id, 1) })
        .collect();
    (assignments, likelihoods)
}

// EM clustering on reads on specified nodes.
// To create the first prediction, we first construct an old model, then extropolate the prediction
fn clustering_on(
    clustered_units: &mut HashSet<usize>,
    genotypes: &HashMap<usize, usize>,
    weights: &mut [Vec<f64>],
    paths: &[(usize, Vec<(usize, usize)>)],
    selected_path: &HashSet<usize>,
) {
    let mut model = Model::create_on(clustered_units, genotypes, weights, paths);
    let mut likelihoods = selected_path
        .iter()
        .map(|&i| model.lk(&paths[i].1))
        .sum::<f64>();
    loop {
        for &index in selected_path.iter() {
            weights[index] = model.weights(&paths[index].1);
        }
        model = Model::create_on(clustered_units, genotypes, weights, paths);
        let new_lk = selected_path
            .iter()
            .map(|&i| model.lk(&paths[i].1))
            .sum::<f64>();
        if new_lk - likelihoods < 0.0001 {
            break;
        } else {
            likelihoods = new_lk;
        }
    }
}

#[derive(Debug, Clone)]
pub struct Model {
    node_num: usize,
    // The probability to observe the n-th node.
    node_fraction: Vec<f64>,
    // The probability to hop in hap1 or hap2 given the node n.
    hap1_node: Vec<f64>,
    hap2_node: Vec<f64>,
    // The probability for observing cluster c given n in the haplotype.
    hap1: Vec<Vec<f64>>,
    hap2: Vec<Vec<f64>>,
}

impl Model {
    const SMALL: f64 = 0.0000001;
    pub fn create(paths: &[Vec<(usize, usize)>], weights: &[Vec<f64>]) -> Self {
        // +1 is needed!
        let node_num = paths
            .iter()
            .filter_map(|x| x.iter().map(|x| x.0 + 1).max())
            .max()
            .unwrap();
        let mut cluster_num = vec![0; node_num];
        for path in paths.iter() {
            for &(n, c) in path.iter() {
                cluster_num[n] = cluster_num[n].max(c + 1);
            }
        }
        let mut hap1: Vec<_> = cluster_num.iter().map(|&c| vec![PSEUDO_COUNT; c]).collect();
        let mut hap2 = hap1.clone();
        for (ws, path) in weights.iter().zip(paths.iter()) {
            for &(n, c) in path.iter() {
                hap1[n][c] += ws[0];
                hap2[n][c] += ws[1];
            }
        }
        Self::new_inner(hap1, hap2, node_num)
    }
    pub fn create_with_id(paths: &[(usize, &[(usize, usize)])], assign: &[(usize, u8)]) -> Self {
        // +1 is needed!
        let node_num = paths
            .iter()
            .filter_map(|(_, x)| x.iter().map(|x| x.0 + 1).max())
            .max()
            .unwrap();
        let mut cluster_num = vec![0; node_num];
        let assign: HashMap<usize, u8> = assign.iter().copied().collect();
        for (_, path) in paths.iter() {
            for &(n, c) in path.iter() {
                cluster_num[n] = cluster_num[n].max(c + 1);
            }
        }
        let mut hap1: Vec<_> = cluster_num.iter().map(|&c| vec![PSEUDO_COUNT; c]).collect();
        let mut hap2 = hap1.clone();
        for (id, path) in paths.iter() {
            let hap = if assign[id] == 0 {
                &mut hap1
            } else {
                &mut hap2
            };
            for &(n, c) in path.iter() {
                hap[n][c] += 1f64;
            }
        }
        Self::new_inner(hap1, hap2, node_num)
    }
    fn new_inner(mut hap1: Vec<Vec<f64>>, mut hap2: Vec<Vec<f64>>, node_num: usize) -> Self {
        let total_weight = hap1.iter().chain(hap2.iter()).flatten().sum::<f64>();
        let mut node_fraction: Vec<_> = hap1
            .iter()
            .zip(hap2.iter())
            .map(|(h1, h2)| h1.iter().chain(h2).sum::<f64>())
            .collect();
        let hap1_node: Vec<_> = hap1
            .iter()
            .zip(node_fraction.iter())
            .map(|(h1, total)| h1.iter().sum::<f64>() / total)
            .collect();
        let hap2_node: Vec<_> = hap2
            .iter()
            .zip(node_fraction.iter())
            .map(|(h2, total)| h2.iter().sum::<f64>() / total)
            .collect();
        node_fraction.iter_mut().for_each(|x| *x /= total_weight);
        hap1.iter_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        hap2.iter_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        assert!((1f64 - node_fraction.iter().sum::<f64>()).abs() < 0.001);
        for (h1, h2) in hap1_node.iter().zip(hap2_node.iter()) {
            assert!((1f64 - h1 - h2).abs() < 0.0001);
        }
        Self {
            node_fraction,
            node_num,
            hap1_node,
            hap2_node,
            hap1,
            hap2,
        }
    }
    pub fn create_dip(paths: &[Vec<(usize, usize)>], weights: &[f64]) -> Self {
        // +1 is needed!
        let node_num = paths
            .iter()
            .filter_map(|x| x.iter().map(|x| x.0 + 1).max())
            .max()
            .unwrap();
        let mut cluster_num = vec![0; node_num];
        for path in paths.iter() {
            for &(n, c) in path.iter() {
                cluster_num[n] = cluster_num[n].max(c + 1);
            }
        }
        let mut hap1: Vec<_> = cluster_num.iter().map(|&c| vec![0f64; c]).collect();
        let mut hap2 = hap1.clone();
        for (&w, path) in weights.iter().zip(paths.iter()) {
            for &(n, c) in path.iter() {
                hap1[n][c] += w;
                hap2[n][c] += 1. - w;
            }
        }
        Self::new_inner(hap1, hap2, node_num)
    }
    pub fn update(&mut self, paths: &[Vec<(usize, usize)>], weights: &[Vec<f64>]) {
        self.hap1
            .iter_mut()
            .for_each(|xs| xs.iter_mut().for_each(|x| *x = PSEUDO_COUNT));
        self.hap2
            .iter_mut()
            .for_each(|xs| xs.iter_mut().for_each(|x| *x = PSEUDO_COUNT));
        for (ws, path) in weights.iter().zip(paths.iter()) {
            for &(n, c) in path.iter() {
                self.hap1[n][c] += ws[0];
                self.hap2[n][c] += ws[1];
            }
        }
        let total_weight = self
            .hap1
            .iter()
            .chain(self.hap2.iter())
            .flatten()
            .sum::<f64>();
        self.node_fraction
            .iter_mut()
            .zip(self.hap1.iter().zip(self.hap2.iter()))
            .for_each(|(w, (h1, h2))| {
                *w = h1.iter().chain(h2).sum::<f64>();
            });
        self.hap1_node
            .iter_mut()
            .zip(self.hap1.iter().zip(self.node_fraction.iter()))
            .for_each(|(w, (h1, total))| {
                *w = h1.iter().sum::<f64>() / total;
            });
        self.hap2_node
            .iter_mut()
            .zip(self.hap2.iter().zip(self.node_fraction.iter()))
            .for_each(|(w, (h2, total))| {
                *w = h2.iter().sum::<f64>() / total;
            });
        self.node_fraction
            .iter_mut()
            .for_each(|x| *x /= total_weight);
        self.hap1.iter_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        self.hap2.iter_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        assert!((1f64 - self.node_fraction.iter().sum::<f64>()).abs() < 0.001);
        for (h1, h2) in self.hap1_node.iter().zip(self.hap2_node.iter()) {
            assert!((1f64 - h1 - h2).abs() < 0.0001);
        }
    }
    pub fn update_dip(&mut self, paths: &[Vec<(usize, usize)>], weights: &[f64]) {
        self.hap1
            .iter_mut()
            .for_each(|xs| xs.iter_mut().for_each(|x| *x = PSEUDO_COUNT));
        self.hap2
            .iter_mut()
            .for_each(|xs| xs.iter_mut().for_each(|x| *x = PSEUDO_COUNT));
        for (&w, path) in weights.iter().zip(paths.iter()) {
            for &(n, c) in path.iter() {
                self.hap1[n][c] += w;
                self.hap2[n][c] += 1. - w;
            }
        }
        let total_weight = self
            .hap1
            .iter()
            .chain(self.hap2.iter())
            .flatten()
            .sum::<f64>();
        self.node_fraction
            .iter_mut()
            .zip(self.hap1.iter().zip(self.hap2.iter()))
            .for_each(|(w, (h1, h2))| {
                *w = h1.iter().chain(h2).sum::<f64>();
            });
        self.hap1_node
            .iter_mut()
            .zip(self.hap1.iter().zip(self.node_fraction.iter()))
            .for_each(|(w, (h1, total))| {
                *w = h1.iter().sum::<f64>() / total;
            });
        self.hap2_node
            .iter_mut()
            .zip(self.hap2.iter().zip(self.node_fraction.iter()))
            .for_each(|(w, (h2, total))| {
                *w = h2.iter().sum::<f64>() / total;
            });
        self.node_fraction
            .iter_mut()
            .for_each(|x| *x /= total_weight);
        self.hap1.iter_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        self.hap2.iter_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        assert!((1f64 - self.node_fraction.iter().sum::<f64>()).abs() < 0.001);
        for (h1, h2) in self.hap1_node.iter().zip(self.hap2_node.iter()) {
            assert!((1f64 - h1 - h2).abs() < 0.0001);
        }
    }
    pub fn create_on(
        units: &HashSet<usize>,
        genotypes: &HashMap<usize, usize>,
        weights: &[Vec<f64>],
        paths: &[(usize, Vec<(usize, usize)>)],
    ) -> Self {
        let node_num = units.iter().max().copied().unwrap_or_default() + 1;
        let mut hap1: Vec<_> = (0..node_num)
            .map(|n| vec![PSEUDO_COUNT; genotypes.get(&n).copied().unwrap_or_default()])
            .collect();
        let mut hap2 = hap1.clone();
        for (weight, (_, path)) in weights.iter().zip(paths.iter()) {
            for &(n, c) in path.iter().filter(|&(n, _)| units.contains(n)) {
                hap1[n][c] += weight[0];
                hap2[n][c] += weight[1];
            }
        }
        Self::new_inner(hap1, hap2, node_num)
    }
    fn weights_dip(&self, path: &[(usize, usize)]) -> f64 {
        let hap1_w = self.hap1_lk(path);
        let hap2_w = self.hap2_lk(path);
        let ln_sum = logsumexp(hap1_w, hap2_w);
        (hap1_w - ln_sum).exp()
    }

    fn weights(&self, path: &[(usize, usize)]) -> Vec<f64> {
        let hap1_w = self.hap1_lk(path);
        let hap2_w = self.hap2_lk(path);
        let ln_sum = logsumexp(hap1_w, hap2_w);
        vec![(hap1_w - ln_sum).exp(), (hap2_w - ln_sum).exp()]
    }
    #[allow(dead_code)]
    fn weights_update(&self, path: &[(usize, usize)], ws: &mut [f64]) {
        assert_eq!(ws.len(), 2);
        let hap1_w = self.hap1_lk(path);
        let hap2_w = self.hap2_lk(path);
        let ln_sum = logsumexp(hap1_w, hap2_w);
        ws[0] = (hap1_w - ln_sum).exp();
        ws[1] = (hap2_w - ln_sum).exp();
    }

    fn hap1_lk(&self, path: &[(usize, usize)]) -> f64 {
        path.iter()
            .map(|&(n, c)| {
                self.node_fraction[n].max(Self::SMALL).ln()
                    + self.hap1_node[n].max(Self::SMALL).ln()
                    + self.hap1[n].get(c).unwrap_or(&Self::SMALL).ln()
            })
            .sum::<f64>()
    }
    fn hap2_lk(&self, path: &[(usize, usize)]) -> f64 {
        path.iter()
            .map(|&(n, c)| {
                self.node_fraction[n].max(Self::SMALL).ln()
                    + self.hap2_node[n].max(Self::SMALL).ln()
                    + self.hap2[n].get(c).unwrap_or(&Self::SMALL).ln()
            })
            .sum::<f64>()
    }
    fn lk(&self, path: &[(usize, usize)]) -> f64 {
        let hap1_w = self.hap1_lk(path);
        let hap2_w = self.hap2_lk(path);
        logsumexp(hap1_w, hap2_w)
    }
    pub fn predict(&self, path: &[(usize, usize)]) -> u8 {
        (self.hap2_lk(path) < self.hap1_lk(path)) as u8
    }
}

fn logsumexp(x: f64, y: f64) -> f64 {
    let (max, min) = if x < y { (y, x) } else { (x, y) };
    max + (1f64 + (min - max).exp()).ln()
}

fn sort_graph(num_of_nodes: usize, edges: &[Vec<usize>]) -> Vec<usize> {
    use std::collections::VecDeque;
    let mut order = vec![];
    let mut is_arrived = vec![false; num_of_nodes];
    let mut queue = VecDeque::new();
    let start_node = edges
        .iter()
        .enumerate()
        .find(|(_, x)| x.len() == 1)
        .map(|x| x.0)
        .unwrap_or(0);
    queue.push_back(start_node);
    is_arrived[start_node] = true;
    while !queue.is_empty() {
        let node = queue.pop_front().unwrap();
        order.push(node);
        for &to in edges[node].iter() {
            if !is_arrived[to] {
                is_arrived[to] = true;
                queue.push_back(to);
            }
        }
    }
    order
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro128StarStar;
    #[test]
    fn it_works() {
        assert!(true);
    }
    fn sim_path<R: Rng>(
        template: &[(usize, usize)],
        rng: &mut R,
        min_len: usize,
        max_len: usize,
    ) -> Vec<(usize, usize)> {
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
    #[test]
    fn em_test_1() {
        let template1: Vec<(usize, usize)> = (0..10).map(|x| (x, 0)).collect();
        let template2: Vec<(usize, usize)> = (0..10).map(|x| (x, 1)).collect();
        let max_len = 6;
        let min_len = 3;
        let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(39420);
        let path_num = 500;
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                }
            })
            .collect();
        let init: Vec<_> = (0..path_num).map(|i| (2 * i / path_num) as u8).collect();
        let (result, _) = em_clustering(&paths, &init);
        for (i, &pred) in result.iter().enumerate() {
            assert_eq!(pred, (2 * i / path_num) as u8);
        }
    }
    fn sim_path_error<R: Rng>(
        template: &[(usize, usize)],
        cluster_num: &HashMap<usize, usize>,
        rng: &mut R,
        min_len: usize,
        max_len: usize,
        error_rate: f64,
    ) -> Vec<(usize, usize)> {
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
    fn em_test_2() {
        let mut cluster_num: HashMap<_, _> = HashMap::new();
        let template1: Vec<_> = (0..100)
            .map(|x| {
                let c = cluster_num.entry(x).or_default();
                *c += 1;
                (x, *c - 1)
            })
            .collect();
        let template2: Vec<_> = (0..100)
            .map(|x| {
                let c = cluster_num.entry(x).or_default();
                if !(20usize..22usize).contains(&x) {
                    *c += 1;
                }
                (x, *c - 1)
            })
            .collect();
        let max_len = 6;
        let min_len = 3;
        let error_rate = 0.05;
        let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(39420);
        let path_num = 1000;
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                let temp = if i < path_num / 2 {
                    &template1
                } else {
                    &template2
                };
                sim_path_error(&temp, &cluster_num, &mut rng, min_len, max_len, error_rate)
            })
            .collect();
        let init: Vec<_> = (0..path_num)
            .map(|i| {
                let is_error = rng.gen_bool(error_rate) as u8;
                let answer = (2 * i / path_num) as u8;
                answer ^ is_error
            })
            .collect();
        let (result, _) = em_clustering(&paths, &init);
        let mut error = 0;
        for (i, &pred) in result.iter().enumerate() {
            error += (pred != (2 * i / path_num) as u8) as u32;
            if pred != (2 * i / path_num) as u8 {
                eprintln!("{}\t{:?}", i, paths[i]);
            }
        }
        assert!((error as f64 / path_num as f64) < error_rate);
    }
    #[test]
    fn em_test_4() {
        let template1: Vec<usize> = vec![
            (0..=90).collect::<Vec<_>>(),
            (80..=90).collect::<Vec<_>>(),
            (80..=90).collect::<Vec<_>>(),
            (91..100).collect::<Vec<_>>(),
            vec![60],
            (100..120).collect::<Vec<_>>(),
        ]
        .concat();
        let template2: Vec<usize> = {
            let mut copied = template1.clone();
            let mut index = 0;
            copied.retain(|_| {
                index += 1;
                !(20..50).contains(&(index - 1))
            });
            copied
        };
        let mut cluster_num: HashMap<_, _> = HashMap::new();
        let template1: Vec<_> = template1
            .iter()
            .map(|&x| {
                let c = cluster_num.entry(x).or_default();
                *c += 1;
                (x, *c - 1)
            })
            .collect();
        let template2: Vec<_> = template2
            .iter()
            .map(|&x| {
                let c = cluster_num.entry(x).or_default();
                *c += 1;
                (x, *c - 1)
            })
            .collect();
        let max_len = 6;
        let min_len = 3;
        let error_rate = 0.05;
        let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(39420);
        let path_num = 1000;
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                let temp = if i < path_num / 2 {
                    &template1
                } else {
                    &template2
                };
                sim_path_error(&temp, &cluster_num, &mut rng, min_len, max_len, error_rate)
            })
            .collect();
        let init: Vec<_> = (0..path_num)
            .map(|i| {
                let is_error = rng.gen_bool(error_rate) as u8;
                let answer = (2 * i / path_num) as u8;
                answer ^ is_error
            })
            .collect();
        let error_prev = init
            .iter()
            .enumerate()
            .filter(|&(i, &x)| (2 * i / path_num) as u8 != x)
            .count() as f64
            / path_num as f64;
        let (result, _) = em_clustering(&paths, &init);
        let mut error = 0;
        for (i, &pred) in result.iter().enumerate() {
            error += (pred != (2 * i / path_num) as u8) as u32;
            if pred != (2 * i / path_num) as u8 {
                eprintln!("{}\t{:?}", i, paths[i]);
            }
        }
        let error = error as f64 / path_num as f64;
        assert!(error < error_rate, "{}->{}", error_prev, error);
    }
    #[test]
    fn em_test_5() {
        let template1: Vec<usize> = (0..10).collect();
        let mut cluster_num: HashMap<_, _> = HashMap::new();
        let template1: Vec<_> = template1
            .iter()
            .map(|&x| {
                let c = cluster_num.entry(x).or_default();
                *c += 2;
                (x, *c - 2)
            })
            .collect();
        let max_len = 6;
        let min_len = 3;
        let error_rate = 0.05;
        let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(39420);
        let path_num = 1000;
        let paths: Vec<_> = (0..path_num)
            .map(|_| {
                let temp = &template1;
                sim_path_error(&temp, &cluster_num, &mut rng, min_len, max_len, error_rate)
            })
            .collect();
        let init: Vec<_> = (0..path_num).map(|_| rng.gen_range(0..=1)).collect();
        let (result, _) = em_clustering(&paths, &init);
        let errors = result.iter().map(|&x| x as u32).sum::<u32>();
        let error = errors as f64 / path_num as f64;
        let error = (1f64 - error).min(error);
        for ((b, a), p) in init.iter().zip(result.iter()).zip(paths.iter()) {
            eprintln!("{}\t{}\t{:?}", b, a, p);
        }
        assert!(error < error_rate, "{}", error);
    }
}
