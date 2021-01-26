use log::debug;
use rand::Rng;
use std::collections::{HashMap, HashSet};
pub fn em_clustering(paths: &[Vec<(usize, usize)>], init_cluster: &[u8]) -> (Vec<u8>, f64) {
    let mut weights: Vec<_> = init_cluster
        .iter()
        .map(|&x| {
            if x == 0 {
                vec![1f64, 0f64]
            } else {
                vec![0f64, 1f64]
            }
        })
        .collect();
    let mut model = Model::create(&paths, &weights);
    debug!("Start EM algorithm.");
    let mut lk = paths.iter().map(|path| model.lk(path)).sum::<f64>();
    for i in 0.. {
        debug!("LK\t{}\t{:?}", i, lk);
        for (ws, path) in weights.iter_mut().zip(paths.iter()) {
            *ws = model.weights(&path);
        }
        model.update(&paths, &weights);
        let next_lk = paths.iter().map(|path| model.lk(path)).sum::<f64>();
        if next_lk - lk < 0.0001 {
            break;
        } else {
            lk = next_lk;
        }
    }
    let assignments: Vec<_> = weights
        .iter()
        .map(|ws| if ws[0] < ws[1] { 1 } else { 0 })
        .collect();
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
        // debug!("{}", likelihoods);
        // for node in model.nodes.iter() {
        //     debug!(
        //         "Model\t{}\t{:?}\t{:?}",
        //         node, model.hap1[node], model.hap2[node]
        //     );
        // }
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
    // debug!("{}", likelihoods);
}

#[derive(Debug, Clone)]
pub struct Model {
    nodes: HashSet<usize>,
    hap1: HashMap<usize, Vec<f64>>,
    hap2: HashMap<usize, Vec<f64>>,
}

impl Model {
    const SMALL: f64 = 0.0000001;
    pub fn create(paths: &[Vec<(usize, usize)>], weights: &[Vec<f64>]) -> Self {
        let mut cluster_num: HashMap<usize, usize> = HashMap::new();
        for path in paths.iter() {
            for &(n, c) in path.iter() {
                cluster_num
                    .entry(n)
                    .and_modify(|x| *x = (*x).max(c + 1))
                    .or_insert(c + 1);
            }
        }
        let mut hap1: HashMap<_, _> = cluster_num
            .iter()
            .map(|(&node, &cluster)| (node, vec![0f64; cluster]))
            .collect();
        let mut hap2 = hap1.clone();
        for (ws, path) in weights.iter().zip(paths.iter()) {
            for &(n, c) in path.iter() {
                hap1.get_mut(&n).unwrap()[c] += ws[0];
                hap2.get_mut(&n).unwrap()[c] += ws[1];
            }
        }
        hap1.values_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        hap2.values_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        let nodes: HashSet<usize> = cluster_num.keys().copied().collect();
        Self { hap1, hap2, nodes }
    }
    pub fn update(&mut self, paths: &[Vec<(usize, usize)>], weights: &[Vec<f64>]) {
        self.hap1
            .values_mut()
            .for_each(|xs| xs.iter_mut().for_each(|x| *x = 0f64));
        self.hap2
            .values_mut()
            .for_each(|xs| xs.iter_mut().for_each(|x| *x = 0f64));
        for (ws, path) in weights.iter().zip(paths.iter()) {
            for &(n, c) in path.iter() {
                self.hap1.get_mut(&n).unwrap()[c] += ws[0];
                self.hap2.get_mut(&n).unwrap()[c] += ws[1];
            }
        }
        self.hap1.values_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        self.hap2.values_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
    }
    pub fn create_on(
        units: &HashSet<usize>,
        genotypes: &HashMap<usize, usize>,
        weights: &[Vec<f64>],
        paths: &[(usize, Vec<(usize, usize)>)],
    ) -> Self {
        let mut hap1: HashMap<_, _> = units
            .iter()
            .map(|&node| (node, vec![0.; genotypes[&node]]))
            .collect();
        let mut hap2 = hap1.clone();
        for (weight, (_, path)) in weights.iter().zip(paths.iter()) {
            for &(n, c) in path.iter().filter(|&(n, _)| units.contains(n)) {
                if let Some(ws) = hap1.get_mut(&n) {
                    ws[c] += weight[0];
                }
                if let Some(ws) = hap2.get_mut(&n) {
                    ws[c] += weight[1];
                }
            }
        }
        hap1.values_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        hap2.values_mut().for_each(|ws| {
            let sum = ws.iter().sum::<f64>();
            if sum > Self::SMALL {
                ws.iter_mut().for_each(|w| *w /= sum);
            }
        });
        let nodes = units.clone();
        Self { hap1, hap2, nodes }
    }
    fn weights(&self, path: &[(usize, usize)]) -> Vec<f64> {
        let hap1_w = self.hap1_lk(path);
        let hap2_w = self.hap2_lk(path);
        let ln_sum = logsumexp(hap1_w, hap2_w);
        vec![(hap1_w - ln_sum).exp(), (hap2_w - ln_sum).exp()]
    }
    fn hap1_lk(&self, path: &[(usize, usize)]) -> f64 {
        path.iter()
            .map(|&(n, c)| match self.hap1.get(&n) {
                Some(ws) => ws[c].max(Self::SMALL).ln(),
                None => 0f64,
            })
            .sum::<f64>()
            + (0.5f64).ln()
    }
    fn hap2_lk(&self, path: &[(usize, usize)]) -> f64 {
        path.iter()
            .map(|&(n, c)| match self.hap2.get(&n) {
                Some(ws) => ws[c].max(Self::SMALL).ln(),
                None => 0.,
            })
            .sum::<f64>()
            + (0.5f64).ln()
    }
    fn lk(&self, path: &[(usize, usize)]) -> f64 {
        let hap1_w = self.hap1_lk(path);
        let hap2_w = self.hap2_lk(path);
        logsumexp(hap1_w, hap2_w)
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

// #[cfg(test)]
// mod test {
//     use super::*;
//     use rand::SeedableRng;
//     use rand_xoshiro::Xoshiro128StarStar;
//     #[test]
//     fn it_works() {
//         assert!(true);
//     }
//     #[test]
//     fn phase_test_1() {
//         let paths = vec![
//             (0, vec![(0, 0), (1, 0), (2, 0), (3, 0)]),
//             (1, vec![(0, 0), (1, 0), (2, 0), (3, 0)]),
//             (2, vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
//             (3, vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
//         ];
//         let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(238390);
//         let (result, lk) = em_progressive_clustering(&paths, &mut rng);
//         assert!(lk > -0.1, "{}", lk);
//         assert_eq!(result[0], result[1]);
//         assert_eq!(result[3], result[2]);
//         assert_ne!(result[0], result[3]);
//     }
//     #[test]
//     fn phase_test_2() {
//         // More complicated example.
//         let path1 = (0, vec![(0, 0), (1, 0), (10, 0), (3, 0), (10, 4)]);
//         let path2 = (1, vec![(1, 0), (10, 0), (3, 0), (10, 4), (110, 1)]);
//         let path3 = (2, vec![(0, 1), (1, 1), (10, 1), (3, 1), (10, 1), (110, 0)]);
//         let path4 = (3, vec![(10, 1), (3, 1), (10, 1)]);
//         let paths = vec![path1, path2, path3, path4];
//         let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(238390);
//         let (result, lk) = em_progressive_clustering(&paths, &mut rng);
//         assert!(lk > -0.1, "{}", lk);
//         assert_eq!(result[0], result[1]);
//         assert_eq!(result[2], result[3]);
//         assert_ne!(result[0], result[3]);
//     }
// }
