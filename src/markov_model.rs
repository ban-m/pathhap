#![allow(unused_variables, dead_code)]
//! A tiny Markov model to clustering paths.
//! For a fixed set of path P={p={(n_i,c_i)|i}}, a Markov model is G=(V,E), where V={(n,c)|n:node, c:cluster} and E={(n,c)->(n',c') | all edge in the paths}.
//! We encode a transition probability on each edge w(e) for all e in E.
//! The likelihood of a path p is computed by traversing the graph G along p, multiplicating the weight of each edge.
//! Note that for a path with length of 1, we set a "iniial probability" for a node (n,c), the fraction of the node in the entire nodeset.
//! We (hopefully) can derive a clustering algorithm for this model, like other mixture models, by first computing the
//! likelihood ratio of the model, then impute the parameters based on the obtained weights.
#[derive(Debug, Clone)]
pub struct MarkovModel {
    num_cluster: Vec<usize>,
    num_nodes: usize,
    // Node->Node->Cluster->Cluster
    // Note that this is a adjacency list in Node->Node
    // and adjacency matrix in Cluster->Cluster,
    // So, to obtain the weight (n,c)->(n',c'),
    // hap1[n].find(n')[c][c'] would suffice.
    hap1: Vec<Vec<Vec<Vec<f64>>>>,
    hap2: Vec<Vec<Vec<Vec<f64>>>>,
    // This is the pseudo count for each edge.
    // In other words, even if there is no edges between (n,c)->(n',c'),
    // we allow to jump from (n,c) to (n',c') with the probability of error_rate/D
    // where D is the degree of
    error_rate: f64,
}

/// Polish inital clustering `init_assignment` by Markov clustering.
/// Panic if `paths.len()!=init_assignment.len()`.
pub fn markov_clustering(
    paths: &[Vec<(usize, usize)>],
    init_assignment: &[u8],
    _error_rate: f64,
) -> Vec<u8> {
    // Weight for haplotype 1.
    let mut weights: Vec<_> = init_assignment
        .iter()
        .map(|&x| if x == 0 { 1. } else { 0. })
        .collect();
    let mut lk = std::f64::NEG_INFINITY;
    let mut model = MarkovModel::new(&paths, &weights);
    loop {
        weights
            .iter_mut()
            .zip(paths)
            .for_each(|(w, p)| *w = model.posterior_hap1(p));
        model.update(&paths, &weights);
        let next_lk = paths.iter().map(|p| model.lk(p)).sum::<f64>();
        if next_lk - lk < 0.00001 {
            break;
        }
        lk = next_lk;
    }
    weights
        .iter()
        .map(|&x| if 0.5 < x { 0 } else { 1 })
        .collect()
}

impl MarkovModel {
    pub fn new(paths: &[Vec<(usize, usize)>], weights: &[f64]) -> Self {
        unimplemented!()
    }
    pub fn update(&mut self, paths: &[Vec<(usize, usize)>], weights: &[f64]) {
        unimplemented!()
    }
    pub fn posterior_hap1(&self, path: &[(usize, usize)]) -> f64 {
        unimplemented!()
    }
    pub fn lk(&self, path: &[(usize, usize)]) -> f64 {
        unimplemented!()
    }
}
