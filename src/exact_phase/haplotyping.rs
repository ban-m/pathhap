use super::graph_traversal::*;
use super::model;
use log::{debug, warn};
use std::collections::{HashMap, HashSet};
const ERROR_FRACTION: f64 = 0.05;
pub fn haplotype_cc(paths: &[Vec<(usize, usize)>], max_len: usize, er: f64) -> (Vec<u8>, f64) {
    debug!("start.");
    let num_of_nodes = number_of_nodes(paths);
    let node_traverse_order = determine_traversal_order(num_of_nodes, paths);
    let split_paths = downsampling_up_to(&node_traverse_order, paths, max_len);
    let boundary_nodes: Vec<_> = get_boundary_nodes_on(&node_traverse_order, &split_paths);
    let beginning_paths = get_beginning_paths(&boundary_nodes, &split_paths);
    let prefetch_nodes = get_prefetch_nodes(&node_traverse_order, &boundary_nodes, &split_paths);
    let result = haplotyping_exact(
        &node_traverse_order,
        &prefetch_nodes,
        &beginning_paths,
        &split_paths,
        er,
    );
    let model = model::Model::new(&split_paths, &result, 0);
    let result: Vec<_> = paths.iter().map(|path| model.predict_path(path)).collect();
    let lk = {
        let paths: Vec<(usize, &[(usize, usize)])> =
            paths.iter().map(|x| x.as_slice()).enumerate().collect();
        let result: Vec<_> = result.iter().copied().enumerate().collect();
        let model = model::Model::new(&paths, &result, 1);
        paths.iter().map(|(_, x)| model.likelihood(x)).sum::<f64>()
    };
    (result, lk)
}

fn get_cluster_num(paths: &[(usize, &[(usize, usize)])]) -> HashMap<usize, usize> {
    let mut cluster_num: HashMap<usize, usize> = HashMap::new();
    for (_, p) in paths.iter() {
        for &(n, c) in p.iter() {
            match cluster_num.get_mut(&n) {
                Some(res) => *res = (*res).max(c + 1),
                None => {
                    cluster_num.insert(n, c + 1);
                }
            }
        }
    }
    cluster_num
}

fn haplotyping_exact(
    order: &[usize],
    prefetch_nodes: &[SortedNodes],
    beginning_paths: &[Paths],
    paths: &[(usize, &[(usize, usize)])],
    error_rate: f64,
) -> Vec<(usize, u8)> {
    debug!("PATH\tID\tPath");
    for (id, path) in paths.iter() {
        debug!("PATH\t{}\t{:?}", id, path);
    }
    // Calculate the number of cluster.
    let cluster_num = get_cluster_num(paths);
    debug!("DUMP\tITER\tNODE\tPrefetches\tLen\tBeginPath");
    for (i, ((node, bn), begpath)) in order
        .iter()
        .zip(prefetch_nodes.iter())
        .zip(beginning_paths.iter())
        .enumerate()
    {
        let len = bn.iter().map(|n| cluster_num[n]).sum::<usize>();
        debug!("\t{}\t{}\t{:?}\t{}\t{:?}", i, node, bn, len, begpath);
    }

    let paths: HashMap<usize, _> = paths.iter().copied().collect();
    // Haplotyping.
    // Return the haplotype information giving maximum value for the next iteration.
    let mut argmax: Vec<Vec<Haplotype>> = vec![];
    let mut prev_score: Vec<f64> = Vec::new();
    for (t, (prefetches, begin_paths)) in prefetch_nodes
        .iter()
        .zip(beginning_paths.iter())
        .enumerate()
    {
        debug!("Start\t{}\t{:?}\t{:?}", t, prefetches, begin_paths);
        for p in begin_paths.iter() {
            debug!("{:?}", paths[p]);
        }
        let hp = Haplotyping::new(prefetches, &cluster_num);
        if t == 0 {
            // Need not to look previous haplotyping.
            let next_converter = hp.calc_converter(&prefetch_nodes[t + 1]);
            prev_score = vec![std::f64::NEG_INFINITY; next_converter.num_pattern()];
            argmax.push(vec![Haplotype(0); next_converter.num_pattern()]);
            for haplotype in hp.clone() {
                let current_score = begin_paths
                    .iter()
                    .map(|id| hp.score(haplotype, &paths[id]))
                    .sum::<f64>();
                let next_pos = next_converter.convert(haplotype).pos();
                debug!("HAP\t{:b}\t{:b}\t{}", haplotype.0, next_pos, current_score);
                if prev_score[next_pos] < current_score {
                    prev_score[next_pos] = current_score;
                    argmax[t][next_pos] = haplotype;
                }
            }
        } else if t + 1 == order.len() {
            // Need not to look next haplotyping.
            let prev_converter = hp.calc_converter(&prefetch_nodes[t]);
            let mut score = vec![std::f64::NEG_INFINITY; hp.num_pattern()];
            argmax.push(vec![Haplotype(0); hp.num_pattern()]);
            for haplotype in hp.clone() {
                let score_on_node = begin_paths
                    .iter()
                    .map(|id| hp.score(haplotype, &paths[id]))
                    .sum::<f64>();
                let prev_position = prev_converter.convert(haplotype).pos();
                let current_score = prev_score[prev_position] + score_on_node;
                let pos = haplotype.pos();
                debug!(
                    "HAP\t{:b}\t{:b}\t{:.3}\t{:.3}",
                    haplotype.0, pos, score_on_node, current_score
                );
                score[pos] = current_score;
                argmax[t][pos] = haplotype;
            }
            prev_score = score;
        } else {
            let next_converter = hp.calc_converter(&prefetch_nodes[t + 1]);
            let prev_converter = hp.calc_converter(&prefetch_nodes[t - 1]);
            let mut score = vec![std::f64::NEG_INFINITY; next_converter.num_pattern()];
            argmax.push(vec![Haplotype(0); next_converter.num_pattern()]);
            for haplotype in hp.clone() {
                let score_on_node = begin_paths
                    .iter()
                    .map(|id| hp.score(haplotype, &paths[id]))
                    .sum::<f64>();
                let prev_position = prev_converter.convert(haplotype).pos();
                let current_score = prev_score[prev_position] + score_on_node;
                let next_pos = next_converter.convert(haplotype).pos();
                debug!(
                    "HAP\t{:b}\t{:b}\t{:b}\t{:.3}\t{:.3}",
                    haplotype.0, prev_position, next_pos, score_on_node, current_score
                );
                if score[next_pos] < current_score {
                    score[next_pos] = current_score;
                    argmax[t][next_pos] = haplotype;
                }
            }
            prev_score = score;
        };
    }
    // Traceback.
    let (lk, mut current_argmax_hap): (f64, Haplotype) = prev_score
        .iter()
        .zip(argmax.pop().unwrap())
        .max_by(|x, y| (x.0).partial_cmp(&y.0).unwrap())
        .map(|(&x, y)| (x, y))
        .unwrap();
    debug!("MAX LK IS {:.3}", lk);
    let mut result = vec![(prefetch_nodes.last().unwrap(), current_argmax_hap)];
    // The last node is removed from argmax already.
    for (prefetches, argmax) in prefetch_nodes.windows(2).zip(argmax.iter()).rev() {
        let hp = Haplotyping::new(&prefetches[1], &cluster_num);
        let converter = hp.calc_converter(&prefetches[0]);
        let pos = converter.convert(current_argmax_hap).pos();
        current_argmax_hap = argmax[pos];
        result.push((&prefetches[0], current_argmax_hap));
    }
    for (nodes, hap) in result.iter() {
        debug!("HAPLOTYPE\t{:?}\t{:b}", nodes, hap.0);
    }
    let model = HapModel::new(&result, &cluster_num);
    paths
        .iter()
        .map(|(&id, path)| (id, model.predict(path, error_rate)))
        .collect()
}

#[derive(Debug, Clone)]
pub struct HapModel {
    // If (n,c) is in hap1, `is_in_hap1.contains(&(n,c)) == true`
    is_in_hap1: HashMap<usize, Vec<bool>>,
}
impl HapModel {
    fn new(result: &[(&SortedNodes, Haplotype)], cluster_num: &HashMap<usize, usize>) -> Self {
        let mut is_in_hap1: HashMap<usize, Vec<bool>> = HashMap::new();
        for &(nodes, hap) in result {
            for (node, clusters) in Haplotyping::recover(hap, nodes, cluster_num) {
                match is_in_hap1.get(&node) {
                    Some(res) => assert_eq!(res, &clusters),
                    None => {
                        is_in_hap1.insert(node, clusters);
                    }
                }
            }
        }
        for (node, &hap) in cluster_num.iter() {
            if !is_in_hap1.contains_key(node) {
                is_in_hap1.insert(*node, vec![false; hap]);
            }
        }
        Self { is_in_hap1 }
    }
    fn predict(&self, path: &[(usize, usize)], error_rate: f64) -> u8 {
        let (mut hap1, mut hap2) = (0f64, 0f64);
        for &(ref node, cluster) in path {
            if let Some(clusters) = self.is_in_hap1.get(node) {
                let hap1_size = clusters.iter().filter(|&&x| x).count();
                let hap2_size = clusters.len() - hap1_size;
                if clusters[cluster] {
                    let hap1score = (1f64 - error_rate * hap2_size as f64) / hap1_size as f64;
                    let hap2score = error_rate;
                    hap1 += hap1score.ln();
                    hap2 += hap2score.ln();
                } else {
                    let hap1score = error_rate;
                    let hap2score = (1f64 - error_rate * hap1_size as f64) / hap2_size as f64;
                    hap1 += hap1score.ln();
                    hap2 += hap2score.ln();
                };
            } else {
                // This node never happened in the data.
                // It should not be, but anyway we leave it alone,
                // just by indicating that this node never appreared in the data.
                warn!("{}-{} is not recorded in the data.", node, cluster);
            }
        }
        if hap2 < hap1 {
            0
        } else {
            1
        }
    }
}

#[derive(Clone, Debug)]
pub struct Haplotyping<'a, 'b> {
    nodes: &'a SortedNodes,
    cluster_num: &'b HashMap<usize, usize>,
    offsets: Vec<usize>,
}

impl<'a, 'b> Haplotyping<'a, 'b> {
    fn recover(
        hap: Haplotype,
        nodes: &SortedNodes,
        cluster_num: &HashMap<usize, usize>,
    ) -> Vec<(usize, Vec<bool>)> {
        let hap = hap.0;
        let mut recovered_hap = vec![];
        let mut current_position = 0;
        for &node in nodes.iter() {
            let cl = cluster_num[&node];
            let hap_on_node: Vec<bool> = (0..cl)
                .map(|i| hap >> (i + current_position) & 0b1 == 0b1)
                .collect();
            recovered_hap.push((node, hap_on_node));
            current_position += cl;
        }
        recovered_hap
    }
    fn new(nodes: &'a SortedNodes, cluster_num: &'b HashMap<usize, usize>) -> Self {
        let mut offsets = vec![];
        let mut current_position = 0;
        for &node in nodes.iter() {
            offsets.push(current_position);
            current_position += cluster_num[&node];
        }
        Self {
            nodes,
            cluster_num,
            offsets,
        }
    }
    fn calc_converter(&self, nodes: &SortedNodes) -> Converter {
        let len = self
            .nodes
            .iter()
            .map(|n| self.cluster_num[n])
            .sum::<usize>();
        let mut pattern = 0;
        // TODO: There should be more efficient algorithm by marge sorting.
        let mut current_position = 0;
        for node in self.nodes.iter() {
            let cl = self.cluster_num[node];
            if nodes.contains(node) {
                // Fill the current_position-th bit to (current_position+cluster_num)-th bit.
                for i in current_position..(current_position + cl) {
                    pattern |= 1 << i;
                }
            }
            current_position += cl;
        }
        Converter { len, pattern }
    }
    fn num_pattern(&self) -> usize {
        let total_len = self
            .nodes
            .iter()
            .map(|n| self.cluster_num[n])
            .sum::<usize>();
        1 << total_len
    }
    fn score_with_error(&self, Haplotype(hap): Haplotype, path: &[(usize, usize)], er: f64) -> f64 {
        // Haplotypes. Starts with (1/2, 1/2).ln();
        // let (mut hap1, mut hap2) = (-2f64.ln(), -2f64.ln());
        let (mut hap1, mut hap2) = (0f64, 0f64);
        for (node, cluster) in path.iter() {
            let idx = match self.nodes.binary_search(node) {
                Ok(res) => res,
                Err(_) => panic!("{:?} not contains {}", self.nodes, node),
            };
            let bit_from = self.offsets[idx];
            let bit_size = self.cluster_num[node];
            let mask = (1 << bit_size) - 1;
            let cl_in_1 = ((hap >> bit_from) & mask).count_ones();
            let cl_in_2 = bit_size as u32 - cl_in_1;
            if (hap >> (bit_from + cluster) & 0b1) == 0b1 {
                hap1 += ((1f64 - cl_in_2 as f64 * er) / cl_in_1 as f64).ln();
                hap2 += er.ln();
            } else {
                hap1 += er.ln();
                hap2 += ((1f64 - cl_in_1 as f64 * er) / cl_in_2 as f64).ln();
            }
        }
        //hap1 + hap2
        assert!(!hap1.is_nan() && !hap2.is_nan());
        if hap1 < hap2 {
            hap2
        } else {
            hap1
        }
    }
    fn score(&self, hap: Haplotype, path: &[(usize, usize)]) -> f64 {
        self.score_with_error(hap, path, ERROR_FRACTION)
    }
}

impl<'a, 'b> std::iter::IntoIterator for Haplotyping<'a, 'b> {
    type Item = Haplotype;
    type IntoIter = HaplotypeEnumerator<'a, 'b>;
    fn into_iter(self) -> Self::IntoIter {
        HaplotypeEnumerator::new(self)
    }
}

#[derive(Debug, Clone)]
pub struct HaplotypeEnumerator<'a, 'b> {
    haplotyping: Haplotyping<'a, 'b>,
    pattern_num: usize,
    counter: usize,
}

impl<'a, 'b> HaplotypeEnumerator<'a, 'b> {
    pub fn new(hp: Haplotyping<'a, 'b>) -> Self {
        Self {
            pattern_num: hp.num_pattern(),
            haplotyping: hp,
            counter: 0,
        }
    }
}

impl<'a, 'b> std::iter::Iterator for HaplotypeEnumerator<'a, 'b> {
    type Item = Haplotype;
    fn next(&mut self) -> Option<Self::Item> {
        if self.counter < self.pattern_num {
            let pat = self.counter as u64;
            self.counter += 1;
            Some(Haplotype(pat))
        } else {
            None
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Haplotype(u64);

impl std::fmt::Debug for Haplotype {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "H({:b})", self.0)
    }
}

impl Haplotype {
    fn pos(&self) -> usize {
        self.0 as usize
    }
}

// Convert from a haplotype into another.
#[derive(Debug, Clone, Copy)]
pub struct Converter {
    // If the i-th bit is true, then
    // the i-th bit should be in the converted haplotyping.
    pattern: u64,
    // The length of the haplotype.
    len: usize,
}

impl Converter {
    fn num_pattern(&self) -> usize {
        let total_size = self.pattern.count_ones() as usize;
        1 << total_size
    }
    fn convert(&self, h: Haplotype) -> Haplotype {
        let hap = h.0;
        let mut result = 0;
        let mut slot = 0;
        for i in 0..self.len {
            if ((self.pattern >> i) & 0b1) == 0b1 {
                result |= ((hap >> i) & 0b1) << slot;
                slot += 1;
            }
        }
        Haplotype(result)
    }
}

fn get_prefetch_nodes(
    order: &[usize],
    boundary_nodes: &[Nodes],
    split_paths: &[(usize, &[(usize, usize)])],
) -> Vec<SortedNodes> {
    let beginning_paths = get_beginning_paths(&boundary_nodes, &split_paths);
    let paths: HashMap<_, _> = split_paths.iter().copied().collect();
    let mut haplotyped_nodes: HashSet<usize> = HashSet::new();
    let mut current_focus_nodes: HashSet<usize> = HashSet::new();
    let mut prefetched_nodes = vec![];
    for (&focal_node, begin_paths) in order.iter().zip(beginning_paths) {
        let prefetch = get_prefetch_nodes_on(&begin_paths, &paths, &haplotyped_nodes);
        current_focus_nodes.extend(prefetch);
        let prefetch = SortedNodes::new(current_focus_nodes.iter().copied().collect::<Vec<_>>());
        prefetched_nodes.push(prefetch);
        haplotyped_nodes.insert(focal_node);
        current_focus_nodes.remove(&focal_node);
    }
    prefetched_nodes
}

fn get_prefetch_nodes_on(
    path_ids: &Paths,
    paths: &HashMap<usize, &[(usize, usize)]>,
    already_phased_nodes: &HashSet<usize>,
) -> HashSet<usize> {
    path_ids
        .iter()
        .flat_map(|id| {
            paths[id]
                .iter()
                .map(|x| x.0)
                .filter(|x| !already_phased_nodes.contains(x))
        })
        .collect()
}

fn downsampling_up_to<'a>(
    order: &[usize],
    paths: &'a [Vec<(usize, usize)>],
    max_len: usize,
) -> Vec<(usize, &'a [(usize, usize)])> {
    let paths: Vec<_> = paths.iter().map(|x| x.as_slice()).enumerate().collect();
    let boundary_nodes = get_boundary_nodes_on(&order, &paths);
    let beginning_paths = get_beginning_paths(&boundary_nodes, &paths);
    let mut paths: HashMap<_, _> = paths.into_iter().collect();
    let mut haplotyped_nodes: HashSet<usize> = HashSet::new();
    let mut current_focus_nodes: HashSet<usize> = HashSet::new();
    for (&focal_node, begin_paths) in order.iter().zip(beginning_paths) {
        loop {
            if begin_paths.is_empty() {
                break;
            }
            let new_nodes: HashSet<_> =
                get_prefetch_nodes_on(&begin_paths, &paths, &haplotyped_nodes);
            let current_haplotyping_size = new_nodes.union(&current_focus_nodes).count();
            if current_haplotyping_size > max_len {
                // Split the longest reads in begin_paths into half.
                let max_path = *begin_paths.iter().max_by_key(|id| paths[id].len()).unwrap();
                let longest_path = paths.remove(&max_path).unwrap();
                // Maybe we can use the remaining path somehow.
                let split_path = split_path_at(longest_path, focal_node);
                paths.insert(max_path, split_path);
            } else {
                break;
            }
        }
        // Register all the nodes into current_focus.
        let new_nodes = get_prefetch_nodes_on(&begin_paths, &paths, &haplotyped_nodes);
        current_focus_nodes.extend(new_nodes);
        current_focus_nodes.remove(&focal_node);
        haplotyped_nodes.insert(focal_node);
    }
    // Sanity check.
    let order: HashSet<_> = order.iter().copied().collect();
    assert_eq!(order, haplotyped_nodes);
    let mut paths: Vec<(usize, &[_])> = paths.into_iter().collect();
    paths.sort();
    paths
}

fn split_path_at(path: &[(usize, usize)], node: usize) -> &[(usize, usize)] {
    let first_half = &path[..path.len() / 2];
    if first_half.iter().any(|&(n, _)| n == node) {
        first_half
    } else {
        let second_half = &path[path.len() / 2..];
        assert!(second_half.iter().any(|&(n, _)| n == node));
        second_half
    }
}

fn get_beginning_paths(
    boundary_nodes: &[Nodes],
    paths: &[(usize, &[(usize, usize)])],
) -> Vec<Paths> {
    let mut paths: Vec<_> = paths.to_vec();
    let mut beginning_paths = vec![];
    for bn in boundary_nodes {
        let mut ids: Vec<usize> = vec![];
        paths.retain(|&(id, path)| {
            if path.iter().any(|&(n, _)| bn.contains(&n)) {
                ids.push(id);
                false
            } else {
                true
            }
        });
        beginning_paths.push(Paths::new(ids));
    }
    beginning_paths
}

#[derive(Debug, Clone)]
pub struct SortedNodes {
    nodes: Vec<usize>,
}

impl SortedNodes {
    fn contains(&self, node: &usize) -> bool {
        self.nodes.binary_search(node).is_ok()
    }
    fn new(mut nodes: Vec<usize>) -> Self {
        nodes.sort_unstable();
        Self { nodes }
    }
    fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.nodes.iter()
    }
    fn binary_search(&self, node: &usize) -> Result<usize, usize> {
        self.nodes.binary_search(node)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Paths {
    // ID. Not index.
    paths: Vec<usize>,
}

impl std::iter::IntoIterator for Paths {
    type Item = usize;
    type IntoIter = std::vec::IntoIter<Self::Item>;
    fn into_iter(self) -> Self::IntoIter {
        self.paths.into_iter()
    }
}

impl Paths {
    pub fn is_empty(&self) -> bool {
        self.paths.is_empty()
    }
    pub fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.paths.iter()
    }
    pub fn new(paths: Vec<usize>) -> Self {
        Self { paths }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro256PlusPlus;
    #[test]
    fn it_works() {}
    #[test]
    fn get_boundary_nodes_test() {
        let order = vec![0, 1, 2, 3, 4, 7, 6, 5, 8, 10, 9, 11];
        let path: Vec<_> = vec![
            vec![0, 1, 2, 3],
            vec![2, 3, 4],
            vec![8, 7, 4, 3],
            vec![8, 10, 11],
            vec![11, 10, 9, 8],
            vec![7, 6, 5],
            vec![5, 1],
            vec![6, 5],
            vec![0],
        ];
        let path: Vec<Vec<_>> = path
            .into_iter()
            .map(|p| p.into_iter().map(|x| (x, 0)).collect())
            .collect();
        let path: Vec<_> = path.iter().map(|x| x.as_slice()).enumerate().collect();
        let boundary_nodes = get_boundary_nodes_on(&order, &path);
        assert_eq!(boundary_nodes.len(), order.len());
        let mut boundary: Vec<Vec<usize>> = boundary_nodes
            .iter()
            .map(|x| x.nodes.iter().copied().collect())
            .collect();
        boundary.iter_mut().for_each(|x| x.sort());
        assert_eq!(boundary[0], vec![0]);
        assert_eq!(boundary[1], vec![1]);
        assert_eq!(boundary[2], vec![1, 2]);
        assert_eq!(boundary[3], vec![1, 3]);
        assert_eq!(boundary[4], vec![1, 4]);
        assert_eq!(boundary[5], vec![1, 7]);
        assert_eq!(boundary[6], vec![1, 6, 7]);
        assert_eq!(boundary[7], vec![5, 7]);
        assert_eq!(boundary[8], vec![8]);
        assert_eq!(boundary[9], vec![8, 10]);
        assert_eq!(boundary[10], vec![9, 10]);
        assert_eq!(boundary[11], vec![11]);
    }
    #[test]
    fn get_boundary_nodes_test_2() {
        let order = vec![0, 2, 3, 5, 4, 6, 1];
        let paths = vec![
            vec![3, 2, 4],
            vec![0, 2, 4, 5, 6],
            vec![1, 2, 3, 5, 6],
            vec![6],
            vec![0, 2],
            vec![3, 5],
        ];
        let paths: Vec<Vec<_>> = paths
            .into_iter()
            .map(|p| p.into_iter().map(|x| (x, 0)).collect::<Vec<_>>())
            .collect();
        let paths: Vec<(usize, &[_])> = paths.iter().map(|x| x.as_slice()).enumerate().collect();
        let boundary_nodes = get_boundary_nodes_on(&order, &paths);
        assert_eq!(boundary_nodes.len(), order.len());
        let boundary: Vec<Vec<_>> = boundary_nodes
            .into_iter()
            .map(|x| {
                let mut xs: Vec<_> = x.nodes.iter().copied().collect();
                xs.sort();
                xs
            })
            .collect();
        assert_eq!(boundary[0], vec![0]);
        assert_eq!(boundary[1], vec![2]);
        assert_eq!(boundary[2], vec![2, 3]);
        assert_eq!(boundary[3], vec![2, 5]);
        assert_eq!(boundary[4], vec![2, 4, 5]);
        assert_eq!(boundary[5], vec![2, 6]);
        assert_eq!(boundary[6], vec![1]);
    }
    #[test]
    fn get_beginning_paths_test() {
        let order = vec![0, 1, 2, 3, 4, 7, 6, 5, 8, 10, 11, 9];
        let path: Vec<_> = vec![
            vec![0, 1, 2, 3],
            vec![2, 3, 4],
            vec![8, 7, 4, 3],
            vec![8, 10, 11],
            vec![11, 10, 9, 8],
            vec![7, 6, 5],
            vec![5, 1],
            vec![6, 5],
            vec![0],
        ];
        let path: Vec<Vec<_>> = path
            .into_iter()
            .map(|p| p.into_iter().map(|x| (x, 0)).collect())
            .collect();
        let path: Vec<_> = path.iter().map(|x| x.as_slice()).enumerate().collect();
        let boundary_nodes: Vec<_> = get_boundary_nodes_on(&order, &path);
        let result = get_beginning_paths(&boundary_nodes, &path);
        let answer: Vec<_> = vec![
            vec![0, 8],
            vec![6],
            vec![1],
            vec![2],
            vec![],
            vec![5],
            vec![7],
            vec![],
            vec![3, 4],
            vec![],
            vec![],
        ]
        .into_iter()
        .map(Paths::new)
        .collect();
        for (i, (r, a)) in result.iter().zip(answer.iter()).enumerate() {
            assert_eq!(r, a, "{}", i);
        }
    }
    #[test]
    fn get_beggining_paths_test_2() {
        let order = vec![0, 2, 3, 5, 4, 6, 1];
        let paths = vec![
            vec![3, 2, 4],
            vec![0, 2, 4, 5, 6],
            vec![1, 2, 3, 5, 6],
            vec![6],
            vec![0, 2],
            vec![3, 5],
        ];
        let paths: Vec<Vec<_>> = paths
            .into_iter()
            .map(|p| p.into_iter().map(|x| (x, 0)).collect())
            .collect::<Vec<_>>();
        let paths: Vec<_> = paths.iter().map(|x| x.as_slice()).enumerate().collect();
        let boundary_nodes: Vec<_> = get_boundary_nodes_on(&order, &paths);
        let result = get_beginning_paths(&boundary_nodes, &paths);
        let answer: Vec<_> = vec![
            vec![1, 4],
            vec![0, 2],
            vec![5],
            vec![],
            vec![],
            vec![3],
            vec![],
        ]
        .into_iter()
        .map(Paths::new)
        .collect();
        for (i, (r, a)) in result.iter().zip(answer.iter()).enumerate() {
            assert_eq!(r, a, "{}", i);
        }
    }
    #[test]
    fn hap_recover_test() {
        let nodes = SortedNodes::new(vec![0, 1, 3]);
        let cluster_num: HashMap<usize, usize> =
            vec![(0, 2), (1, 1), (2, 3), (3, 3)].into_iter().collect();
        let hap = Haplotype(0b000_0_00);
        let recovered = Haplotyping::recover(hap, &nodes, &cluster_num);
        let answer = vec![
            (0, vec![false, false]),
            (1, vec![false]),
            (3, vec![false, false, false]),
        ];
        assert_eq!(recovered, answer);
        let hap = Haplotype(0b011_0_01);
        let recovered = Haplotyping::recover(hap, &nodes, &cluster_num);
        let answer = vec![
            (0, vec![true, false]),
            (1, vec![false]),
            (3, vec![true, true, false]),
        ];
        assert_eq!(recovered, answer);
    }
    #[test]
    fn hap_new_test() {
        let nodes = SortedNodes::new(vec![0, 1, 3]);
        let cluster_num: HashMap<usize, usize> =
            vec![(0, 2), (1, 1), (2, 3), (3, 3)].into_iter().collect();
        let nodes2 = SortedNodes::new(vec![0, 1, 2]);
        // in total 2+1+3 = 6 choises.
        let hp = Haplotyping::new(&nodes, &cluster_num);
        assert_eq!(hp.num_pattern(), 1 << 6);
        let hp = Haplotyping::new(&nodes2, &cluster_num);
        assert_eq!(hp.num_pattern(), 1 << 6);
    }
    #[test]
    fn hap_converter_test() {
        let nodes = SortedNodes::new(vec![0, 1, 3]);
        let cluster_num: HashMap<usize, usize> =
            vec![(0, 2), (1, 1), (2, 3), (3, 3)].into_iter().collect();
        let nodes2 = SortedNodes::new(vec![0, 1, 2]);
        let hp = Haplotyping::new(&nodes, &cluster_num);
        let converter = hp.calc_converter(&nodes2);
        let pat = converter.num_pattern();
        // 2 cluster for 0, 1 cluster for 1, thus in total 3 choises, 2^3 = 8 pattern.
        assert_eq!(pat, 8);
        let hap = Haplotype(0b000_0_00);
        assert_eq!(converter.convert(hap), Haplotype(0b0_00));
        let hap = Haplotype(0b111_0_00);
        assert_eq!(converter.convert(hap), Haplotype(0b0_00));
        let hap = Haplotype(0b000_1_00);
        assert_eq!(converter.convert(hap), Haplotype(0b1_00));
        let hap = Haplotype(0b000_0_01);
        assert_eq!(hap, Haplotype(0b_01));
    }
    #[test]
    fn hap_converter_test_2() {
        let nodes = SortedNodes::new(vec![0, 1, 2, 3, 4]);
        let cluster_num: HashMap<usize, usize> = vec![(0, 2), (1, 2), (2, 3), (3, 2), (4, 2)]
            .into_iter()
            .collect();
        let nodes2 = SortedNodes::new(vec![1, 2, 3, 4]);
        let hp = Haplotyping::new(&nodes, &cluster_num);
        let converter = hp.calc_converter(&nodes2);
        assert_eq!(converter.len, 11);
        let ans_pat = 0b11_11_111_11_00;
        assert_eq!(
            converter.pattern, ans_pat,
            "{:b}-{:b}",
            converter.pattern, ans_pat
        );
        let pat = converter.num_pattern();
        assert_eq!(pat, 1 << 9);
        let hap = Haplotype(0b10_01_101_01_01);
        assert_eq!(converter.convert(hap), Haplotype(0b10_01_101_01));
    }
    #[test]
    fn hap_score_test() {
        let nodes = SortedNodes::new(vec![0, 1, 3]);
        let cluster_num: HashMap<usize, usize> =
            vec![(0, 2), (1, 1), (2, 3), (3, 3)].into_iter().collect();
        let hp = Haplotyping::new(&nodes, &cluster_num);
        let path = &[(0, 1), (1, 0), (3, 0)];
        let er = 0.00001;
        let hap = Haplotype(0b111_1_11);
        let score = hp.score_with_error(hap, path, er);
        let opt = (1f64 / 3f64).ln() + 1f64.ln() + (1f64 / 2f64).ln();
        assert!((opt - score).abs() < 0.0001, "{},{}", score, opt);
        let hap = Haplotype(0b000_0_00);
        let score = hp.score_with_error(hap, path, er);
        assert!((opt - score).abs() < 0.0001, "{},{}", score, opt);
        let hap = Haplotype(0b001_1_10);
        let score = hp.score_with_error(hap, path, er);
        assert!(score.abs() < 0.001, "{}", score);
        let hap = Haplotype(0b110_0_01);
        let score = hp.score_with_error(hap, path, er);
        assert!(score.abs() < 0.001, "{}", score);
    }
    #[test]
    fn hap_test_1() {
        let paths = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1)],
        ];
        let (result, _) = haplotype_cc(&paths, 14, ERROR_FRACTION);
        eprintln!("{:?}", result);
        let is_ok = result == vec![0, 0, 1, 1] || result == vec![1, 1, 0, 0];
        assert!(is_ok, "{:?}", result);
    }
    #[test]
    fn hap_test_2() {
        // More complicated example.
        let path1 = vec![(0, 0), (1, 0), (2, 0), (3, 0), (2, 2)];
        let path2 = vec![(1, 0), (2, 0), (3, 0), (2, 2), (4, 1)];
        let path3 = vec![(0, 1), (1, 1), (2, 1), (3, 1), (2, 1), (4, 0)];
        let path4 = vec![(2, 1), (3, 1), (2, 1)];
        let paths = vec![path1, path2, path3, path4];
        let (result, _) = haplotype_cc(&paths, 14, ERROR_FRACTION);
        let is_ok = vec![0, 0, 1, 1] == result || vec![1, 1, 0, 0] == result;
        assert!(is_ok, "{:?}", result);
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
    fn hap_test_random_linear() {
        let template1: Vec<_> = (0..10).map(|x| (x, 0)).collect();
        let template2: Vec<_> = (0..10).map(|x| (x, 1)).collect();
        let min_len = 3;
        let max_len = 6;
        let path_num = 10;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let paths: Vec<Vec<(usize, usize)>> = (0..path_num)
            .map(|i| {
                if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                }
            })
            .collect();
        // paths.shuffle(&mut rng);
        let (result, _) = haplotype_cc(&paths, 20, ERROR_FRACTION);
        let answer: Vec<_> = (0..path_num)
            .map(|i| match i < path_num / 2 {
                true => 0,
                false => 1,
            })
            .collect();
        let answer2: Vec<_> = answer
            .iter()
            .map(|&i| match i {
                0 => 1,
                1 => 0,
                _ => panic!(),
            })
            .collect();
        let is_ok = result == answer || result == answer2;
        eprintln!("{:?}, {:?}", answer, answer2);
        assert!(is_ok, "{:?}", result);
    }
    #[test]
    fn hap_test_random_loop() {
        let mut count: HashMap<_, usize> = HashMap::new();
        let template1: Vec<(usize, usize)> = {
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
        let template2: Vec<(usize, usize)> = {
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
        // eprintln!("{:?}", template2);
        let min_len = 3;
        let max_len = 6;
        let path_num = 10;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                }
            })
            .collect();
        // paths.shuffle(&mut rng);
        let (result, _) = haplotype_cc(&paths, 20, ERROR_FRACTION);
        let answer: Vec<_> = (0..path_num)
            .map(|i| match i < path_num / 2 {
                true => 0,
                false => 1,
            })
            .collect();
        let answer2: Vec<_> = answer
            .iter()
            .map(|&i| match i {
                0 => 1,
                1 => 0,
                _ => panic!(),
            })
            .collect();
        let is_ok = result == answer || result == answer2;
        assert!(is_ok, "{:?}", result);
    }
    #[test]
    fn hap_test_random_branch() {
        let mut count: HashMap<_, usize> = HashMap::new();
        let template1: Vec<(usize, usize)> = vec![0, 1, 2, 3, 4, 5, 8, 9, 10]
            .into_iter()
            .map(|n| {
                let count = count.entry(n).or_default();
                *count += 1;
                (n, *count - 1)
            })
            .collect();
        // eprintln!("{:?}", template1);
        let template2: Vec<(usize, usize)> = vec![0, 1, 2, 3, 6, 7, 8, 9, 10]
            .into_iter()
            .map(|n| {
                let count = count.entry(n).or_default();
                *count += 1;
                (n, *count - 1)
            })
            .collect();
        // eprintln!("{:?}", template2);
        let min_len = 3;
        let max_len = 6;
        let path_num = 10;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                }
            })
            .collect();
        let (result, _) = haplotype_cc(&paths, 20, ERROR_FRACTION);
        let answer: Vec<_> = (0..path_num)
            .map(|i| match i < path_num / 2 {
                true => 0,
                false => 1,
            })
            .collect();
        let answer2: Vec<_> = answer
            .iter()
            .map(|&i| match i {
                0 => 1,
                1 => 0,
                _ => panic!(),
            })
            .collect();
        let is_ok = result == answer || result == answer2;
        assert!(is_ok, "{:?}", result);
    }
    #[test]
    fn hap_test_random_loop_branch() {
        let mut count: HashMap<_, usize> = HashMap::new();
        let template1: Vec<(usize, usize)> = {
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
        let template2: Vec<(usize, usize)> = {
            let nodes = vec![0, 1, 2, 5, 4, 3, 2, 6, 9, 10, 11, 12, 13, 11, 12, 14];
            nodes
                .iter()
                .map(|&node| {
                    let count = count.entry(node).or_default();
                    *count += 1;
                    (node, *count - 1)
                })
                .collect()
        };
        let min_len = 3;
        let max_len = 6;
        let path_num = 10;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(34089);
        let paths: Vec<_> = (0..path_num)
            .map(|i| {
                if i < path_num / 2 {
                    sim_path(&template1, &mut rng, min_len, max_len)
                } else {
                    sim_path(&template2, &mut rng, min_len, max_len)
                }
            })
            .collect();
        let (result, _) = haplotype_cc(&paths, 20, ERROR_FRACTION);
        let answer: Vec<_> = (0..path_num)
            .map(|i| match i < path_num / 2 {
                true => 0,
                false => 1,
            })
            .collect();
        let answer2: Vec<_> = answer
            .iter()
            .map(|&i| match i {
                0 => 1,
                1 => 0,
                _ => panic!(),
            })
            .collect();
        let is_ok = result == answer || result == answer2;
        assert!(is_ok, "{:?}", result);
    }
}
