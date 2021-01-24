//! phasing module by GraphWhatsHap.
use super::graph_traversal::*;
use super::model;
use log::{debug, error};
use std::collections::{HashMap, HashSet};
// Phaing connected components.
pub fn phase_cc(paths: &[Vec<(usize, usize)>], max_occ: usize) -> (Vec<u8>, f64) {
    // Convert paths into paths with ID.
    debug!("Start");
    // Construct a graph.
    let num_of_nodes = number_of_nodes(paths);
    // BFS to determine the
    // CAUSION: ARRAYS SHOULD BE ALIGNED W.R.T BFS ORDER!
    let node_traverse_order = determine_traversal_order(num_of_nodes, paths);
    let (path_to_be_used, path_unused) = downsampling_up_to(&node_traverse_order, paths, max_occ);
    let node_paths: Vec<_> = node_traverse_order
        .iter()
        .map(|&n| get_paths_on(&path_to_be_used, &Nodes::new(&[n])))
        .collect();
    let boundary_nodes: Vec<_> = get_boundary_nodes_on(&node_traverse_order, &path_to_be_used);
    let boundary_paths: Vec<_> = boundary_nodes
        .iter()
        .map(|bn| get_paths_on(&path_to_be_used, bn))
        .collect();
    debug!(
        "Removed:Remaines={}:{}",
        path_unused.len(),
        path_to_be_used.len()
    );
    (0..num_of_nodes).for_each(|i| assert!(boundary_nodes[i].contains(&node_traverse_order[i])));
    let mut result = exact_phase(
        num_of_nodes,
        &path_to_be_used,
        &node_traverse_order,
        &node_paths,
        &boundary_nodes,
        &boundary_paths,
    );
    // PSEUOD count 1.
    let model = model::Model::new(&path_to_be_used, &result, 1);
    debug!("Finished");
    let mut paths_with_id = path_to_be_used;
    for (idx, path) in path_unused {
        result.push((idx, model.predict_path(&path)));
        paths_with_id.push((idx, path));
    }
    result.sort_by_key(|x| x.0);
    paths_with_id.sort_by_key(|x| x.0);
    {
        paths_with_id
            .iter()
            .zip(result.iter())
            .for_each(|(x, y)| assert_eq!(x.0, y.0));
    }
    let model = model::Model::new(&paths_with_id, &result, 0);
    let likelihoods = paths.iter().map(|x| model.likelihood(x)).sum::<f64>();
    debug!("LK:{:.4}", likelihoods);
    let result: Vec<u8> = result.into_iter().map(|x| x.1).collect();
    (result, likelihoods)
}

type PathWithID<'a> = (usize, &'a [(usize, usize)]);
fn downsampling_up_to<'a>(
    order: &[usize],
    paths: &'a [Vec<(usize, usize)>],
    max_occ: usize,
) -> (Vec<PathWithID<'a>>, Vec<PathWithID<'a>>) {
    let mut reads_to_be_used: Vec<_> = paths.iter().map(|p| p.as_slice()).enumerate().collect();
    reads_to_be_used.sort_by_key(|x| x.1.len());
    let boundary_nodes = get_boundary_nodes_on(&order, &reads_to_be_used);
    for (i, ns) in boundary_nodes.iter().enumerate() {
        debug!("D\t{}\t{}", i, ns.len());
    }
    let boundary_path_ids: Vec<Vec<usize>> = boundary_nodes
        .iter()
        .map(|bn| get_paths_on(&reads_to_be_used, bn))
        .map(|indices| indices.iter().map(|&i| reads_to_be_used[i].0).collect())
        .collect();
    let mut reads_unused = vec![];
    let mut node_counts: HashMap<usize, u32> = HashMap::new();
    for path in paths.iter() {
        for &(n, _) in path {
            *node_counts.entry(n).or_default() += 1;
        }
    }
    {
        let mut node_counts: Vec<_> = node_counts.iter().collect();
        node_counts.sort();
        for (n, val) in node_counts {
            debug!("BEFORE\t{}\t{}", n, val);
        }
    }
    // ReadID => Boundary ID.
    let reverse_index = {
        let mut reverse_index = vec![vec![]; reads_to_be_used.len()];
        for (i, ids) in boundary_path_ids.iter().enumerate() {
            for &id in ids.iter() {
                reverse_index[id].push(i);
            }
        }
        reverse_index
    };
    let mut boundary_path_number: Vec<_> = boundary_path_ids.iter().map(|x| x.len()).collect();
    loop {
        let next_remove_probe = {
            let mut excess_boundary: Vec<_> = boundary_path_number
                .iter()
                .enumerate()
                .filter(|&(_, &occ)| occ > max_occ)
                .collect();
            excess_boundary.sort_by_key(|x| x.1);
            // IS THIS...?
            excess_boundary.reverse();
            excess_boundary
                .into_iter()
                .filter_map(|(idx, _)| {
                    let rm_id = reads_to_be_used
                        .iter()
                        .map(|&(id, path)| {
                            let min_occ_along = path
                                .iter()
                                .filter_map(|&(n, _)| node_counts.get(&n))
                                .min()
                                .unwrap();
                            let is_ok = reverse_index[id].contains(&idx);
                            (min_occ_along, is_ok)
                        })
                        .enumerate()
                        .filter(|&(_, (&min_occ, is_ok))| min_occ > 1 && is_ok)
                        .max_by_key(|(_, x)| x.0)
                        .map(|x| x.0);
                    match rm_id {
                        Some(idx) => Some(reads_to_be_used.remove(idx)),
                        None => None,
                    }
                })
                .next()
        };
        let (removed_path_id, removed_path) = match next_remove_probe {
            Some(res) => res,
            None => break,
        };
        for &b in reverse_index[removed_path_id].iter() {
            boundary_path_number[b] -= 1;
        }
        for (n, _) in removed_path.iter() {
            *node_counts.get_mut(n).unwrap() -= 1;
        }
        reads_unused.push((removed_path_id, removed_path));
    }
    // Dump metrics from here to....
    {
        for (i, n) in boundary_path_number.iter().enumerate() {
            debug!("{}\t{}", i, n);
        }
        let mut node_counts: Vec<_> = node_counts.iter().collect();
        node_counts.sort();
        for (n, val) in node_counts {
            debug!("AFTER\t{}\t{}", n, val);
        }
    }
    // Here.
    (reads_to_be_used, reads_unused)
}

fn get_paths_on(paths: &[(usize, &[(usize, usize)])], nodes: &Nodes) -> PathSet {
    let paths: Vec<_> = paths
        .iter()
        .enumerate()
        .filter(|(_, (_, path))| path.iter().any(|(n, _)| nodes.contains(n)))
        .map(|(index, _)| index)
        .collect();
    PathSet::new(paths)
}

// Exact phasing.
fn exact_phase(
    num_of_nodes: usize,
    paths: &[(usize, &[(usize, usize)])],
    order: &[usize],
    node_paths: &[PathSet],
    boundary_nodes: &[Nodes],
    boundary_paths: &[PathSet],
) -> Vec<(usize, u8)> {
    // Strip ids from paths.
    debug!("Begin");
    for (idx, (id, p)) in paths.iter().enumerate() {
        debug!("Path:{}\t{}\t{:?}", idx, id, p);
    }
    for (i, (((br, bn), nr), on)) in boundary_paths
        .iter()
        .zip(boundary_nodes.iter())
        .zip(node_paths.iter())
        .zip(order.iter())
        .enumerate()
    {
        debug!("{}\t{:?}\t{:?}\t{:?}\t{}", i, br, bn, nr, on);
    }
    // calculate each bipartition on each nodes.
    let ls_node: Vec<Vec<_>> = order
        .iter()
        .zip(node_paths.iter())
        .map(|(&node, paths_on_node)| enumerate_all_bipartition(&paths, paths_on_node, node))
        .collect();
    assert!(ls_node.iter().all(|xs| xs.iter().all(|&x| x < 0.000001)));
    let (mut ls_hat, mut ls_hat_arg) = {
        debug!("Summarizing 0-th node");
        let convert_current_to_next =
            boundary_paths[0].get_intersection_pattern(&boundary_paths[1]);
        let intersection_size = boundary_paths[0].count_intersection(&boundary_paths[1]) as u32;
        let mut ls_hat = vec![std::f64::NEG_INFINITY; 2usize.pow(intersection_size)];
        let mut ls_hat_i_arg = vec![0; 2usize.pow(intersection_size)];
        // partition on R(D_0) and partition on R(D_0) and R(D_1)
        let mut pattern = 0;
        for partition in 0..(1 << boundary_paths[0].len() as usize) {
            if partition != 0 {
                let flip_bit = (((partition - 1) ^ partition) as usize).trailing_ones();
                pattern = convert_current_to_next.flip_from(pattern, flip_bit);
            }
            let update = ls_node[0][partition];
            if ls_hat[pattern] < update {
                ls_hat[pattern] = update;
                ls_hat_i_arg[pattern] = partition;
            }
        }
        (ls_hat, vec![ls_hat_i_arg])
    };
    for i in 1..num_of_nodes - 1 {
        debug!("Computing \\hat{{L}}[{}] from previous \\hat{{L}}", i);
        let convert_pattern_current_to_prev =
            boundary_paths[i].get_intersection_pattern(&boundary_paths[i - 1]);
        let convert_pattern_current_to_next =
            boundary_paths[i].get_intersection_pattern(&boundary_paths[i + 1]);
        let convert_pattern_current_to_node =
            boundary_paths[i].get_intersection_pattern(&node_paths[i]);
        let intersection_size = boundary_paths[i].count_intersection(&boundary_paths[i + 1]) as u32;
        let mut ls_hat_next = vec![std::f64::NEG_INFINITY; 2usize.pow(intersection_size)];
        let mut ls_hat_arg_next = vec![0; 2usize.pow(intersection_size)];
        let (mut prev_pattern, mut next_pattern, mut node_pattern) = (0, 0, 0);
        for partition in 0..(1 << boundary_paths[i].len()) as usize {
            if partition != 0 {
                let flip_bit = (((partition - 1) ^ partition) as usize).trailing_ones();
                prev_pattern = convert_pattern_current_to_prev.flip_from(prev_pattern, flip_bit);
                next_pattern = convert_pattern_current_to_next.flip_from(next_pattern, flip_bit);
                node_pattern = convert_pattern_current_to_node.flip_from(node_pattern, flip_bit);
            }
            let update = ls_hat[prev_pattern] + ls_node[i][node_pattern];
            if ls_hat_next[next_pattern] < update {
                ls_hat_next[next_pattern] = update;
                ls_hat_arg_next[next_pattern] = partition;
            }
        }
        ls_hat_arg.push(ls_hat_arg_next);
        ls_hat = ls_hat_next;
    }
    let (argmax, ls_max) = {
        let i = num_of_nodes - 1;
        let convert_pattern_current_to_prev =
            boundary_paths[i].get_intersection_pattern(&boundary_paths[i - 1]);
        let convert_pattern_current_to_node =
            boundary_paths[i].get_intersection_pattern(&node_paths[i]);
        let (mut prev_pattern, mut node_pattern) = (0, 0);
        (0..(1 << boundary_paths[i].len() as u64) as usize)
            .map(|partition| {
                if partition != 0 {
                    let flip_bit = (((partition - 1) ^ partition) as usize).trailing_ones();
                    prev_pattern =
                        convert_pattern_current_to_prev.flip_from(prev_pattern, flip_bit);
                    node_pattern =
                        convert_pattern_current_to_node.flip_from(node_pattern, flip_bit);
                }
                assert!(boundary_nodes[i].contains(&order[i]));
                (partition, ls_hat[prev_pattern] + ls_node[i][node_pattern])
            })
            .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap())
            .unwrap()
    };
    // Traceback.
    debug!("Max LK is {:?}, {:?}", ls_max, argmax);
    let (mut hap1, mut hap2) = (HashSet::new(), HashSet::new());
    // partition of 0 is the unique partition for the empty set, ((),()).
    // Current partition is the argmax patition for the R(D(`current_index`)).
    let (mut current_index, mut current_partition) = (num_of_nodes - 1, argmax);
    for (i, &r) in boundary_paths[current_index].iter().enumerate() {
        match (current_partition >> i) & 1 == 0 {
            true => hap1.insert(r),
            false => hap2.insert(r),
        };
    }
    while 0 < current_index {
        // Let's get the bipartition of R(D_{current_index-1}).
        let converter = boundary_paths[current_index]
            .get_intersection_pattern(&boundary_paths[current_index - 1]);
        let prev_argmax = ls_hat_arg[current_index - 1][converter.convert(current_partition)];
        debug!(
            "Best R(D_{}) = [{:?}] partition = {:b}",
            current_index - 1,
            boundary_paths[current_index - 1],
            prev_argmax,
        );
        // Decompose this partition.
        for (i, &r) in boundary_paths[current_index - 1].iter().enumerate() {
            match (prev_argmax >> i) & 1 == 0 {
                true => hap1.insert(r),
                false => hap2.insert(r),
            };
        }
        current_partition = prev_argmax;
        current_index -= 1;
    }
    debug!("HAP\tHAPID\tReadID\tRead");
    for &r in hap1.iter() {
        debug!("HAP\t0\t{}\t{:?}", r, paths[r]);
    }
    for &r in hap2.iter() {
        debug!("HAP\t1\t{}\t{:?}", r, paths[r]);
    }
    assert!(hap1.is_disjoint(&hap2));
    let mut phasing = vec![];
    phasing.extend(hap1.into_iter().map(|r| (paths[r].0, 0)));
    phasing.extend(hap2.into_iter().map(|r| (paths[r].0, 1)));
    phasing
}

fn enumerate_all_bipartition(
    paths: &[(usize, &[(usize, usize)])],
    path_indices: &PathSet,
    node: usize,
) -> Vec<f64> {
    let path_number = path_indices.len();
    if paths
        .iter()
        .all(|(_, path)| path.iter().all(|&(n, _)| n != node))
    {
        for (_, p) in paths.iter() {
            error!("{:?}", p);
        }
        panic!("Node {} never appears in the dataset.", node);
    }
    let cluster_num = *paths
        .iter()
        .flat_map(|(_, path)| path.iter().filter(|&&(n, _)| n == node))
        .map(|(_, c)| c)
        .max()
        .unwrap()
        + 1;
    (0..2usize.pow(path_number as u32))
        .map(|path_pattern| {
            // For read in S and T.
            let mut bi_counts = [vec![0u32; cluster_num], vec![0u32; cluster_num]];
            for (i, &path_index) in path_indices.iter().enumerate() {
                // If path should be in S if the i-th bit is 0, otherwise T.
                let bucket = (path_pattern >> i) & 1;
                for &(_, c) in paths[path_index].1.iter().filter(|&&(n, _)| n == node) {
                    bi_counts[bucket][c] += 1;
                }
            }
            // Note: sum_i u_i*ln(u_i/T) = sum_i u_i*ln(u_i) - u_i*ln(T) = sum_i u_i*ln(u_i) - T*ln(T)
            let log_likelihood = bi_counts
                .iter()
                .map(|xs| xs.iter().copied().map(xlnx).sum::<f64>())
                .sum::<f64>();
            let s_num = bi_counts[0].iter().sum::<u32>();
            let t_num = bi_counts[1].iter().sum::<u32>();
            let lk = log_likelihood - xlnx(s_num) - xlnx(t_num);
            assert!(
                lk < 0.0000001,
                "{:.3},{:.3},{:.3},{:.3}",
                lk,
                log_likelihood,
                s_num,
                t_num
            );
            lk
        })
        .collect()
}

// Return xlog_e(x). If x==0, return 0.
fn xlnx(x: u32) -> f64 {
    if x == 0 {
        0.
    } else {
        let x = x as f64;
        x * x.ln()
    }
}

#[derive(Clone, Default, PartialEq, Eq)]
struct PathSet {
    // The index of the path.
    // NOT THE ID OF THE PATH
    paths: Vec<usize>,
}

impl std::fmt::Debug for PathSet {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "(R:{:?})", self.paths)
    }
}

impl PathSet {
    fn new(mut paths: Vec<usize>) -> Self {
        paths.sort_unstable();
        Self { paths }
    }
    fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.paths.iter()
    }
    fn len(&self) -> usize {
        self.paths.len()
    }
    #[allow(dead_code)]
    fn is_empty(&self) -> bool {
        self.paths.is_empty()
    }
    fn count_intersection(&self, other: &Self) -> usize {
        // TODO: faster method is available(by merge sort like procedure.)
        let set = self.paths.iter().copied().collect::<HashSet<_>>();
        other.paths.iter().filter(|x| set.contains(x)).count()
    }
    // Create a convert pattern. It can used to convert a bipartition pattern
    // from `self` to a bipartition patten on the intersection set of `to` and `self`
    fn get_intersection_pattern(&self, to: &Self) -> IntersectPattern {
        // TODO: Faster method is available.
        let mut bitpattern = 0;
        let backset: HashSet<_> = to.iter().copied().collect();
        for (i, r) in self.iter().enumerate() {
            if backset.contains(r) {
                bitpattern |= 1 << i;
            }
        }
        IntersectPattern::new(bitpattern, self.paths.len())
    }
}

/// Bit Convert struct.
/// It behaves like `mask-squash`. In other words, if this struct contains a bit pattern x,
/// then, it convert n into m by picking every i-th bit having 1 in x.
#[derive(Debug, Clone, Default)]
pub struct IntersectPattern {
    len: usize,
    pattern: usize,
    // If the i-th bit of `pattern` is 1,
    // then, the i-th read is in the `fat_pattern[i]`-th bit in the
    // converted pattern.
    fat_pattern: Vec<usize>,
}

impl IntersectPattern {
    fn new(pattern: usize, len: usize) -> Self {
        let (mut fat_pattern, mut pointer) = (vec![], 0);
        for i in 0..len {
            fat_pattern.push(pointer);
            if ((pattern >> i) & 0b1) == 0b1 {
                pointer += 1;
            }
        }
        Self {
            pattern,
            len,
            fat_pattern,
        }
    }
    #[allow(dead_code)]
    fn update(&self, prev_pattern: usize, flip_bit: usize) -> usize {
        if ((self.pattern >> flip_bit) & 0b1) == 0b1 {
            // The flipped bit is in the intersection.
            let flip_bit = self.fat_pattern[flip_bit];
            prev_pattern ^ (0b1 << flip_bit)
        } else {
            // The flipped bit is irrelavant to the intersection.
            prev_pattern
        }
    }
    #[inline]
    fn flip_from(&self, pattern: usize, flip_bit: u32) -> usize {
        let mut pattern = pattern;
        if flip_bit == 1 {
            pattern ^= self.pattern & 0b1;
        } else {
            for i in 0..flip_bit {
                let probe = (self.pattern >> i) & 0b1;
                let mask = (1 << i) - 1;
                let location = (self.pattern & mask).count_ones();
                pattern ^= probe << location;
            }
        }
        pattern
    }
    fn convert(&self, pattern: usize) -> usize {
        assert!((pattern >> self.len).count_ones() == 0);
        let mut result = 0;
        let mut slot = 0;
        for i in 0..self.len {
            if ((self.pattern >> i) & 0b1) == 1 {
                result |= ((pattern >> i) & 1) << slot;
                slot += 1;
            }
        }
        result
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    fn flipped_bit(x: usize) -> usize {
        (x ^ (x + 1)).trailing_ones() as usize
    }
    fn increment(x: usize, bit: usize) -> usize {
        let mut x = x;
        for i in 0..bit {
            x ^= 0b1 << i;
        }
        x
    }
    #[test]
    fn bit_flip_test() {
        let x = 12918043;
        let bit = flipped_bit(x);
        assert_eq!(x + 1, increment(x, bit));
    }
    #[test]
    fn number_of_nodes_test() {
        let paths = vec![vec![(0, 1), (1, 0), (2, 0)], vec![(0, 1), (1, 1), (3, 0)]];
        let num_of_nodes = number_of_nodes(&paths);
        assert_eq!(num_of_nodes, 4);
    }
    #[test]
    fn get_boundary_paths_test() {
        let paths: Vec<Vec<(usize, usize)>> = vec![
            vec![1, 2, 6, 0],    // 0
            vec![2, 6, 0, 7],    // 1
            vec![6, 0, 7, 9],    // 2
            vec![8, 9, 7, 0],    // 3
            vec![1, 2, 4, 5],    // 4
            vec![4, 5, 3, 9, 8], // 5
            vec![3, 5, 4, 2, 1], // 6
        ]
        .into_iter()
        .map(|ps| ps.into_iter().map(|x| (x, 0)).collect())
        .collect();
        let paths: Vec<_> = paths.iter().map(|x| x.as_slice()).enumerate().collect();
        let order = vec![1, 2, 4, 6, 5, 0, 3, 7, 9, 8];
        let mut edges: HashSet<(usize, usize)> = HashSet::new();
        for (_, ps) in paths.iter() {
            for w in ps.windows(2) {
                edges.insert((w[0].0, w[1].0));
            }
        }
        assert_eq!(
            get_boundary_nodes(&paths, &order, 0, &edges),
            Nodes::new(&vec![1])
        );
        assert_eq!(
            get_boundary_nodes(&paths, &order, 4, &edges),
            Nodes::new(&vec![6, 5])
        );
        assert_eq!(
            get_boundary_nodes(&paths, &order, 7, &edges),
            Nodes::new(&vec![3, 7])
        );
        let nodes = get_boundary_nodes(&paths, &order, 0, &edges);
        assert_eq!(get_paths_on(&paths, &nodes), PathSet::new(vec![0, 4, 6]));
        let nodes = get_boundary_nodes(&paths, &order, 3, &edges);
        assert_eq!(nodes, Nodes::new(&vec![4, 6]));
        assert_eq!(
            get_paths_on(&paths, &nodes),
            PathSet::new(vec![0, 1, 2, 4, 5, 6])
        );
        let nodes = get_boundary_nodes(&paths, &order, 6, &edges);
        assert_eq!(
            get_paths_on(&paths, &nodes),
            PathSet::new(vec![0, 1, 2, 3, 5, 6])
        );
        let nodes = get_boundary_nodes(&paths, &order, 8, &edges);
        assert_eq!(get_paths_on(&paths, &nodes), PathSet::new(vec![2, 3, 5]));
        let nodes = get_boundary_nodes(&paths, &order, 2, &edges);
        assert_eq!(nodes, Nodes::new(&vec![2, 4]));
        assert_eq!(
            get_paths_on(&paths, &nodes),
            PathSet::new(vec![0, 1, 4, 5, 6])
        );
    }
    #[test]
    fn convet_test() {}
}
