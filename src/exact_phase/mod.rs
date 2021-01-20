// TODO: We eliminate the every first read from the bipartition, as
// we can assume that the first read always belongs to the haplotype S.
// TODO: We should not compute the T haplotype, as
// all reads should be eigher S or T haplotype and we can reconstruct the haplotype S.
use log::debug;
// use rayon::prelude::*;
use std::collections::{HashMap, HashSet, VecDeque};
mod model;

fn number_of_nodes(paths: &[Vec<(usize, usize)>]) -> usize {
    paths
        .iter()
        .filter_map(|x| x.iter().map(|x| x.0).max())
        .max()
        .unwrap()
        + 1
}

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
        .map(|&n| get_paths_on(&path_to_be_used, &Nodes::new(&vec![n])))
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

// BFS and return the order.
// We choose the un-traversed smallest children at each step.
fn determine_traversal_order(num_of_nodes: usize, paths: &[Vec<(usize, usize)>]) -> Vec<usize> {
    let mut edges = vec![HashSet::new(); num_of_nodes];
    for path in paths.iter() {
        for w in path.windows(2) {
            let (from, to) = (w[0].0, w[1].0);
            edges[from].insert(to);
            edges[to].insert(from);
        }
    }
    let edges: Vec<Vec<_>> = edges
        .into_iter()
        .map(|eds| {
            let mut eds: Vec<_> = eds.into_iter().collect();
            eds.sort();
            eds
        })
        .collect();
    bfs(num_of_nodes, &edges)
}

// Breadth first search and returns the arriving order.
fn bfs(num_of_nodes: usize, edges: &[Vec<usize>]) -> Vec<usize> {
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
    assert!(is_arrived.iter().all(|&x| x));
    order
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
    // ID, not the indices.
    let boundary_path_ids: Vec<Vec<usize>> = boundary_nodes
        .iter()
        .map(|bn| get_paths_on(&reads_to_be_used, bn))
        .map(|indices| indices.iter().map(|&i| reads_to_be_used[i].0).collect())
        .collect();
    let mut reads_unused = vec![];
    let mut edge_counts: HashMap<(usize, usize), u32> = {
        let mut count = HashMap::new();
        for path in paths.iter() {
            for w in path.windows(2) {
                let (f, t) = ((w[0].0).min(w[1].0), (w[0].0).max(w[0].1));
                *count.entry((f, t)).or_default() += 1;
            }
        }
        count
    };
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
        let (idx, &max) = boundary_path_number
            .iter()
            .enumerate()
            .max_by_key(|x| x.1)
            .unwrap();
        if max <= max_occ {
            break;
        }
        let (removed_path_id, removed_path) = {
            let rm_id = reads_to_be_used.iter().position(|&(id, path)| {
                let on_max_boundary = reverse_index[id].contains(&idx);
                let ok_to_remove = path.windows(2).all(|w| {
                    let (f, t) = ((w[0].0).min(w[1].0), (w[0].0).max(w[0].1));
                    *edge_counts.get(&(f, t)).unwrap() > 0
                });
                on_max_boundary && ok_to_remove
            });
            match rm_id {
                Some(idx) => reads_to_be_used.remove(idx),
                None => break,
            }
        };
        for &b in reverse_index[removed_path_id].iter() {
            boundary_path_number[b] -= 1;
        }
        for w in removed_path.windows(2) {
            let (f, t) = ((w[0].0).min(w[1].0), (w[0].0).max(w[0].1));
            *edge_counts.entry((f, t)).or_default() -= 1;
        }
        reads_unused.push((removed_path_id, removed_path));
    }
    (reads_to_be_used, reads_unused)
}

fn get_boundary_nodes_on(order: &[usize], paths: &[(usize, &[(usize, usize)])]) -> Vec<Nodes> {
    let edges: HashSet<(usize, usize)> = {
        let mut edges = HashSet::new();
        for (_, path) in paths.iter() {
            for w in path.windows(2) {
                edges.insert((w[0].0, w[1].0));
            }
        }
        edges
    };
    (0..order.len())
        .map(|i| get_boundary_nodes(&paths, &order, i, &edges))
        .collect()
}

// Take P={p_1,..,p_n}, V={v_1, ...,v_m}, i, and edges E, then
// let U = {v_1,..., v_i} and return
// D(U) = {v_i} + {u \in U | there is some node w not in U such that (w,u) \in E }
fn get_boundary_nodes(
    _paths: &[(usize, &[(usize, usize)])],
    order: &[usize],
    iteration: usize,
    edges: &HashSet<(usize, usize)>,
) -> Nodes {
    let mut nodes: HashSet<_> = order
        .iter()
        .take(iteration + 1)
        .filter(|&&node| {
            let is_connected_to_outer = order
                .iter()
                .skip(iteration + 1)
                .any(|&m| edges.contains(&(node, m)) || edges.contains(&(m, node)));
            is_connected_to_outer
        })
        .copied()
        .collect();
    nodes.insert(order[iteration]);
    Nodes { nodes }
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
    debug!("{:?}", order);
    debug!("DUMP\tORDERE\tNODE\tLK");
    for (i, (v, xs)) in order.iter().zip(ls_node.iter()).enumerate() {
        let xs: Vec<_> = xs.iter().map(|x| format!("{:.2}", x)).collect();
        debug!("DUMP\t{}\t{}\t[{}]", i, v, xs.join(","));
    }
    let (mut ls_hat, mut ls_hat_arg) = {
        debug!("Summarizing 0-th node");
        let convert_pattern_boundary =
            boundary_paths[0].get_intersection_pattern(&boundary_paths[1]);
        let intersection_size = boundary_paths[0].count_intersection(&boundary_paths[1]) as u32;
        let mut ls_hat = vec![std::f64::NEG_INFINITY; 2usize.pow(intersection_size)];
        let mut ls_hat_i_arg = vec![0; 2usize.pow(intersection_size)];
        for partition in 0..2usize.pow(boundary_paths[0].len() as u32) {
            let pattern = convert_pattern_boundary.convert(partition);
            let update = ls_node[0][partition];
            if ls_hat[pattern] < update {
                ls_hat[pattern] = update;
                ls_hat_i_arg[pattern] = partition;
            }
        }
        (ls_hat, vec![ls_hat_i_arg])
    };
    // let (mut ls, mut ls_hat, mut ls_hat_arg) = (ls_node[0].clone(), vec![], vec![]);
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
        for partition in 0..2usize.pow(boundary_paths[i].len() as u32) {
            let prev_pattern = convert_pattern_current_to_prev.convert(partition);
            let next_pattern = convert_pattern_current_to_next.convert(partition);
            let node_pattern = convert_pattern_current_to_node.convert(partition);
            assert!(boundary_nodes[i].contains(&order[i]));
            let update = ls_hat[prev_pattern] + ls_node[i][node_pattern];
            if ls_hat_next[next_pattern] < update {
                ls_hat_next[next_pattern] = update;
                ls_hat_arg_next[next_pattern] = partition;
            }
        }
        ls_hat_arg.push(ls_hat_arg_next);
        ls_hat = ls_hat_next;
    }
    // debug!("Summarizing {}-th node", i - 1);
    // let convert_pattern_boundary =
    //     boundary_paths[i - 1].get_intersection_pattern(&boundary_paths[i]);
    // let intersection_size = boundary_paths[i - 1].count_intersection(&boundary_paths[i]) as u32;
    // ls_hat.clear();
    // ls_hat
    //     .extend(std::iter::repeat(std::f64::NEG_INFINITY).take(2usize.pow(intersection_size)));
    // let mut ls_hat_i_arg = vec![0; 2usize.pow(intersection_size)];
    // for partition in 0..2usize.pow(boundary_paths[i - 1].len() as u32) {
    //     let pattern = convert_pattern_boundary.convert(partition);
    //     assert!(boundary_nodes[i].contains(&order[i]));
    //     let update = ls[partition];
    //     if ls_hat[pattern] < update {
    //         ls_hat[pattern] = update;
    //         ls_hat_i_arg[pattern] = partition;
    //     }
    // }
    // ls_hat_arg.push(ls_hat_i_arg);
    // debug!("partitioning {}-th node", i);
    // let convert_pattern_node = boundary_paths[i].get_intersection_pattern(&node_paths[i]);
    // let convert_pattern_boundary =
    //     boundary_paths[i].get_intersection_pattern(&boundary_paths[i - 1]);
    // ls.clear();
    // ls.extend(std::iter::repeat(0.).take(2usize.pow(boundary_paths[i].len() as u32)));
    // for partition in 0..2usize.pow(boundary_paths[i].len() as u32) {
    //     let converted_partition = convert_pattern_boundary.convert(partition);
    //     ls[partition] =
    //         ls_hat[converted_partition] + ls_node[i][convert_pattern_node.convert(partition)];
    // }
    // }
    let (argmax, ls_max) = {
        let i = num_of_nodes - 1;
        let convert_pattern_current_to_prev =
            boundary_paths[i - 1].get_intersection_pattern(&boundary_paths[i]);
        let convert_pattern_current_to_node =
            boundary_paths[i].get_intersection_pattern(&node_paths[i]);
        (0..2usize.pow(boundary_paths[i].len() as u32))
            .map(|partition| {
                let prev_pattern = convert_pattern_current_to_prev.convert(partition);
                let node_pattern = convert_pattern_current_to_node.convert(partition);
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
        paths.sort();
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
/// #Example
/// ```rust
/// let pattern = 0b10001;
#[derive(Debug, Clone, Default)]
pub struct IntersectPattern {
    len: usize,
    pattern: usize,
}

impl IntersectPattern {
    fn new(pattern: usize, len: usize) -> Self {
        Self { pattern, len }
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

#[derive(Debug, Clone, Default)]
struct Nodes {
    nodes: HashSet<usize>,
}

impl Nodes {
    pub fn new(node: &[usize]) -> Self {
        let nodes: HashSet<_> = node.iter().copied().collect();
        Nodes { nodes }
    }
    pub fn contains(&self, c: &usize) -> bool {
        self.nodes.contains(&c)
    }
    pub fn len(&self) -> usize {
        self.nodes.len()
    }
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

impl std::cmp::PartialEq for Nodes {
    fn eq(&self, other: &Self) -> bool {
        self.nodes.intersection(&other.nodes).count() == self.len()
            && other.len() == self.nodes.len()
    }
}

impl std::cmp::Eq for Nodes {}

#[cfg(test)]
pub mod test {
    use super::*;
    use std::collections::HashMap;
    #[test]
    fn number_of_nodes_test() {
        let paths = vec![vec![(0, 1), (1, 0), (2, 0)], vec![(0, 1), (1, 1), (3, 0)]];
        let num_of_nodes = number_of_nodes(&paths);
        assert_eq!(num_of_nodes, 4);
    }
    #[test]
    fn traversal_order_test() {
        let paths: Vec<Vec<(usize, usize)>> = vec![vec![0, 1, 2, 3, 2, 4]]
            .into_iter()
            .map(|xs| xs.into_iter().map(|x| (x, 0)).collect())
            .collect();
        let order = determine_traversal_order(5, &paths);
        assert_eq!(order, vec![0, 1, 2, 3, 4]);
        let paths: Vec<Vec<(usize, usize)>> = vec![
            vec![1, 0],
            vec![0, 3],
            vec![5, 3],
            vec![3, 4],
            vec![4, 2],
            vec![2, 5],
            vec![4, 2, 6],
        ]
        .into_iter()
        .map(|xs| xs.into_iter().map(|x| (x, 0)).collect())
        .collect();
        let order = determine_traversal_order(7, &paths);
        assert_eq!(order, vec![1, 0, 3, 4, 5, 2, 6]);
    }
    #[test]
    fn bfs_test() {
        let edges = vec![vec![1], vec![0, 2, 3, 4], vec![1, 4], vec![1], vec![1, 2]];
        let order = bfs(5, &edges);
        assert_eq!(order, vec![0, 1, 2, 3, 4]);
        let edges = vec![
            vec![1, 3],
            vec![0],
            vec![4, 5, 6],
            vec![0, 4, 5],
            vec![2, 3],
            vec![2, 3],
            vec![2],
        ];
        let order = bfs(7, &edges);
        assert_eq!(order, vec![1, 0, 3, 4, 5, 2, 6]);
    }
    #[test]
    fn bfs_test_2() {
        let num_nodes = 7;
        let edges: Vec<Vec<usize>> = vec![
            vec![1],
            vec![1, 2, 3],
            vec![1],
            vec![1, 4, 5, 6],
            vec![3],
            vec![3],
            vec![3],
        ]
        .into_iter()
        .map(|x| x.into_iter().collect())
        .collect();
        let orders = bfs(num_nodes, &edges);
        // eprintln!("Order:{:?}", orders);
        let mut arrived_nodes: Vec<usize> = vec![];
        let mut count: HashMap<usize, u32> = (0..num_nodes).map(|x| (x, 0)).collect();
        for i in orders {
            *count.get_mut(&i).unwrap() += 1;
            if !arrived_nodes.is_empty() {
                let is_ok = arrived_nodes.iter().any(|&f| edges[f].contains(&i));
                assert!(is_ok);
            }
            arrived_nodes.push(i);
        }
        assert!(count.values().all(|&x| x == 1));
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
