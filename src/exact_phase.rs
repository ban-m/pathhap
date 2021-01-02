// TODO: We eliminate the every first read from the bipartition, as
// we can assume that the first read always belongs to the haplotype S.
// TODO: We should not compute the T haplotype, as
// all reads should be eigher S or T haplotype and we can reconstruct the haplotype S.
use log::debug;
// use rayon::prelude::*;
use std::collections::{HashSet, VecDeque};
// Exact phasing.
pub fn exact_phase(paths_w_id: &[(usize, Vec<(usize, usize)>)]) -> (Vec<(usize, u8)>, Model) {
    // Strip ids from paths.
    debug!("Begin");
    let paths: Vec<&[(usize, usize)]> = paths_w_id.iter().map(|x| x.1.as_slice()).collect();
    for (idx, p) in paths.iter().enumerate() {
        debug!("Path:{}:{:?}", idx, p);
    }
    let num_of_nodes = paths
        .iter()
        .flat_map(|x| x.iter().map(|x| x.0))
        .max()
        .unwrap()
        + 1;
    // BFS to determine the
    // CAUSION: ARRAYS SHOULD BE ALIGNED W.R.T BFS ORDER!
    let (boundary_reads, boundary_nodes, node_reads, ordered_nodes) = breadth_first_search(&paths);
    for (i, (((br, bn), nr), on)) in boundary_reads
        .iter()
        .zip(boundary_nodes.iter())
        .zip(node_reads.iter())
        .zip(ordered_nodes.iter())
        .enumerate()
    {
        debug!("{}\t{:?}\t{:?}\t{:?}\t{}", i, br, bn, nr, on);
    }
    // calculate each bipartition on each nodes.
    let ls_node: Vec<Vec<_>> = ordered_nodes
        .iter()
        .zip(node_reads.iter())
        .map(|(node, reads_on_node)| enumerate_all_bipartition(&paths, reads_on_node, *node))
        .collect();
    assert!(ls_node.iter().all(|xs| xs.iter().all(|&x| x < 0.000001)));
    debug!("{:?}", ordered_nodes);
    debug!("DUMP\tORDERE\tNODE\tLK");
    for (i, (v, xs)) in ordered_nodes.iter().zip(ls_node.iter()).enumerate() {
        let xs: Vec<_> = xs.iter().map(|x| format!("{:.2}", x)).collect();
        debug!("DUMP\t{}\t{}\t[{}]", i, v, xs.join(","));
    }
    // Calculate scores of each bipartition of the reads
    // Occured both in R(v_i) and R(D_{i-1}).
    // `arg` return the maximum partition of *R(v_i)*.
    let (ls_node_star, ls_node_star_arg): (Vec<_>, Vec<_>) = (0..num_of_nodes)
        .map(|i| match boundary_nodes[i].contains(ordered_nodes[i]) {
            true => (Vec::new(), Vec::new()),
            false => {
                assert!(i > 0);
                let intersection_size =
                    boundary_reads[i - 1].count_intersection(&node_reads[i]) as u32;
                let mut ls_node_star_i =
                    vec![std::f64::NEG_INFINITY; 2usize.pow(intersection_size)];
                let mut ls_node_arg_i = vec![0; 2usize.pow(intersection_size)];
                let intersection_pattern =
                    node_reads[i].get_intersection_pattern(&boundary_reads[i - 1]);
                for (read_pattern, &lk) in ls_node[i].iter().enumerate() {
                    let pattern = intersection_pattern.convert(read_pattern);
                    if ls_node_star_i[pattern] < lk {
                        ls_node_star_i[pattern] = lk;
                        ls_node_arg_i[pattern] = read_pattern;
                    };
                }
                (ls_node_star_i, ls_node_arg_i)
            }
        })
        .unzip();
    debug!("ARGMAX\tORDER\tNODE\tLKS");
    for (i, (star, arg)) in ls_node_star.iter().zip(ls_node_star_arg.iter()).enumerate() {
        if !boundary_nodes[i].contains(ordered_nodes[i]) {
            let star: Vec<_> = star
                .iter()
                .zip(arg)
                .map(|(lk, arg)| format!("{:.2}-{:b}", lk, arg))
                .collect();
            debug!("ARGMAX\t{}\t{:?}\t{}", i, ordered_nodes[i], star.join(","));
        }
    }
    let (mut ls, mut ls_hat, mut ls_hat_arg) = {
        // Fill hat{Ls}[0]
        // Intersection pattern of R(D_0)) => R(D_0) and R(D_1)
        let convert_pattern_boundary =
            boundary_reads[0].get_intersection_pattern(&boundary_reads[1]);
        // Intersection pattern of R(D_0)) => R(D_0) and R(v_1)
        let convert_pattern_node = boundary_reads[0].get_intersection_pattern(&node_reads[1]);
        let intersection_size = boundary_reads[0].count_intersection(&boundary_reads[1]) as u32;
        let mut ls_hat = vec![std::f64::NEG_INFINITY; 2usize.pow(intersection_size)];
        let mut ls_hat_arg = vec![0; 2usize.pow(intersection_size)];
        let ls = ls_node[0].clone();
        for partition in 0..2usize.pow(boundary_reads[0].len() as u32) {
            let pattern = convert_pattern_boundary.convert(partition);
            let update = if boundary_nodes[1].contains(ordered_nodes[1]) {
                ls[partition]
            } else {
                ls[partition] + ls_node_star[1][convert_pattern_node.convert(partition)]
            };
            if ls_hat[pattern] < update {
                ls_hat_arg[pattern] = partition;
                ls_hat[pattern] = update;
            }
        }
        debug!("====The 0th node====:{:?}", boundary_reads[0]);
        for (pat, lk) in ls.iter().enumerate() {
            debug!("{}\t{:b}", lk, pat);
        }
        debug!("====The 0th node<->the 1st node:{:?}", boundary_reads[1]);
        for (pat, lk) in ls_hat.iter().enumerate() {
            debug!("{}\t{:b}", lk, pat);
        }
        // Nobug.
        (vec![ls], vec![ls_hat], vec![ls_hat_arg])
    };
    // This is for the next node, node 1.
    for i in 1..num_of_nodes {
        debug!("partitioning {}-th node", i);
        let convert_pattern_node = boundary_reads[i].get_intersection_pattern(&node_reads[i]);
        let convert_pattern_boundary =
            boundary_reads[i].get_intersection_pattern(&boundary_reads[i - 1]);
        let mut ls_i = vec![0.; 2usize.pow(boundary_reads[i].len() as u32)];
        for partition in 0..2usize.pow(boundary_reads[i].len() as u32) {
            let converted_partition = convert_pattern_boundary.convert(partition);
            debug!(
                "Pattern:{:b}->{:b}({:?}->{:?})",
                partition,
                converted_partition,
                boundary_reads[i],
                boundary_reads[i - 1]
            );
            debug!(
                "Pattern:{:b}->{:b}({:?}->{:?})",
                partition,
                convert_pattern_node.convert(partition),
                boundary_reads[i],
                node_reads[i],
            );
            ls_i[partition] = if boundary_nodes[i].contains(ordered_nodes[i]) {
                ls_hat[i - 1][converted_partition]
                    + ls_node[i][convert_pattern_node.convert(partition)]
            } else {
                ls_hat[i - 1][converted_partition]
            };
        }
        debug!("====The {}th node====:{:?}", i, boundary_reads[i]);
        for (pat, lk) in ls_i.iter().enumerate() {
            debug!("{}\t{:b}", lk, pat);
        }
        ls.push(ls_i);
        // Update for the next iteration.
        if i + 1 < num_of_nodes {
            let convert_pattern_boundary =
                boundary_reads[i].get_intersection_pattern(&boundary_reads[i + 1]);
            let convert_pattern_node =
                boundary_reads[i].get_intersection_pattern(&node_reads[i + 1]);
            let intersection_size =
                boundary_reads[i].count_intersection(&boundary_reads[i + 1]) as u32;
            let mut ls_hat_i = vec![std::f64::NEG_INFINITY; 2usize.pow(intersection_size)];
            let mut ls_hat_i_arg = vec![0; 2usize.pow(intersection_size)];
            for partition in 0..2usize.pow(boundary_reads[i].len() as u32) {
                let pattern = convert_pattern_boundary.convert(partition);
                let update = if boundary_nodes[i + 1].contains(ordered_nodes[i + 1]) {
                    ls[i][partition]
                } else {
                    let node_pattern = convert_pattern_node.convert(partition);
                    ls[i][partition] + ls_node_star[i + 1][node_pattern]
                };
                if ls_hat_i[pattern] < update {
                    ls_hat_i[pattern] = update;
                    ls_hat_i_arg[pattern] = partition;
                }
            }
            debug!("The {}th node.", i + 1);
            for (pattern, (lk, pat)) in ls_hat_i.iter().zip(ls_hat_i_arg.iter()).enumerate() {
                debug!("{}\t{}\t{:b}", pattern, lk, pat);
            }
            ls_hat.push(ls_hat_i);
            ls_hat_arg.push(ls_hat_i_arg);
        }
    }
    // Traceback.
    assert_eq!(ls.last().unwrap().len(), 1);
    assert_eq!(ls_hat.last().unwrap().len(), 1);
    assert_eq!(ls.len(), num_of_nodes);
    assert_eq!(ls_hat.len(), num_of_nodes - 1);
    // This score is given by ls_hat[num_of_nodes-2].
    debug!("Max LK is {:?}", ls.last().unwrap());
    let (mut hap1, mut hap2) = (HashSet::new(), HashSet::new());
    // partition of 0 is the unique partition for the empty set, ((),()).
    // Current partition is the argmax patition for the R(D(`current_index`)).
    let (mut current_index, mut current_partition) = (num_of_nodes - 1, 0);
    while 0 < current_index {
        // Let's get the bipartition of R(D_{current_index-1}).
        let converter = boundary_reads[current_index]
            .get_intersection_pattern(&boundary_reads[current_index - 1]);
        let prev_argmax = ls_hat_arg[current_index - 1][converter.convert(current_partition)];
        debug!(
            "Best R(D_{}) partition = {:b}",
            current_index - 1,
            prev_argmax
        );
        // Decompose this partition.
        for (i, &r) in boundary_reads[current_index - 1].iter().enumerate() {
            match (prev_argmax >> i) & 1 == 0 {
                true => hap1.insert(r),
                false => hap2.insert(r),
            };
        }
        // Let's get the bipartition of R(v_current_index) if this vertex is not treated in the current_index-th iteration.
        if !boundary_nodes[current_index].contains(ordered_nodes[current_index]) {
            // `prev_argmax` is a bipartition on R(D_{current_index-1}).
            let converter = boundary_reads[current_index - 1]
                .get_intersection_pattern(&node_reads[current_index]);
            let partition_on_reads =
                ls_node_star_arg[current_index][converter.convert(prev_argmax)];
            for (i, &r) in node_reads[current_index].iter().enumerate() {
                match (partition_on_reads >> i) & 1 == 0 {
                    true => hap1.insert(r),
                    false => hap2.insert(r),
                };
            }
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
    let model = Model::new(&paths, &hap1, &hap2);
    let mut phasing = vec![];
    phasing.extend(hap1.into_iter().map(|r| (paths_w_id[r].0, 0)));
    phasing.extend(hap2.into_iter().map(|r| (paths_w_id[r].0, 1)));
    (phasing, model)
}

fn enumerate_all_bipartition(
    paths: &[&[(usize, usize)]],
    path_indices: &PathSet,
    node: usize,
) -> Vec<f64> {
    let path_number = path_indices.len();
    let cluster_num = *paths
        .iter()
        .flat_map(|path| path.iter().filter(|&&(n, _)| n == node))
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
                for &(_, c) in paths[path_index].iter().filter(|&&(n, _)| n == node) {
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

// i-th read set(R(D_i)), i-th boundary nodes D_i, and the i-th nodes v_i.
fn breadth_first_search(
    paths: &[&[(usize, usize)]],
) -> (Vec<PathSet>, Vec<Nodes>, Vec<PathSet>, Vec<usize>) {
    let num_nodes = paths
        .iter()
        .flat_map(|x| x.iter().map(|x| x.0))
        .max()
        .unwrap()
        + 1;
    let mut edges = vec![HashSet::new(); num_nodes];
    for path in paths.iter() {
        for w in path.windows(2) {
            let (from, to) = (w[0].0, w[1].0);
            edges[from].insert(to);
            edges[to].insert(from);
        }
    }
    let (boundary_nodes, bfs_order) = breadth_first_search_inner(&edges, num_nodes);
    let boundary_reads: Vec<_> = boundary_nodes
        .iter()
        .map(|nodes| {
            paths
                .iter()
                .enumerate()
                .filter_map(|(idx, path)| {
                    if path.iter().any(|(n, _)| nodes.contains(*n)) {
                        Some(idx)
                    } else {
                        None
                    }
                })
                .collect()
        })
        .map(PathSet::new)
        .collect();
    let node_reads: Vec<_> = bfs_order
        .iter()
        .map(|node| {
            paths
                .iter()
                .enumerate()
                .filter_map(|(idx, path)| {
                    if path.iter().any(|(n, _)| n == node) {
                        Some(idx)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        .map(PathSet::new)
        .collect();
    (boundary_reads, boundary_nodes, node_reads, bfs_order)
}

fn breadth_first_search_inner(
    edges: &[HashSet<usize>],
    num_of_nodes: usize,
) -> (Vec<Nodes>, Vec<usize>) {
    let order = bfs(edges, num_of_nodes);
    let mut current_set = HashSet::new();
    let boundary_nodes: Vec<_> = order
        .iter()
        .map(|&node| {
            current_set.insert(node);
            let boundary: HashSet<_> = current_set
                .iter()
                .filter(|&&n| {
                    // Check if there is v, not in `current_set`.
                    edges[n].iter().any(|to| !current_set.contains(to))
                })
                .copied()
                .collect();
            Nodes { nodes: boundary }
        })
        .collect();
    (boundary_nodes, order)
}

// Breadth first search and returns the arriving order.
fn bfs(edges: &[HashSet<usize>], num_of_nodes: usize) -> Vec<usize> {
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

#[derive(Debug, Clone, Default)]
pub struct Model {
    // haplotype -> position -> distribution.
    haps: Vec<Vec<Vec<f64>>>,
}

impl Model {
    const PLOIDY: usize = 2;
    const PSEUDO_COUNT: usize = 1;
    // Predict short path
    pub fn predict_path(&self, path: &[(usize, usize)]) -> u8 {
        self.haps
            .iter()
            .map(|hap| path.iter().map(|&(n, c)| hap[n][c].ln()).sum::<f64>())
            .enumerate()
            .map(|(x, y)| (x as u8, y))
            .max_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
            .unwrap()
            .0
    }
    pub fn new(paths: &[&[(usize, usize)]], hap1: &HashSet<usize>, hap2: &HashSet<usize>) -> Self {
        let nodes_and_clusters = {
            let num_of_nodes = paths
                .iter()
                .flat_map(|x| x.iter().map(|&(n, _)| n))
                .max()
                .unwrap()
                + 1;
            let mut clusters = vec![Self::PSEUDO_COUNT; num_of_nodes];
            for path in paths.iter() {
                for &(n, c) in path.iter() {
                    clusters[n] = clusters[n].max(c + 1);
                }
            }
            clusters
        };

        let mut hap_count: Vec<Vec<Vec<u32>>> = (0..Self::PLOIDY)
            .map(|_| nodes_and_clusters.iter().map(|&cl| vec![0; cl]).collect())
            .collect();
        for (r, path) in paths.iter().enumerate() {
            assert!(hap1.contains(&r) || hap2.contains(&r));
            let index = if hap1.contains(&r) { 0 } else { 1 };
            for &(n, c) in path.iter() {
                hap_count[index][n][c] += 1;
            }
        }
        let haps: Vec<Vec<Vec<f64>>> = hap_count
            .iter()
            .map(|xss| {
                xss.iter()
                    .map(|xs| {
                        let sum = xs.iter().sum::<u32>();
                        xs.iter().map(|&x| x as f64 / sum as f64).collect()
                    })
                    .collect()
            })
            .collect();
        Self { haps }
    }
}

#[derive(Clone, Default)]
struct PathSet {
    // The index of the path.
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
    pub fn contains(&self, c: usize) -> bool {
        self.nodes.contains(&c)
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use std::collections::HashMap;
    #[test]
    fn bfs_test() {
        let num_nodes = 7;
        let edges: Vec<HashSet<usize>> = vec![
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
        let orders = bfs(&edges, num_nodes);
        eprintln!("Order:{:?}", orders);
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
    fn convet_test() {}
}
