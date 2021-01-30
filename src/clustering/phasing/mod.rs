//! phasing module by GraphWhatsHap.
use super::graph_traversal::*;
use super::model;
use log::{debug, error, trace};
use rayon::prelude::*;
use std::collections::HashSet;
mod xlogx;
// use std::collections::BTreeMap;
// use std::collections::HashMap;
// Phaing connected components.
pub fn phase_cc(
    paths: &[Vec<(usize, usize)>],
    max_occ: usize,
    subsample_size: Option<usize>,
) -> Vec<u8> {
    let min_occ = 6;
    // Convert paths into paths with ID.
    debug!("Start");
    // Construct a graph.
    let num_of_nodes = number_of_nodes(paths);
    // BFS to determine the
    // CAUSION: ARRAYS SHOULD BE ALIGNED W.R.T BFS ORDER!
    let node_traverse_order = determine_traversal_order(num_of_nodes, paths);
    let start = std::time::Instant::now();
    let (path_to_be_used, path_unused) =
        downsampling_up_to(&node_traverse_order, paths, max_occ, min_occ);
    let end = std::time::Instant::now();
    debug!(
        "Removed:Remains={}:{}({}ms)",
        path_unused.len(),
        path_to_be_used.len(),
        (end - start).as_millis()
    );
    let node_paths: Vec<_> = node_traverse_order
        .iter()
        .map(|&n| get_paths_on(&path_to_be_used, &Nodes::new(&[n])))
        .collect();
    let boundary_nodes: Vec<_> = get_boundary_nodes_on(&node_traverse_order, &path_to_be_used);
    let boundary_paths: Vec<_> = boundary_nodes
        .iter()
        .map(|bn| get_paths_on(&path_to_be_used, bn))
        .collect();
    (0..num_of_nodes).for_each(|i| assert!(boundary_nodes[i].contains(&node_traverse_order[i])));
    for bp in boundary_paths.iter() {
        if bp.len() <= 60 {
            error!("{:?},{}", bp, bp.len());
            error!("The boundary is too large.");
            error!("Please check input and clean up your paths a little bit.");
            panic!();
        }
    }
    let result = match subsample_size {
        None => phase_paths(
            num_of_nodes,
            &path_to_be_used,
            &node_traverse_order,
            &node_paths,
            &boundary_paths,
        ),
        Some(s) => phase_paths_fast(
            num_of_nodes,
            &path_to_be_used,
            &node_traverse_order,
            &node_paths,
            &boundary_paths,
            s,
        ),
    };
    let model = model::Model::new(&path_to_be_used, &result, 0.1);
    debug!("Finished");
    let mut paths_with_id = path_to_be_used;
    paths_with_id.extend(path_unused);
    // for (idx, path) in path_unused {
    //     result.push((idx, model.predict_path(&path)));
    //     paths_with_id.push((idx, path));
    // }
    let mut result: Vec<_> = paths_with_id
        .iter()
        .map(|&(idx, path)| (idx, model.predict_path(&path)))
        .collect();
    result.sort_by_key(|x| x.0);
    paths_with_id.sort_by_key(|x| x.0);
    {
        paths_with_id
            .iter()
            .zip(result.iter())
            .for_each(|(x, y)| assert_eq!(x.0, y.0));
    }
    result.into_iter().map(|x| x.1).collect()
}

type PathWithID<'a> = (usize, &'a [(usize, usize)]);
fn downsampling_up_to<'a>(
    order: &[usize],
    paths: &'a [Vec<(usize, usize)>],
    max_occ: usize,
    min_occ: usize,
) -> (Vec<PathWithID<'a>>, Vec<PathWithID<'a>>) {
    debug!("Sorting reads.");
    let num_of_nodes = super::graph_traversal::number_of_nodes(&paths);
    let mut reads_unused: Vec<_> = paths.iter().map(|p| p.as_slice()).enumerate().collect();
    let mut reads_to_be_used = vec![];
    reads_unused.sort_by(|(_, x), (_, y)| x.len().cmp(&y.len()).reverse());
    let boundary_nodes = get_boundary_nodes_on(&order, &reads_unused);
    for (n, bn) in boundary_nodes.iter().enumerate() {
        debug!("{}\t{}", n, bn.len());
    }
    debug!("Start Downsampling");
    debug!("Target coverage:{}", max_occ);
    let boundary_path_ids: Vec<Vec<usize>> = boundary_nodes
        .iter()
        .map(|bn| get_paths_on(&reads_unused, bn))
        .map(|indices| indices.iter().map(|&i| reads_unused[i].0).collect())
        .collect();
    let mut node_counts: Vec<_> = vec![0; num_of_nodes];
    // ReadID => Boundary ID.
    let reverse_index = {
        let mut reverse_index = vec![vec![]; reads_unused.len()];
        for (i, ids) in boundary_path_ids.iter().enumerate() {
            for &id in ids.iter() {
                reverse_index[id].push(i);
            }
        }
        reverse_index
    };
    let mut boundary_path_number: Vec<_> = vec![0; boundary_nodes.len()];
    loop {
        let next_cand = get_next_candidate(
            &boundary_path_number,
            &node_counts,
            &reads_unused,
            &reverse_index,
            max_occ,
        );
        let (removed_path_id, removed_path) = match next_cand {
            Some(idx) => reads_unused.remove(idx),
            None => break,
        };
        for &b in reverse_index[removed_path_id].iter() {
            boundary_path_number[b] += 1;
        }
        for (n, _) in removed_path.iter() {
            node_counts[*n] += 1;
        }
        reads_to_be_used.push((removed_path_id, removed_path))
    }
    // Rescue reads.
    reads_unused.reverse();
    while let Some((node, pos)) = get_rescued_path(&node_counts, &reads_unused, min_occ) {
        let (removed_path_id, removed_path) = reads_unused.remove(pos);
        for &b in reverse_index[removed_path_id].iter() {
            boundary_path_number[b] += 1;
        }
        for (n, _) in removed_path.iter() {
            node_counts[*n] += 1;
        }
        let removed_path = cut_path_contains(node, &removed_path);
        reads_to_be_used.push((removed_path_id, removed_path))
    }
    for (i, n) in boundary_path_number.iter().enumerate() {
        debug!("BOUND\t{}\t{}", i, n);
    }
    for (n, val) in node_counts.iter().enumerate() {
        debug!("NODE\t{}\t{}", n, val);
    }
    for (node, _) in node_counts.iter().enumerate() {
        if reads_to_be_used
            .iter()
            .all(|(_, path)| path.iter().all(|&(n, _)| n != node))
        {
            for (_, p) in reads_to_be_used.iter() {
                error!("{:?}", p);
            }
            error!("Node {} never appears in the dataset.", node);
            error!("Please set more large value for maximum coverage.");
            error!("Current:{}", max_occ);
            panic!()
        }
    }
    (reads_to_be_used, reads_unused)
}

fn get_next_candidate(
    boundary_path_number: &[usize],
    node_counts: &[usize],
    paths: &[(usize, &[(usize, usize)])],
    reverse_index: &[Vec<usize>],
    max_occ: usize,
) -> Option<usize> {
    let mut node_counts: Vec<_> = node_counts.iter().enumerate().collect();
    node_counts.sort_by_key(|x| x.1);
    node_counts.iter().find_map(|&(node, _)| {
        paths.iter().position(|&(id, path)| {
            path.iter().any(|&(n, _)| n == node)
                && reverse_index[id]
                    .iter()
                    .all(|&b_id| boundary_path_number[b_id] < max_occ)
        })
    })
}

fn get_rescued_path(
    node_counts: &[usize],
    paths: &[(usize, &[(usize, usize)])],
    min_occ: usize,
) -> Option<(usize, usize)> {
    let mut node_counts: Vec<_> = node_counts
        .iter()
        .enumerate()
        .filter(|&(_, &occ)| occ < min_occ)
        .collect();
    node_counts.sort_by_key(|x| x.1);
    node_counts.reverse();
    node_counts
        .iter()
        .filter(|&(_, &occ)| occ < min_occ)
        .find_map(|&(node, _)| {
            paths
                .iter()
                .position(|&(_, path)| path.iter().any(|&(n, _)| n == node))
                .map(|p| (node, p))
        })
}

fn cut_path_contains(node: usize, path: &[(usize, usize)]) -> &[(usize, usize)] {
    if path.len() <= 3 {
        path
    } else {
        let position = path.iter().position(|&(n, _)| n == node).unwrap();
        if position == 0 {
            &path[..3]
        } else if position == path.len() - 1 {
            &path[path.len() - 2..]
        } else {
            &path[position - 1..position + 1]
        }
    }
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

fn phase_paths_fast(
    num_of_nodes: usize,
    paths: &[(usize, &[(usize, usize)])],
    order: &[usize],
    node_paths: &[PathSet],
    boundary_paths: &[PathSet],
    sub: usize,
) -> Vec<(usize, u8)> {
    assert_eq!(num_of_nodes, boundary_paths.len());
    // Strip ids from paths.
    debug!("Begin");
    // calculate each bipartition on each nodes.
    let start = std::time::Instant::now();
    let ls_node: Vec<Vec<_>> = order
        .par_iter()
        .zip(node_paths.par_iter())
        .map(|(&node, paths_on_node)| enumerate_all_bipartition(&paths, paths_on_node, node))
        .collect();
    assert!(ls_node.iter().all(|xs| xs.iter().all(|&x| x < 0.000001)));
    let end = std::time::Instant::now();
    debug!("Precompute partitions:{}ms", (end - start).as_millis());
    let mut ls_hat = vec![];
    let start = std::time::Instant::now();
    debug!("Marging matrices.");
    for (i, ls_node) in ls_node.iter().enumerate() {
        let last_ls = ls_hat.last();
        // let start = std::time::Instant::now();
        let next_ls = fill_next_ls_sub(&boundary_paths, i, &node_paths[i], ls_node, last_ls, sub);
        ls_hat.push(next_ls);
        // let end = std::time::Instant::now();
        // eprintln!("{}\t{}", i, (end - start).as_millis());
    }
    debug!("GHap:{}", (std::time::Instant::now() - start).as_millis());
    // Traceback.
    //let (argmax, (_, max_lk)) = ls_hat
    let (argmax, _, max_lk) = ls_hat
        .pop()
        .unwrap()
        .into_iter()
        .max_by(|a, b| (a.1.partial_cmp(&b.1).unwrap()))
        //.max_by(|(_, a), (_, b)| (a.0.partial_cmp(&b.0).unwrap()))
        .unwrap();
    debug!("Max LK is {:?}, {:?}", max_lk, argmax);
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
        let converter = boundary_paths[current_index]
            .get_intersection_pattern(&boundary_paths[current_index - 1]);
        let prev_partition = converter.convert(current_partition);
        let &(_, _, prev_argmax) = ls_hat[current_index - 1]
            .iter()
            .find(|&x| x.0 == prev_partition)
            .unwrap();
        // let prev_argmax = ls_hat[current_index - 1][&converter.convert(current_partition)].1;
        for (i, &r) in boundary_paths[current_index - 1].iter().enumerate() {
            match (prev_argmax >> i) & 1 == 0 {
                true => hap1.insert(r),
                false => hap2.insert(r),
            };
        }
        current_partition = prev_argmax;
        current_index -= 1;
    }
    if log::log_enabled!(log::Level::Trace) {
        trace!("HAP\tHAPID\tReadID\tRead");
        for &r in hap1.iter() {
            trace!("HAP\t0\t{}\t{:?}", r, paths[r]);
        }
        for &r in hap2.iter() {
            trace!("HAP\t1\t{}\t{:?}", r, paths[r]);
        }
    }
    assert!(hap1.is_disjoint(&hap2));
    hap1.into_iter()
        .map(|r| (paths[r].0, 0))
        .chain(hap2.into_iter().map(|r| (paths[r].0, 1)))
        .collect()
}

// Exact phasing.
fn phase_paths(
    num_of_nodes: usize,
    paths: &[(usize, &[(usize, usize)])],
    order: &[usize],
    node_paths: &[PathSet],
    boundary_paths: &[PathSet],
) -> Vec<(usize, u8)> {
    assert_eq!(num_of_nodes, boundary_paths.len());
    // Strip ids from paths.
    use histgram_viz::Histgram;
    debug!("Begin");
    if log::log_enabled!(log::Level::Debug) {
        let lens: Vec<_> = paths.iter().map(|x| x.1.len()).collect();
        let hist = Histgram::new(&lens);
        eprintln!("{}", hist.format(20, 20));
    }
    // for (idx, (id, p)) in paths.iter().enumerate() {
    //     debug!("Path:{}\t{}\t{:?}", idx, id, p);
    // }
    // for (i, (((br, bn), nr), on)) in boundary_paths
    //     .iter()
    //     .zip(boundary_nodes.iter())
    //     .zip(node_paths.iter())
    //     .zip(order.iter())
    //     .enumerate()
    // {
    //     debug!("{}\t{:?}\t{:?}\t{:?}\t{}", i, br, bn, nr, on);
    // }
    // calculate each bipartition on each nodes.
    let start = std::time::Instant::now();
    let ls_node: Vec<Vec<_>> = order
        .par_iter()
        .zip(node_paths.par_iter())
        .map(|(&node, paths_on_node)| enumerate_all_bipartition(&paths, paths_on_node, node))
        .collect();
    assert!(ls_node.iter().all(|xs| xs.iter().all(|&x| x < 0.000001)));
    let end = std::time::Instant::now();
    debug!("Precompute partitions:{}ms", (end - start).as_millis());
    let (mut ls_hat, mut ls_hat_arg) = (vec![], vec![]);
    let start = std::time::Instant::now();
    debug!("Marging matrices.");
    for (i, ls_node) in ls_node.iter().enumerate() {
        // debug!("Computing \\hat{{L}}[{}] from previous \\hat{{L}}", i);
        let (ls_hat_next, ls_hat_arg_next) =
            fill_next_ls_par(&boundary_paths, i, &node_paths[i], ls_node, &ls_hat);
        ls_hat = ls_hat_next;
        ls_hat_arg.push(ls_hat_arg_next);
    }
    debug!("GHap:{}", (std::time::Instant::now() - start).as_millis());
    // Traceback.
    let (argmax, ls_max) = ls_hat_arg
        .pop()
        .unwrap()
        .into_iter()
        .zip(ls_hat)
        .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap())
        .unwrap();
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
        // debug!(
        //     "Best R(D_{}) = [{:?}] partition = {:b}",
        //     current_index - 1,
        //     boundary_paths[current_index - 1],
        //     prev_argmax,
        // );
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
    if log::log_enabled!(log::Level::Trace) {
        trace!("HAP\tHAPID\tReadID\tRead");
        for &r in hap1.iter() {
            trace!("HAP\t0\t{}\t{:?}", r, paths[r]);
        }
        for &r in hap2.iter() {
            trace!("HAP\t1\t{}\t{:?}", r, paths[r]);
        }
    }
    assert!(hap1.is_disjoint(&hap2));
    hap1.into_iter()
        .map(|r| (paths[r].0, 0))
        .chain(hap2.into_iter().map(|r| (paths[r].0, 1)))
        .collect()
}

#[allow(dead_code)]
fn fill_next_ls(
    boundary_paths: &[PathSet],
    i: usize,
    node_paths: &PathSet,
    ls_node: &[f64],
    ls_hat: &[f64],
) -> (Vec<f64>, Vec<usize>) {
    if i == 0 {
        let convert_current_to_next =
            boundary_paths[0].get_intersection_pattern(&boundary_paths[1]);
        let intersection_size = boundary_paths[0].count_intersection(&boundary_paths[1]) as u32;
        let mut ls_hat = vec![std::f64::NEG_INFINITY; 2usize.pow(intersection_size)];
        let mut ls_hat_i_arg = vec![0; 2usize.pow(intersection_size)];
        // partition on R(D_0) and partition on R(D_0) and R(D_1)
        let mut pattern = 0;
        for partition in 0usize..(1 << boundary_paths[0].len() as usize) {
            if partition != 0 {
                let flip_bit = ((partition - 1) ^ partition).trailing_ones();
                pattern = convert_current_to_next.flip_from_fast(pattern, flip_bit);
            }
            let update = ls_node[partition];
            if ls_hat[pattern] < update {
                ls_hat[pattern] = update;
                ls_hat_i_arg[pattern] = partition;
            }
        }
        (ls_hat, ls_hat_i_arg)
    } else if i == boundary_paths.len() - 1 {
        let convert_pattern_current_to_prev =
            boundary_paths[i].get_intersection_pattern(&boundary_paths[i - 1]);
        let convert_pattern_current_to_node =
            boundary_paths[i].get_intersection_pattern(&node_paths);
        let (mut prev_pattern, mut node_pattern) = (0, 0);
        let pattern_num = 1 << boundary_paths[i].len();
        (0..(pattern_num as usize))
            .map(|partition| {
                if partition != 0 {
                    let flip_bit = ((partition - 1) ^ partition).trailing_ones();
                    prev_pattern =
                        convert_pattern_current_to_prev.flip_from_fast(prev_pattern, flip_bit);
                    node_pattern =
                        convert_pattern_current_to_node.flip_from_fast(node_pattern, flip_bit);
                }
                (ls_hat[prev_pattern] + ls_node[node_pattern], partition)
            })
            .unzip()
    } else {
        let convert_pattern_current_to_prev =
            boundary_paths[i].get_intersection_pattern(&boundary_paths[i - 1]);
        let convert_pattern_current_to_next =
            boundary_paths[i].get_intersection_pattern(&boundary_paths[i + 1]);
        let convert_pattern_current_to_node =
            boundary_paths[i].get_intersection_pattern(node_paths);
        let intersection_size = boundary_paths[i].count_intersection(&boundary_paths[i + 1]) as u32;
        let mut ls_hat_next = vec![std::f64::NEG_INFINITY; 2usize.pow(intersection_size)];
        let mut ls_hat_arg_next = vec![0; 2usize.pow(intersection_size)];
        let (mut prev_pattern, mut next_pattern, mut node_pattern) = (0, 0, 0);
        for partition in 0..(1 << boundary_paths[i].len()) as usize {
            if partition != 0 {
                let flip_bit = ((partition - 1) ^ partition).trailing_ones();
                prev_pattern =
                    convert_pattern_current_to_prev.flip_from_fast(prev_pattern, flip_bit);
                next_pattern =
                    convert_pattern_current_to_next.flip_from_fast(next_pattern, flip_bit);
                node_pattern =
                    convert_pattern_current_to_node.flip_from_fast(node_pattern, flip_bit);
            }
            let update = ls_hat[prev_pattern] + ls_node[node_pattern];
            if ls_hat_next[next_pattern] < update {
                ls_hat_next[next_pattern] = update;
                ls_hat_arg_next[next_pattern] = partition;
            }
        }
        (ls_hat_next, ls_hat_arg_next)
    }
}

fn fill_next_ls_par(
    boundary_paths: &[PathSet],
    i: usize,
    node_paths: &PathSet,
    ls_node: &[f64],
    ls_hat: &[f64],
) -> (Vec<f64>, Vec<usize>) {
    if i == 0 {
        let (current, next) = match &boundary_paths[i..=i + 1] {
            [x, y] => (x, y),
            _ => panic!(),
        };
        let (free_position, fat_pattern) = get_select_array(current, next);
        let patterns = get_expanded_next_patterns(&fat_pattern);
        let flip_pattern = get_flip_pattern(&free_position);
        let flip_pattern_node = get_node_flip_pattern(current, next, node_paths);
        let convert_pattern_current_to_node = current.get_intersection_pattern(&node_paths);
        patterns
            .par_iter()
            .flat_map(|&free_pattern| {
                let mut pattern = free_pattern;
                let mut node_pattern = convert_pattern_current_to_node.convert(pattern);
                (0..(1 << free_position.len()))
                    .map(|subst: usize| {
                        if subst != 0 {
                            let flipped_bit = ((subst - 1) ^ subst).count_ones() as usize;
                            pattern ^= flip_pattern[flipped_bit];
                            node_pattern ^= flip_pattern_node[flipped_bit];
                        }
                        (ls_node[node_pattern], pattern)
                    })
                    .max_by(|a, b| (a.0).partial_cmp(&b.0).unwrap())
            })
            .unzip()
    } else if i == boundary_paths.len() - 1 {
        let (prev, current) = match &boundary_paths[i - 1..=i] {
            [x, y] => (x, y),
            _ => panic!(),
        };
        // Usual computation is OK.
        let convert_pattern_current_to_prev = current.get_intersection_pattern(prev);
        let convert_pattern_current_to_node = current.get_intersection_pattern(&node_paths);
        let (mut prev_pattern, mut node_pattern) = (0, 0);
        let pattern_num = 1 << current.len();
        (0..(pattern_num as usize))
            .map(|partition| {
                if partition != 0 {
                    let flip_bit = ((partition - 1) ^ partition).trailing_ones();
                    prev_pattern =
                        convert_pattern_current_to_prev.flip_from_fast(prev_pattern, flip_bit);
                    node_pattern =
                        convert_pattern_current_to_node.flip_from_fast(node_pattern, flip_bit);
                }
                (ls_hat[prev_pattern] + ls_node[node_pattern], partition)
            })
            .unzip()
    } else {
        let (prev, current, next) = match &boundary_paths[i - 1..=i + 1] {
            [x, y, z] => (x, y, z),
            _ => panic!(),
        };
        let (free_position, fat_pattern) = get_select_array(current, next);
        let convert_pattern_current_to_prev = current.get_intersection_pattern(prev);
        let convert_pattern_current_to_node = current.get_intersection_pattern(&node_paths);
        // The position of the path not appears in the next paths.
        let patterns = get_expanded_next_patterns(&fat_pattern);
        let flip_pattern = get_flip_pattern(&free_position);
        let flip_pattern_node = get_node_flip_pattern(current, next, node_paths);
        let flip_pattern_prev = get_boundary_flip_pattern(current, next, prev);
        patterns
            .par_iter()
            .map(|&free_pattern| {
                let mut pattern = free_pattern;
                let mut prev_pattern = convert_pattern_current_to_prev.convert(pattern);
                let mut node_pattern = convert_pattern_current_to_node.convert(pattern);
                (0..(1 << free_position.len()))
                    .map(|subst: usize| {
                        if subst != 0 {
                            let flipped_bit = ((subst - 1) ^ subst).count_ones() as usize;
                            pattern ^= flip_pattern[flipped_bit];
                            prev_pattern ^= flip_pattern_prev[flipped_bit];
                            node_pattern ^= flip_pattern_node[flipped_bit];
                        }
                        (ls_hat[prev_pattern] + ls_node[node_pattern], pattern)
                    })
                    .max_by(|a, b| (a.0).partial_cmp(&b.0).unwrap())
                    .unwrap()
            })
            .unzip()
    }
}

fn fill_next_ls_sub(
    boundary_paths: &[PathSet],
    i: usize,
    node_paths: &PathSet,
    ls_node: &[f64],
    ls_hat: Option<&Vec<(usize, f64, usize)>>,
    subsample_size: usize,
) -> Vec<(usize, f64, usize)> {
    if i == 0 {
        let convert_current_to_next =
            boundary_paths[0].get_intersection_pattern(&boundary_paths[1]);
        // partition on R(D_0) and partition on R(D_0) and R(D_1)
        let mut pattern = 0;
        let mut ls_hat: Vec<(usize, f64, usize)> = Vec::new();
        for partition in 0usize..(1 << boundary_paths[0].len() as usize) {
            if partition != 0 {
                let flip_bit = ((partition - 1) ^ partition).trailing_ones();
                pattern = convert_current_to_next.flip_from_fast(pattern, flip_bit);
            }
            let update = ls_node[partition];
            ls_hat.push((pattern, update, partition));
        }
        dedup(&mut ls_hat);
        ls_hat
    } else if i == boundary_paths.len() - 1 {
        let (prev, current) = match &boundary_paths[i - 1..=i] {
            [x, y] => (x, y),
            _ => panic!(),
        };
        let convert_pattern_current_to_node = current.get_intersection_pattern(&node_paths);
        let (free_position, fat_pattern) = get_select_array(current, prev);
        let flip_pattern_node = get_node_flip_pattern(current, prev, node_paths);
        let flip_pattern = get_flip_pattern(&free_position);
        let mut ls_hat_next = vec![];
        ls_hat.unwrap().into_iter().for_each(|&(pattern, lk, _)| {
            // Convert `pattern` (i-1 and i) intersection pattern
            // into gappy i-th pattern.
            let mut pattern_current = get_expanded_next_pattern(&fat_pattern, pattern);
            let mut node_pattern = convert_pattern_current_to_node.convert(pattern_current);
            (0..((1 << free_position.len()) as usize)).for_each(|subst| {
                if subst != 0 {
                    let flipped_bit = ((subst - 1) ^ subst).count_ones() as usize;
                    pattern_current ^= flip_pattern[flipped_bit];
                    node_pattern ^= flip_pattern_node[flipped_bit];
                }
                let lk = lk + ls_node[node_pattern];
                ls_hat_next.push((pattern_current, lk, pattern_current));
            });
        });
        dedup(&mut ls_hat_next);
        ls_hat_next
    } else {
        let (prev, current, next) = match &boundary_paths[i - 1..=i + 1] {
            [x, y, z] => (x, y, z),
            _ => panic!(),
        };
        let (free_position, fat_pattern) = get_select_array(current, prev);
        let flip_pattern = get_flip_pattern(&free_position);
        let flip_pattern_nega = get_flip_pattern(&fat_pattern);
        let convert_pattern_current_to_node = current.get_intersection_pattern(&node_paths);
        let flip_node = get_node_flip_pattern(current, prev, node_paths);
        let flip_node_nega = get_node_flip_pattern_nega(current, prev, node_paths);
        let flip_next = get_boundary_flip_pattern(current, prev, next);
        let flip_next_nega = get_boundary_flip_pattern_nega(current, prev, next);
        // This is the (i-1)-th pattern and lk.
        let convert_pattern_current_to_next = boundary_paths[i].get_intersection_pattern(next);
        let (mut current_pattern, mut node_pattern, mut next_pattern) = (0, 0, 0);
        let mut prev_pattern = std::usize::MAX - 1;
        let mut ls_hat_next = Vec::with_capacity(subsample_size);
        for &(pattern, lk, _) in ls_hat.unwrap() {
            // Convert `pattern` (i-1 and i) intersection pattern
            // into gappy i-th pattern.
            if prev_pattern + 1 != pattern {
                current_pattern = get_expanded_next_pattern(&fat_pattern, pattern);
                node_pattern = convert_pattern_current_to_node.convert(current_pattern);
                next_pattern = convert_pattern_current_to_next.convert(current_pattern);
            } else {
                // Update by flipping...!
                // First, erasing all the 1's in the free position.
                current_pattern ^= flip_pattern[free_position.len()];
                node_pattern ^= flip_node[free_position.len()];
                next_pattern ^= flip_next[free_position.len()];
                // Next, flip bit it should when the i-th bit of the intersection flop.
                let flipped_bit = ((pattern - 1) ^ pattern).count_ones() as usize;
                current_pattern ^= flip_pattern_nega[flipped_bit];
                node_pattern ^= flip_node_nega[flipped_bit];
                next_pattern ^= flip_next_nega[flipped_bit];
            }
            prev_pattern = pattern;
            for subst in 0..((1 << free_position.len()) as usize) {
                if subst != 0 {
                    let flipped_bit = ((subst - 1) ^ subst).count_ones() as usize;
                    current_pattern ^= flip_pattern[flipped_bit];
                    node_pattern ^= flip_node[flipped_bit];
                    next_pattern ^= flip_next[flipped_bit];
                }
                let lk = lk + ls_node[node_pattern];
                ls_hat_next.push((next_pattern, lk, current_pattern));
                if ls_hat_next.len() > 4 * subsample_size {
                    discard_small_elements(&mut ls_hat_next, subsample_size);
                }
            }
        }
        discard_small_elements(&mut ls_hat_next, subsample_size);
        dedup(&mut ls_hat_next);
        ls_hat_next
    }
}

fn discard_small_elements(xs: &mut Vec<(usize, f64, usize)>, num: usize) {
    if num < xs.len() {
        xs.select_nth_unstable_by(num, |a, b| (a.1).partial_cmp(&b.1).unwrap().reverse());
        xs.truncate(num);
    }
}

// Sort given vector by the first argument,
// merge entries with the same first argument into the maximum second value.
fn dedup(xs: &mut Vec<(usize, f64, usize)>) {
    if 1 < xs.len() {
        xs.par_sort_unstable_by_key(|x| x.0);
        let mut pointer = 0;
        let mut current = xs[0];
        let len = xs.len();
        for i in 1..len {
            let x = xs[i];
            if current.0 != x.0 {
                xs[pointer] = current;
                pointer += 1;
                current = x;
            } else if current.0 == x.0 && current.1 < x.1 {
                current = x;
            }
        }
        xs[pointer] = current;
        xs.truncate(pointer + 1);
    }
}

#[derive(Debug, Clone)]
struct HeapElm(f64);

impl std::cmp::PartialEq for HeapElm {
    fn eq(&self, other: &Self) -> bool {
        self.0 <= other.0 && other.0 <= self.0
    }
}

impl std::cmp::Eq for HeapElm {}

impl std::cmp::PartialOrd for HeapElm {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::cmp::Ord for HeapElm {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.0 < other.0 {
            std::cmp::Ordering::Greater
        } else if other.0 < self.0 {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

// Get flipping pattern.
// In other words, the i-th element of the return value would be
// what we want to apply by XOR(^) operation to the node-partition pattern when the `i`-th read in the
// *non*-intersection element is flipped.
fn get_node_flip_pattern(current: &PathSet, next: &PathSet, paths_on_node: &PathSet) -> Vec<usize> {
    let mut flip_pattern = vec![0];
    let mut filled_bit = 0;
    for path_index in current
        .iter()
        .filter(|&&path_index| !next.contains(path_index))
    {
        if let Ok(pointer) = paths_on_node.paths.binary_search(path_index) {
            filled_bit |= 1 << pointer;
        }
        flip_pattern.push(filled_bit);
    }
    flip_pattern
}

// Get flipping pattern.
// In other words, the i-th element of the return value would be
// what we want to apply by XOR(^) operation to the node-partition pattern when the `i`-th read in the
// intersection element is flipped.
fn get_node_flip_pattern_nega(
    current: &PathSet,
    next: &PathSet,
    paths_on_node: &PathSet,
) -> Vec<usize> {
    let mut flip_pattern = vec![0];
    let mut filled_bit = 0;
    for path_index in current
        .iter()
        .filter(|&&path_index| next.contains(path_index))
    {
        if let Ok(pointer) = paths_on_node.paths.binary_search(path_index) {
            filled_bit |= 1 << pointer;
        }
        flip_pattern.push(filled_bit);
    }
    flip_pattern
}

// The i-th element of the return value is what we XOR-ing
// to convert a partition on `flip`
// when the i-th element of the `from` value not in the `to` is flipped.
fn get_boundary_flip_pattern(from: &PathSet, to: &PathSet, flip: &PathSet) -> Vec<usize> {
    let mut flip_pattern = vec![0];
    let mut filled_bit = 0;
    let flip: Vec<_> = flip.iter().copied().filter(|&n| from.contains(n)).collect();
    for path_index in from.iter().filter(|&&path_index| !to.contains(path_index)) {
        if let Ok(pointer) = flip.binary_search(path_index) {
            filled_bit |= 1 << pointer;
        }
        flip_pattern.push(filled_bit);
    }
    flip_pattern
}

// The i-th element of the return value is what we XOR-ing
// to convert a partition on `flip`
// when the i-th element of the intersection of `from` and `to` is flipped.
fn get_boundary_flip_pattern_nega(from: &PathSet, to: &PathSet, flip: &PathSet) -> Vec<usize> {
    let mut flip_pattern = vec![0];
    let mut filled_bit = 0;
    let flip: Vec<_> = flip.iter().copied().filter(|&n| from.contains(n)).collect();
    for path_index in from.iter().filter(|&&path_index| to.contains(path_index)) {
        if let Ok(pointer) = flip.binary_search(path_index) {
            filled_bit |= 1 << pointer;
        }
        flip_pattern.push(filled_bit);
    }
    flip_pattern
}

// If the i-th bit of the absent position read is flipped,
// XOR with the i-th value of the return value.
fn get_flip_pattern(position: &[usize]) -> Vec<usize> {
    let mut flip_pattern = vec![0];
    let mut filled_bit = 0;
    for &pointer in position {
        filled_bit |= 1 << pointer;
        flip_pattern.push(filled_bit);
    }
    flip_pattern
}

// Return value:
// The sum of the length of the two vector is equel to the length of `current`.
// The i-th element of the 1st vector is the position of the `i`-th read in the
// `current` set. Here, I mean the `i`-th is the rank in the intersection between `current` and `next`.
// The i-th elment of the 2nd vector is the position of the `i`-th not-in-the-intersection element in the `current` vector.
fn get_select_array(current: &PathSet, next: &PathSet) -> (Vec<usize>, Vec<usize>) {
    let (mut free_position, mut fat_pattern) = (vec![], vec![]);
    for (idx, &path_index) in current.iter().enumerate() {
        if next.contains(path_index) {
            fat_pattern.push(idx);
        } else {
            free_position.push(idx);
        }
    }
    (free_position, fat_pattern)
}

// Expand the i-th bit of the `from` into `corresponding_location[i]`-th bit.
fn get_expanded_next_pattern(corresponding_location: &[usize], from: usize) -> usize {
    corresponding_location
        .iter()
        .enumerate()
        .fold(0, |acc, (i, location)| {
            acc | (((from >> i) & 0b1) << location)
        })
}

// `corresponding_location`:The i-th element is the location of the i-th bit of the next_pattern in the current pattern.
fn get_expanded_next_patterns(corresponding_location: &[usize]) -> Vec<usize> {
    let mut pattern: usize = 0;
    (0..(1usize << corresponding_location.len()))
        .map(|next_pattern| {
            if next_pattern != 0 {
                let fliped_bit = ((next_pattern - 1) ^ next_pattern).trailing_ones() as usize;
                for &bit in corresponding_location[..fliped_bit].iter() {
                    pattern ^= 1 << bit;
                }
            }
            pattern
        })
        .collect()
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
    let paths: Vec<Vec<_>> = path_indices
        .iter()
        .map(|&path_index| {
            paths[path_index]
                .1
                .iter()
                .filter_map(|&(n, c)| if n == node { Some(c) } else { None })
                .collect()
        })
        .collect();
    let cluster_num = paths.iter().filter_map(|cs| cs.iter().max()).max().unwrap() + 1;
    // let mut bi_counts = [vec![0u32; cluster_num], vec![0u32; cluster_num]];
    let mut bi_counts = vec![0usize; cluster_num * 2];
    for &c in paths.iter().flat_map(|x| x) {
        bi_counts[c] += 1;
    }
    // Note: sum_i u_i*ln(u_i/T) = sum_i u_i*ln(u_i) - u_i*ln(T) = sum_i u_i*ln(u_i) - T*ln(T)
    let mut counts = [bi_counts.iter().sum::<usize>(), 0];
    (0..1 << path_number)
        .map(|pattern: usize| {
            if pattern != 0 {
                let flip_path = ((pattern - 1) ^ pattern).count_ones();
                for (i, path) in paths.iter().take(flip_path as usize).enumerate() {
                    let next_bucket = (pattern >> i) & 0b1;
                    let prev_bucket = next_bucket ^ 0b1;
                    counts[prev_bucket] -= path.len();
                    counts[next_bucket] += path.len();
                    let prev_bucket = prev_bucket * cluster_num;
                    let next_bucket = next_bucket * cluster_num;
                    // Remove from previous bucket, add to next bucket.
                    for &c in path {
                        bi_counts[prev_bucket + c] -= 1;
                        bi_counts[next_bucket + c] += 1;
                    }
                }
            }
            let log_likelihood = bi_counts.iter().map(|&x| xlogx::XLOGX[x]).sum::<f64>();
            log_likelihood - xlogx::XLOGX[counts[0]] - xlogx::XLOGX[counts[1]]
        })
        .collect()
}

// Return xlog_e(x). If x==0, return 0.
#[allow(dead_code)]
fn xlnx(&x: &usize) -> f64 {
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
    fn contains(&self, path_index: usize) -> bool {
        self.paths.binary_search(&path_index).is_ok()
    }
}

/// Bit Convert struct.
/// It behaves like `mask-squash`. In other words, if this struct contains a bit pattern x,
/// then, it convert n into m by picking every i-th bit having 1 in x.
#[derive(Debug, Clone, Default)]
pub struct IntersectPattern {
    len: usize,
    pattern: usize,
    // Return the flipping bits if the i-th bit of the original pattern was flipped.
    flip_pattern: Vec<usize>,
}

impl IntersectPattern {
    fn new(pattern: usize, len: usize) -> Self {
        // First put zero. Mock value.
        let mut flip_pattern = vec![0];
        let (mut pointer, mut filled_bit) = (1, 0);
        for flip_bit in 0..len {
            if (pattern >> flip_bit) & 0b1 == 0b1 {
                filled_bit |= pointer;
                pointer = pointer << 1;
            }
            flip_pattern.push(filled_bit);
        }
        Self {
            pattern,
            len,
            flip_pattern,
        }
    }
    #[allow(dead_code)]
    #[inline]
    fn flip_from(&self, pattern: usize, flip_bit: u32) -> usize {
        (0..flip_bit)
            .fold((pattern, 0), |(pattern, location), i| {
                let probe = (self.pattern >> i) & 0b1;
                (pattern ^ (probe << location), location + probe)
            })
            .0
    }
    #[inline]
    fn flip_from_fast(&self, pattern: usize, flip_bit: u32) -> usize {
        pattern ^ self.flip_pattern[flip_bit as usize]
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
    fn flip_test() {
        let pattern = 0b0111000;
        let conv = IntersectPattern::new(pattern, 7);
        let current = 0b000;
        assert_eq!(conv.flip_from(current, 1), current);
        assert_eq!(conv.flip_from(current, 2), current);
        assert_eq!(conv.flip_from(current, 3), current);
        assert_eq!(conv.flip_from(current, 4), 0b001);
        assert_eq!(conv.flip_from(current, 7), 0b111);
        assert_eq!(conv.flip_from_fast(current, 1), current);
        assert_eq!(conv.flip_from_fast(current, 2), current);
        assert_eq!(conv.flip_from_fast(current, 3), current);
        assert_eq!(conv.flip_from_fast(current, 4), 0b001);
        assert_eq!(conv.flip_from_fast(current, 7), 0b111);
        assert_eq!(conv.flip_from_fast(0b111, 7), 0b000);
    }
    #[test]
    fn get_flip_pattern_test() {
        let position = vec![0, 1, 2, 3, 4];
        assert_eq!(
            get_flip_pattern(&position),
            vec![0, 0b1, 0b11, 0b111, 0b1111, 0b11111]
        );
        let position = vec![0, 1, 3, 5, 6];
        let answer = vec![
            0,
            0b0_000_001,
            0b0_000_011,
            0b0_001_011,
            0b0_101_011,
            0b1_101_011,
        ];
        assert_eq!(get_flip_pattern(&position), answer);
    }
    #[test]
    fn get_select_array_test() {
        let current = PathSet::new(vec![0, 1, 2, 3]);
        let next = PathSet::new(vec![1, 2, 3, 7]);
        let (free, fat) = get_select_array(&current, &next);
        assert_eq!(free, vec![0]);
        assert_eq!(fat, vec![1, 2, 3]);
        let current = PathSet::new(vec![1, 2, 4, 6, 7]);
        let next = PathSet::new(vec![2, 3, 4, 8]);
        let (free, fat) = get_select_array(&current, &next);
        assert_eq!(free, vec![0, 3, 4]);
        assert_eq!(fat, vec![1, 2]);
    }
    #[test]
    fn get_expanded_next_patterns_test() {
        let current = PathSet::new(vec![0, 1, 2, 3]);
        let next = PathSet::new(vec![1, 2, 3, 4]);
        let (_, fat) = get_select_array(&current, &next);
        let pats = get_expanded_next_patterns(&fat);
        let answer = vec![
            0b0000, 0b0010, 0b0100, 0b0110, 0b1000, 0b1010, 0b1100, 0b1110,
        ];
        assert_eq!(pats, answer);
        let current = PathSet::new(vec![1, 2, 4, 6, 7]);
        let next = PathSet::new(vec![2, 6]);
        let (_, fat) = get_select_array(&current, &next);
        let pats = get_expanded_next_patterns(&fat);
        let answer = vec![0b00_000, 0b00_010, 0b01_000, 0b01_010];
        assert_eq!(pats, answer);
    }
    #[test]
    fn get_node_flip_pattern_test() {
        let current = PathSet::new(vec![0, 1, 2, 3, 4]);
        let next = PathSet::new(vec![1, 2, 3, 4, 5]);
        let node = PathSet::new(vec![3, 4]);
        let pattern = get_node_flip_pattern(&current, &next, &node);
        let answer = vec![0, 0b0];
        assert_eq!(pattern, answer);
        let current = PathSet::new(vec![0, 1, 2, 3, 4]);
        let next = PathSet::new(vec![1, 2, 3, 4, 5]);
        let node = PathSet::new(vec![0, 2, 4]);
        let pattern = get_node_flip_pattern(&current, &next, &node);
        let answer = vec![0, 0b001];
        assert_eq!(pattern, answer);
        let current = PathSet::new(vec![0, 1, 3, 5, 6, 7]);
        let next = PathSet::new(vec![1, 3]);
        let node = PathSet::new(vec![0, 1, 5]);
        let pattern = get_node_flip_pattern(&current, &next, &node);
        let answer = vec![0, 0b001, 0b101, 0b101, 0b101];
        assert_eq!(pattern, answer);
    }
    #[test]
    fn get_prev_flip_pattern_test() {
        let prev = PathSet::new(vec![0, 1, 2, 3]);
        let current = PathSet::new(vec![1, 2, 3, 4]);
        let next = PathSet::new(vec![4, 5, 6]);
        let pattern = get_boundary_flip_pattern(&current, &next, &prev);
        let answer = vec![0, 0b001, 0b011, 0b111];
        assert_eq!(answer, pattern);
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
    fn dedup_test() {
        let mut res = vec![(0, 1., 1), (0, 1.1, 2), (0, 0.9, 4)];
        dedup(&mut res);
        assert_eq!(res[0].0, 0);
        assert_eq!(res[0].2, 2);
        assert_eq!(res.len(), 1);
        let mut res = vec![
            (1, 0.9, 2),
            (2, 0.8, 3),
            (2, 0.7, 4),
            (3, 0.9, 4),
            (3, 10., 5),
            (4, 10., 4),
        ];
        dedup(&mut res);
        assert_eq!(res.len(), 4);
        let res: Vec<_> = res.into_iter().map(|(x, _, y)| (x, y)).collect();
        assert_eq!(res, vec![(1, 2), (2, 3), (3, 5), (4, 4)]);
    }
}
