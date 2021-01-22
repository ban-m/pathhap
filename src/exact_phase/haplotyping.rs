use super::graph_traversal::*;
pub fn haplotype_cc(paths: &[Vec<(usize, usize)>], max_len: usize) -> (Vec<u8>, f64) {
    debug!("start.");
    let num_of_nodes = number_of_nodes(paths);
    let node_traverse_order = determine_traversal_order(num_of_nodes, paths);
    let split_paths = downsampling_up_to(&node_traverse_order, paths, max_occ);
    let boundary_nodes = get_boundary_nodes_on(&node_traverse_order, &split_paths);
    let haplotyping_regions = ;
    (vec![], 0.)
}
