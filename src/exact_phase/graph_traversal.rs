use log::debug;
use std::collections::{HashSet, VecDeque};
pub fn number_of_nodes(paths: &[Vec<(usize, usize)>]) -> usize {
    paths
        .iter()
        .filter_map(|x| x.iter().map(|x| x.0).max())
        .max()
        .unwrap()
        + 1
}

// BFS and return the order.
// We choose the un-traversed smallest children at each step.
pub fn determine_traversal_order(num_of_nodes: usize, paths: &[Vec<(usize, usize)>]) -> Vec<usize> {
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
            eds.sort_unstable();
            eds
        })
        .collect();
    traverse_graph(num_of_nodes, &edges)
    //bfs(num_of_nodes, &edges)
}

// Traverse graph so that the maximum number of the boundary would be as small as possible.
// To this end, we first
fn traverse_graph(num_of_nodes: usize, edges: &[Vec<usize>]) -> Vec<usize> {
    let bridges: Vec<Vec<bool>> = enumerate_bridges(num_of_nodes, edges);
    // DFS. Edge would be selected, first non-bridge edge with shortest return path,
    // second any bridging edge.
    let mut is_used = vec![false; num_of_nodes];
    let mut stack = vec![get_start_node(num_of_nodes, edges)];
    let mut order = stack.clone();
    'dfs: while !stack.is_empty() {
        let last = *stack.last().unwrap();
        is_used[last] = true;
        let next_non_bridge_node = edges[last]
            .iter()
            .zip(bridges[last].iter())
            .filter(|&(&to, b)| !b && !is_used[to]) // non-bridge non-traversed.
            .map(|(&to, _)| (to, min_return_time(&edges, &is_used, to, last)))
            .min_by(|(x, x_ret), (y, y_ret)| match x_ret.cmp(y_ret) {
                std::cmp::Ordering::Equal => x.cmp(y),
                other => other,
            });
        if let Some((next, _)) = next_non_bridge_node {
            order.push(next);
            stack.push(next);
            continue 'dfs;
        }
        let next_bridge_node = edges[last]
            .iter()
            .zip(bridges[last].iter())
            .find(|&(&to, &b)| b && !is_used[to]); // bridge non-arived node
        if let Some((&next, _)) = next_bridge_node {
            stack.push(next);
            order.push(next);
            continue 'dfs;
        }
        stack.pop().unwrap();
    }
    // Sanity check.
    let mut is_used = vec![false; num_of_nodes];
    debug!("{:?}", order);
    for &n in order.iter() {
        assert!(!is_used[n]);
        is_used[n] = true;
    }
    assert!(is_used.iter().all(|&b| b));
    order
}

// starting from non-used node `start`, traverse the graph and return the distance from
fn min_return_time(edges: &[Vec<usize>], is_used: &[bool], start: usize, parent: usize) -> u32 {
    let mut is_arrived = is_used.to_vec();
    let mut queue: Vec<_> = edges[start]
        .iter()
        .filter(|&&to| to != parent)
        .copied()
        .collect();
    let mut distance = 0;
    while !queue.is_empty() {
        let mut next_queue = vec![];
        for &to in queue.iter() {
            if is_used[to] {
                // Returned to the used node.
                break;
            } else if !is_arrived[to] {
                is_arrived[to] = true;
                next_queue.push(to);
            }
        }
        queue = next_queue;
        distance += 1;
    }
    distance
}

fn get_start_node(_: usize, edges: &[Vec<usize>]) -> usize {
    edges
        .iter()
        .enumerate()
        .find(|(_, x)| x.len() == 1)
        .map(|x| x.0)
        .unwrap_or(0)
}

fn enumerate_bridges(num_of_nodes: usize, edges: &[Vec<usize>]) -> Vec<Vec<bool>> {
    let mut is_used = vec![false; num_of_nodes];
    let mut is_edge_used: Vec<_> = edges.iter().map(|x| vec![false; x.len()]).collect();
    let mut stack = vec![0];
    let mut order = vec![0; num_of_nodes];
    let mut low_link = vec![0; num_of_nodes];
    let mut iteration = 0;
    let mut parent = vec![0; num_of_nodes];
    'dfs: while !stack.is_empty() {
        let last = *stack.last().unwrap();
        is_used[last] = true;
        for (i, &to) in edges[last].iter().enumerate() {
            if !is_used[to] {
                stack.push(to);
                is_edge_used[last][i] = true;
                iteration += 1;
                order[to] = iteration;
                parent[to] = last;
                continue 'dfs;
            }
        }
        let last = stack.pop().unwrap();
        let ll = edges[last]
            .iter()
            .zip(is_edge_used[last].iter())
            .filter(|&(&to, _)| to != parent[last])
            .map(|(&to, &is_used)| if is_used { low_link[to] } else { order[to] })
            .min();
        low_link[last] = match ll {
            Some(x) => x.min(order[last]),
            None => order[last],
        };
    }
    edges
        .iter()
        .enumerate()
        .map(|(from, eds)| {
            eds.iter()
                .map(|&to| (order[from] < low_link[to]) | (order[to] < low_link[from]))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}

// Breadth first search and returns the arriving order.
#[allow(dead_code)]
fn bfs(num_of_nodes: usize, edges: &[Vec<usize>]) -> Vec<usize> {
    let mut order = vec![];
    let mut is_arrived = vec![false; num_of_nodes];
    let mut queue = VecDeque::new();
    let start_node = get_start_node(num_of_nodes, edges);
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

pub fn get_boundary_nodes_on(order: &[usize], paths: &[(usize, &[(usize, usize)])]) -> Vec<Nodes> {
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
pub fn get_boundary_nodes(
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
#[derive(Debug, Clone, Default)]
pub struct Nodes {
    pub nodes: HashSet<usize>,
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
mod test {
    use super::*;
    use std::collections::HashMap;
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
        assert!(order == vec![1, 0, 3, 4, 2, 5, 6]);
    }
    #[test]
    fn bridge_test() {
        let num_of_nodes = 4;
        let edges = vec![vec![1], vec![0, 2, 3], vec![1], vec![1]];
        let bridges = enumerate_bridges(num_of_nodes, &edges);
        let answer = vec![vec![true], vec![true, true, true], vec![true], vec![true]];
        assert_eq!(bridges, answer);
    }
    #[test]
    fn bridge_test_2() {
        let num_of_nodes = 6;
        let edges = vec![
            vec![1],
            vec![0, 2, 3],
            vec![1, 4],
            vec![1, 4],
            vec![2, 3, 5],
            vec![4],
        ];
        let bridges = enumerate_bridges(num_of_nodes, &edges);
        let answer = vec![
            vec![true],
            vec![true, false, false],
            vec![false, false],
            vec![false, false],
            vec![false, false, true],
            vec![true],
        ];
        assert_eq!(bridges, answer);
    }
    #[test]
    fn bridge_test_3() {
        let num_of_nodes = 8;
        let edges = vec![
            vec![1, 2],
            vec![0, 5],
            vec![0, 3, 4, 5],
            vec![2, 4],
            vec![2, 3],
            vec![1, 2, 6, 7],
            vec![5],
            vec![5],
        ];
        let bridges = enumerate_bridges(num_of_nodes, &edges);
        let answer = vec![
            vec![false, false],
            vec![false, false],
            vec![false, false, false, false],
            vec![false, false],
            vec![false, false],
            vec![false, false, true, true],
            vec![true],
            vec![true],
        ];
        assert_eq!(bridges, answer);
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
}
