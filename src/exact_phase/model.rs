use log::debug;
use std::collections::HashMap;
/// A statistical model for haplotyping. This struct is only to calculate the entire likelihood,
/// or predict short reads.
#[derive(Debug, Clone, Default)]
pub struct Model {
    // haplotype -> position -> distribution.
    haps: Vec<Vec<Vec<f64>>>,
}

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

const BIG_SMALL: f64 = -100000000000000000000000000.;
impl Model {
    const PLOIDY: usize = 2;
    pub fn likelihood(&self, path: &[(usize, usize)]) -> f64 {
        let lks: Vec<_> = self
            .haps
            .iter()
            .map(|hap| {
                path.iter()
                    .map(|&(n, c)| {
                        if hap[n][c] < 0.000001 {
                            BIG_SMALL
                        } else {
                            hap[n][c].ln()
                        }
                    })
                    .sum::<f64>()
            })
            .collect();
        if lks.iter().any(|x| x.is_nan()) {
            eprintln!("{:?}\t{:?}", lks, path);
            eprintln!("{:?}", self);
        }
        logsumexp(&lks)
    }
    // Predict short path
    pub fn predict_path(&self, path: &[(usize, usize)]) -> u8 {
        self.haps
            .iter()
            .map(|hap| path.iter().map(|&(n, c)| hap[n][c].ln()).sum::<f64>())
            .enumerate()
            .map(|(x, y)| (x as u8, y))
            .max_by(|x, y| match x.1.partial_cmp(&(y.1)) {
                Some(x) => x,
                None => {
                    for (i, hap) in self.haps.iter().enumerate() {
                        for (pos, frac) in hap.iter().enumerate() {
                            let frac: Vec<_> = frac.iter().map(|x| format!("{:.3}", x)).collect();
                            debug!("{}\t{}\t[{}]", i, pos, frac.join(","));
                        }
                    }
                    panic!("COMP\t({},{})-({},{})", x.0, x.1, y.0, y.1)
                }
            })
            .unwrap()
            .0
    }
    pub fn new(
        paths: &[(usize, Vec<(usize, usize)>)],
        haplotypes: &[(usize, u8)],
        pseudo_count: u32,
    ) -> Self {
        let haplotypes: HashMap<usize, u8> = haplotypes.iter().copied().collect();
        let nodes_and_clusters = {
            let num_of_nodes = paths
                .iter()
                .flat_map(|(_, x)| x.iter().map(|&(n, _)| n))
                .max()
                .unwrap()
                + 1;
            let mut clusters = vec![0; num_of_nodes];
            for (_, path) in paths.iter() {
                for &(n, c) in path.iter() {
                    clusters[n] = clusters[n].max(c + 1);
                }
            }
            clusters
        };
        let mut hap_count: Vec<Vec<Vec<u32>>> = (0..Self::PLOIDY)
            .map(|_| {
                nodes_and_clusters
                    .iter()
                    .map(|&cl| vec![pseudo_count; cl])
                    .collect()
            })
            .collect();
        for (r, path) in paths.iter() {
            assert!(haplotypes.contains_key(r));
            let index = haplotypes[r] as usize;
            assert!(index < 2);
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
                        if sum > 0 {
                            xs.iter().map(|&x| x as f64 / sum as f64).collect()
                        } else {
                            vec![0.; xs.len()]
                        }
                    })
                    .collect()
            })
            .collect();
        Self { haps }
    }
}
