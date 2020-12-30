/// Phasing API.
pub fn phase(paths: &[(String, Vec<(u64, u64)>)], _max_occ: usize) -> Vec<(String, u8)> {
    paths.iter().map(|(id, _)| (id.clone(), 0)).collect()
}
