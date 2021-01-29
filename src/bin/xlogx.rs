fn main() {
    let xs: Vec<_> = (0..50)
        .map(|x| {
            if x == 0 {
                0.
            } else {
                let x = x as f64;
                x * x.ln()
            }
        })
        .collect();
    let xs: Vec<_> = xs.iter().map(|x| format!("{}", x)).collect();
    println!("const XLOGX:[f64;50] = [{}];", xs.join(","));
}
