use ark_ff::FftField;

pub fn compute_roots_of_unity<F: FftField>(order: u64) -> Result<Vec<F>, String> {
    if !order.is_power_of_two() {
        return Err("Order must be a power of two".to_string());
    }

    let root = F::get_root_of_unity(order).unwrap();

    let mut roots = Vec::with_capacity(order as usize);
    let mut current = F::one();
    for _ in 0..order {
        roots.push(current);
        current *= root;
    }
    Ok(roots)
}
