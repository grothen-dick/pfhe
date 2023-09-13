pub fn modular_inverse(num1: u32, num2: u32) -> i64 {
    let (mut x0, mut x1) = if num1 > num2 {
        (num1, num2)
    } else {
        (num2, num1)
    };
    let (mut y0, mut y1): (i64, i64) = (0, 1);
    while x1 > 0 {
        let q = x0 / x1;
        (x0, x1) = (x1, x0 - q * x1);
        (y0, y1) = (y1, y0 - i64::from(q) * y1);
    }
    y0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_finds_modular_inverse() {
        let inv_x = modular_inverse(37, 5);
        assert_eq!(inv_x, 15);
    }

    #[test]
    fn it_finds_modular_inverse_wrong_order() {
        let inv_x = modular_inverse(374533, 52343);
        assert_eq!(inv_x, 142177);
    }
}
