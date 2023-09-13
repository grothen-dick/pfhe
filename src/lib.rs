pub fn extended_euclidean_algorithm(num1: i32, num2: i32) -> (i32, i32) {
    let (mut x0, mut x1) = if num1 > num2 {
        (num1, num2)
    } else {
        (num2, num1)
    };
    let (mut y0, mut y1) = (0, 1);
    while x1 > 0 {
        let q = x0 / x1;
        let mut temp = x0;
        x0 = x1;
        x1 = temp - q * x1;
        temp = y0;
        y0 = y1;
        y1 = temp - q * y1;
    }
    (y0, x0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_finds_modular_inverse() {
        let (inv_x, gcd) = extended_euclidean_algorithm(37, 5);
        assert_eq!(inv_x, 15);
        assert_eq!(gcd, 1);
    }

    #[test]
    fn it_finds_modular_inverse_wrong_order() {
        let (inv_x, gcd) = extended_euclidean_algorithm(374533, 52343);
        assert_eq!(inv_x, 142177);
        assert_eq!(gcd, 1);
    }
}
