use dict::{Felt252Dict, Felt252DictTrait};
use nullable::{NullableTrait, nullable_from_box, match_nullable, FromNullableResult};
use yas_core::numbers::signed_integer::{i128::{i128, u128Intoi128}, integer_trait::IntegerTrait};

#[derive(Destruct)]
struct Matrix {
    data: Felt252Dict<Nullable<i128>>,
    rows: u8,
    cols: u8,
}

mod errors {
    const INVALID_DIMENSION: felt252 = 'Matrix: invalid dimension';
    const INVALID_MATRIX_INVERSION: felt252 = 'Matrix: matrix not invertible';
}

trait MatrixTrait {
    fn new(rows: u8, cols: u8) -> Matrix;

    fn get(ref self: Matrix, row: u8, col: u8) -> i128;

    fn set(ref self: Matrix, row: u8, col: u8, value: i128) -> ();

    fn transpose(ref self: Matrix) -> Matrix;

    fn minor(ref self: Matrix, exclude_row: u8, exclude_col: u8) -> Matrix;

    fn det(ref self: Matrix) -> i128;

    fn inv(ref self: Matrix) -> Matrix;
}

impl MatrixImpl of MatrixTrait {
    fn new(rows: u8, cols: u8) -> Matrix {
        Matrix { data: Default::default(), rows, cols }
    }

    fn get(ref self: Matrix, row: u8, col: u8) -> i128 {
        let key: u8 = row * self.cols + col;
        match match_nullable(self.data.get(key.into())) {
            FromNullableResult::Null => 0_u128.into(),
            FromNullableResult::NotNull(value) => value.unbox(),
        }
    }

    fn set(ref self: Matrix, row: u8, col: u8, value: i128) {
        let key: u8 = row * self.cols.into() + col;
        self.data.insert(key.into(), nullable_from_box(BoxTrait::new(value)));
    }

    fn transpose(ref self: Matrix) -> Matrix {
        let mut result = MatrixTrait::new(self.cols, self.rows);
        let max_index = self.rows * self.cols;
        let mut index: u8 = 0;
        loop {
            if index == max_index {
                break;
            }
            let row = index / self.cols;
            let col = index % self.cols;
            result.set(col, row, self.get(row, col));
            index += 1;
        };
        result
    }

    fn minor(ref self: Matrix, exclude_row: u8, exclude_col: u8) -> Matrix {
        let mut minor_matrix = MatrixTrait::new(self.rows - 1, self.cols - 1);

        let mut index: u8 = 0;
        let max_index: u8 = self.rows * self.cols;
        loop {
            if index >= max_index {
                break;
            };

            let row = index / self.cols;
            let col = index % self.cols;

            if row != exclude_row && col != exclude_col {
                let row_offset = if row > exclude_row {
                    1
                } else {
                    0
                };
                let col_offset = if col > exclude_col {
                    1
                } else {
                    0
                };

                let value = self.get(row, col);
                minor_matrix.set(row - row_offset, col - col_offset, value);
            };

            index += 1;
        };

        minor_matrix
    }

    fn det(ref self: Matrix) -> i128 {
        // [Check] Matrix is square
        assert(self.rows == self.cols, errors::INVALID_DIMENSION);
        if self.rows == 1 {
            return self.get(0, 0);
        }

        if self.rows == 2 {
            return (self.get(0, 0) * self.get(1, 1)) - (self.get(0, 1) * self.get(1, 0));
        }

        let mut det: i128 = 0_u128.into();
        let mut col: u8 = 0;
        loop {
            if col >= self.cols {
                break;
            }

            let coef = self.get(0, col);
            let mut minor = self.minor(0, col);
            if col % 2 == 0 {
                det += coef * minor.det();
            } else {
                det -= coef * minor.det();
            };

            col += 1;
        };

        return det;
    }

    fn inv(ref self: Matrix) -> Matrix {
        let determinant = self.det();
        assert(determinant != 0_u128.into(), errors::INVALID_MATRIX_INVERSION);

        let mut inverse_matrix = MatrixTrait::new(self.rows, self.cols);
        let max_index = self.rows * self.cols;
        let mut index: u8 = 0;
        loop {
            if index == max_index {
                break;
            }
            let row = index / self.cols;
            let col = index % self.cols;

            // Cofactor computation
            let sign = if (row + col) % 2 == 0 {
                1_u128.into()
            } else {
                -1
            };
            let mut minor = self.minor(row, col);
            let cofactor = if (row + col) % 2 == 0 {
                minor.det()
            } else {
                -minor.det()
            };
            inverse_matrix.set(row, col, cofactor / determinant);

            index += 1;
        };
        inverse_matrix
    }
}

impl MatrixAdd of Add<Matrix> {
    fn add(mut lhs: Matrix, mut rhs: Matrix) -> Matrix {
        // [Check] Dimesions are compatible
        assert(lhs.rows == rhs.rows && lhs.cols == rhs.cols, errors::INVALID_DIMENSION);
        let mut result = MatrixTrait::new(lhs.rows, lhs.cols);
        let max_index = lhs.rows * lhs.cols;
        let mut index = 0;
        loop {
            if index == max_index {
                break;
            }
            let row = index / lhs.cols;
            let col = index % lhs.cols;
            let value = lhs.get(row, col) + rhs.get(row, col);
            result.set(row, col, value);
            index += 1;
        };
        result
    }
}

impl MatrixSub of Sub<Matrix> {
    fn sub(mut lhs: Matrix, mut rhs: Matrix) -> Matrix {
        // [Check] Dimesions are compatible
        assert(lhs.rows == rhs.rows && lhs.cols == rhs.cols, errors::INVALID_DIMENSION);
        let mut result = MatrixTrait::new(lhs.rows, lhs.cols);
        let max_index = lhs.rows * lhs.cols;
        let mut index = 0;
        loop {
            if index == max_index {
                break;
            }
            let row = index / lhs.cols;
            let col = index % lhs.cols;
            let value = lhs.get(row, col) - rhs.get(row, col);
            result.set(row, col, value);
            index += 1;
        };
        result
    }
}

impl MatrixMul of Mul<Matrix> {
    fn mul(mut lhs: Matrix, mut rhs: Matrix) -> Matrix {
        // [Check] Dimesions are compatible
        assert(lhs.cols == rhs.rows, errors::INVALID_DIMENSION);
        let mut result = MatrixTrait::new(lhs.rows, rhs.cols);
        let max_index = lhs.rows * rhs.cols;
        let mut index: u8 = 0;
        loop {
            if index == max_index {
                break;
            }

            let row = index / rhs.cols;
            let col = index % rhs.cols;

            let mut sum: i128 = 0_u128.into();
            let mut k: u8 = 0;
            loop {
                if k == lhs.cols {
                    break;
                }

                sum += lhs.get(row, k) * rhs.get(k, col);
                k += 1;
            };

            result.set(row, col, sum);
            index += 1;
        };

        result
    }
}


#[cfg(test)]
mod tests {
    use super::{Matrix, MatrixTrait, errors};
    use debug::PrintTrait;

    #[test]
    #[available_gas(1_000_000)]
    fn test_matrix_get_set() {
        let rows: u8 = 3;
        let cols: u8 = 4;
        let mut matrix: Matrix = MatrixTrait::new(rows, cols);
        matrix.set(0, 1, 100_u128.into());
        matrix.set(2, 3, 50_u128.into());
        assert(matrix.get(0, 1) == 100_u128.into(), 'Matrix: get or set failed');
        assert(matrix.get(2, 3) == 50_u128.into(), 'Matrix: get or set failed');
    }

    #[test]
    #[available_gas(1_000_000)]
    fn test_matrix_transpose() {
        let rows: u8 = 2;
        let cols: u8 = 3;
        let mut matrix = MatrixTrait::new(rows, cols);
        matrix.set(0, 2, 100_u128.into());
        matrix.set(1, 1, 50_u128.into());
        let mut transposed = matrix.transpose();
        assert(transposed.get(2, 0) == 100_u128.into(), 'Matrix: transpose failed');
        assert(transposed.get(1, 1) == 50_u128.into(), 'Matrix: transpose failed');
    }

    #[test]
    #[available_gas(1_000_000)]
    fn test_matrix_addition() {
        let rows: u8 = 2;
        let cols: u8 = 3;
        let mut matrix1 = MatrixTrait::new(rows, cols);
        let mut matrix2 = MatrixTrait::new(rows, cols);
        matrix1.set(0, 0, 10_u128.into());
        matrix1.set(1, 1, 20_u128.into());
        matrix2.set(0, 0, 5_u128.into());
        matrix2.set(1, 1, 15_u128.into());
        let mut result = matrix1 + matrix2;
        assert(result.get(0, 0) == 15_u128.into(), 'Matrix: addition failed');
        assert(result.get(1, 1) == 35_u128.into(), 'Matrix: addition failed');
    }

    #[test]
    #[available_gas(1_000_000)]
    fn test_matrix_subtraction() {
        let rows: u8 = 2;
        let cols: u8 = 3;
        let mut matrix1 = MatrixTrait::new(rows, cols);
        let mut matrix2 = MatrixTrait::new(rows, cols);
        matrix1.set(0, 0, 20_u128.into());
        matrix1.set(1, 1, 30_u128.into());
        matrix2.set(0, 0, 5_u128.into());
        matrix2.set(1, 1, 10_u128.into());
        let mut result = matrix1 - matrix2;
        assert(result.get(0, 0) == 15_u128.into(), 'Matrix: subtraction failed');
        assert(result.get(1, 1) == 20_u128.into(), 'Matrix: subtraction failed');
    }

    #[test]
    #[available_gas(10_000_000)]
    fn test_matrix_square_multiplication() {
        let size: u8 = 2;
        let mut matrix1 = MatrixTrait::new(size, size);
        let mut matrix2 = MatrixTrait::new(size, size);
        matrix1.set(0, 0, 1_u128.into());
        matrix1.set(0, 1, 2_u128.into());
        matrix1.set(1, 0, 3_u128.into());
        matrix1.set(1, 1, 4_u128.into());
        matrix2.set(0, 0, 2_u128.into());
        matrix2.set(0, 1, 0_u128.into());
        matrix2.set(1, 0, 1_u128.into());
        matrix2.set(1, 1, 3_u128.into());
        let mut result = matrix1 * matrix2;
        assert(result.get(0, 0) == 4_u128.into(), 'Matrix: multiplication failed');
        assert(result.get(0, 1) == 6_u128.into(), 'Matrix: multiplication failed');
        assert(result.get(1, 0) == 10_u128.into(), 'Matrix: multiplication failed');
        assert(result.get(1, 1) == 12_u128.into(), 'Matrix: multiplication failed');
    }

    #[test]
    #[available_gas(10_000_000)]
    fn test_matrix_rectangle_multiplication() {
        let mut matrix1 = MatrixTrait::new(2, 3);
        let mut matrix2 = MatrixTrait::new(3, 2);
        matrix1.set(0, 0, 1_u128.into());
        matrix1.set(0, 1, 2_u128.into());
        matrix1.set(0, 2, 3_u128.into());
        matrix1.set(1, 0, 4_u128.into());
        matrix1.set(1, 1, 5_u128.into());
        matrix1.set(1, 2, 6_u128.into());
        matrix2.set(0, 0, 7_u128.into());
        matrix2.set(0, 1, 8_u128.into());
        matrix2.set(1, 0, 9_u128.into());
        matrix2.set(1, 1, 10_u128.into());
        matrix2.set(2, 0, 11_u128.into());
        matrix2.set(2, 1, 12_u128.into());
        let mut result = matrix1 * matrix2;
        assert(result.get(0, 0) == 58_u128.into(), 'Matrix: multiplication failed');
        assert(result.get(0, 1) == 64_u128.into(), 'Matrix: multiplication failed');
        assert(result.get(1, 0) == 139_u128.into(), 'Matrix: multiplication failed');
        assert(result.get(1, 1) == 154_u128.into(), 'Matrix: multiplication failed');
    }

    #[test]
    #[available_gas(5_000_000)]
    fn test_matrix_determinant_2x2() {
        let mut matrix = MatrixTrait::new(2, 2);
        matrix.set(0, 0, 4_u128.into());
        matrix.set(0, 1, 3_u128.into());
        matrix.set(1, 0, 1_u128.into());
        matrix.set(1, 1, 2_u128.into());
        assert(matrix.det() == 5_u128.into(), 'Matrix: det computation failed');
    }

    #[test]
    #[available_gas(10_000_000)]
    fn test_matrix_determinant_3x3() {
        let mut matrix = MatrixTrait::new(3, 3);
        matrix.set(0, 0, 6_u128.into());
        matrix.set(0, 1, 1_u128.into());
        matrix.set(0, 2, 1_u128.into());
        matrix.set(1, 0, 4_u128.into());
        matrix.set(1, 1, -2_u128.into());
        matrix.set(1, 2, 5_u128.into());
        matrix.set(2, 0, 2_u128.into());
        matrix.set(2, 1, 8_u128.into());
        matrix.set(2, 2, 7_u128.into());
        let det = matrix.det();
        assert(matrix.det() == -306_u128.into(), 'Matrix: det computation failed');
    }

    #[test]
    #[available_gas(10_000_000)]
    fn test_matrix_inverse_2x2() {
        let mut matrix = MatrixTrait::new(2, 2);
        matrix.set(0, 0, 1_u128.into());
        matrix.set(0, 1, 2_u128.into());
        matrix.set(1, 0, 3_u128.into());
        matrix.set(1, 1, 4_u128.into());
        let mut inverse = matrix.inv();
        assert(inverse.get(0, 0) == -2_u128.into(), 'Matrix: inversion failed');
        assert(inverse.get(0, 1) == 1_u128.into(), 'Matrix: inversion failed');
        assert(inverse.get(1, 0) == 1_u128.into(), 'Matrix: inversion failed');
        assert(inverse.get(1, 1) == 0_u128.into(), 'Matrix: inversion failed');
    }
}
