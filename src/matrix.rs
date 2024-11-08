use std::{
    fmt,
    ops::{Index, IndexMut},
};

use crate::{
    letter::Letter,
    score::{score_digit_count, Score},
};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AlignmentMatrix {
    buf: Vec<Score>,
    width: usize,
}

impl AlignmentMatrix {
    pub fn zeroed(height: usize, width: usize) -> Self {
        Self { buf: vec![0; height * width], width }
    }

    pub fn height(&self) -> usize {
        self.buf.len() / self.width()
    }

    pub fn width(&self) -> usize {
        self.width
    }

    fn pack_index(&self, i: usize, j: usize) -> Option<usize> {
        if j >= self.width {
            None
        } else {
            Some(i * self.width + j)
        }
    }

    fn unpack_index(&self, index: usize) -> (usize, usize) {
        (index / self.width, index % self.width)
    }

    pub fn get_ref(&self, i: usize, j: usize) -> Option<&Score> {
        let packed_index = self.pack_index(i, j)?;
        self.buf.get(packed_index)
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> Option<&mut Score> {
        let packed_index = self.pack_index(i, j)?;
        self.buf.get_mut(packed_index)
    }

    pub fn get(&self, i: usize, j: usize) -> Option<Score> {
        self.get_ref(i, j).copied()
    }

    #[must_use]
    pub fn set(&mut self, i: usize, j: usize, score: Score) -> bool {
        if let Some(ref_mut) = self.get_mut(i, j) {
            *ref_mut = score;
            true
        } else {
            false
        }
    }

    pub fn max(&self) -> Option<Score> {
        self.buf.iter().copied().max()
    }

    pub fn min(&self) -> Option<Score> {
        self.buf.iter().copied().min()
    }

    pub fn argmax(&self) -> Option<(usize, usize)> {
        self.buf
            .iter()
            .copied()
            .enumerate()
            .max_by_key(|(_, value)| *value)
            .map(|(k, _)| self.unpack_index(k))
    }

    pub fn argmin(&self) -> Option<(usize, usize)> {
        self.buf
            .iter()
            .copied()
            .enumerate()
            .min_by_key(|(_, value)| *value)
            .map(|(k, _)| self.unpack_index(k))
    }

    pub fn argmax_many(&self) -> Vec<(usize, usize)> {
        let maybe_max = self.max();
        if let Some(max) = maybe_max {
            self.buf
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, value)| *value == max)
                .map(|(k, _)| self.unpack_index(k))
                .collect()
        } else {
            Vec::new()
        }
    }
}

impl Index<(usize, usize)> for AlignmentMatrix {
    type Output = Score;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        let height = self.height();
        let width = self.width();
        self.get_ref(i, j).unwrap_or_else(|| invalid_index(i, j, height, width))
    }
}

impl IndexMut<(usize, usize)> for AlignmentMatrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        let height = self.height();
        let width = self.width();
        self.get_mut(i, j).unwrap_or_else(|| invalid_index(i, j, height, width))
    }
}

impl Index<[usize; 2]> for AlignmentMatrix {
    type Output = Score;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self[(index[0], index[1])]
    }
}

impl IndexMut<[usize; 2]> for AlignmentMatrix {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self[(index[0], index[1])]
    }
}

#[cold]
#[inline(never)]
fn invalid_index(i: usize, j: usize, height: usize, width: usize) -> ! {
    panic!(
        "invalid indices [{i}, {j}] for matrix dimensions [{height}, {width}]",
    )
}

fn index_digit_count(k: usize) -> u32 {
    if k > 0 {
        k.ilog10() + 1
    } else {
        1
    }
}

/**
 * Example:
```text
matrix 5x4
 |0 |1 |2 |3 |
-|==+==+==+==|
0| 0|-2|-4|-6|
-|--+--+--+--|
1|-2| 1|-1|-3|
-|--+--+--+--|
2|-4|-1| 2| 0|
-|--+--+--+--|
3|-6|-3| 0| 1|
-|--+--+--+--|
4|-8|-5|-2|-1|
-|==+==+==+==|
```
 */
#[derive(Debug, Clone, Copy)]
pub struct PrettyPrint<'a>(pub &'a AlignmentMatrix);

impl fmt::Display for PrettyPrint<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self(matrix) = self;
        write!(f, "matrix {}x{}\n", matrix.height(), matrix.width())?;
        let Some(min_score) = matrix.min() else {
            return Ok(());
        };
        let Some(max_score) = matrix.max() else {
            return Ok(());
        };
        let score_max_digits =
            score_digit_count(min_score).max(score_digit_count(max_score));
        let height_max_digits = index_digit_count(matrix.height());
        let width_max_digits = index_digit_count(matrix.width());
        let max_digits = score_max_digits.max(width_max_digits);
        for _ in 0 .. height_max_digits {
            write!(f, " ")?;
        }
        write!(f, "|")?;
        for j in 0 .. matrix.width() {
            write!(f, "{:<textwidth$}|", j, textwidth = max_digits as usize)?;
        }
        write!(f, "\n")?;
        for i in 0 .. matrix.height() {
            for _ in 0 .. height_max_digits {
                write!(f, "-")?;
            }
            for w in 0 .. matrix.width() {
                if w == 0 {
                    write!(f, "|")?;
                } else {
                    write!(f, "+")?;
                }
                for _ in 0 .. max_digits {
                    if i == 0 {
                        write!(f, "=")?;
                    } else {
                        write!(f, "-")?;
                    }
                }
            }
            write!(f, "|\n")?;
            write!(
                f,
                "{:<textwidth$}|",
                i,
                textwidth = height_max_digits as usize
            )?;
            for j in 0 .. matrix.width() {
                write!(
                    f,
                    "{:>textwidth$}|",
                    matrix[[i, j]],
                    textwidth = max_digits as usize
                )?;
            }
            write!(f, "\n")?;
        }
        for _ in 0 .. height_max_digits {
            write!(f, "-")?;
        }
        for w in 0 .. matrix.width() {
            if w == 0 {
                write!(f, "|")?;
            } else {
                write!(f, "+")?;
            }
            for _ in 0 .. max_digits {
                write!(f, "=")?;
            }
        }
        write!(f, "|\n")?;
        Ok(())
    }
}

/**
 * Example:
```text
matrix 5x4
   |0 |1 |2 |3 |
   |--+--+--+--|
   |  |W |H |Y |
-+-|==+==+==+==|
0| | 0|-2|-4|-6|
-+-|--+--+--+--|
1|W|-2| 1|-1|-3|
-+-|--+--+--+--|
2|H|-4|-1| 2| 0|
-+-|--+--+--+--|
3|A|-6|-3| 0| 1|
-+-|--+--+--+--|
4|T|-8|-5|-2|-1|
-+-|==+==+==+==|
```
 */
#[derive(Debug, Clone, Copy)]
pub struct LabeledPrettyPrint<'a>(
    pub &'a AlignmentMatrix,
    pub &'a [Letter],
    pub &'a [Letter],
);

impl fmt::Display for LabeledPrettyPrint<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self(matrix, row_seq, col_seq) = self;
        write!(f, "matrix {}x{}\n", matrix.height(), matrix.width())?;
        let Some(min_score) = matrix.min() else {
            return Ok(());
        };
        let Some(max_score) = matrix.max() else {
            return Ok(());
        };
        let score_max_digits =
            score_digit_count(min_score).max(score_digit_count(max_score));
        let height_max_digits = index_digit_count(matrix.height());
        let width_max_digits = index_digit_count(matrix.width());
        let max_digits = score_max_digits.max(width_max_digits);
        for _ in 0 .. height_max_digits + 2 {
            write!(f, " ")?;
        }
        write!(f, "|")?;
        for j in 0 .. matrix.width() {
            write!(f, "{:<textwidth$}|", j, textwidth = max_digits as usize)?;
        }
        write!(f, "\n")?;
        for _ in 0 .. height_max_digits + 2 {
            write!(f, " ")?;
        }
        for w in 0 .. matrix.width() {
            if w == 0 {
                write!(f, "|")?;
            } else {
                write!(f, "+")?;
            }
            for _ in 0 .. max_digits {
                write!(f, "-")?;
            }
        }
        write!(f, "|\n")?;
        for _ in 0 .. height_max_digits + 2 {
            write!(f, " ")?;
        }
        write!(f, "|")?;
        for j in 0 .. matrix.width() {
            write!(
                f,
                "{:<textwidth$}|",
                j.checked_sub(matrix.width().saturating_sub(col_seq.len()))
                    .and_then(|adjusted_j| { col_seq.get(adjusted_j).copied() })
                    .unwrap_or(' '),
                textwidth = max_digits as usize
            )?;
        }
        write!(f, "\n")?;
        for i in 0 .. matrix.height() {
            for _ in 0 .. height_max_digits {
                write!(f, "-")?;
            }
            write!(f, "+-")?;
            for w in 0 .. matrix.width() {
                if w == 0 {
                    write!(f, "|")?;
                } else {
                    write!(f, "+")?;
                }
                for _ in 0 .. max_digits {
                    if i == 0 {
                        write!(f, "=")?;
                    } else {
                        write!(f, "-")?;
                    }
                }
            }
            write!(f, "|\n")?;
            write!(
                f,
                "{:<textwidth$}|{}|",
                i,
                i.checked_sub(matrix.height().saturating_sub(row_seq.len()))
                    .and_then(|adjusted_i| { row_seq.get(adjusted_i).copied() })
                    .unwrap_or(' '),
                textwidth = height_max_digits as usize
            )?;
            for j in 0 .. matrix.width() {
                write!(
                    f,
                    "{:>textwidth$}|",
                    matrix[[i, j]],
                    textwidth = max_digits as usize
                )?;
            }
            write!(f, "\n")?;
        }
        for _ in 0 .. height_max_digits {
            write!(f, "-")?;
        }
        write!(f, "+-")?;
        for w in 0 .. matrix.width() {
            if w == 0 {
                write!(f, "|")?;
            } else {
                write!(f, "+")?;
            }
            for _ in 0 .. max_digits {
                write!(f, "=")?;
            }
        }
        write!(f, "|\n")?;
        Ok(())
    }
}
