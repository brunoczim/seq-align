use std::fmt;

use crate::{
    letter::{Letter, NormalizeLetter, GAP},
    matrix::AlignmentMatrix,
    score::Score,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GlobalAlignmentConfig {
    pub match_penalty: Score,
    pub mismatch_penalty: Score,
    pub gap_penalty: Score,
}

impl Default for GlobalAlignmentConfig {
    fn default() -> Self {
        Self { match_penalty: 1, mismatch_penalty: -1, gap_penalty: -2 }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GlobalAlignmentResult {
    pub aligned_row_seq: Vec<Letter>,
    pub aligned_column_seq: Vec<Letter>,
    pub score: Score,
    pub identity_numer: u32,
    pub identity_denom: u32,
}

impl GlobalAlignmentResult {
    pub fn identity(&self) -> f64 {
        f64::from(self.identity_numer) / f64::from(self.identity_denom)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TracebackStep {
    TopLeft,
    Top,
    Left,
}

pub fn needleman_wunsch(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
) -> GlobalAlignmentResult {
    let matrix = compute_nw_matrix(row_seq, column_seq, config);
    traceback_nw_best_alignment(row_seq, column_seq, config, &matrix)
}

pub fn traceback_nw_best_alignment(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
    matrix: &AlignmentMatrix,
) -> GlobalAlignmentResult {
    let mut current_i = matrix.height() - 1;
    let mut current_j = matrix.width() - 1;

    let initial_capacity = row_seq.len() + column_seq.len();
    let mut result = GlobalAlignmentResult {
        aligned_row_seq: Vec::with_capacity(initial_capacity),
        aligned_column_seq: Vec::with_capacity(initial_capacity),
        score: matrix[[current_i, current_j]],
        identity_numer: 0,
        identity_denom: 0,
    };

    while current_i > 0 || current_j > 0 {
        let current_score = matrix[[current_i, current_j]];
        let mut maybe_step = None;
        if current_i > 0 {
            let previous_score = matrix[[current_i - 1, current_j]];
            let penalty = config.gap_penalty;
            if current_score == previous_score + penalty {
                maybe_step = Some(TracebackStep::Top);
            }
        }
        if maybe_step.is_none() && current_j > 0 {
            let previous_score = matrix[[current_i, current_j - 1]];
            let penalty = config.gap_penalty;
            if current_score == previous_score + penalty {
                maybe_step = Some(TracebackStep::Left);
            }
        }
        let step = maybe_step.unwrap_or(TracebackStep::TopLeft);

        match step {
            TracebackStep::TopLeft => {
                current_i -= 1;
                current_j -= 1;
                traceback_nw_top_left(
                    row_seq,
                    column_seq,
                    &mut result,
                    current_i,
                    current_j,
                );
            },
            TracebackStep::Top => {
                current_i -= 1;
                traceback_nw_top(row_seq, &mut result, current_i);
            },
            TracebackStep::Left => {
                current_j -= 1;
                traceback_nw_left(column_seq, &mut result, current_j);
            },
        }
    }

    result.aligned_row_seq.shrink_to_fit();
    result.aligned_column_seq.shrink_to_fit();
    result.aligned_row_seq.reverse();
    result.aligned_column_seq.reverse();
    result.identity_denom = result.identity_denom.max(1);
    result
}

pub fn compute_nw_matrix(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
) -> AlignmentMatrix {
    let row_count = row_seq.len() + 1;
    let column_count = column_seq.len() + 1;
    let mut matrix = AlignmentMatrix::zeroed(row_count, column_count);
    fill_nw_matrix_base(row_seq, column_seq, config, &mut matrix);
    fill_nw_matrix_content(row_seq, column_seq, config, &mut matrix);
    matrix
}

fn fill_nw_matrix_base(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
    matrix: &mut AlignmentMatrix,
) {
    for j in 1 ..= column_seq.len() {
        let score = (j as Score) * config.gap_penalty;
        matrix[[0, j]] = score;
    }
    for i in 1 ..= row_seq.len() {
        let score = (i as Score) * config.gap_penalty;
        matrix[[i, 0]] = score;
    }
}

fn fill_nw_matrix_content(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
    matrix: &mut AlignmentMatrix,
) {
    let mut base_i = 0;
    let mut base_j = 0;
    loop {
        if base_j >= column_seq.len() {
            break;
        }
        for j in base_j .. column_seq.len() {
            compute_nw_matrix_cell(
                row_seq, column_seq, config, matrix, base_i, j,
            );
        }
        base_i += 1;

        if base_i >= row_seq.len() {
            break;
        }
        for i in base_i .. row_seq.len() {
            compute_nw_matrix_cell(
                row_seq, column_seq, config, matrix, i, base_j,
            );
        }
        base_j += 1;
    }
}

fn compute_nw_matrix_cell(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
    matrix: &mut AlignmentMatrix,
    pred_i: usize,
    pred_j: usize,
) {
    let top_left = matrix[[pred_i, pred_j]];
    let top = matrix[[pred_i, pred_j + 1]];
    let left = matrix[[pred_i + 1, pred_j]];

    let row_letter = row_seq.get(pred_i).normalize_letter();
    let column_letter = column_seq.get(pred_j).normalize_letter();
    let no_gap_penalty = if row_letter == column_letter {
        config.match_penalty
    } else {
        config.mismatch_penalty
    };
    let no_gap_score = top_left + no_gap_penalty;

    let best_gap_neighbor = top.max(left);
    let best_gap_score = best_gap_neighbor + config.gap_penalty;

    matrix[[pred_i + 1, pred_j + 1]] = best_gap_score.max(no_gap_score);
}

fn traceback_nw_top_left(
    row_seq: &[Letter],
    column_seq: &[Letter],
    result: &mut GlobalAlignmentResult,
    current_i: usize,
    current_j: usize,
) {
    let row_letter = row_seq.get(current_i).normalize_letter();
    let column_letter = column_seq.get(current_j).normalize_letter();
    result.aligned_row_seq.push(row_letter);
    result.aligned_column_seq.push(column_letter);
    result.identity_denom += 1;
    if row_letter == column_letter && row_letter != GAP {
        result.identity_numer += 1;
    }
}

fn traceback_nw_top(
    row_seq: &[Letter],
    result: &mut GlobalAlignmentResult,
    current_i: usize,
) {
    let row_letter = row_seq.get(current_i).normalize_letter();
    result.aligned_row_seq.push(row_letter);
    result.aligned_column_seq.push(GAP);
}

fn traceback_nw_left(
    column_seq: &[Letter],
    result: &mut GlobalAlignmentResult,
    current_j: usize,
) {
    let column_letter = column_seq.get(current_j).normalize_letter();
    result.aligned_row_seq.push(GAP);
    result.aligned_column_seq.push(column_letter);
}

#[derive(Debug, Clone, Copy)]
pub struct PrettyPrint<'a> {
    pub row_seq_name: &'a str,
    pub column_seq_name: &'a str,
    pub result: &'a GlobalAlignmentResult,
    pub max_width: usize,
}

impl<'a> fmt::Display for PrettyPrint<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let identity = (100_000.0 * self.result.identity()).round() / 1000.0;
        write!(f, "# sequence above : {}\n", self.row_seq_name)?;
        write!(f, "# sequence below : {}\n", self.column_seq_name)?;
        write!(f, "# identity       : {}%\n", identity)?;
        write!(f, "# score          : {}\n", self.result.score)?;
        write!(f, "\n")?;

        let length = self
            .result
            .aligned_row_seq
            .len()
            .max(self.result.aligned_column_seq.len());
        let mut i = 0;
        while i < length {
            let block_start = i;
            let block_end = length.min(block_start + self.max_width);
            write!(f, "# block : {block_start}..{block_end}\n")?;
            for k in block_start .. block_end {
                write!(
                    f,
                    "{}",
                    self.result.aligned_row_seq.get(k).normalize_letter()
                )?;
            }
            write!(f, "\n")?;
            for k in block_start .. block_end {
                write!(
                    f,
                    "{}",
                    self.result.aligned_column_seq.get(k).normalize_letter()
                )?;
            }
            write!(f, "\n")?;

            let row_block =
                &self.result.aligned_row_seq[block_start .. block_end];
            let column_block =
                &self.result.aligned_column_seq[block_start .. block_end];
            let mut identity_iter = row_block.iter().zip(column_block);
            while let Some(k) =
                (&mut identity_iter).position(|(row_letter, column_letter)| {
                    row_letter == column_letter
                })
            {
                for _ in 0 .. k {
                    write!(f, " ")?;
                }
                write!(f, "*")?;
            }
            write!(f, "\n\n")?;
            i = block_end;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::global::GlobalAlignmentResult;

    use super::{needleman_wunsch, GlobalAlignmentConfig};

    #[test]
    fn simple_what_why_with_gap() {
        let input_row_seq = ['W', 'H', 'A', 'T'];
        let input_column_seq = ['W', 'H', 'Y'];
        let input_config = GlobalAlignmentConfig {
            match_penalty: 1,
            mismatch_penalty: -1,
            gap_penalty: -2,
        };

        let expected_result = GlobalAlignmentResult {
            aligned_row_seq: vec!['W', 'H', 'A', 'T'],
            aligned_column_seq: vec!['W', 'H', 'Y', '-'],
            score: -1,
            identity_numer: 2,
            identity_denom: 3,
        };

        let actual_result = needleman_wunsch(
            &input_row_seq[..],
            &input_column_seq[..],
            input_config,
        );

        assert_eq!(actual_result, expected_result);
    }

    #[test]
    fn multiple_inner_gaps() {
        let input_row_seq = ['G', 'C', 'A', 'T', 'G', 'C', 'G'];
        let input_column_seq = ['G', 'A', 'T', 'T', 'A', 'C', 'A'];
        let input_config = GlobalAlignmentConfig {
            match_penalty: 1,
            mismatch_penalty: -1,
            gap_penalty: -1,
        };

        let expected_result = GlobalAlignmentResult {
            aligned_row_seq: vec!['G', 'C', 'A', 'T', 'G', '-', 'C', 'G'],
            aligned_column_seq: vec!['G', '-', 'A', 'T', 'T', 'A', 'C', 'A'],
            score: 0,
            identity_numer: 4,
            identity_denom: 6,
        };

        let actual_result = needleman_wunsch(
            &input_row_seq[..],
            &input_column_seq[..],
            input_config,
        );

        assert_eq!(actual_result, expected_result);
    }
}
