use crate::{
    letter::{Letter, NormalizeLetter, GAP},
    matrix::AlignmentMatrix,
    score::Score,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocalAlignmentConfig {
    pub match_weight: Score,
    pub mismatch_weight: Score,
    pub gap_weight: Score,
}

impl Default for LocalAlignmentConfig {
    fn default() -> Self {
        Self { match_weight: 1, mismatch_weight: -1, gap_weight: -2 }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocallyAlignedSeq {
    pub start: usize,
    pub end: usize,
    pub data: Vec<Letter>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalAlignmentResult {
    pub aligned_row_seq: LocallyAlignedSeq,
    pub aligned_column_seq: LocallyAlignedSeq,
    pub score: Score,
}

pub fn best_smith_waterman(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: LocalAlignmentConfig,
) -> LocalAlignmentResult {
    let matrix = compute_sw_matrix(row_seq, column_seq, config);
    traceback_best_sw_alignment(row_seq, column_seq, &matrix)
}

pub fn traceback_best_sw_alignment(
    row_seq: &[Letter],
    column_seq: &[Letter],
    matrix: &AlignmentMatrix,
) -> LocalAlignmentResult {
    let (max_i, max_j) = matrix
        .argmax_rev()
        .unwrap_or((matrix.height() - 1, matrix.width() - 1));
    let mut current_i = max_i;
    let mut current_j = max_j;

    let initial_capacity = row_seq.len() + column_seq.len();
    let mut result = LocalAlignmentResult {
        aligned_row_seq: LocallyAlignedSeq {
            start: max_i,
            end: max_i,
            data: Vec::with_capacity(initial_capacity),
        },
        aligned_column_seq: LocallyAlignedSeq {
            start: max_j,
            end: max_j,
            data: Vec::with_capacity(initial_capacity),
        },
        score: matrix[[max_i, max_j]],
    };

    while matrix[[current_i, current_j]] != 0 {
        let top_left = matrix[[current_i - 1, current_j - 1]];
        let top = matrix[[current_i - 1, current_j]];
        let left = matrix[[current_i, current_j - 1]];
        let maximum = top_left.max(top).max(left);

        if current_i > 0
            && current_j > 0
            && (top_left == maximum || top_left == 0)
        {
            current_i -= 1;
            current_j -= 1;
            traceback_sw_top_left(
                row_seq,
                column_seq,
                &mut result,
                current_i,
                current_j,
            );
        } else if current_i > 0 && (top == maximum || top == 0) {
            current_i -= 1;
            traceback_sw_top(row_seq, &mut result, current_i);
        } else {
            current_j -= 1;
            traceback_sw_left(column_seq, &mut result, current_j);
        }
    }

    result.aligned_row_seq.data.reverse();
    result.aligned_column_seq.data.reverse();

    result
}

pub fn compute_sw_matrix(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: LocalAlignmentConfig,
) -> AlignmentMatrix {
    let row_count = row_seq.len() + 1;
    let column_count = column_seq.len() + 1;
    let mut matrix = AlignmentMatrix::zeroed(row_count, column_count);
    fill_sw_matrix_content(row_seq, column_seq, config, &mut matrix);
    matrix
}

fn fill_sw_matrix_content(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: LocalAlignmentConfig,
    matrix: &mut AlignmentMatrix,
) {
    let mut base_i = 0;
    let mut base_j = 0;
    loop {
        if base_j >= column_seq.len() {
            break;
        }
        for j in base_j .. column_seq.len() {
            compute_sw_matrix_cell(
                row_seq, column_seq, config, matrix, base_i, j,
            );
        }
        base_i += 1;

        if base_i >= row_seq.len() {
            break;
        }
        for i in base_i .. row_seq.len() {
            compute_sw_matrix_cell(
                row_seq, column_seq, config, matrix, i, base_j,
            );
        }
        base_j += 1;
    }
}

fn compute_sw_matrix_cell(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: LocalAlignmentConfig,
    matrix: &mut AlignmentMatrix,
    pred_i: usize,
    pred_j: usize,
) {
    let top_left = matrix[[pred_i, pred_j]];
    let top = matrix[[pred_i, pred_j + 1]];
    let left = matrix[[pred_i + 1, pred_j]];

    let row_letter = row_seq.get(pred_i).normalize_letter();
    let column_letter = column_seq.get(pred_j).normalize_letter();
    let no_gap_weight = if row_letter == column_letter {
        config.match_weight
    } else {
        config.mismatch_weight
    };
    let no_gap_score = top_left + no_gap_weight;

    let best_gap_neighbor = top.max(left);
    let best_gap_score = best_gap_neighbor + config.gap_weight;

    matrix[[pred_i + 1, pred_j + 1]] = best_gap_score.max(no_gap_score).max(0);
}

fn traceback_sw_top_left(
    row_seq: &[Letter],
    column_seq: &[Letter],
    result: &mut LocalAlignmentResult,
    current_i: usize,
    current_j: usize,
) {
    let row_letter = row_seq.get(current_i).normalize_letter();
    let column_letter = column_seq.get(current_j).normalize_letter();
    result.aligned_row_seq.start -= 1;
    result.aligned_row_seq.data.push(row_letter);
    result.aligned_column_seq.start -= 1;
    result.aligned_column_seq.data.push(column_letter);
}

fn traceback_sw_top(
    row_seq: &[Letter],
    result: &mut LocalAlignmentResult,
    current_i: usize,
) {
    let row_letter = row_seq.get(current_i).normalize_letter();
    result.aligned_row_seq.start -= 1;
    result.aligned_row_seq.data.push(row_letter);
    result.aligned_column_seq.data.push(GAP);
}

fn traceback_sw_left(
    column_seq: &[Letter],
    result: &mut LocalAlignmentResult,
    current_j: usize,
) {
    let column_letter = column_seq.get(current_j).normalize_letter();
    result.aligned_row_seq.data.push(GAP);
    result.aligned_column_seq.start -= 1;
    result.aligned_column_seq.data.push(column_letter);
}

#[cfg(test)]
mod test {
    use super::{
        best_smith_waterman,
        LocalAlignmentConfig,
        LocalAlignmentResult,
        LocallyAlignedSeq,
    };

    #[test]
    fn easy_case() {
        let input_row_seq = ['G', 'G', 'T', 'T', 'G', 'A', 'C', 'T', 'A'];
        let input_column_seq = ['T', 'G', 'T', 'T', 'A', 'C', 'G', 'G'];
        let input_config = LocalAlignmentConfig {
            match_weight: 3,
            mismatch_weight: -3,
            gap_weight: -2,
        };

        let expected_result = LocalAlignmentResult {
            aligned_row_seq: LocallyAlignedSeq {
                start: 1,
                end: 7,
                data: vec!['G', 'T', 'T', 'G', 'A', 'C'],
            },
            aligned_column_seq: LocallyAlignedSeq {
                start: 1,
                end: 6,
                data: vec!['G', 'T', 'T', '-', 'A', 'C'],
            },
            score: 13,
        };

        let actual_result = best_smith_waterman(
            &input_row_seq[..],
            &input_column_seq[..],
            input_config,
        );

        assert_eq!(actual_result, expected_result);
    }
}
