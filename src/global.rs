use crate::{
    matrix::AlignmentMatrix,
    values::{Letter, NormalizeLetter, Score, GAP},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct GlobalAlignmentConfig {
    pub match_weight: Score,
    pub mismatch_weight: Score,
    pub gap_weight: Score,
}

impl Default for GlobalAlignmentConfig {
    fn default() -> Self {
        Self { match_weight: 1, mismatch_weight: -1, gap_weight: -2 }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GlobalAlignmentResult {
    pub aligned_row_seq: Vec<Letter>,
    pub aligned_column_seq: Vec<Letter>,
    pub score: Score,
}

pub fn needleman_wunsch(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
) -> GlobalAlignmentResult {
    let matrix = compute_nw_matrix(row_seq, column_seq, config);
    mount_best_nw_alignment(row_seq, column_seq, config, &matrix)
}

pub fn mount_best_nw_alignment(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
    matrix: &AlignmentMatrix,
) -> GlobalAlignmentResult {
    let initial_capacity = row_seq.len() + column_seq.len();
    let mut result = GlobalAlignmentResult {
        aligned_row_seq: Vec::with_capacity(initial_capacity),
        aligned_column_seq: Vec::with_capacity(initial_capacity),
        score: 0,
    };

    let mut current_i = matrix.height() - 1;
    let mut current_j = matrix.width() - 1;
    while current_i > 0 && current_j > 0 {
        let top_left = matrix[[current_i - 1, current_j - 1]];
        let top = matrix[[current_i - 1, current_j]];
        let left = matrix[[current_i, current_j - 1]];
        let maximum = top_left.max(top).max(left);

        if current_i > 0 && current_j > 0 && top_left == maximum {
            current_i -= 1;
            current_j -= 1;
            let row_letter = row_seq.get(current_i).normalize_letter();
            let column_letter = column_seq.get(current_j).normalize_letter();

            result.aligned_row_seq.push(row_letter);
            result.aligned_column_seq.push(column_letter);

            result.score += if row_letter == GAP || column_letter == GAP {
                config.gap_weight
            } else if row_letter == column_letter {
                config.match_weight
            } else {
                config.mismatch_weight
            };
        } else if current_i > 0 && top == maximum {
            current_i -= 1;
            let row_letter = row_seq.get(current_i).normalize_letter();

            result.aligned_row_seq.push(row_letter);
            result.aligned_column_seq.push(GAP);

            result.score += config.gap_weight;
        } else {
            current_j -= 1;
            let column_letter = column_seq.get(current_j).normalize_letter();

            result.aligned_row_seq.push(GAP);
            result.aligned_column_seq.push(column_letter);

            result.score += config.gap_weight;
        }
    }
    result.aligned_row_seq.shrink_to_fit();
    result.aligned_column_seq.shrink_to_fit();
    result.aligned_row_seq.reverse();
    result.aligned_column_seq.reverse();
    result
}

pub fn compute_nw_matrix(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: GlobalAlignmentConfig,
) -> AlignmentMatrix {
    let row_count = row_seq.len() + 1;
    let column_count = column_seq.len() + 1;
    let mut matrix = AlignmentMatrix::new(row_count, column_count);
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
        let score = (j as Score) * config.gap_weight;
        matrix[[0, j]] = score;
    }
    for i in 1 ..= row_seq.len() {
        let score = (i as Score) * config.gap_weight;
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
    let no_gap_weight = if row_letter == column_letter {
        config.match_weight
    } else {
        config.mismatch_weight
    };
    let no_gap_score = top_left + no_gap_weight;

    let best_gap_neighbor = top.max(left);
    let best_gap_score = best_gap_neighbor + config.gap_weight;

    matrix[[pred_i + 1, pred_j + 1]] = best_gap_score.max(no_gap_score);
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
            match_weight: 1,
            mismatch_weight: -1,
            gap_weight: -2,
        };

        let expected_result = GlobalAlignmentResult {
            aligned_row_seq: vec!['W', 'H', 'A', 'T'],
            aligned_column_seq: vec!['W', 'H', 'Y', '-'],
            score: -1,
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
            match_weight: 1,
            mismatch_weight: -1,
            gap_weight: -1,
        };

        let expected_result = GlobalAlignmentResult {
            aligned_row_seq: vec!['G', 'C', 'A', 'T', '-', 'G', 'C', 'G'],
            aligned_column_seq: vec!['G', '-', 'A', 'T', 'T', 'A', 'C', 'A'],
            score: 0,
        };

        let actual_result = needleman_wunsch(
            &input_row_seq[..],
            &input_column_seq[..],
            input_config,
        );

        assert_eq!(actual_result, expected_result);
    }
}
