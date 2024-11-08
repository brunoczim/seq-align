use std::fmt;

use crate::{
    letter::{Letter, NormalizeLetter, GAP},
    matrix::AlignmentMatrix,
    score::Score,
};

/// Penalty/base score system of a global alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocalAlignmentConfig {
    /// Added when letters match.
    pub match_penalty: Score,
    /// Added when letters do not match, but it is not a gap.
    pub mismatch_penalty: Score,
    /// Added when there's a gap.
    pub gap_penalty: Score,
}

impl Default for LocalAlignmentConfig {
    fn default() -> Self {
        Self { match_penalty: 1, mismatch_penalty: -1, gap_penalty: -2 }
    }
}

/// An aligned sequence, used in local alignment results.
///
/// Corresponds to a slice of an input sequence, possibly with gaps inserted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocallyAlignedSeq {
    /// Position in the input sequence that delimits where the local alignment
    /// starts.
    pub start: usize,
    /// Position in the input sequence that delimits where the local alignment
    /// ends.
    pub end: usize,
    /// The aligned slice of the input sequence, with potential gaps.
    pub data: Vec<Letter>,
}

/// A local alignment, computed by Smith-Waterman.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalAlignmentResult {
    /// The aligned version of the input sequence that was associated with a
    /// "row" display in the matrix. It is aligned with the sequence displayed
    /// as a "column".
    pub aligned_row_seq: LocallyAlignedSeq,
    /// The aligned version of the input sequence that was associated with a
    /// "column" display in the matrix. It is aligned with the sequence
    /// displayed as a "row".
    pub aligned_column_seq: LocallyAlignedSeq,
    /// Total score of the alignment.
    pub score: Score,
    /// Numerator of the identity fraction (32-bit).
    pub identity_numer: u32,
    /// Denominator of the identity fraction (32-bit).
    pub identity_denom: u32,
}

impl LocalAlignmentResult {
    /// Computes the identity as a percentage.
    pub fn identity(&self) -> f64 {
        f64::from(self.identity_numer) / f64::from(self.identity_denom)
    }
}

/// Possible directions during traceback phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TracebackStep {
    /// Towards i - 1, j -1
    TopLeft,
    /// Towards i - 1, j
    Top,
    /// Towards i, j - 1
    Left,
}

/// Computes the Smith-Waterman algorithm, and returns all the local alignments
/// with the best score.
/// `row_seq` and `column_seq` are the sequences to be aligned.
/// `row_seq` will be displayed as a row in the matrix, while `column_seq` will
/// be displayed as a column in the matrix.
pub fn best_smith_waterman(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: LocalAlignmentConfig,
) -> Vec<LocalAlignmentResult> {
    let matrix = compute_sw_matrix(row_seq, column_seq, config);
    traceback_best_sw_alignment(row_seq, column_seq, config, &matrix)
}

/// Given Smit-Waterman input and a score matrix already populated, this
/// function computes the alignment.
pub fn traceback_best_sw_alignment(
    row_seq: &[Letter],
    column_seq: &[Letter],
    config: LocalAlignmentConfig,
    matrix: &AlignmentMatrix,
) -> Vec<LocalAlignmentResult> {
    let mut results = Vec::new();
    for (max_i, max_j) in matrix.argmax_many() {
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
            identity_numer: 0,
            identity_denom: 0,
        };

        while matrix[[current_i, current_j]] != 0 {
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
                    traceback_sw_top_left(
                        row_seq,
                        column_seq,
                        &mut result,
                        current_i,
                        current_j,
                    );
                },
                TracebackStep::Top => {
                    current_i -= 1;
                    traceback_sw_top(row_seq, &mut result, current_i);
                },
                TracebackStep::Left => {
                    current_j -= 1;
                    traceback_sw_left(column_seq, &mut result, current_j);
                },
            }

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
        result.identity_denom = result.identity_denom.max(1);

        results.push(result);
    }
    results
}

/// This function fills a Smith-Waterman score matrix.
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

/// This function fills the scores of a Smith-Waterman matrix.
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

/// Computes the score of an individual cell of a Smith-Waterman matrix,
/// assuming that their top-left (pred_i, pred_j), top (pred_i, pred_j + 1) and
/// left (pred_i + 1, pred_j) cells are already computed and filled.
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
    let no_gap_penalty = if row_letter == column_letter {
        config.match_penalty
    } else {
        config.mismatch_penalty
    };
    let no_gap_score = top_left + no_gap_penalty;

    let best_gap_neighbor = top.max(left);
    let best_gap_score = best_gap_neighbor + config.gap_penalty;

    matrix[[pred_i + 1, pred_j + 1]] = best_gap_score.max(no_gap_score).max(0);
}

/// Registers result of a traceback going to a previous top-left cell in a
/// Smith-Waterman local alignment.
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
    result.identity_denom += 1;
    if row_letter == column_letter && row_letter != GAP {
        result.identity_numer += 1;
    }
}

/// Registers result of a traceback going to a previous top cell in a
/// Smith-Waterman local alignment.
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

/// Registers result of a traceback going to a previous left cell in a
/// Smith-Waterman local alignment.
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

/// Pretty print formatting of _one_ local alignment, as in a report.
#[derive(Debug, Clone, Copy)]
pub struct PrettyPrintOne<'a> {
    /// Print name of the sequence that was associated with a row display.
    pub row_seq_name: &'a str,
    /// Print name of the sequence that was associated with a column display.
    pub column_seq_name: &'a str,
    /// An already finished local alignment result.
    pub result: &'a LocalAlignmentResult,
    /// Maximum width in terms of characters.
    pub max_width: usize,
}

impl<'a> fmt::Display for PrettyPrintOne<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let identity = (100_000.0 * self.result.identity()).round() / 1000.0;
        write!(f, "# sequence above : {}\n", self.row_seq_name)?;
        write!(f, "# sequence below : {}\n", self.column_seq_name)?;
        write!(
            f,
            "# range above    : {}..{}\n",
            self.result.aligned_row_seq.start, self.result.aligned_row_seq.end
        )?;
        write!(
            f,
            "# range below    : {}..{}\n",
            self.result.aligned_column_seq.start,
            self.result.aligned_column_seq.end
        )?;
        write!(f, "# identity       : {}%\n", identity)?;
        write!(f, "# score          : {}\n", self.result.score)?;
        write!(f, "\n")?;

        let length = self
            .result
            .aligned_row_seq
            .data
            .len()
            .max(self.result.aligned_column_seq.data.len());
        let mut i = 0;
        while i < length {
            let block_start = i;
            let block_end = length.min(block_start + self.max_width);
            write!(f, "# block : {block_start}..{block_end}\n")?;
            for k in block_start .. block_end {
                write!(
                    f,
                    "{}",
                    self.result.aligned_row_seq.data.get(k).normalize_letter()
                )?;
            }
            write!(f, "\n")?;
            for k in block_start .. block_end {
                write!(
                    f,
                    "{}",
                    self.result
                        .aligned_column_seq
                        .data
                        .get(k)
                        .normalize_letter()
                )?;
            }
            write!(f, "\n")?;

            let row_block =
                &self.result.aligned_row_seq.data[block_start .. block_end];
            let column_block =
                &self.result.aligned_column_seq.data[block_start .. block_end];
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

/// Pretty print in report formatting of all local alignment in a list of
/// results.
#[derive(Debug, Clone, Copy)]
pub struct PrettyPrintMany<'a> {
    /// Print name of the sequence that was associated with a row display.
    pub row_seq_name: &'a str,
    /// Print name of the sequence that was associated with a column display.
    pub column_seq_name: &'a str,
    /// A list of already finished local alignments with best scores.
    pub results: &'a [LocalAlignmentResult],
    /// Maximum width in terms of characters.
    pub max_width: usize,
}

impl fmt::Display for PrettyPrintMany<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.results.is_empty() {
            write!(f, "No local alignment found.")?;
        }
        for (i, result) in self.results.iter().enumerate() {
            write!(f, "#### #### #### #### #### #### #### ####\n")?;
            write!(f, "Best local alignment #{i}\n")?;
            write!(f, "#### #### #### #### #### #### #### ####\n")?;
            write!(f, "\n")?;
            let pretty_print_one = PrettyPrintOne {
                result,
                row_seq_name: self.row_seq_name,
                column_seq_name: self.column_seq_name,
                max_width: self.max_width,
            };
            write!(f, "{}\n", pretty_print_one)?;
        }
        Ok(())
    }
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
            match_penalty: 3,
            mismatch_penalty: -3,
            gap_penalty: -2,
        };

        let expected_result = vec![LocalAlignmentResult {
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
            identity_numer: 5,
            identity_denom: 5,
        }];

        let actual_result = best_smith_waterman(
            &input_row_seq[..],
            &input_column_seq[..],
            input_config,
        );

        assert_eq!(actual_result, expected_result);
    }
}
