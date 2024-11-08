use seq_align::{
    letter::Letter,
    local::{best_smith_waterman, LocalAlignmentConfig, PrettyPrintMany},
};

const ROW_SEQUENCE: &[Letter] = &[
    'M', 'T', 'E', 'N', 'S', 'T', 'S', 'T', 'P', 'A', 'A', 'K', 'P', 'K', 'R',
    'A', 'K', 'A', 'S', 'K', 'K', 'S', 'T', 'D', 'H', 'P', 'K', 'Y', 'S', 'D',
    'M', 'I', 'V', 'A', 'A', 'I', 'Q', 'A', 'E', 'K', 'N', 'R', 'A', 'G', 'S',
    'S', 'R', 'Q', 'S', 'I', 'Q', 'K', 'Y', 'I', 'K', 'S', 'H', 'Y', 'K', 'V',
    'G', 'E', 'N', 'A', 'D', 'S', 'Q', 'I', 'K', 'L', 'S', 'I', 'K', 'R', 'L',
    'V', 'T', 'T', 'G', 'V', 'L', 'K', 'Q', 'T', 'K', 'G', 'V', 'G', 'A', 'S',
    'G', 'S', 'F', 'R', 'L', 'A', 'K', 'S', 'D', 'E', 'P', 'K', 'R', 'S', 'V',
    'A', 'F', 'K', 'K', 'T', 'K', 'K', 'E', 'V', 'K', 'K', 'V', 'A', 'T', 'P',
    'K', 'K', 'A', 'A', 'K', 'P', 'K', 'K', 'A', 'A', 'S', 'K', 'A', 'P', 'S',
    'K', 'K', 'P', 'K', 'A', 'T', 'P', 'V', 'K', 'K', 'A', 'K', 'K', 'K', 'P',
    'A', 'A', 'T', 'P', 'K', 'K', 'T', 'K', 'K', 'P', 'K', 'T', 'V', 'K', 'A',
    'K', 'P', 'V', 'K', 'A', 'S', 'K', 'P', 'K', 'K', 'T', 'K', 'P', 'V', 'K',
    'P', 'K', 'A', 'K', 'S', 'S', 'A', 'K', 'R', 'T', 'G', 'K', 'K', 'K',
];

const COLUMN_SEQUENCE: &[Letter] = &[
    'M', 'S', 'E', 'T', 'A', 'P', 'V', 'P', 'Q', 'P', 'A', 'S', 'V', 'A', 'P',
    'E', 'K', 'P', 'A', 'A', 'T', 'K', 'K', 'T', 'R', 'K', 'P', 'A', 'K', 'A',
    'A', 'V', 'P', 'R', 'K', 'K', 'P', 'A', 'G', 'P', 'S', 'V', 'S', 'E', 'L',
    'I', 'V', 'Q', 'A', 'V', 'S', 'S', 'S', 'K', 'E', 'R', 'S', 'G', 'V', 'S',
    'L', 'A', 'A', 'L', 'K', 'K', 'S', 'L', 'A', 'A', 'A', 'G', 'Y', 'D', 'V',
    'E', 'K', 'N', 'N', 'S', 'R', 'I', 'K', 'L', 'G', 'L', 'K', 'S', 'L', 'V',
    'N', 'K', 'G', 'T', 'L', 'V', 'Q', 'T', 'K', 'G', 'T', 'G', 'A', 'A', 'G',
    'S', 'F', 'K', 'L', 'N', 'K', 'K', 'A', 'E', 'S', 'K', 'A', 'S', 'T', 'T',
    'K', 'V', 'T', 'V', 'K', 'A', 'K', 'A', 'S', 'G', 'A', 'A', 'K', 'K', 'P',
    'K', 'K', 'T', 'A', 'G', 'A', 'A', 'A', 'K', 'K', 'T', 'V', 'K', 'T', 'P',
    'K', 'K', 'P', 'K', 'K', 'P', 'A', 'V', 'S', 'K', 'K', 'T', 'S', 'S', 'K',
    'S', 'P', 'K', 'K', 'P', 'K', 'V', 'V', 'K', 'A', 'K', 'K', 'V', 'A', 'K',
    'S', 'P', 'A', 'K', 'A', 'K', 'A', 'V', 'K', 'P', 'K', 'A', 'A', 'K', 'V',
    'K', 'V', 'T', 'K', 'P', 'K', 'T', 'P', 'A', 'K', 'P', 'K', 'K', 'A', 'A',
    'P', 'K', 'K', 'K',
];

const CONFIG: LocalAlignmentConfig = LocalAlignmentConfig {
    gap_penalty: -2,
    match_penalty: 1,
    mismatch_penalty: -1,
};

fn main() {
    let results = best_smith_waterman(ROW_SEQUENCE, COLUMN_SEQUENCE, CONFIG);
    println!(
        "{}",
        PrettyPrintMany {
            row_seq_name: "<row sequence>",
            column_seq_name: "<column sequence>",
            max_width: 80,
            results: &results,
        }
    );
}
