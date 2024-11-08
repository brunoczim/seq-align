use seq_align::{
    global::{needleman_wunsch, GlobalAlignmentConfig, PrettyPrint},
    letter::Letter,
};

const ROW_SEQUENCE: &[Letter] = &['G', 'C', 'C', 'G', 'C', 'C', 'G', 'G', 'C'];

const COLUMN_SEQUENCE: &[Letter] = &['C', 'C', 'C', 'C'];

const CONFIG: GlobalAlignmentConfig = GlobalAlignmentConfig {
    gap_penalty: -4,
    match_penalty: 7,
    mismatch_penalty: -3,
};

fn main() {
    let result = needleman_wunsch(ROW_SEQUENCE, COLUMN_SEQUENCE, CONFIG);
    println!(
        "{}",
        PrettyPrint {
            row_seq_name: "<row sequence>",
            column_seq_name: "<column sequence>",
            max_width: 80,
            result: &result,
        }
    );
}
