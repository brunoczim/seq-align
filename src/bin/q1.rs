use seq_align::{
    global::{needleman_wunsch, GlobalAlignmentConfig, PrettyPrint},
    letter::Letter,
};

fn main() {
    let candidates = [
        ("Equus Caballus", EQUUS_CABALLUS),
        ("Odocoileus Virginianus", ODOCOILEUS_VIRGINIANUS),
        ("Bos Taurus", BOS_TAURUS),
        ("Sus Scrofa", SUS_SCROFA),
        ("Chrysocyon Brachyurus", CHRYSOCYON_BRACHYURUS),
        ("Gallus Gallus", GALLUS_GALLUS),
        ("Oncorhynchus Mykiss", ONCORHYNCHUS_MYKISS),
    ];

    let human_name = "Homo Sapiens";
    let human_sequence = HOMO_SAPIENS;

    for (candidate_name, candidate_sequence) in candidates {
        let result =
            needleman_wunsch(human_sequence, candidate_sequence, CONFIG);
        println!(
            "{}",
            PrettyPrint {
                row_seq_name: human_name,
                column_seq_name: candidate_name,
                max_width: 80,
                result: &result,
            }
        );
    }
}

const CONFIG: GlobalAlignmentConfig = GlobalAlignmentConfig {
    gap_penalty: -2,
    match_penalty: 1,
    mismatch_penalty: -1,
};

const HOMO_SAPIENS: &[Letter] = &[
    'V', 'L', 'S', 'P', 'A', 'D', 'K', 'T', 'N', 'V', 'K', 'A', 'A', 'W', 'G',
    'K', 'V', 'G', 'A', 'H', 'A', 'G', 'E', 'Y', 'G', 'A', 'E', 'A', 'L', 'E',
    'R', 'M', 'F', 'L', 'S', 'F', 'P', 'T', 'T', 'K', 'T', 'Y', 'F', 'P', 'H',
    'F', 'D', 'L', 'S', 'H', 'G', 'S', 'A', 'Q', 'V', 'K', 'G', 'H', 'G', 'K',
    'K', 'V', 'A', 'D', 'A', 'L', 'T', 'N', 'A', 'V', 'A', 'H', 'V', 'D', 'D',
    'M', 'P', 'N', 'A', 'L', 'S', 'A', 'L', 'S', 'D', 'L', 'H', 'A', 'H', 'K',
    'L', 'R', 'V', 'D', 'P', 'V', 'N', 'F', 'K', 'L', 'L', 'S', 'H', 'C', 'L',
    'L', 'V', 'T', 'L', 'A', 'A', 'H', 'L', 'P', 'A', 'E', 'F', 'T', 'P', 'A',
    'V', 'H', 'A', 'S', 'L', 'D', 'K', 'F', 'L', 'A', 'S', 'V', 'S', 'T', 'V',
    'L', 'T', 'S', 'K', 'Y',
];

const EQUUS_CABALLUS: &[Letter] = &[
    'V', 'L', 'S', 'A', 'A', 'D', 'K', 'T', 'N', 'V', 'K', 'A', 'A', 'W', 'S',
    'K', 'V', 'G', 'G', 'H', 'A', 'G', 'E', 'Y', 'G', 'A', 'E', 'A', 'L', 'E',
    'R', 'M', 'F', 'L', 'G', 'F', 'P', 'T', 'T', 'K', 'T', 'Y', 'F', 'P', 'H',
    'F', 'D', 'L', 'S', 'H', 'G', 'S', 'A', 'Q', 'V', 'K', 'A', 'H', 'G', 'K',
    'K', 'V', 'G', 'D', 'A', 'L', 'T', 'L', 'A', 'V', 'G', 'H', 'L', 'D', 'D',
    'L', 'P', 'G', 'A', 'L', 'S', 'N', 'L', 'S', 'D', 'L', 'H', 'A', 'H', 'K',
    'L', 'R', 'V', 'D', 'P', 'V', 'N', 'F', 'K', 'L', 'L', 'S', 'H', 'C', 'L',
    'L', 'S', 'T', 'L', 'A', 'V', 'H', 'L', 'P', 'N', 'D', 'F', 'T', 'P', 'A',
    'V', 'H', 'A', 'S', 'L', 'D', 'K', 'F', 'L', 'S', 'S', 'V', 'S', 'T', 'V',
    'L', 'T', 'S', 'K', 'Y', 'R',
];

const ODOCOILEUS_VIRGINIANUS: &[Letter] = &[
    'V', 'L', 'S', 'A', 'A', 'N', 'K', 'S', 'N', 'V', 'K', 'A', 'A', 'W', 'G',
    'K', 'V', 'G', 'G', 'N', 'A', 'P', 'A', 'Y', 'G', 'A', 'Q', 'A', 'L', 'Q',
    'R', 'M', 'F', 'L', 'S', 'F', 'P', 'T', 'T', 'K', 'T', 'Y', 'F', 'P', 'H',
    'F', 'D', 'L', 'S', 'H', 'G', 'S', 'A', 'Q', 'Q', 'K', 'A', 'H', 'G', 'Q',
    'K', 'V', 'A', 'N', 'A', 'L', 'T', 'K', 'A', 'Q', 'G', 'H', 'L', 'N', 'D',
    'L', 'P', 'G', 'T', 'L', 'S', 'N', 'L', 'S', 'N', 'L', 'H', 'A', 'H', 'K',
    'L', 'R', 'V', 'N', 'P', 'V', 'N', 'F', 'K', 'L', 'L', 'S', 'H', 'S', 'L',
    'L', 'V', 'T', 'L', 'A', 'S', 'H', 'L', 'P', 'T', 'N', 'F', 'T', 'P', 'A',
    'V', 'H', 'A', 'N', 'L', 'N', 'K', 'F', 'L', 'A', 'N', 'D', 'S', 'T', 'V',
    'L', 'T', 'S', 'K', 'Y', 'R',
];

const BOS_TAURUS: &[Letter] = &[
    'V', 'L', 'S', 'A', 'A', 'D', 'K', 'G', 'N', 'V', 'K', 'A', 'A', 'W', 'G',
    'K', 'V', 'G', 'G', 'H', 'A', 'A', 'E', 'Y', 'G', 'A', 'E', 'A', 'L', 'E',
    'R', 'M', 'F', 'L', 'S', 'F', 'P', 'T', 'T', 'K', 'T', 'Y', 'F', 'P', 'H',
    'F', 'D', 'L', 'S', 'H', 'G', 'S', 'A', 'Q', 'V', 'K', 'G', 'H', 'G', 'A',
    'K', 'V', 'A', 'A', 'A', 'L', 'T', 'K', 'A', 'V', 'E', 'H', 'L', 'D', 'D',
    'L', 'P', 'G', 'A', 'L', 'S', 'E', 'L', 'S', 'D', 'L', 'H', 'A', 'H', 'K',
    'L', 'R', 'V', 'D', 'P', 'V', 'N', 'F', 'K', 'L', 'L', 'S', 'H', 'S', 'L',
    'L', 'V', 'T', 'L', 'A', 'S', 'H', 'L', 'P', 'S', 'D', 'F', 'T', 'P', 'A',
    'V', 'H', 'A', 'S', 'L', 'D', 'K', 'F', 'L', 'A', 'N', 'V', 'S', 'T', 'V',
    'L', 'T', 'S', 'K', 'Y', 'R',
];

const SUS_SCROFA: &[Letter] = &[
    'V', 'L', 'S', 'A', 'A', 'D', 'K', 'A', 'N', 'V', 'K', 'A', 'A', 'W', 'G',
    'K', 'V', 'G', 'G', 'Q', 'A', 'G', 'A', 'H', 'G', 'A', 'E', 'A', 'L', 'E',
    'R', 'M', 'F', 'L', 'G', 'F', 'P', 'T', 'T', 'K', 'T', 'Y', 'F', 'P', 'H',
    'F', 'N', 'L', 'S', 'H', 'G', 'S', 'D', 'Q', 'V', 'K', 'A', 'H', 'G', 'Q',
    'K', 'V', 'A', 'D', 'A', 'L', 'T', 'K', 'A', 'V', 'G', 'H', 'L', 'D', 'D',
    'L', 'P', 'G', 'A', 'L', 'S', 'A', 'L', 'S', 'D', 'L', 'H', 'A', 'H', 'K',
    'L', 'R', 'V', 'D', 'P', 'V', 'N', 'F', 'K', 'L', 'L', 'S', 'H', 'C', 'L',
    'L', 'V', 'T', 'L', 'A', 'A', 'H', 'H', 'P', 'D', 'D', 'F', 'N', 'P', 'S',
    'V', 'H', 'A', 'S', 'L', 'D', 'K', 'F', 'L', 'A', 'N', 'V', 'S', 'T', 'V',
    'L', 'T', 'S', 'K', 'Y', 'R',
];

const CHRYSOCYON_BRACHYURUS: &[Letter] = &[
    'V', 'L', 'S', 'P', 'A', 'D', 'K', 'T', 'N', 'I', 'K', 'S', 'T', 'W', 'D',
    'K', 'I', 'G', 'G', 'H', 'A', 'G', 'D', 'Y', 'G', 'G', 'E', 'A', 'L', 'D',
    'R', 'T', 'F', 'Q', 'S', 'F', 'P', 'T', 'T', 'K', 'T', 'Y', 'F', 'P', 'H',
    'F', 'D', 'L', 'S', 'P', 'G', 'S', 'A', 'Q', 'V', 'K', 'A', 'H', 'G', 'K',
    'K', 'V', 'A', 'D', 'A', 'L', 'T', 'T', 'A', 'V', 'A', 'H', 'L', 'D', 'D',
    'L', 'P', 'G', 'A', 'L', 'S', 'A', 'L', 'S', 'D', 'L', 'H', 'A', 'Y', 'K',
    'L', 'R', 'V', 'D', 'P', 'V', 'N', 'F', 'K', 'L', 'L', 'S', 'H', 'C', 'L',
    'L', 'V', 'T', 'L', 'A', 'C', 'H', 'H', 'P', 'T', 'E', 'F', 'T', 'P', 'A',
    'V', 'H', 'A', 'S', 'L', 'D', 'K', 'F', 'F', 'T', 'A', 'V', 'S', 'T', 'V',
    'L', 'T', 'S', 'K', 'Y', 'R',
];

const GALLUS_GALLUS: &[Letter] = &[
    'M', 'L', 'T', 'A', 'E', 'D', 'K', 'K', 'L', 'I', 'Q', 'Q', 'A', 'W', 'E',
    'K', 'A', 'A', 'S', 'H', 'Q', 'E', 'E', 'F', 'G', 'A', 'E', 'A', 'L', 'T',
    'R', 'M', 'F', 'T', 'T', 'Y', 'P', 'Q', 'T', 'K', 'T', 'Y', 'F', 'P', 'H',
    'F', 'D', 'L', 'S', 'P', 'G', 'S', 'D', 'Q', 'V', 'R', 'G', 'H', 'G', 'K',
    'K', 'V', 'L', 'G', 'A', 'L', 'G', 'N', 'A', 'V', 'K', 'N', 'V', 'D', 'N',
    'L', 'S', 'Q', 'A', 'M', 'A', 'E', 'L', 'S', 'N', 'L', 'H', 'A', 'Y', 'N',
    'L', 'R', 'V', 'D', 'P', 'V', 'N', 'F', 'K', 'L', 'L', 'S', 'Q', 'C', 'I',
    'Q', 'V', 'V', 'L', 'A', 'V', 'H', 'M', 'G', 'K', 'D', 'Y', 'T', 'P', 'E',
    'V', 'H', 'A', 'A', 'F', 'D', 'K', 'F', 'L', 'S', 'A', 'V', 'S', 'A', 'V',
    'L', 'A', 'E', 'K', 'Y', 'R',
];

const ONCORHYNCHUS_MYKISS: &[Letter] = &[
    'X', 'S', 'L', 'T', 'A', 'K', 'D', 'K', 'S', 'V', 'V', 'K', 'A', 'F', 'W',
    'G', 'K', 'I', 'S', 'G', 'K', 'A', 'D', 'V', 'V', 'G', 'A', 'E', 'A', 'L',
    'G', 'R', 'M', 'L', 'T', 'A', 'Y', 'P', 'Q', 'T', 'K', 'T', 'Y', 'F', 'S',
    'H', 'W', 'A', 'D', 'L', 'S', 'P', 'G', 'S', 'G', 'P', 'V', 'K', 'K', 'H',
    'G', 'G', 'I', 'I', 'M', 'G', 'A', 'I', 'G', 'K', 'A', 'V', 'G', 'L', 'M',
    'D', 'D', 'L', 'V', 'G', 'G', 'M', 'S', 'A', 'L', 'S', 'D', 'L', 'H', 'A',
    'F', 'K', 'L', 'R', 'V', 'D', 'P', 'G', 'N', 'F', 'K', 'I', 'L', 'S', 'H',
    'N', 'I', 'L', 'V', 'T', 'L', 'A', 'I', 'H', 'F', 'P', 'S', 'D', 'F', 'T',
    'P', 'E', 'V', 'H', 'I', 'A', 'V', 'D', 'K', 'F', 'L', 'A', 'A', 'V', 'S',
    'A', 'A', 'L', 'A', 'D', 'K', 'Y', 'R',
];
