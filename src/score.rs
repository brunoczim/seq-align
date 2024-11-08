/// Score is an 64-bit signed integer (allows negative values).
pub type Score = i64;

// Counts how many decimal digits a score needs to be rendered.
pub fn score_digit_count(score: Score) -> u32 {
    if score > 0 {
        score.ilog10() + 1
    } else if score < 0 {
        (-score).ilog10() + 2
    } else {
        1
    }
}
