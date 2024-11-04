pub type Score = i64;

pub fn score_digit_count(score: Score) -> u32 {
    if score > 0 {
        score.ilog10() + 1
    } else if score < 0 {
        (-score).ilog10() + 2
    } else {
        1
    }
}
