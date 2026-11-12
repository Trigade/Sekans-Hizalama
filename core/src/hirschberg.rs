use crate::nw::needleman_wunsch;
use crate::scoring::{validate_sequence, AlignError, AlignmentResult, AlignmentStats, Scoring};

pub fn hirschberg(a: &[u8], b: &[u8], s: Scoring) -> Result<AlignmentResult, AlignError> {
    s.validate()?;
    s.require_linear()?;
    validate_sequence(a, 'A')?;
    validate_sequence(b, 'B')?;

    let mut out_a = Vec::with_capacity(a.len() + b.len());
    let mut out_b = Vec::with_capacity(a.len() + b.len());
    let mut fwd = vec![0i32; b.len() + 1];
    let mut rev = vec![0i32; b.len() + 1];

    solve(a, b, s, &mut out_a, &mut out_b, &mut fwd, &mut rev)?;

    let score = out_a
        .iter()
        .zip(&out_b)
        .map(|(&x, &y)| if x == b'-' || y == b'-' { s.gap() } else { s.substitution(x, y) })
        .sum();

    let mut path = Vec::with_capacity(out_a.len() + 1);
    let (mut i, mut j) = (0usize, 0usize);
    path.push((i, j));
    for (&x, &y) in out_a.iter().zip(&out_b) {
        if x != b'-' {
            i += 1;
        }
        if y != b'-' {
            j += 1;
        }
        path.push((i, j));
    }

    let stats = AlignmentStats::from_aligned(&out_a, &out_b);

    Ok(AlignmentResult {
        score,
        aligned_a: out_a.into_iter().map(char::from).collect(),
        aligned_b: out_b.into_iter().map(char::from).collect(),
        rows: a.len() + 1,
        cols: b.len() + 1,
        matrix: Vec::new(),
        dirs: Vec::new(),
        path,
        start: (0, 0),
        end: (a.len(), b.len()),
        stats,
        algorithm: "hirschberg",
        layers: None,
    })
}

fn solve(
    a: &[u8],
    b: &[u8],
    s: Scoring,
    out_a: &mut Vec<u8>,
    out_b: &mut Vec<u8>,
    fwd: &mut [i32],
    rev: &mut [i32],
) -> Result<(), AlignError> {
    if a.is_empty() {
        out_a.extend(std::iter::repeat_n(b'-', b.len()));
        out_b.extend_from_slice(b);
        return Ok(());
    }
    if b.is_empty() {
        out_a.extend_from_slice(a);
        out_b.extend(std::iter::repeat_n(b'-', a.len()));
        return Ok(());
    }
    if a.len() == 1 || b.len() == 1 {
        let r = needleman_wunsch(a, b, s)?;
        out_a.extend_from_slice(r.aligned_a.as_bytes());
        out_b.extend_from_slice(r.aligned_b.as_bytes());
        return Ok(());
    }

    let mid = a.len() / 2;
    let m = b.len();

    last_row_forward(&a[..mid], b, s, &mut fwd[..=m]);
    last_row_reverse(&a[mid..], b, s, &mut rev[..=m]);

    let mut split = 0;
    let mut best = i32::MIN;
    for k in 0..=m {
        let v = fwd[k] + rev[m - k];
        if v > best {
            best = v;
            split = k;
        }
    }

    solve(&a[..mid], &b[..split], s, out_a, out_b, fwd, rev)?;
    solve(&a[mid..], &b[split..], s, out_a, out_b, fwd, rev)
}

fn last_row_forward(a: &[u8], b: &[u8], s: Scoring, row: &mut [i32]) {
    let gap = s.gap();
    for (j, cell) in row.iter_mut().enumerate() {
        *cell = j as i32 * gap;
    }
    for (i, &ca) in a.iter().enumerate() {
        let mut diag = row[0];
        row[0] = (i as i32 + 1) * gap;
        for j in 1..row.len() {
            let up = row[j];
            let v = (diag + s.substitution(ca, b[j - 1])).max(up + gap).max(row[j - 1] + gap);
            diag = up;
            row[j] = v;
        }
    }
}

fn last_row_reverse(a: &[u8], b: &[u8], s: Scoring, row: &mut [i32]) {
    let gap = s.gap();
    let m = b.len();
    for (j, cell) in row.iter_mut().enumerate() {
        *cell = j as i32 * gap;
    }
    for (i, &ca) in a.iter().rev().enumerate() {
        let mut diag = row[0];
        row[0] = (i as i32 + 1) * gap;
        for j in 1..=m {
            let up = row[j];
            let v = (diag + s.substitution(ca, b[m - j])).max(up + gap).max(row[j - 1] + gap);
            diag = up;
            row[j] = v;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nw::needleman_wunsch_score;

    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *seed >> 33
    }

    fn random_seq(len: usize, alphabet: &[u8], seed: &mut u64) -> Vec<u8> {
        (0..len).map(|_| alphabet[(lcg(seed) as usize) % alphabet.len()]).collect()
    }

    #[test]
    fn rastgele_dizilerde_nw_ile_ayni_skor() {
        let mut seed = 42u64;
        let schemes = [Scoring::new(1, -1, -1), Scoring::new(2, -1, -2), Scoring::blosum62(0, -4)];
        for scheme in schemes {
            for _ in 0..60 {
                let la = (lcg(&mut seed) % 40) as usize;
                let lb = (lcg(&mut seed) % 40) as usize;
                let alpha: &[u8] = if matches!(scheme.subst, crate::scoring::Substitution::Blosum62) {
                    b"ARNDCQEGHILKMFPSTWYV"
                } else {
                    b"ACGT"
                };
                let a = random_seq(la, alpha, &mut seed);
                let b = random_seq(lb, alpha, &mut seed);
                let h = hirschberg(&a, &b, scheme).unwrap();
                let expected = needleman_wunsch_score(&a, &b, scheme).unwrap();
                assert_eq!(h.score, expected);
                let clean_a: Vec<u8> = h.aligned_a.bytes().filter(|&c| c != b'-').collect();
                let clean_b: Vec<u8> = h.aligned_b.bytes().filter(|&c| c != b'-').collect();
                assert_eq!(clean_a, a);
                assert_eq!(clean_b, b);
            }
        }
    }

    #[test]
    fn yol_koseden_koseye() {
        let r = hirschberg(b"GATTACA", b"GCATGCU", Scoring::new(1, -1, -1)).unwrap();
        assert_eq!(r.score, 0);
        assert_eq!(r.path.first(), Some(&(0, 0)));
        assert_eq!(r.path.last(), Some(&(7, 7)));
        assert!(r.matrix.is_empty());
    }

    #[test]
    fn uzun_dizi() {
        let mut seed = 7u64;
        let a = random_seq(3000, b"ACGT", &mut seed);
        let b = random_seq(2500, b"ACGT", &mut seed);
        let s = Scoring::new(1, -1, -1);
        let h = hirschberg(&a, &b, s).unwrap();
        assert_eq!(h.score, needleman_wunsch_score(&a, &b, s).unwrap());
    }

    #[test]
    fn afin_reddedilir() {
        assert!(hirschberg(b"A", b"A", Scoring::blosum62(-11, -1)).is_err());
    }
}