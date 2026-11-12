use crate::scoring::{
    validate_sequence, AlignError, AlignmentResult, AlignmentStats, Dir, Scoring,
};

pub fn smith_waterman(a: &[u8], b: &[u8], s: Scoring) -> Result<AlignmentResult, AlignError> {
    s.validate()?;
    s.require_linear()?;
    validate_sequence(a, 'A')?;
    validate_sequence(b, 'B')?;

    let rows = a.len() + 1;
    let cols = b.len() + 1;

    let mut matrix = vec![0i32; rows * cols];
    let mut dirs = vec![Dir::None; rows * cols];

    let mut best_score = 0i32;
    let mut best_cell = (0usize, 0usize);

    for i in 1..rows {
        let ca = a[i - 1];
        let row = i * cols;
        let prev = row - cols;

        for j in 1..cols {
            let diag = matrix[prev + j - 1] + s.substitution(ca, b[j - 1]);
            let up = matrix[prev + j] + s.gap();
            let left = matrix[row + j - 1] + s.gap();

            let (cell, dir) = if diag >= up && diag >= left {
                (diag, Dir::Diag)
            } else if up >= left {
                (up, Dir::Up)
            } else {
                (left, Dir::Left)
            };

            let (cell, dir) = if cell > 0 { (cell, dir) } else { (0, Dir::None) };

            matrix[row + j] = cell;
            dirs[row + j] = dir;

            if cell > best_score {
                best_score = cell;
                best_cell = (i, j);
            }
        }
    }

    let cap = rows + cols;
    let mut out_a: Vec<u8> = Vec::with_capacity(cap);
    let mut out_b: Vec<u8> = Vec::with_capacity(cap);
    let mut path: Vec<(usize, usize)> = Vec::with_capacity(cap);

    let (mut i, mut j) = best_cell;

    if best_score > 0 {
        path.push((i, j));
        loop {
            match dirs[i * cols + j] {
                Dir::Diag => {
                    out_a.push(a[i - 1]);
                    out_b.push(b[j - 1]);
                    i -= 1;
                    j -= 1;
                }
                Dir::Up => {
                    out_a.push(a[i - 1]);
                    out_b.push(b'-');
                    i -= 1;
                }
                Dir::Left => {
                    out_a.push(b'-');
                    out_b.push(b[j - 1]);
                    j -= 1;
                }
                Dir::None => break,
            }
            path.push((i, j));
        }
    }

    out_a.reverse();
    out_b.reverse();
    path.reverse();

    let stats = AlignmentStats::from_aligned(&out_a, &out_b);
    let start = path.first().copied().unwrap_or(best_cell);

    Ok(AlignmentResult {
        score: best_score,
        aligned_a: out_a.into_iter().map(char::from).collect(),
        aligned_b: out_b.into_iter().map(char::from).collect(),
        rows,
        cols,
        matrix,
        dirs,
        path,
        start,
        end: best_cell,
        stats,
        algorithm: "smith-waterman",
        layers: None,
    })
}

pub fn smith_waterman_score(a: &[u8], b: &[u8], s: Scoring) -> Result<i32, AlignError> {
    s.validate()?;
    s.require_linear()?;
    validate_sequence(a, 'A')?;
    validate_sequence(b, 'B')?;

    let cols = b.len() + 1;
    let mut prev = vec![0i32; cols];
    let mut curr = vec![0i32; cols];
    let mut best = 0i32;

    for &ca in a {
        curr[0] = 0;
        for j in 1..cols {
            let diag = prev[j - 1] + s.substitution(ca, b[j - 1]);
            let up = prev[j] + s.gap();
            let left = curr[j - 1] + s.gap();
            let v = diag.max(up).max(left).max(0);
            curr[j] = v;
            best = best.max(v);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    Ok(best)
}

#[cfg(test)]
mod tests {
    use super::*;

    const S: Scoring = Scoring::new(1, -1, -1);

    #[test]
    fn wikipedia_ornegi() {
        let s = Scoring::new(3, -3, -2);
        let r = smith_waterman(b"TGTTACGG", b"GGTTGACTA", s).unwrap();
        assert_eq!(r.score, 13);
        assert_eq!(r.aligned_a, "GTT-AC");
        assert_eq!(r.aligned_b, "GTTGAC");
    }

    #[test]
    fn alt_dizi_bulunur() {
        let r = smith_waterman(b"TTTTACGTTTTT", b"GGACGGG", S).unwrap();
        assert_eq!(r.score, 3);
        assert_eq!(r.aligned_a, "ACG");
        assert_eq!(r.aligned_b, "ACG");
        assert_eq!(r.start, (4, 2));
        assert_eq!(r.end, (7, 5));
    }

    #[test]
    fn ortak_harf_yoksa_bos_sonuc() {
        let r = smith_waterman(b"AAAA", b"TTTT", S).unwrap();
        assert_eq!(r.score, 0);
        assert!(r.aligned_a.is_empty());
        assert!(r.path.is_empty());
    }

    #[test]
    fn bos_girdi() {
        let r = smith_waterman(b"", b"ACGT", S).unwrap();
        assert_eq!(r.score, 0);
        assert!(r.path.is_empty());
    }

    #[test]
    fn matriste_negatif_deger_yok() {
        let r = smith_waterman(b"GATTACA", b"GCATGCU", S).unwrap();
        assert!(r.matrix.iter().all(|&v| v >= 0));
    }

    #[test]
    fn yol_sifirdan_baslar_maksimumda_biter() {
        let r = smith_waterman(b"TGTTACGG", b"GGTTGACTA", Scoring::new(3, -3, -2)).unwrap();
        let (si, sj) = *r.path.first().unwrap();
        let (ei, ej) = *r.path.last().unwrap();
        assert_eq!(r.cell(si, sj), 0);
        assert_eq!(r.cell(ei, ej), r.score);
        assert_eq!((ei, ej), r.end);
    }

    #[test]
    fn hizalanan_skor_matris_skoruna_esit() {
        let s = Scoring::new(3, -3, -2);
        let r = smith_waterman(b"TGTTACGG", b"GGTTGACTA", s).unwrap();
        let rescored: i32 = r
            .aligned_a
            .bytes()
            .zip(r.aligned_b.bytes())
            .map(|(x, y)| if x == b'-' || y == b'-' { s.gap() } else { s.substitution(x, y) })
            .sum();
        assert_eq!(rescored, r.score);
    }

    #[test]
    fn lineer_bellek_surumu_ayni_skoru_verir() {
        let cases: &[(&[u8], &[u8], Scoring)] = &[
            (b"TGTTACGG", b"GGTTGACTA", Scoring::new(3, -3, -2)),
            (b"GATTACA", b"GCATGCU", S),
            (b"HEAGAWGHEE", b"PAWHEAE", S),
            (b"AAAA", b"TTTT", S),
        ];
        for &(a, b, s) in cases {
            let full = smith_waterman(a, b, s).unwrap().score;
            let lin = smith_waterman_score(a, b, s).unwrap();
            assert_eq!(full, lin);
        }
    }
}