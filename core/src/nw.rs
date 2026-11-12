use crate::scoring::{
    validate_sequence, AlignError, AlignmentResult, AlignmentStats, Dir, Scoring,
};

pub fn needleman_wunsch(a: &[u8], b: &[u8], s: Scoring) -> Result<AlignmentResult, AlignError> {
    s.validate()?;
    s.require_linear()?;
    validate_sequence(a, 'A')?;
    validate_sequence(b, 'B')?;

    let rows = a.len() + 1;
    let cols = b.len() + 1;

    let mut matrix = vec![0i32; rows * cols];
    let mut dirs = vec![Dir::None; rows * cols];

    for i in 1..rows {
        matrix[i * cols] = i as i32 * s.gap();
        dirs[i * cols] = Dir::Up;
    }
    for j in 1..cols {
        matrix[j] = j as i32 * s.gap();
        dirs[j] = Dir::Left;
    }

    for i in 1..rows {
        let ca = a[i - 1];
        let row = i * cols;
        let prev = row - cols;

        for j in 1..cols {
            let diag = matrix[prev + j - 1] + s.substitution(ca, b[j - 1]);
            let up = matrix[prev + j] + s.gap();
            let left = matrix[row + j - 1] + s.gap();

            let (best, dir) = if diag >= up && diag >= left {
                (diag, Dir::Diag)
            } else if up >= left {
                (up, Dir::Up)
            } else {
                (left, Dir::Left)
            };

            matrix[row + j] = best;
            dirs[row + j] = dir;
        }
    }

    let cap = rows + cols;
    let mut out_a: Vec<u8> = Vec::with_capacity(cap);
    let mut out_b: Vec<u8> = Vec::with_capacity(cap);
    let mut path: Vec<(usize, usize)> = Vec::with_capacity(cap);

    let (mut i, mut j) = (rows - 1, cols - 1);
    path.push((i, j));

    while i > 0 || j > 0 {
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
            Dir::None => unreachable!("traceback ({i},{j}) hücresinde yönsüz kaldı"),
        }
        path.push((i, j));
    }

    out_a.reverse();
    out_b.reverse();
    path.reverse();

    let stats = AlignmentStats::from_aligned(&out_a, &out_b);

    Ok(AlignmentResult {
        score: matrix[rows * cols - 1],
        aligned_a: out_a.into_iter().map(char::from).collect(),
        aligned_b: out_b.into_iter().map(char::from).collect(),
        rows,
        cols,
        matrix,
        dirs,
        path,
        start: (0, 0),
        end: (rows - 1, cols - 1),
        stats,
        algorithm: "needleman-wunsch",
        layers: None,
    })
}

pub fn needleman_wunsch_score(a: &[u8], b: &[u8], s: Scoring) -> Result<i32, AlignError> {
    s.validate()?;
    s.require_linear()?;
    validate_sequence(a, 'A')?;
    validate_sequence(b, 'B')?;

    let cols = b.len() + 1;
    let mut prev: Vec<i32> = (0..cols).map(|j| j as i32 * s.gap()).collect();
    let mut curr = vec![0i32; cols];

    for (i, &ca) in a.iter().enumerate() {
        curr[0] = (i as i32 + 1) * s.gap();
        for j in 1..cols {
            let diag = prev[j - 1] + s.substitution(ca, b[j - 1]);
            let up = prev[j] + s.gap();
            let left = curr[j - 1] + s.gap();
            curr[j] = diag.max(up).max(left);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    Ok(prev[cols - 1])
}

#[cfg(test)]
mod tests {
    use super::*;

    const S: Scoring = Scoring::new(1, -1, -1);

    #[test]
    fn wikipedia_ornegi() {
        let r = needleman_wunsch(b"GATTACA", b"GCATGCU", S).unwrap();
        assert_eq!(r.score, 0);
        assert_eq!(r.aligned_a.len(), r.aligned_b.len());
    }

    #[test]
    fn ozdes_diziler() {
        let r = needleman_wunsch(b"ACGT", b"ACGT", S).unwrap();
        assert_eq!(r.score, 4);
        assert_eq!(r.aligned_a, "ACGT");
        assert_eq!(r.aligned_b, "ACGT");
        assert_eq!(r.stats.identity, 100.0);
    }

    #[test]
    fn bos_dizi_tamamen_bosluk() {
        let r = needleman_wunsch(b"", b"ACG", S).unwrap();
        assert_eq!(r.score, -3);
        assert_eq!(r.aligned_a, "---");
        assert_eq!(r.aligned_b, "ACG");
    }

    #[test]
    fn iki_bos_dizi() {
        let r = needleman_wunsch(b"", b"", S).unwrap();
        assert_eq!(r.score, 0);
        assert_eq!(r.path, vec![(0, 0)]);
    }

    #[test]
    fn yol_koseden_koseye() {
        let r = needleman_wunsch(b"GATTACA", b"GCATGCU", S).unwrap();
        assert_eq!(r.path.first(), Some(&(0, 0)));
        assert_eq!(r.path.last(), Some(&(r.rows - 1, r.cols - 1)));
    }

    #[test]
    fn hizalanan_skor_matris_skoruna_esit() {
        let r = needleman_wunsch(b"GATTACA", b"GCATGCU", S).unwrap();
        let rescored: i32 = r
            .aligned_a
            .bytes()
            .zip(r.aligned_b.bytes())
            .map(|(x, y)| if x == b'-' || y == b'-' { S.gap() } else { S.substitution(x, y) })
            .sum();
        assert_eq!(rescored, r.score);
    }

    #[test]
    fn lineer_bellek_surumu_ayni_skoru_verir() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"GATTACA", b"GCATGCU"),
            (b"ACGTACGT", b"TGCATGCA"),
            (b"", b"AAA"),
            (b"HEAGAWGHEE", b"PAWHEAE"),
        ];
        for &(a, b) in cases {
            let full = needleman_wunsch(a, b, S).unwrap().score;
            let lin = needleman_wunsch_score(a, b, S).unwrap();
            assert_eq!(full, lin, "{:?} vs {:?}", a, b);
        }
    }

    #[test]
    fn gecersiz_karakter_reddedilir() {
        let err = needleman_wunsch(b"ACgT", b"ACGT", S).unwrap_err();
        assert!(matches!(err, AlignError::InvalidChar { seq: 'A', pos: 2, .. }));
    }

    #[test]
    fn gecersiz_skorlama_reddedilir() {
        assert!(needleman_wunsch(b"A", b"A", Scoring::new(1, -1, 2)).is_err());
        assert!(needleman_wunsch(b"A", b"A", Scoring::new(1, 1, -1)).is_err());
        assert!(needleman_wunsch(b"A", b"A", Scoring::new(0, -1, -1)).is_err());
    }
}