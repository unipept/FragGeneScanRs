use crate::dna::{get_rc_dna, nt2int, trinucleotide};
use crate::{gene, hmm};

pub fn viterbi(
    global: &hmm::Global,
    local: &hmm::Local,
    head: Vec<u8>,
    seq: Vec<u8>,
    whole_genome: bool,
    formatted: bool,
) -> Vec<gene::Gene> {
    let gene_len = if whole_genome { 120 } else { 60 }; // minimum length to be output

    let mut alpha: Vec<[f64; hmm::NUM_STATE]> = vec![];
    let mut path: Vec<[hmm::State; hmm::NUM_STATE]> = vec![];
    let mut temp_i: [usize; 6] = [0; 6];
    let mut temp_i_1: [usize; 6] = [0; 6];

    for _ in 0..seq.len() {
        alpha.push([0.0; hmm::NUM_STATE]);
        path.push([hmm::State::N; hmm::NUM_STATE]);
    }
    alpha[0].copy_from_slice(&global.pi);
    for i in &mut alpha[0] {
        *i *= -1.0
    } // TODO broadcast operation

    // If the sequence starts with a stop codon
    if seq[0] == b'T'
        && ((seq[1] == b'A' && seq[2] == b'A')
            || (seq[1] == b'A' && seq[2] == b'G')
            || (seq[1] == b'G' && seq[2] == b'A'))
    {
        alpha[0][hmm::State::E] = f64::INFINITY;
        alpha[1][hmm::State::E] = f64::INFINITY;
        path[1][hmm::State::E] = hmm::State::E;
        path[2][hmm::State::E] = hmm::State::E;

        alpha[2][hmm::State::M6] = f64::INFINITY;
        alpha[1][hmm::State::M5] = f64::INFINITY;
        alpha[0][hmm::State::M4] = f64::INFINITY;
        alpha[2][hmm::State::M3] = f64::INFINITY;
        alpha[1][hmm::State::M2] = f64::INFINITY;
        alpha[0][hmm::State::M1] = f64::INFINITY;

        alpha[2][hmm::State::E] -= if seq[1] == b'A' && seq[2] == b'A' {
            0.53_f64.ln()
        } else if seq[1] == b'A' && seq[2] == b'G' {
            0.16_f64.ln()
        } else {
            0.30_f64.ln()
        }
    }

    // If the sequence starts with a reverse stop codon
    if seq[2] == b'A'
        && ((seq[1] == b'T' && seq[0] == b'T')
            || (seq[1] == b'T' && seq[0] == b'C')
            || (seq[1] == b'C' && seq[0] == b'T'))
    {
        alpha[0][hmm::State::Sr] = f64::INFINITY;
        alpha[1][hmm::State::Sr] = f64::INFINITY;
        alpha[2][hmm::State::Sr] = alpha[0][hmm::State::S];
        path[1][hmm::State::Sr] = hmm::State::Sr;
        path[2][hmm::State::Sr] = hmm::State::Sr;

        alpha[2][hmm::State::M3r] = f64::INFINITY;
        alpha[2][hmm::State::M6r] = f64::INFINITY;

        alpha[2][hmm::State::Sr] = if seq[1] == b'T' && seq[0] == b'T' {
            0.53_f64.ln()
        } else if seq[1] == b'T' && seq[0] == b'C' {
            0.16_f64.ln()
        } else {
            0.30_f64.ln()
        }
    }

    let mut num_noncoding = 0; // number of invalid nts in sequence
    for t in 1..seq.len() {
        let from = nt2int(seq[t - 1]).unwrap_or(2);
        let from0 = if t > 1 {
            nt2int(seq[t - 2]).unwrap_or(2)
        } else {
            2
        };
        let to = nt2int(seq[t]).unwrap_or_else(|| {
            num_noncoding += 1;
            2
        });
        let from2 = from0 * 4 + from;

        // M state

        for (k, i) in hmm::State::M1.up_through(hmm::State::M6).enumerate() {
            if alpha[t][i].is_finite() {
                if i == hmm::State::M1 {
                    // from M state
                    alpha[t][i] = alpha[t - 1][hmm::State::M6]
                        - global.tr.gg
                        - global.tr.mm
                        - local.e_m[0][from2][to];
                    path[t][i] = hmm::State::M6;

                    // from D state
                    if !whole_genome {
                        for j in hmm::State::M1.up_through(hmm::State::M5).rev() {
                            let num_d = if j >= i {
                                i as i32 + 6 - j as i32
                            } else if j.next() < i {
                                i as i32 - j as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j]
                                    - global.tr.md
                                    - local.e_m[0][from2][to]
                                    - 0.25_f64.ln() * (num_d - 1) as f64
                                    - global.tr.dd * (num_d - 2) as f64
                                    - global.tr.dm;
                                if temp_alpha < alpha[t][i] {
                                    alpha[t][i] = temp_alpha;
                                    path[t][i] = j;
                                }
                            }
                        }
                    }

                    // from Start state
                    let temp_alpha = alpha[t - 1][hmm::State::S] - local.e_m[0][from2][to];
                    if temp_alpha < alpha[t][i] {
                        alpha[t][i] = temp_alpha;
                        path[t][i] = hmm::State::S;
                    }
                } else {
                    // from M state
                    alpha[t][i] =
                        alpha[t - 1][i.previous()] - global.tr.mm - local.e_m[k][from2][to];
                    path[t][i] = i.previous();

                    // from D state
                    if !whole_genome {
                        for j in hmm::State::M1.up_through(hmm::State::M6).rev() {
                            let num_d = if j >= i {
                                i as i32 + 6 - j as i32
                            } else if j.next() < i {
                                i as i32 - j as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j]
                                    - global.tr.md
                                    - local.e_m[k][from2][to]
                                    - 0.25_f64.ln() * (num_d - 1) as f64
                                    - global.tr.dd * (num_d - 2) as f64
                                    - global.tr.dm;
                                if temp_alpha < alpha[t][i] {
                                    alpha[t][i] = temp_alpha;
                                    path[t][i] = j;
                                }
                            }
                        }
                    }
                }

                // from I state (l251)
                let j = if i == hmm::State::M1 {
                    hmm::State::I6
                } else {
                    hmm::State::I1.up_through(hmm::State::N).nth(k - 1).unwrap()
                    // TODO
                };

                // to avoid stop codon
                if t < 2 {
                } else if (i == hmm::State::M2 || i == hmm::State::M5)
                    && t + 1 < seq.len()
                    && seq[temp_i[j as usize - hmm::State::I1 as usize]] == b'T'
                    && ((seq[t] == b'A' && seq[t + 1] == b'A')
                        || (seq[t] == b'A' && seq[t + 1] == b'G')
                        || (seq[t] == b'G' && seq[t + 1] == b'A'))
                {
                } else if (i == hmm::State::M3 || i == hmm::State::M6)
                    && temp_i[j as usize - hmm::State::I1 as usize] > 0
                    && (seq[temp_i[j as usize - hmm::State::I1 as usize] - 1] == b'T')
                    && ((seq[temp_i[j as usize - hmm::State::I1 as usize]] == b'A'
                        && seq[t] == b'A')
                        || (seq[temp_i[j as usize - hmm::State::I1 as usize]] == b'A'
                            && seq[t] == b'G')
                        || (seq[temp_i[j as usize - hmm::State::I1 as usize]] == b'G'
                            && seq[t] == b'A'))
                {
                } else {
                    let temp_alpha = alpha[t - 1][j] - global.tr.im - 0.25_f64.ln();
                    if temp_alpha < alpha[t][i] {
                        alpha[t][i] = temp_alpha;
                        path[t][i] = j;
                    }
                }
            }
        }

        // I state
        for (i, m) in hmm::State::I1
            .up_through(hmm::State::I6)
            .zip(hmm::State::M1.up_through(hmm::State::M6))
        {
            // from I state
            alpha[t][i] = alpha[t - 1][i] - global.tr.ii - global.tr_ii[from][to];
            path[t][i] = i;

            // from M state
            let temp_alpha = alpha[t - 1][m]
                - global.tr.mi
                - global.tr_mi[from][to]
                - if i == hmm::State::I6 {
                    global.tr.gg
                } else {
                    0.0
                };
            if temp_alpha < alpha[t][i] {
                alpha[t][i] = temp_alpha;
                path[t][i] = m;
                temp_i[i as usize - hmm::State::I1 as usize] = t - 1;
            }
        }

        // M' state
        for (k, i) in hmm::State::M1r.up_through(hmm::State::M6r).enumerate() {
            if (i == hmm::State::M1r || i == hmm::State::M4r)
                && t >= 3
                && seq[t - 1] == b'A'
                && ((seq[t - 2] == b'T' && seq[t - 3] == b'T')
                    || (seq[t - 2] == b'T' && seq[t - 3] == b'C')
                    || (seq[t - 2] == b'C' && seq[t - 3] == b'T'))
            {
                // from Start state since this is actually a stop codon in minus strand
                alpha[t][i] = alpha[t - 1][hmm::State::Sr] - local.e_m1[k][from2][to];
                path[t][i] = hmm::State::Sr;
            } else {
                if i == hmm::State::M1r {
                    // from M state
                    alpha[t][i] = alpha[t - 1][hmm::State::M6r]
                        - global.tr.gg
                        - global.tr.mm
                        - local.e_m1[0][from2][to];
                    path[t][i] = hmm::State::M6r;

                    // from D state
                    if !whole_genome {
                        for j in hmm::State::M1r.up_through(hmm::State::M5r).rev() {
                            let num_d = if j >= i {
                                i as i32 + 6 - j as i32
                            } else if j.next() < i {
                                i as i32 - j as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j] - global.tr.md
                                               - local.e_m1[0][from2][to] // TODO different from forward but merge?
                                               - 0.25_f64.ln() * (num_d - 1) as f64
                                               - global.tr.dd * (num_d - 2) as f64
                                               - global.tr.dm;
                                if temp_alpha < alpha[t][i] {
                                    alpha[t][i] = temp_alpha;
                                    path[t][i] = j;
                                }
                            }
                        }
                    }
                } else {
                    // from M state
                    alpha[t][i] =
                        alpha[t - 1][i.previous()] - global.tr.mm - local.e_m1[k][from2][to];
                    path[t][i] = i.previous();

                    // from D state
                    if !whole_genome {
                        for j in hmm::State::M1r.up_through(hmm::State::M6r).rev() {
                            let num_d = if j >= i {
                                i as i32 + 6 - j as i32
                            } else if j.next() < i {
                                i as i32 - j as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j]
                                    - global.tr.md
                                    - local.e_m1[k][from2][to]
                                    - 0.25_f64.ln() * (num_d - 1) as f64
                                    - global.tr.dd * (num_d - 2) as f64
                                    - global.tr.dm;
                                if temp_alpha < alpha[t][i] {
                                    alpha[t][i] = temp_alpha;
                                    path[t][i] = j;
                                }
                            }
                        }
                    }
                }

                // from I state
                let j = if i == hmm::State::M1r {
                    hmm::State::I6r
                } else {
                    hmm::State::I1r
                        .up_through(hmm::State::N)
                        .nth(k - 1)
                        .unwrap() // TODO
                };

                // to avoid stop codon
                if t < 2 {
                } else if (i == hmm::State::M2r || i == hmm::State::M5r)
                    && t + 1 < seq.len()
                    && seq[t + 1] == b'A'
                    && ((seq[t] == b'T'
                        && seq[temp_i_1[j as usize - hmm::State::I1r as usize]] == b'T')
                        || (seq[t] == b'T'
                            && seq[temp_i_1[j as usize - hmm::State::I1r as usize]] == b'C')
                        || (seq[t] == b'A'
                            && seq[temp_i_1[j as usize - hmm::State::I1r as usize]] == b'T'))
                {
                } else if (i == hmm::State::M3r || i == hmm::State::M6r)
                    && seq[t] == b'A'
                    && temp_i_1[j as usize - hmm::State::I1r as usize] > 1
                    && ((seq[temp_i_1[j as usize - hmm::State::I1r as usize]] == b'T'
                        && seq[temp_i_1[j as usize - hmm::State::I1r as usize] - 1] == b'T')
                        || (seq[temp_i_1[j as usize - hmm::State::I1r as usize]] == b'T'
                            && seq[temp_i_1[j as usize - hmm::State::I1r as usize] - 1] == b'C')
                        || (seq[temp_i_1[j as usize - hmm::State::I1r as usize]] == b'C'
                            && seq[temp_i_1[j as usize - hmm::State::I1r as usize] - 1] == b'T'))
                {
                } else {
                    let temp_alpha = alpha[t - 1][j] - global.tr.im - 0.25_f64.ln();
                    if temp_alpha < alpha[t][i] {
                        alpha[t][i] = temp_alpha;
                        path[t][i] = j;
                    }
                }
            }
        }

        // I' state
        for (i, m) in hmm::State::I1r
            .up_through(hmm::State::I6r)
            .zip(hmm::State::M1r.up_through(hmm::State::M6r))
        {
            // from I state
            alpha[t][i] = alpha[t - 1][i] - global.tr.ii - global.tr_ii[from][to];
            path[t][i] = i;

            // from M state
            if (t >= 3 && path[t - 3][hmm::State::Sr] != hmm::State::R)
                && (t >= 4 && path[t - 4][hmm::State::Sr] != hmm::State::R)
                && (t >= 5 && path[t - 5][hmm::State::Sr] != hmm::State::R)
            {
                let temp_alpha = alpha[t - 1][m]
                    - global.tr.mi
                    - global.tr_mi[from][to]
                    - if i == hmm::State::I6r {
                        global.tr.gg
                    } else {
                        0.0
                    };
                if temp_alpha < alpha[t][i] {
                    alpha[t][i] = temp_alpha;
                    path[t][i] = m;
                    temp_i_1[i as usize - hmm::State::I1r as usize] = t - 1;
                }
            }
        }

        // non_coding state
        // TODO just a minimum of three
        alpha[t][hmm::State::R] =
            alpha[t - 1][hmm::State::R] - local.tr_rr[from][to] - global.tr.rr;
        path[t][hmm::State::R] = hmm::State::R;

        let temp_alpha = alpha[t - 1][hmm::State::E] - global.tr.er;
        if temp_alpha < alpha[t][hmm::State::R] {
            alpha[t][hmm::State::R] = temp_alpha;
            path[t][hmm::State::R] = hmm::State::E;
        }

        let temp_alpha = alpha[t - 1][hmm::State::Er] - global.tr.er;
        if temp_alpha < alpha[t][hmm::State::R] {
            alpha[t][hmm::State::R] = temp_alpha;
            path[t][hmm::State::R] = hmm::State::Er;
        }

        alpha[t][hmm::State::R] -= 0.95_f64.ln();

        // end state
        if alpha[t][hmm::State::E] == 0.0 {
            alpha[t][hmm::State::E] = f64::INFINITY;
            path[t][hmm::State::E] = hmm::State::N;

            if t < seq.len() - 2
                && seq[t] == b'T'
                && ((seq[t + 1] == b'A' && seq[t + 2] == b'A')
                    || (seq[t + 1] == b'A' && seq[t + 2] == b'G')
                    || (seq[t + 1] == b'G' && seq[t + 2] == b'A'))
            {
                alpha[t + 2][hmm::State::E] = f64::INFINITY;

                // transition from frame4, frame5 and frame6
                let temp_alpha = alpha[t - 1][hmm::State::M6] - global.tr.ge;
                if temp_alpha < alpha[t + 2][hmm::State::E] {
                    alpha[t + 2][hmm::State::E] = temp_alpha;
                    path[t][hmm::State::E] = hmm::State::M6;
                }

                // transition from frame1, frame2 and frame3
                let temp_alpha = alpha[t - 1][hmm::State::M3] - global.tr.ge;
                if temp_alpha < alpha[t + 2][hmm::State::E] {
                    alpha[t + 2][hmm::State::E] = temp_alpha;
                    path[t][hmm::State::E] = hmm::State::M3;
                }

                alpha[t][hmm::State::E] = f64::INFINITY;
                alpha[t + 1][hmm::State::E] = f64::INFINITY;
                path[t + 1][hmm::State::E] = hmm::State::E;
                path[t + 2][hmm::State::E] = hmm::State::E;

                alpha[t + 2][hmm::State::M6] = f64::INFINITY;
                alpha[t + 1][hmm::State::M5] = f64::INFINITY;
                alpha[t][hmm::State::M4] = f64::INFINITY;
                alpha[t + 2][hmm::State::M3] = f64::INFINITY;
                alpha[t + 1][hmm::State::M2] = f64::INFINITY;
                alpha[t][hmm::State::M1] = f64::INFINITY;

                alpha[t + 2][hmm::State::E] -= if seq[t + 1] == b'A' && seq[t + 2] == b'A' {
                    0.54_f64.ln()
                } else if seq[t + 1] == b'A' && seq[t + 2] == b'G' {
                    0.16_f64.ln()
                } else {
                    0.30_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                let mut sub_sum = 0.0;

                if t >= 60 {
                    // bug reported by Yu-Wei ------ TODO 60 is the incomplete minimum length? can be merged?
                    for i in (t - 60)..=(t - 3) {
                        if i + 2 < seq.len() {
                            start_freq -= local.tr_e[i + 60 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2]).unwrap_or(0)];
                        }
                    }
                } else if t > 3 {
                    for i in 0..=(t - 3) {
                        if i + 2 < seq.len() {
                            sub_sum += local.tr_e[i + 60 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2]).unwrap_or(0)];
                        }
                    }
                    sub_sum = sub_sum * 58.0 / (t - 3 + 1) as f64;
                    start_freq -= sub_sum;
                }

                let h_kd = local.dist_e[2]
                    * (-1.0 * (start_freq - local.dist_e[1]).powi(2)
                        / (local.dist_e[0]).powi(2)
                        / 2.0)
                        .exp();
                let r_kd = local.dist_e[5]
                    * (-1.0 * (start_freq - local.dist_e[4]).powi(2)
                        / (local.dist_e[3]).powi(2)
                        / 2.0)
                        .exp();
                alpha[t + 2][hmm::State::E] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // start' state
        // originally stop codon of genes in - strand
        if alpha[t][hmm::State::Sr] == 0.0 {
            alpha[t][hmm::State::Sr] = f64::INFINITY;
            path[t][hmm::State::Sr] = hmm::State::N;

            if t < seq.len() - 2
                && seq[t + 2] == b'A'
                && ((seq[t + 1] == b'T' && seq[t] == b'T')
                    || (seq[t + 1] == b'T' && seq[t] == b'C')
                    || (seq[t + 1] == b'C' && seq[t] == b'T'))
            {
                alpha[t][hmm::State::Sr] = f64::INFINITY;
                alpha[t + 1][hmm::State::Sr] = f64::INFINITY;
                alpha[t + 2][hmm::State::Sr] = alpha[t - 1][hmm::State::R] - global.tr.rs;
                path[t][hmm::State::Sr] = hmm::State::R;
                path[t + 1][hmm::State::Sr] = hmm::State::Sr;
                path[t + 2][hmm::State::Sr] = hmm::State::Sr;

                let temp_alpha = alpha[t - 1][hmm::State::Er] - global.tr.es;
                if temp_alpha < alpha[t + 2][hmm::State::Sr] {
                    alpha[t + 2][hmm::State::Sr] = temp_alpha;
                    path[t][hmm::State::Sr] = hmm::State::Er;
                }

                let temp_alpha = alpha[t - 1][hmm::State::E] - global.tr.es1;
                if temp_alpha < alpha[t + 2][hmm::State::Sr] {
                    alpha[t + 2][hmm::State::Sr] = temp_alpha;
                    path[t][hmm::State::Sr] = hmm::State::E;
                }

                alpha[t + 2][hmm::State::M3r] = f64::INFINITY;
                alpha[t + 2][hmm::State::M6r] = f64::INFINITY;

                alpha[t + 2][hmm::State::Sr] -= if seq[t + 1] == b'T' && seq[t] == b'T' {
                    0.54_f64.ln()
                } else if seq[t + 1] == b'T' && seq[t] == b'C' {
                    0.16_f64.ln()
                } else {
                    0.30_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;

                // TODO needs same 60-edgecase as above?
                for i in 3..=60 {
                    if t + i + 2 < seq.len() {
                        start_freq += local.tr_s1[i - 3][trinucleotide(
                            seq[t + i],
                            seq[t + i + 1],
                            seq[t + i + 2],
                        )
                        .unwrap_or(0)];
                    }
                }

                let h_kd = local.dist_s1[2]
                    * (-1.0 * (start_freq - local.dist_s1[1]).powi(2)
                        / (local.dist_s1[0]).powi(2)
                        / 2.0)
                        .exp();
                let r_kd = local.dist_s1[5]
                    * (-1.0 * (start_freq - local.dist_s1[4]).powi(2)
                        / (local.dist_s1[3]).powi(2)
                        / 2.0)
                        .exp();
                alpha[t + 2][hmm::State::Sr] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // start state
        if alpha[t][hmm::State::S] == 0.0 {
            alpha[t][hmm::State::S] = f64::INFINITY;
            path[t][hmm::State::S] = hmm::State::N;

            if t < seq.len() - 2
                && (seq[t] == b'A' || seq[t] == b'G' || seq[t] == b'T')
                && seq[t + 1] == b'A'
                && seq[t + 2] == b'G'
            {
                alpha[t][hmm::State::S] = f64::INFINITY;
                alpha[t + 1][hmm::State::S] = f64::INFINITY;
                alpha[t + 2][hmm::State::S] = alpha[t - 1][hmm::State::R] - global.tr.rs;
                path[t][hmm::State::S] = hmm::State::R;
                path[t + 1][hmm::State::S] = hmm::State::S;
                path[t + 2][hmm::State::S] = hmm::State::S;

                let temp_alpha = alpha[t - 1][hmm::State::E] - global.tr.es;
                if temp_alpha < alpha[t + 2][hmm::State::S] {
                    alpha[t + 2][hmm::State::S] = temp_alpha;
                    path[t][hmm::State::S] = hmm::State::E;
                }

                let temp_alpha = alpha[t - 1][hmm::State::Er] - global.tr.es1;
                if temp_alpha < alpha[t + 2][hmm::State::S] {
                    alpha[t + 2][hmm::State::S] = temp_alpha;
                    path[t][hmm::State::S] = hmm::State::Er;
                }

                alpha[t + 2][hmm::State::S] -= if seq[t] == b'A' {
                    0.83_f64.ln()
                } else if seq[t] == b'G' {
                    0.10_f64.ln()
                } else {
                    0.07_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                let mut sub_sum = 0.0;

                if t >= 30 {
                    // TODO why 30?
                    for i in (t - 30)..=(t + 30) {
                        if i + 2 < seq.len() {
                            start_freq -= local.tr_s[i + 30 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2]).unwrap_or(0)];
                        }
                    }
                } else {
                    for i in 0..=(t + 30) {
                        if i + 2 < seq.len() {
                            sub_sum += local.tr_s[i + 30 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2]).unwrap_or(0)];
                        }
                    }
                    sub_sum *= 61.0 / (t + 30 + 1) as f64;
                    start_freq -= sub_sum;
                }

                let h_kd = local.dist_s[2]
                    * (-1.0 * (start_freq - local.dist_s[1]).powi(2)
                        / (local.dist_s[0]).powi(2)
                        / 2.0)
                        .exp();
                let r_kd = local.dist_s[5]
                    * (-1.0 * (start_freq - local.dist_s[4]).powi(2)
                        / (local.dist_s[3]).powi(2)
                        / 2.0)
                        .exp();
                alpha[t + 2][hmm::State::S] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // end' state
        // originally start codon of genes in - strand
        if alpha[t][hmm::State::Er] == 0.0 {
            alpha[t][hmm::State::Er] = f64::INFINITY;
            path[t][hmm::State::Er] = hmm::State::N;

            if t < seq.len() - 2
                && seq[t] == b'C'
                && seq[t + 1] == b'A'
                && (seq[t + 2] == b'T' || seq[t + 2] == b'C' || seq[t + 2] == b'A')
            {
                // transition from frame6
                alpha[t][hmm::State::Er] = f64::INFINITY;
                alpha[t + 1][hmm::State::Er] = f64::INFINITY;
                alpha[t + 2][hmm::State::Er] = alpha[t - 1][hmm::State::M6r] - global.tr.ge;
                path[t][hmm::State::Er] = hmm::State::M6r;
                path[t + 1][hmm::State::Er] = hmm::State::Er;
                path[t + 2][hmm::State::Er] = hmm::State::Er;

                alpha[t + 2][hmm::State::Er] -= if seq[t + 2] == b'T' {
                    0.83_f64.ln()
                } else if seq[t + 2] == b'C' {
                    0.10_f64.ln()
                } else {
                    0.07_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                let mut sub_sum = 0.0;

                if t >= 30 {
                    for i in (t - 30)..=(t + 30) {
                        if i + 2 < seq.len() {
                            start_freq -= local.tr_e1[i + 30 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2]).unwrap_or(0)];
                        }
                    }
                } else {
                    for i in 0..=(t + 30) {
                        if i + 2 < seq.len() {
                            sub_sum += local.tr_e1[i + 30 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2]).unwrap_or(0)];
                        }
                    }
                    sub_sum *= 61.0 / (t + 30 + 1) as f64;
                    start_freq -= sub_sum;
                }

                let h_kd = local.dist_e1[2]
                    * (-1.0 * (start_freq - local.dist_e1[1]).powi(2)
                        / (local.dist_e1[0]).powi(2)
                        / 2.0)
                        .exp();
                let r_kd = local.dist_e1[5]
                    * (-1.0 * (start_freq - local.dist_e1[4]).powi(2)
                        / (local.dist_e1[3]).powi(2)
                        / 2.0)
                        .exp();
                alpha[t + 2][hmm::State::Er] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        if num_noncoding > 9 {
            for i in hmm::states() {
                if i != hmm::State::R {
                    alpha[t][i] = f64::INFINITY;
                    path[t][i] = hmm::State::R;
                }
            }
        }
    }

    let vpath = backtrack(&alpha, path);
    output(
        &local,
        head,
        seq,
        whole_genome,
        formatted,
        vpath,
        gene_len,
        alpha,
    )
}

fn backtrack(
    alpha: &Vec<[f64; hmm::NUM_STATE]>,
    path: Vec<[hmm::State; hmm::NUM_STATE]>,
) -> Vec<hmm::State> {
    // backtrack array to find the optimal path
    let mut vpath: Vec<hmm::State> = vec![hmm::State::S]; // or null
    let mut prob = f64::INFINITY;
    for (&prob_, i) in alpha.last().expect("empty seq").iter().zip(hmm::states()) {
        if prob_ < prob {
            vpath[0] = i;
            prob = prob_;
        }
    }

    // backtrack the optimal path
    for t in (0..=path.len() - 2).rev() {
        vpath.push(path[t + 1][*vpath.last().unwrap()]);
    }
    vpath.reverse();
    vpath
}

fn output(
    local: &hmm::Local,
    head: Vec<u8>,
    seq: Vec<u8>,
    whole_genome: bool,
    formatted: bool,
    vpath: Vec<hmm::State>,
    gene_len: usize,
    alpha: Vec<[f64; hmm::NUM_STATE]>,
) -> Vec<gene::Gene> {
    let mut genes = vec![];
    let mut codon_start = 0; // ternaire boolean?
    let mut start_t: isize = -1;
    let mut dna_start_t_withstop: usize = 0;
    let mut dna_start_t: usize = 0;

    let mut dna = vec![];
    let mut dna_f = vec![];
    let mut insert = vec![];
    let mut delete = vec![];

    let mut prev_match = hmm::State::S; // or no state

    let mut start_orf = 0; // initialize?

    for t in 0..seq.len() {
        if codon_start == 0
            && start_t < 0
            && ((vpath[t] >= hmm::State::M1 && vpath[t] <= hmm::State::M6)
                || (vpath[t] >= hmm::State::M1r && vpath[t] <= hmm::State::M6r)
                || vpath[t] == hmm::State::S
                || vpath[t] == hmm::State::Sr)
        {
            dna_start_t_withstop = t + 1;
            dna_start_t = t + 1;
            start_t = (t + 1) as isize;
            // introduce dna_start_t_withstop YY July 2018
        }

        if codon_start == 0
            && (vpath[t] == hmm::State::M1
                || vpath[t] == hmm::State::M4
                || vpath[t] == hmm::State::M1r
                || vpath[t] == hmm::State::M4r)
        {
            dna.clear();
            dna_f.clear();
            insert.clear();
            delete.clear();

            dna.push(seq[t]);
            dna_f.push(seq[t]);
            dna_start_t_withstop = t + 1;
            dna_start_t = t + 1;
            if vpath[t] == hmm::State::M1 || vpath[t] == hmm::State::M4r {
                if t > 2 {
                    dna_start_t_withstop = t - 2;
                }
            }

            start_orf = t + 1;
            prev_match = vpath[t];

            codon_start = if vpath[t] < hmm::State::M6 { 1 } else { -1 }
        } else if codon_start != 0
            && (vpath[t] == hmm::State::E || vpath[t] == hmm::State::Er || t == seq.len() - 1)
        {
            let mut end_t;
            if vpath[t] == hmm::State::E || vpath[t] == hmm::State::Er {
                end_t = t + 3
            } else {
                let mut temp_t = t;
                while vpath[temp_t] != hmm::State::M1
                    && vpath[temp_t] != hmm::State::M4
                    && vpath[temp_t] != hmm::State::M1r
                    && vpath[temp_t] != hmm::State::M4r
                {
                    dna.pop();
                    dna_f.pop();
                    temp_t -= 1;
                }
                end_t = temp_t;
            }

            if dna.len() > gene_len {
                let final_score = (alpha[end_t - 4][vpath[end_t - 4]]
                    - alpha[start_t as usize + 2][vpath[start_t as usize + 2]])
                    / (end_t - start_t as usize - 5) as f64;
                let mut frame = start_orf % 3;
                if frame == 0 {
                    frame = 3
                }

                if codon_start == 1 {
                    if start_t == dna_start_t as isize - 3 {
                        dna_start_t -= 3;
                    }

                    if whole_genome {
                        // add refinement of the start codons here
                        let start_old = start_t as usize;
                        let mut codon = &seq[start_old..start_old + 3];
                        let mut s = 0;
                        // find the optimal start codon within 30bp up- and downstream of start codon
                        let mut e_save = 0.0;
                        let mut s_save = 0;
                        while !(codon != b"TAA" || codon != b"TAG" || codon != b"TGA")
                            && start_old >= 1 + s + 35
                        {
                            if codon != b"ATG" || codon != b"GTG" || codon != b"TTG" {
                                let utr = &seq[start_old - 1 - s - 30..start_old - 1 - s - 30 + 63];
                                let mut freq_sum = 0.0;
                                for j in 0..utr.len() - 2 {
                                    freq_sum -= local.tr_s[j][trinucleotide(
                                        utr[j],
                                        utr[j + 1],
                                        utr[j + 2],
                                    )
                                    .unwrap_or(0)];
                                }
                                if s == 0 {
                                    e_save = freq_sum;
                                    s_save = 0;
                                } else if freq_sum < e_save {
                                    e_save = freq_sum;
                                    s_save = -(s as isize); // posivite chain, upstream s_save = -1 * 3
                                }
                            }
                            s += 3;
                            codon = &seq[start_old - 1 - s..start_old - 1 - s + 3];

                            dna_start_t += s_save as usize;
                        }
                    }

                    genes.push(gene::Gene {
                        head: head.clone(),
                        start: dna_start_t,
                        metastart: dna_start_t,
                        end: end_t,
                        frame: frame,
                        score: final_score,
                        dna: dna.clone(),
                        dna_ffn: if formatted {
                            dna_f.clone()
                        } else {
                            dna.clone()
                        },
                        forward_strand: true,
                        inserted: insert.clone(),
                        deleted: delete.clone(),
                    });
                } else if codon_start == -1 {
                    if whole_genome {
                        // add refinement of the start codons here
                        // l 946
                        let end_old = end_t;
                        let mut codon = &seq[end_t - 1 - 2..end_t];
                        let mut s = 0;
                        // find the optimal start codon within 30bp up- and downstream of start codon
                        let mut e_save = 0.0;
                        let mut s_save = 0;
                        while !(codon != b"TAA" || codon != b"TAG" || codon != b"TGA")
                            && end_old - 2 + s + 35 < seq.len()
                        {
                            if codon != b"ATG" || codon != b"GTG" || codon != b"TTG" {
                                let utr = &seq[end_old - 1 - 2 + s - 30..end_old + s + 30];
                                let mut freq_sum = 0.0;
                                for j in 0..utr.len() - 2 {
                                    freq_sum -= local.tr_e1[j][trinucleotide(
                                        utr[j],
                                        utr[j + 1],
                                        utr[j + 2],
                                    )
                                    .unwrap_or(0)]; // TODO stop1? (their note)
                                }
                                if s == 0 || freq_sum < e_save {
                                    e_save = freq_sum;
                                    s_save = s; // negative chain, s_save = s
                                }
                            }
                            s += 3;
                            codon = &seq[end_old - 1 - 2 + s..end_old];
                        }

                        end_t = end_old + s_save;
                    }

                    genes.push(gene::Gene {
                        head: head.clone(),
                        start: dna_start_t_withstop,
                        metastart: dna_start_t,
                        end: end_t,
                        frame: frame,
                        score: final_score,
                        dna: dna.clone(),
                        dna_ffn: get_rc_dna(if formatted { &dna_f } else { &dna }),
                        forward_strand: false,
                        inserted: insert.clone(),
                        deleted: delete.clone(),
                    });
                }
            }

            codon_start = 0;
            start_t = -1;
        } else if codon_start != 0
            && ((vpath[t] >= hmm::State::M1
                && vpath[t] <= hmm::State::M6
                && prev_match >= hmm::State::M1)
                || (vpath[t] >= hmm::State::M1r
                    && vpath[t] <= hmm::State::M6r
                    && prev_match >= hmm::State::M1r))
        {
            let out_nt = if vpath[t] < prev_match {
                vpath[t] as usize + 6 - prev_match as usize
            } else {
                vpath[t] as usize - prev_match as usize
            };
            for kk in 0..out_nt {
                // for deleted nt in reads
                dna.push(b'N');
                dna_f.push(b'x');
                if kk > 0 {
                    delete.push(t + 1);
                }
            }
            dna.pop();
            dna_f.pop();
            dna.push(seq[t]);
            dna_f.push(seq[t]);
            prev_match = vpath[t];
        } else if codon_start != 0
            && ((vpath[t] >= hmm::State::I1 && vpath[t] <= hmm::State::I6)
                || (vpath[t] >= hmm::State::I1r && vpath[t] <= hmm::State::I6r))
        {
            dna_f.push(seq[t].to_ascii_lowercase());
            insert.push(t + 1);
        } else if codon_start != 0 && vpath[t] == hmm::State::R {
            // for long NNNNNNNN, pretend R state
            codon_start = 0;
            start_t = -1;
        }
    }

    genes
}
