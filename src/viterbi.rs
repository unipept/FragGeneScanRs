use crate::dna::Nuc::{A, C, G, T};
use crate::dna::{trinucleotide, Nuc};
use crate::{gene, hmm};

pub fn viterbi(
    global: &hmm::Global,
    local: &hmm::Local,
    head: Vec<u8>,
    seq: Vec<Nuc>,
    whole_genome: bool,
) -> Vec<gene::Gene> {
    let gene_len = if whole_genome { 120 } else { 60 }; // minimum length to be output

    let mut alpha: Vec<[f64; hmm::NUM_STATE]> = vec![];
    let mut path: Vec<[usize; hmm::NUM_STATE]> = vec![];
    let mut temp_i: [usize; 6] = [0; 6];
    let mut temp_i_1: [usize; 6] = [0; 6];

    for _ in 0..seq.len() {
        alpha.push([0.0; hmm::NUM_STATE]);
        path.push([hmm::NOSTATE; hmm::NUM_STATE]);
    }
    alpha[0].copy_from_slice(&global.pi);
    for i in &mut alpha[0] {
        *i *= -1.0
    } // TODO broadcast operation

    // If the sequence starts with a stop codon
    if seq[0] == T
        && ((seq[1] == A && seq[2] == A)
            || (seq[1] == A && seq[2] == G)
            || (seq[1] == G && seq[2] == A))
    {
        alpha[0][hmm::E_STATE] = f64::INFINITY;
        alpha[1][hmm::E_STATE] = f64::INFINITY;
        path[1][hmm::E_STATE] = hmm::E_STATE;
        path[2][hmm::E_STATE] = hmm::E_STATE;

        alpha[2][hmm::M6_STATE] = f64::INFINITY;
        alpha[1][hmm::M5_STATE] = f64::INFINITY;
        alpha[0][hmm::M4_STATE] = f64::INFINITY;
        alpha[2][hmm::M3_STATE] = f64::INFINITY;
        alpha[1][hmm::M2_STATE] = f64::INFINITY;
        alpha[0][hmm::M1_STATE] = f64::INFINITY;

        alpha[2][hmm::E_STATE] -= if seq[1] == A && seq[2] == A {
            0.53_f64.ln()
        } else if seq[1] == A && seq[2] == G {
            0.16_f64.ln()
        } else {
            0.30_f64.ln()
        }
    }

    // If the sequence starts with a reverse stop codon
    if seq[2] == A
        && ((seq[1] == T && seq[0] == T)
            || (seq[1] == T && seq[0] == C)
            || (seq[1] == C && seq[0] == T))
    {
        alpha[0][hmm::S_STATE_1] = f64::INFINITY;
        alpha[1][hmm::S_STATE_1] = f64::INFINITY;
        alpha[2][hmm::S_STATE_1] = alpha[0][hmm::S_STATE];
        path[1][hmm::S_STATE_1] = hmm::S_STATE_1;
        path[2][hmm::S_STATE_1] = hmm::S_STATE_1;

        alpha[2][hmm::M3_STATE_1] = f64::INFINITY;
        alpha[2][hmm::M6_STATE_1] = f64::INFINITY;

        alpha[2][hmm::S_STATE_1] = if seq[1] == T && seq[0] == T {
            0.53_f64.ln()
        } else if seq[1] == T && seq[0] == C {
            0.16_f64.ln()
        } else {
            0.30_f64.ln()
        }
    }

    let mut num_noncoding = 0; // number of invalid nts in sequence
    for t in 1..seq.len() {
        let from = (seq[t - 1]).to_int().unwrap_or(2);
        let from0 = if t > 1 {
            (seq[t - 2]).to_int().unwrap_or(2)
        } else {
            2
        };
        let to = (seq[t]).to_int().unwrap_or_else(|| {
            num_noncoding += 1;
            2
        });
        let from2 = from0 * 4 + from;

        // M state

        for i in hmm::M1_STATE..=hmm::M6_STATE {
            if alpha[t][i].is_finite() {
                if i == hmm::M1_STATE {
                    // from M state
                    alpha[t][i] = alpha[t - 1][hmm::M6_STATE]
                        - global.tr.gg
                        - global.tr.mm
                        - local.e_m[0][from2][to];
                    path[t][i] = hmm::M6_STATE;

                    // from D state
                    if !whole_genome {
                        for j in (hmm::M1_STATE..=hmm::M5_STATE).rev() {
                            let num_d = if j >= i {
                                (i + 6 - j) as i32
                            } else if j + 1 < i {
                                (i - j) as i32
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
                    let temp_alpha = alpha[t - 1][hmm::S_STATE] - local.e_m[0][from2][to];
                    if temp_alpha < alpha[t][i] {
                        alpha[t][i] = temp_alpha;
                        path[t][i] = hmm::S_STATE;
                    }
                } else {
                    // from M state
                    alpha[t][i] = alpha[t - 1][i - 1]
                        - global.tr.mm
                        - local.e_m[i - hmm::M1_STATE][from2][to];
                    path[t][i] = i - 1;

                    // from D state
                    if !whole_genome {
                        for j in (hmm::M1_STATE..=hmm::M6_STATE).rev() {
                            let num_d = if j >= i {
                                (i + 6 - j) as i32
                            } else if j + 1 < i {
                                (i - j) as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j]
                                    - global.tr.md
                                    - local.e_m[i - hmm::M1_STATE][from2][to]
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
                let j = if i == hmm::M1_STATE {
                    hmm::I6_STATE
                } else {
                    hmm::I1_STATE + (i - hmm::M1_STATE - 1)
                };

                // to avoid stop codon
                if t < 2 {
                } else if (i == hmm::M2_STATE || i == hmm::M5_STATE)
                    && t + 1 < seq.len()
                    && seq[temp_i[j - hmm::I1_STATE]] == T
                    && ((seq[t] == A && seq[t + 1] == A)
                        || (seq[t] == A && seq[t + 1] == G)
                        || (seq[t] == G && seq[t + 1] == A))
                {
                } else if (i == hmm::M3_STATE || i == hmm::M6_STATE)
                    && temp_i[j - hmm::I1_STATE] > 0
                    && (seq[temp_i[j - hmm::I1_STATE] - 1] == T)
                    && ((seq[temp_i[j - hmm::I1_STATE]] == A && seq[t] == A)
                        || (seq[temp_i[j - hmm::I1_STATE]] == A && seq[t] == G)
                        || (seq[temp_i[j - hmm::I1_STATE]] == G && seq[t] == A))
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
        for i in hmm::I1_STATE..=hmm::I6_STATE {
            // from I state
            alpha[t][i] = alpha[t - 1][i] - global.tr.ii - global.tr_ii[from][to];
            path[t][i] = i;

            // from M state
            let temp_alpha = alpha[t - 1][i - hmm::I1_STATE + hmm::M1_STATE]
                - global.tr.mi
                - global.tr_mi[from][to]
                - if i == hmm::I6_STATE {
                    global.tr.gg
                } else {
                    0.0
                };
            if temp_alpha < alpha[t][i] {
                alpha[t][i] = temp_alpha;
                path[t][i] = i - hmm::I1_STATE + hmm::M1_STATE;
                temp_i[i - hmm::I1_STATE] = t - 1;
            }
        }

        // M' state
        for i in hmm::M1_STATE_1..=hmm::M6_STATE_1 {
            if (i == hmm::M1_STATE_1 || i == hmm::M4_STATE_1)
                && t >= 3
                && seq[t - 1] == A
                && ((seq[t - 2] == T && seq[t - 3] == T)
                    || (seq[t - 2] == T && seq[t - 3] == C)
                    || (seq[t - 2] == C && seq[t - 3] == T))
            {
                // from Start state since this is actually a stop codon in minus strand
                alpha[t][i] =
                    alpha[t - 1][hmm::S_STATE_1] - local.e_m1[i - hmm::M1_STATE_1][from2][to];
                path[t][i] = hmm::S_STATE_1;
            } else {
                if i == hmm::M1_STATE_1 {
                    // from M state
                    alpha[t][i] = alpha[t - 1][hmm::M6_STATE_1]
                        - global.tr.gg
                        - global.tr.mm
                        - local.e_m1[0][from2][to];
                    path[t][i] = hmm::M6_STATE_1;

                    // from D state
                    if !whole_genome {
                        for j in (hmm::M1_STATE_1..=hmm::M5_STATE_1).rev() {
                            let num_d = if j >= i {
                                (i + 6 - j) as i32
                            } else if j + 1 < i {
                                (i - j) as i32
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
                    alpha[t][i] = alpha[t - 1][i - 1]
                        - global.tr.mm
                        - local.e_m1[i - hmm::M1_STATE_1][from2][to];
                    path[t][i] = i - 1;

                    // from D state
                    if !whole_genome {
                        for j in (hmm::M1_STATE_1..=hmm::M6_STATE_1).rev() {
                            let num_d = if j >= i {
                                (i + 6 - j) as i32
                            } else if j + 1 < i {
                                (i - j) as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j]
                                    - global.tr.md
                                    - local.e_m1[i - hmm::M1_STATE_1][from2][to]
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
                let j = if i == hmm::M1_STATE_1 {
                    hmm::I6_STATE_1
                } else {
                    hmm::I1_STATE_1 + i - hmm::M1_STATE_1 - 1
                };

                // to avoid stop codon
                if t < 2 {
                } else if (i == hmm::M2_STATE_1 || i == hmm::M5_STATE_1)
                    && t + 1 < seq.len()
                    && seq[t + 1] == A
                    && ((seq[t] == T && seq[temp_i_1[j - hmm::I1_STATE_1]] == T)
                        || (seq[t] == T && seq[temp_i_1[j - hmm::I1_STATE_1]] == C)
                        || (seq[t] == A && seq[temp_i_1[j - hmm::I1_STATE_1]] == T))
                {
                } else if (i == hmm::M3_STATE_1 || i == hmm::M6_STATE_1)
                    && seq[t] == A
                    && temp_i_1[j - hmm::I1_STATE_1] > 1
                    && ((seq[temp_i_1[j - hmm::I1_STATE_1]] == T
                        && seq[temp_i_1[j - hmm::I1_STATE_1] - 1] == T)
                        || (seq[temp_i_1[j - hmm::I1_STATE_1]] == T
                            && seq[temp_i_1[j - hmm::I1_STATE_1] - 1] == C)
                        || (seq[temp_i_1[j - hmm::I1_STATE_1]] == C
                            && seq[temp_i_1[j - hmm::I1_STATE_1] - 1] == T))
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
        for i in hmm::I1_STATE_1..=hmm::I6_STATE_1 {
            // from I state
            alpha[t][i] = alpha[t - 1][i] - global.tr.ii - global.tr_ii[from][to];
            path[t][i] = i;

            // from M state
            if (t >= 3 && path[t - 3][hmm::S_STATE_1] != hmm::R_STATE)
                && (t >= 4 && path[t - 4][hmm::S_STATE_1] != hmm::R_STATE)
                && (t >= 5 && path[t - 5][hmm::S_STATE_1] != hmm::R_STATE)
            {
                let temp_alpha = alpha[t - 1][i - hmm::I1_STATE_1 + hmm::M1_STATE_1]
                    - global.tr.mi
                    - global.tr_mi[from][to]
                    - if i == hmm::I6_STATE_1 {
                        global.tr.gg
                    } else {
                        0.0
                    };
                if temp_alpha < alpha[t][i] {
                    alpha[t][i] = temp_alpha;
                    path[t][i] = i - hmm::I1_STATE_1 + hmm::M1_STATE_1;
                    temp_i_1[i - hmm::I1_STATE_1] = t - 1;
                }
            }
        }

        // non_coding state
        // TODO just a minimum of three
        alpha[t][hmm::R_STATE] = alpha[t - 1][hmm::R_STATE] - local.tr_rr[from][to] - global.tr.rr;
        path[t][hmm::R_STATE] = hmm::R_STATE;

        let temp_alpha = alpha[t - 1][hmm::E_STATE] - global.tr.er;
        if temp_alpha < alpha[t][hmm::R_STATE] {
            alpha[t][hmm::R_STATE] = temp_alpha;
            path[t][hmm::R_STATE] = hmm::E_STATE;
        }

        let temp_alpha = alpha[t - 1][hmm::E_STATE_1] - global.tr.er;
        if temp_alpha < alpha[t][hmm::R_STATE] {
            alpha[t][hmm::R_STATE] = temp_alpha;
            path[t][hmm::R_STATE] = hmm::E_STATE_1;
        }

        alpha[t][hmm::R_STATE] -= 0.95_f64.ln();

        // end state
        if alpha[t][hmm::E_STATE] == 0.0 {
            alpha[t][hmm::E_STATE] = f64::INFINITY;
            path[t][hmm::E_STATE] = hmm::NOSTATE;

            if t < seq.len() - 2
                && seq[t] == T
                && ((seq[t + 1] == A && seq[t + 2] == A)
                    || (seq[t + 1] == A && seq[t + 2] == G)
                    || (seq[t + 1] == G && seq[t + 2] == A))
            {
                alpha[t + 2][hmm::E_STATE] = f64::INFINITY;

                // transition from frame4, frame5 and frame6
                let temp_alpha = alpha[t - 1][hmm::M6_STATE] - global.tr.ge;
                if temp_alpha < alpha[t + 2][hmm::E_STATE] {
                    alpha[t + 2][hmm::E_STATE] = temp_alpha;
                    path[t][hmm::E_STATE] = hmm::M6_STATE;
                }

                // transition from frame1, frame2 and frame3
                let temp_alpha = alpha[t - 1][hmm::M3_STATE] - global.tr.ge;
                if temp_alpha < alpha[t + 2][hmm::E_STATE] {
                    alpha[t + 2][hmm::E_STATE] = temp_alpha;
                    path[t][hmm::E_STATE] = hmm::M3_STATE;
                }

                alpha[t][hmm::E_STATE] = f64::INFINITY;
                alpha[t + 1][hmm::E_STATE] = f64::INFINITY;
                path[t + 1][hmm::E_STATE] = hmm::E_STATE;
                path[t + 2][hmm::E_STATE] = hmm::E_STATE;

                alpha[t + 2][hmm::M6_STATE] = f64::INFINITY;
                alpha[t + 1][hmm::M5_STATE] = f64::INFINITY;
                alpha[t][hmm::M4_STATE] = f64::INFINITY;
                alpha[t + 2][hmm::M3_STATE] = f64::INFINITY;
                alpha[t + 1][hmm::M2_STATE] = f64::INFINITY;
                alpha[t][hmm::M1_STATE] = f64::INFINITY;

                alpha[t + 2][hmm::E_STATE] -= if seq[t + 1] == A && seq[t + 2] == A {
                    0.54_f64.ln()
                } else if seq[t + 1] == A && seq[t + 2] == G {
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
                alpha[t + 2][hmm::E_STATE] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // start' state
        // originally stop codon of genes in - strand
        if alpha[t][hmm::S_STATE_1] == 0.0 {
            alpha[t][hmm::S_STATE_1] = f64::INFINITY;
            path[t][hmm::S_STATE_1] = hmm::NOSTATE;

            if t < seq.len() - 2
                && seq[t + 2] == A
                && ((seq[t + 1] == T && seq[t] == T)
                    || (seq[t + 1] == T && seq[t] == C)
                    || (seq[t + 1] == C && seq[t] == T))
            {
                alpha[t][hmm::S_STATE_1] = f64::INFINITY;
                alpha[t + 1][hmm::S_STATE_1] = f64::INFINITY;
                alpha[t + 2][hmm::S_STATE_1] = alpha[t - 1][hmm::R_STATE] - global.tr.rs;
                path[t][hmm::S_STATE_1] = hmm::R_STATE;
                path[t + 1][hmm::S_STATE_1] = hmm::S_STATE_1;
                path[t + 2][hmm::S_STATE_1] = hmm::S_STATE_1;

                let temp_alpha = alpha[t - 1][hmm::E_STATE_1] - global.tr.es;
                if temp_alpha < alpha[t + 2][hmm::S_STATE_1] {
                    alpha[t + 2][hmm::S_STATE_1] = temp_alpha;
                    path[t][hmm::S_STATE_1] = hmm::E_STATE_1;
                }

                let temp_alpha = alpha[t - 1][hmm::E_STATE] - global.tr.es1;
                if temp_alpha < alpha[t + 2][hmm::S_STATE_1] {
                    alpha[t + 2][hmm::S_STATE_1] = temp_alpha;
                    path[t][hmm::S_STATE_1] = hmm::E_STATE;
                }

                alpha[t + 2][hmm::M3_STATE_1] = f64::INFINITY;
                alpha[t + 2][hmm::M6_STATE_1] = f64::INFINITY;

                alpha[t + 2][hmm::S_STATE_1] -= if seq[t + 1] == T && seq[t] == T {
                    0.54_f64.ln()
                } else if seq[t + 1] == T && seq[t] == C {
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
                alpha[t + 2][hmm::S_STATE_1] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // start state
        if alpha[t][hmm::S_STATE] == 0.0 {
            alpha[t][hmm::S_STATE] = f64::INFINITY;
            path[t][hmm::S_STATE] = hmm::NOSTATE;

            if t < seq.len() - 2
                && (seq[t] == A || seq[t] == G || seq[t] == T)
                && seq[t + 1] == A
                && seq[t + 2] == G
            {
                alpha[t][hmm::S_STATE] = f64::INFINITY;
                alpha[t + 1][hmm::S_STATE] = f64::INFINITY;
                alpha[t + 2][hmm::S_STATE] = alpha[t - 1][hmm::R_STATE] - global.tr.rs;
                path[t][hmm::S_STATE] = hmm::R_STATE;
                path[t + 1][hmm::S_STATE] = hmm::S_STATE;
                path[t + 2][hmm::S_STATE] = hmm::S_STATE;

                let temp_alpha = alpha[t - 1][hmm::E_STATE] - global.tr.es;
                if temp_alpha < alpha[t + 2][hmm::S_STATE] {
                    alpha[t + 2][hmm::S_STATE] = temp_alpha;
                    path[t][hmm::S_STATE] = hmm::E_STATE;
                }

                let temp_alpha = alpha[t - 1][hmm::E_STATE_1] - global.tr.es1;
                if temp_alpha < alpha[t + 2][hmm::S_STATE] {
                    alpha[t + 2][hmm::S_STATE] = temp_alpha;
                    path[t][hmm::S_STATE] = hmm::E_STATE_1;
                }

                alpha[t + 2][hmm::S_STATE] -= if seq[t] == A {
                    0.83_f64.ln()
                } else if seq[t] == G {
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
                alpha[t + 2][hmm::S_STATE] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // end' state
        // originally start codon of genes in - strand
        if alpha[t][hmm::E_STATE_1] == 0.0 {
            alpha[t][hmm::E_STATE_1] = f64::INFINITY;
            path[t][hmm::E_STATE_1] = hmm::NOSTATE;

            if t < seq.len() - 2
                && seq[t] == C
                && seq[t + 1] == A
                && (seq[t + 2] == T || seq[t + 2] == C || seq[t + 2] == A)
            {
                // transition from frame6
                alpha[t][hmm::E_STATE_1] = f64::INFINITY;
                alpha[t + 1][hmm::E_STATE_1] = f64::INFINITY;
                alpha[t + 2][hmm::E_STATE_1] = alpha[t - 1][hmm::M6_STATE_1] - global.tr.ge;
                path[t][hmm::E_STATE_1] = hmm::M6_STATE_1;
                path[t + 1][hmm::E_STATE_1] = hmm::E_STATE_1;
                path[t + 2][hmm::E_STATE_1] = hmm::E_STATE_1;

                alpha[t + 2][hmm::E_STATE_1] -= if seq[t + 2] == T {
                    0.83_f64.ln()
                } else if seq[t + 2] == C {
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
                alpha[t + 2][hmm::E_STATE_1] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        if num_noncoding > 9 {
            for i in 0..hmm::NUM_STATE {
                if i != hmm::R_STATE {
                    alpha[t][i] = f64::INFINITY;
                    path[t][i] = hmm::R_STATE;
                }
            }
        }
    }

    let vpath = backtrack(&alpha, path);
    output(&local, head, seq, whole_genome, vpath, gene_len, alpha)
}

fn backtrack(alpha: &Vec<[f64; hmm::NUM_STATE]>, path: Vec<[usize; hmm::NUM_STATE]>) -> Vec<usize> {
    // backtrack array to find the optimal path
    let mut vpath: Vec<usize> = vec![0];
    let mut prob = f64::INFINITY;
    for (i, &prob_) in alpha.last().expect("empty seq").iter().enumerate() {
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
    seq: Vec<Nuc>,
    whole_genome: bool,
    vpath: Vec<usize>,
    gene_len: usize,
    alpha: Vec<[f64; hmm::NUM_STATE]>,
) -> Vec<gene::Gene> {
    let mut genes = vec![];
    let mut codon_start = 0; // ternaire boolean?
    let mut start_t: isize = -1;
    let mut dna_start_t_withstop: usize = 0;
    let mut dna_start_t: usize = 0;

    let mut dna: Vec<Nuc> = vec![];
    let mut insert = vec![];
    let mut delete = vec![];

    let mut prev_match = 0;

    let mut start_orf = 0; // initialize?

    for t in 0..seq.len() {
        if codon_start == 0
            && start_t < 0
            && ((vpath[t] >= hmm::M1_STATE && vpath[t] <= hmm::M6_STATE)
                || (vpath[t] >= hmm::M1_STATE_1 && vpath[t] <= hmm::M6_STATE_1)
                || vpath[t] == hmm::S_STATE
                || vpath[t] == hmm::S_STATE_1)
        {
            dna_start_t_withstop = t + 1;
            dna_start_t = t + 1;
            start_t = (t + 1) as isize;
            // introduce dna_start_t_withstop YY July 2018
        }

        if codon_start == 0
            && (vpath[t] == hmm::M1_STATE
                || vpath[t] == hmm::M4_STATE
                || vpath[t] == hmm::M1_STATE_1
                || vpath[t] == hmm::M4_STATE_1)
        {
            dna.clear();
            insert.clear();
            delete.clear();

            dna.push(seq[t]);
            dna_start_t_withstop = t + 1;
            dna_start_t = t + 1;
            if vpath[t] == hmm::M1_STATE || vpath[t] == hmm::M4_STATE_1 {
                if t > 2 {
                    dna_start_t_withstop = t - 2;
                }
            }

            start_orf = t + 1;
            prev_match = vpath[t];

            codon_start = if vpath[t] < hmm::M6_STATE { 1 } else { -1 }
        } else if codon_start != 0
            && (vpath[t] == hmm::E_STATE || vpath[t] == hmm::E_STATE_1 || t == seq.len() - 1)
        {
            let mut end_t;
            if vpath[t] == hmm::E_STATE || vpath[t] == hmm::E_STATE_1 {
                end_t = t + 3
            } else {
                let mut temp_t = t;
                while vpath[temp_t] != hmm::M1_STATE
                    && vpath[temp_t] != hmm::M4_STATE
                    && vpath[temp_t] != hmm::M1_STATE_1
                    && vpath[temp_t] != hmm::M4_STATE_1
                {
                    dna.pop();
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
                        while !(codon != [T, A, A] || codon != [T, A, G] || codon != [T, G, A])
                            && start_old >= 1 + s + 35
                        {
                            if codon != [A, T, G] || codon != [G, T, G] || codon != [T, T, G] {
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
                        while !(codon != [T, A, A] || codon != [T, A, G] || codon != [T, G, A])
                            && end_old - 2 + s + 35 < seq.len()
                        {
                            if codon != [A, T, G] || codon != [G, T, G] || codon != [T, T, G] {
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
                        forward_strand: false,
                        inserted: insert.clone(),
                        deleted: delete.clone(),
                    });
                }
            }

            codon_start = 0;
            start_t = -1;
        } else if codon_start != 0
            && ((vpath[t] >= hmm::M1_STATE && vpath[t] <= hmm::M6_STATE)
                || (vpath[t] >= hmm::M1_STATE_1 && vpath[t] <= hmm::M6_STATE_1))
            && vpath[t] < prev_match + 6
        {
            let out_nt = if vpath[t] < prev_match {
                vpath[t] + 6 - prev_match
            } else {
                vpath[t] - prev_match
            };
            for kk in 0..out_nt {
                // for deleted nt in reads
                dna.push(Nuc::Xi);
                if kk > 0 {
                    delete.push(t + 1);
                }
            }
            dna.pop();
            dna.push(seq[t]);
            prev_match = vpath[t];
        } else if codon_start != 0
            && ((vpath[t] >= hmm::I1_STATE && vpath[t] <= hmm::I6_STATE)
                || (vpath[t] >= hmm::I1_STATE_1 && vpath[t] <= hmm::I6_STATE_1))
        {
            dna.push(seq[t].to_lower());
            insert.push(t + 1);
        } else if codon_start != 0 && vpath[t] == hmm::R_STATE {
            // for long NNNNNNNN, pretend R state
            codon_start = 0;
            start_t = -1;
        }
    }

    genes
}
