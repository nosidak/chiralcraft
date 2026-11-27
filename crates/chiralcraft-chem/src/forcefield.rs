//! Pure-Rust 3D Molecular Mechanics Geometry Optimizer (UFF/MMFF-lite)
//! and 3D Bond Distance / Angle / Dihedral Measurement Engine.
//!
//! Minimizes an analytical molecular mechanics potential in Ångström space:
//! - Harmonic bond stretching (`k_b (r - r0)^2`)
//! - Harmonic valence angle bending (`k_a (theta - theta0)^2`)
//! - Out-of-plane planarity for sp²/aromatic centers
//! - Soft-core non-bonded steric repulsion (`1–4` and `>=1–5` pairs)

use crate::{BondDisplay, CcAromaticityInfo, CcVal, ChemAtom, ChemMolecule};

#[derive(Debug, Clone, PartialEq)]
pub struct ExplicitHResult {
    pub molecule: ChemMolecule,
    /// Maps each atom index in `molecule` to its parent heavy-atom index in the source molecule.
    /// Heavy atoms map to `Some(i)`; newly added explicit hydrogens map to `Some(parent_heavy_idx)`.
    pub parent_heavy_atom: Vec<Option<usize>>,
    pub added_h_count: usize,
}

/// Computes net formal charge across all atoms in `mol`.
pub fn total_formal_charge(mol: &ChemMolecule) -> i32 {
    mol.atoms.iter().map(|a| a.charge).sum()
}

fn v_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn v_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn v_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn v_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn v_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn v_norm(a: [f64; 3]) -> f64 {
    v_dot(a, a).sqrt()
}

fn v_unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = v_norm(a);
    if n > 1e-8 {
        Some([a[0] / n, a[1] / n, a[2] / n])
    } else {
        None
    }
}

fn perp_unit(u: [f64; 3]) -> [f64; 3] {
    let ref_v = if u[2].abs() < 0.85 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    v_unit(v_cross(u, ref_v)).unwrap_or([1.0, 0.0, 0.0])
}

/// Converts implicit hydrogens (via `CcVal`) into explicit 3D hydrogen atoms with chemically
/// reasonable tetrahedral / trigonal / linear geometry suitable for UFF, xTB, ORCA, and Gaussian.
pub fn add_explicit_hydrogens_3d(mol: &ChemMolecule) -> ExplicitHResult {
    let n = mol.atoms.len();
    let mut out = mol.clone();
    let mut parent_heavy_atom: Vec<Option<usize>> = (0..n).map(Some).collect();
    if n == 0 {
        return ExplicitHResult {
            molecule: out,
            parent_heavy_atom,
            added_h_count: 0,
        };
    }

    let val = CcVal::compute(mol);
    let arom = CcAromaticityInfo::perceive(mol);
    let mut adj: Vec<Vec<(usize, i32)>> = vec![Vec::new(); n];
    for b in &mol.bonds {
        if b.begin_atom < n && b.end_atom < n && b.begin_atom != b.end_atom {
            adj[b.begin_atom].push((b.end_atom, b.order));
            adj[b.end_atom].push((b.begin_atom, b.order));
        }
    }

    let mut next_atom_id = mol.atoms.iter().map(|a| a.id).max().unwrap_or(0).saturating_add(1);
    let mut next_bond_id = mol.bonds.iter().map(|b| b.id).max().unwrap_or(0).saturating_add(1);
    let mut added_h_count = 0usize;

    for i in 0..n {
        let imp_h = val.implicit_h(i) as usize;
        if imp_h == 0 {
            continue;
        }
        let atom = &mol.atoms[i];
        let p0 = [atom.x, atom.y, atom.z_coord];
        let xh_len_ang = match atom.z {
            6 => 1.09,
            7 => 1.01,
            8 => 0.96,
            15 => 1.42,
            16 => 1.34,
            14 => 1.48,
            5 => 1.19,
            _ => ideal_bond_len_ang(atom.z, 1, 1, false),
        };
        let xh_len_pt = xh_len_ang * PT_PER_ANG;

        let has_triple = adj[i].iter().any(|&(_, o)| o >= 3);
        let has_double = adj[i].iter().any(|&(_, o)| o == 2) || arom.is_atom_aromatic(i);
        let dirs: Vec<[f64; 3]> = adj[i]
            .iter()
            .filter_map(|&(nb, _)| {
                let nb_a = &mol.atoms[nb];
                v_unit(v_sub([nb_a.x, nb_a.y, nb_a.z_coord], p0))
            })
            .collect();

        let h_dirs: Vec<[f64; 3]> = match dirs.len() {
            0 => {
                // Isolated atom: CH4 / NH3 / H2O / HF
                let tet = [
                    [1.0, 1.0, 1.0],
                    [1.0, -1.0, -1.0],
                    [-1.0, 1.0, -1.0],
                    [-1.0, -1.0, 1.0],
                ];
                let mut out_d = Vec::with_capacity(imp_h);
                for k in 0..imp_h {
                    let raw = if k < 4 {
                        tet[k]
                    } else {
                        let ang = (k as f64) * 1.25;
                        [ang.cos(), ang.sin(), 0.25]
                    };
                    out_d.push(v_unit(raw).unwrap_or([1.0, 0.0, 0.0]));
                }
                out_d
            }
            1 => {
                let u = dirs[0];
                let opp = v_scale(u, -1.0);
                let e1 = perp_unit(u);
                let e2 = v_unit(v_cross(u, e1)).unwrap_or([0.0, 1.0, 0.0]);
                if imp_h == 1 {
                    if has_triple {
                        vec![opp]
                    } else {
                        // ~115° (O-H / sp2 C-H)
                        let theta = if has_double { 120.0_f64 } else { 109.5_f64 }.to_radians();
                        let ca = (-theta.cos()).clamp(0.0, 1.0);
                        let sa = theta.sin();
                        vec![v_unit(v_add(v_scale(opp, ca), v_scale(e1, sa))).unwrap_or(opp)]
                    }
                } else if imp_h == 2 {
                    if has_double {
                        // Terminal =CH2 planar 120°
                        let ca = 0.5_f64;
                        let sa = (3.0_f64).sqrt() * 0.5;
                        vec![
                            v_unit(v_add(v_scale(opp, ca), v_scale(e1, sa))).unwrap_or(opp),
                            v_unit(v_add(v_scale(opp, ca), v_scale(e1, -sa))).unwrap_or(opp),
                        ]
                    } else {
                        // Pyramidal -NH2 / -PH2 (~109.5°)
                        let theta = 109.47_f64.to_radians();
                        let ca = -theta.cos();
                        let sa = theta.sin();
                        let d1 = v_add(v_scale(opp, ca), v_add(v_scale(e1, sa * 0.866), v_scale(e2, sa * 0.5)));
                        let d2 = v_add(v_scale(opp, ca), v_add(v_scale(e1, -sa * 0.866), v_scale(e2, sa * 0.5)));
                        vec![v_unit(d1).unwrap_or(opp), v_unit(d2).unwrap_or(opp)]
                    }
                } else {
                    // Methyl -CH3 / -NH3+ tetrahedral cone around opp
                    let theta = 109.47_f64.to_radians();
                    let ca = -theta.cos();
                    let sa = theta.sin();
                    (0..imp_h)
                        .map(|k| {
                            let phi = (k as f64) * std::f64::consts::TAU / (imp_h as f64);
                            let radial = v_add(v_scale(e1, phi.cos()), v_scale(e2, phi.sin()));
                            v_unit(v_add(v_scale(opp, ca), v_scale(radial, sa))).unwrap_or(opp)
                        })
                        .collect()
                }
            }
            2 => {
                let sum = v_add(dirs[0], dirs[1]);
                let bis = v_unit(v_scale(sum, -1.0)).unwrap_or_else(|| perp_unit(dirs[0]));
                let norm = v_unit(v_cross(dirs[0], dirs[1])).unwrap_or_else(|| perp_unit(bis));
                if imp_h == 1 {
                    if has_double {
                        vec![bis]
                    } else {
                        // Slightly out-of-plane for sp3 secondary amine / methine
                        let d = v_add(v_scale(bis, 0.82), v_scale(norm, 0.57));
                        vec![v_unit(d).unwrap_or(bis)]
                    }
                } else {
                    // Methylene -CH2- tetrahedral pair
                    let half = (109.47_f64 * 0.5).to_radians();
                    let ch = half.cos();
                    let sh = half.sin();
                    let d1 = v_add(v_scale(bis, ch), v_scale(norm, sh));
                    let d2 = v_add(v_scale(bis, ch), v_scale(norm, -sh));
                    let mut out_d = vec![v_unit(d1).unwrap_or(bis), v_unit(d2).unwrap_or(bis)];
                    while out_d.len() < imp_h {
                        out_d.push(bis);
                    }
                    out_d
                }
            }
            _ => {
                // 3+ neighbors (methine C-H): place H opposite vector sum of neighbors
                let mut sum = [0.0; 3];
                for d in &dirs {
                    sum = v_add(sum, *d);
                }
                let dir = if v_norm(sum) > 0.15 {
                    v_unit(v_scale(sum, -1.0)).unwrap_or([0.0, 0.0, 1.0])
                } else {
                    // Neighbors nearly coplanar: pick normal to first two bonds
                    let n_vec = v_unit(v_cross(dirs[0], dirs[1])).unwrap_or([0.0, 0.0, 1.0]);
                    if p0[2] >= 0.0 { n_vec } else { v_scale(n_vec, -1.0) }
                };
                vec![dir; imp_h]
            }
        };

        for d in h_dirs {
            let mut h_atom = ChemAtom::new(next_atom_id, 1);
            next_atom_id += 1;
            h_atom.x = p0[0] + d[0] * xh_len_pt;
            h_atom.y = p0[1] + d[1] * xh_len_pt;
            h_atom.z_coord = p0[2] + d[2] * xh_len_pt;
            let h_idx = out.add_atom(h_atom);
            out.add_bond(next_bond_id, i, h_idx, 1, BondDisplay::Solid);
            next_bond_id += 1;
            parent_heavy_atom.push(Some(i));
            added_h_count += 1;
        }
        // Keep explicit_num_hydrogens neutral since hydrogens are now materialized as graph nodes.
        out.atoms[i].explicit_num_hydrogens = None;
    }

    out.has_3d_coords = true;
    ExplicitHResult {
        molecule: out,
        parent_heavy_atom,
        added_h_count,
    }
}

/// Prepares a 2D or 3D `ChemMolecule` for quantum chemistry (xTB / ORCA / Gaussian):
/// optionally materializes explicit 3D hydrogens and runs UFF-lite geometry optimization.
pub fn prepare_molecule_for_quantum(
    mol: &ChemMolecule,
    add_explicit_h: bool,
    pre_optimize_uff: bool,
) -> ExplicitHResult {
    let mut res = if add_explicit_h {
        add_explicit_hydrogens_3d(mol)
    } else {
        ExplicitHResult {
            molecule: mol.clone(),
            parent_heavy_atom: (0..mol.atoms.len()).map(Some).collect(),
            added_h_count: 0,
        }
    };
    if pre_optimize_uff && res.molecule.atoms.len() >= 2 {
        let (opt, _) = optimize_geometry_3d(&res.molecule, 180);
        res.molecule = opt;
    }
    res
}

#[derive(Debug, Clone, PartialEq)]
pub struct ForceFieldReport {
    pub initial_energy_kcal: f64,
    pub final_energy_kcal: f64,
    pub iterations: usize,
    pub converged: bool,
}

fn cov_rad(z: i32) -> f64 {
    match z {
        1 => 0.31,
        5 => 0.84,
        6 => 0.76,
        7 => 0.71,
        8 => 0.66,
        9 => 0.57,
        14 => 1.11,
        15 => 1.07,
        16 => 1.05,
        17 => 1.02,
        35 => 1.20,
        53 => 1.39,
        _ => 1.25,
    }
}

fn ideal_bond_len_ang(z1: i32, z2: i32, order: i32, arom: bool) -> f64 {
    let r0 = cov_rad(z1) + cov_rad(z2);
    if arom {
        r0 * 0.915
    } else if order >= 3 {
        r0 * 0.79
    } else if order == 2 {
        r0 * 0.88
    } else {
        r0
    }
}

/// Converts ChiralCraft screen-pt coordinates (`14.4 pt ~ 1.50 Å`) into Ångströms and back.
const PT_PER_ANG: f64 = 14.40 / 1.50;

pub fn measure_distance_angstrom(mol: &ChemMolecule, i: usize, j: usize) -> Option<f64> {
    let a = mol.atoms.get(i)?;
    let b = mol.atoms.get(j)?;
    let dx = (b.x - a.x) / PT_PER_ANG;
    let dy = (b.y - a.y) / PT_PER_ANG;
    let dz = (b.z_coord - a.z_coord) / PT_PER_ANG;
    Some((dx * dx + dy * dy + dz * dz).sqrt())
}

pub fn measure_angle_deg(mol: &ChemMolecule, i: usize, j: usize, k: usize) -> Option<f64> {
    let a = mol.atoms.get(i)?;
    let b = mol.atoms.get(j)?;
    let c = mol.atoms.get(k)?;
    let v1 = (a.x - b.x, a.y - b.y, a.z_coord - b.z_coord);
    let v2 = (c.x - b.x, c.y - b.y, c.z_coord - b.z_coord);
    let n1 = (v1.0 * v1.0 + v1.1 * v1.1 + v1.2 * v1.2).sqrt();
    let n2 = (v2.0 * v2.0 + v2.1 * v2.1 + v2.2 * v2.2).sqrt();
    if n1 < 1e-6 || n2 < 1e-6 {
        return None;
    }
    let cos_t = ((v1.0 * v2.0 + v1.1 * v2.1 + v1.2 * v2.2) / (n1 * n2)).clamp(-1.0, 1.0);
    Some(cos_t.acos().to_degrees())
}

pub fn measure_dihedral_deg(mol: &ChemMolecule, i: usize, j: usize, k: usize, l: usize) -> Option<f64> {
    let p0 = mol.atoms.get(i)?;
    let p1 = mol.atoms.get(j)?;
    let p2 = mol.atoms.get(k)?;
    let p3 = mol.atoms.get(l)?;
    let sub = |a: &crate::ChemAtom, b: &crate::ChemAtom| (a.x - b.x, a.y - b.y, a.z_coord - b.z_coord);
    let cross = |u: (f64, f64, f64), v: (f64, f64, f64)| (u.1 * v.2 - u.2 * v.1, u.2 * v.0 - u.0 * v.2, u.0 * v.1 - u.1 * v.0);
    let dot = |u: (f64, f64, f64), v: (f64, f64, f64)| u.0 * v.0 + u.1 * v.1 + u.2 * v.2;
    let norm = |u: (f64, f64, f64)| dot(u, u).sqrt();

    let b1 = sub(p1, p0);
    let b2 = sub(p2, p1);
    let b3 = sub(p3, p2);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let b2n = norm(b2);
    if norm(n1) < 1e-6 || norm(n2) < 1e-6 || b2n < 1e-6 {
        return None;
    }
    let b2u = (b2.0 / b2n, b2.1 / b2n, b2.2 / b2n);
    let m1 = cross(n1, b2u);
    let x = dot(n1, n2);
    let y = dot(m1, n2);
    Some((-y).atan2(x).to_degrees())
}

/// Runs pure-Rust UFF-lite 3D molecular mechanics optimization with adaptive steepest descent.
pub fn optimize_geometry_3d(mol: &ChemMolecule, max_iters: usize) -> (ChemMolecule, ForceFieldReport) {
    let n = mol.atoms.len();
    if n < 2 {
        return (
            mol.clone(),
            ForceFieldReport {
                initial_energy_kcal: 0.0,
                final_energy_kcal: 0.0,
                iterations: 0,
                converged: true,
            },
        );
    }

    let arom = CcAromaticityInfo::perceive(mol);
    let mut adj: Vec<Vec<(usize, i32)>> = vec![Vec::new(); n];
    let mut bonded = vec![vec![false; n]; n];
    let mut bonds_info: Vec<(usize, usize, f64, f64)> = Vec::new();
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.begin_atom < n && b.end_atom < n && b.begin_atom != b.end_atom {
            let (u, v) = (b.begin_atom, b.end_atom);
            adj[u].push((v, b.order));
            adj[v].push((u, b.order));
            bonded[u][v] = true;
            bonded[v][u] = true;
            let r0 = ideal_bond_len_ang(mol.atoms[u].z, mol.atoms[v].z, b.order, arom.is_bond_aromatic(bi));
            let kb = if b.order >= 2 { 520.0 } else { 320.0 };
            bonds_info.push((u, v, r0, kb));
        }
    }

    // Build angle triplets (i - j - k) and sp2 planarity quadruplets
    let mut angles: Vec<(usize, usize, usize, f64, f64)> = Vec::new();
    let mut planar_centers: Vec<(usize, usize, usize, usize)> = Vec::new();
    let mut pair_13 = vec![vec![false; n]; n];
    for j in 0..n {
        let deg = adj[j].len();
        let has_triple = adj[j].iter().any(|&(_, o)| o >= 3);
        let has_double = adj[j].iter().any(|&(_, o)| o == 2) || arom.is_atom_aromatic(j);
        let theta0 = if has_triple {
            std::f64::consts::PI
        } else if has_double {
            120.0_f64.to_radians()
        } else {
            109.47_f64.to_radians()
        };
        for a in 0..deg {
            for b in (a + 1)..deg {
                let i = adj[j][a].0;
                let k = adj[j][b].0;
                pair_13[i][k] = true;
                pair_13[k][i] = true;
                angles.push((i, j, k, theta0, 55.0));
            }
        }
        if deg == 3 && has_double {
            planar_centers.push((j, adj[j][0].0, adj[j][1].0, adj[j][2].0));
        }
    }

    // Center & convert to Ångströms, adding deterministic tiny Z perturbation if completely flat
    let nf = n as f64;
    let cx = mol.atoms.iter().map(|a| a.x).sum::<f64>() / nf;
    let cy = mol.atoms.iter().map(|a| a.y).sum::<f64>() / nf;
    let cz = mol.atoms.iter().map(|a| a.z_coord).sum::<f64>() / nf;
    let is_flat = mol.atoms.iter().all(|a| (a.z_coord - cz).abs() < 1e-3);
    let mut pos: Vec<[f64; 3]> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let z_init = if is_flat && !arom.is_atom_aromatic(i) && !adj[i].iter().any(|&(_, o)| o >= 2) {
                if i % 2 == 0 { 0.32 } else { -0.32 }
            } else {
                (a.z_coord - cz) / PT_PER_ANG
            };
            [(a.x - cx) / PT_PER_ANG, (a.y - cy) / PT_PER_ANG, z_init]
        })
        .collect();

    let eval_energy_and_grad = |p: &[[f64; 3]], grad: &mut [[f64; 3]]| -> f64 {
        for g in grad.iter_mut() {
            *g = [0.0; 3];
        }
        let mut e = 0.0_f64;

        // 1. Harmonic bond stretching
        for &(u, v, r0, kb) in &bonds_info {
            let dx = p[v][0] - p[u][0];
            let dy = p[v][1] - p[u][1];
            let dz = p[v][2] - p[u][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-5);
            let dr = r - r0;
            e += 0.5 * kb * dr * dr;
            let f = kb * dr / r;
            grad[u][0] -= f * dx; grad[u][1] -= f * dy; grad[u][2] -= f * dz;
            grad[v][0] += f * dx; grad[v][1] += f * dy; grad[v][2] += f * dz;
        }

        // 2. Cosine-harmonic angle bending: 0.5 * ka * (cos(theta) - cos(theta0))^2
        for &(i, j, k, t0, ka) in &angles {
            let u = [p[i][0] - p[j][0], p[i][1] - p[j][1], p[i][2] - p[j][2]];
            let v = [p[k][0] - p[j][0], p[k][1] - p[j][1], p[k][2] - p[j][2]];
            let ru = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt().max(1e-5);
            let rv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-5);
            let dot = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            let c = (dot / (ru * rv)).clamp(-1.0, 1.0);
            let dc = c - t0.cos();
            e += 0.5 * ka * dc * dc;
            let coeff = ka * dc;
            for d in 0..3 {
                let dci = (v[d] / (ru * rv)) - c * (u[d] / (ru * ru));
                let dck = (u[d] / (ru * rv)) - c * (v[d] / (rv * rv));
                grad[i][d] += coeff * dci;
                grad[k][d] += coeff * dck;
                grad[j][d] -= coeff * (dci + dck);
            }
        }

        // 3. Out-of-plane planarity for sp2/aromatic 3-coordinate centers (improper pyramid volume^2)
        for &(c, a, b, d) in &planar_centers {
            let v1 = [p[a][0] - p[c][0], p[a][1] - p[c][1], p[a][2] - p[c][2]];
            let v2 = [p[b][0] - p[c][0], p[b][1] - p[c][1], p[b][2] - p[c][2]];
            let v3 = [p[d][0] - p[c][0], p[d][1] - p[c][1], p[d][2] - p[c][2]];
            let n23 = [v2[1] * v3[2] - v2[2] * v3[1], v2[2] * v3[0] - v2[0] * v3[2], v2[0] * v3[1] - v2[1] * v3[0]];
            let vol = v1[0] * n23[0] + v1[1] * n23[1] + v1[2] * n23[2];
            let kp = 28.0;
            e += 0.5 * kp * vol * vol;
            let n31 = [v3[1] * v1[2] - v3[2] * v1[1], v3[2] * v1[0] - v3[0] * v1[2], v3[0] * v1[1] - v3[1] * v1[0]];
            let n12 = [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]];
            let s = kp * vol;
            for m in 0..3 {
                grad[a][m] += s * n23[m];
                grad[b][m] += s * n31[m];
                grad[d][m] += s * n12[m];
                grad[c][m] -= s * (n23[m] + n31[m] + n12[m]);
            }
        }

        // 4. Soft-core steric repulsion for non-bonded pairs (distance < r_min)
        for i in 0..n {
            for j in (i + 1)..n {
                if bonded[i][j] || pair_13[i][j] {
                    continue;
                }
                let dx = p[j][0] - p[i][0];
                let dy = p[j][1] - p[i][1];
                let dz = p[j][2] - p[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt().max(0.25);
                let r_min = match (mol.atoms[i].z == 1, mol.atoms[j].z == 1) {
                    (true, true) => 1.55,
                    (true, false) | (false, true) => 1.85,
                    (false, false) => 2.45,
                };
                if r < r_min {
                    let dr = r - r_min;
                    let krep = 14.0;
                    e += 0.5 * krep * dr * dr;
                    let f = krep * dr / r;
                    grad[i][0] -= f * dx; grad[i][1] -= f * dy; grad[i][2] -= f * dz;
                    grad[j][0] += f * dx; grad[j][1] += f * dy; grad[j][2] += f * dz;
                }
            }
        }
        e
    };

    let mut grad = vec![[0.0; 3]; n];
    let init_e = eval_energy_and_grad(&pos, &mut grad);
    let mut cur_e = init_e;
    let mut step = 0.0025_f64;
    let mut iters = 0usize;
    let mut converged = false;

    for it in 0..max_iters.clamp(10, 400) {
        iters = it + 1;
        let max_g = grad
            .iter()
            .map(|g| (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt())
            .fold(0.0_f64, f64::max);
        if max_g < 0.08 {
            converged = true;
            break;
        }
        let cap_scale = if max_g * step > 0.12 { 0.12 / (max_g * step) } else { 1.0 };
        let mut trial = pos.clone();
        for i in 0..n {
            for d in 0..3 {
                trial[i][d] -= step * cap_scale * grad[i][d];
            }
        }
        let mut trial_grad = vec![[0.0; 3]; n];
        let trial_e = eval_energy_and_grad(&trial, &mut trial_grad);
        if trial_e <= cur_e {
            pos = trial;
            grad = trial_grad;
            cur_e = trial_e;
            step = (step * 1.12).min(0.015);
        } else {
            step *= 0.45;
            if step < 1e-6 {
                break;
            }
        }
    }

    let mut out = mol.clone();
    for (i, a) in out.atoms.iter_mut().enumerate() {
        a.x = cx + pos[i][0] * PT_PER_ANG;
        a.y = cy + pos[i][1] * PT_PER_ANG;
        a.z_coord = pos[i][2] * PT_PER_ANG;
    }
    out.has_3d_coords = true;
    (
        out,
        ForceFieldReport {
            initial_energy_kcal: init_e,
            final_energy_kcal: cur_e,
            iterations: iters,
            converged,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondDisplay, ChemAtom};

    #[test]
    fn test_uff_lite_3d_optimizer_reduces_energy_and_measures_geometry() {
        // Distorted propane: C-C stretched way too far!
        let mut m = ChemMolecule::new();
        let mut a0 = ChemAtom::new(1, 6); a0.x = 0.0; a0.y = 0.0;
        let mut a1 = ChemAtom::new(2, 6); a1.x = 30.0; a1.y = 0.0; // ~3.1 Å (stretched!)
        let mut a2 = ChemAtom::new(3, 6); a2.x = 55.0; a2.y = 4.0;
        m.add_atom(a0); m.add_atom(a1); m.add_atom(a2);
        m.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        m.add_bond(11, 1, 2, 1, BondDisplay::Solid);

        let (opt, rep) = optimize_geometry_3d(&m, 160);
        assert!(rep.final_energy_kcal < rep.initial_energy_kcal * 0.1, "Energy must drop sharply: {:?}", rep);
        let d01 = measure_distance_angstrom(&opt, 0, 1).unwrap();
        assert!((d01 - 1.52).abs() < 0.08, "C-C bond should relax to ~1.52 Å, got {d01:.3}");
        let ang = measure_angle_deg(&opt, 0, 1, 2).unwrap();
        assert!((ang - 109.5).abs() < 4.0, "C-C-C sp3 angle should relax to ~109.5°, got {ang:.1}°");
    }

    #[test]
    fn test_explicit_h_3d_builder_and_quantum_prep() {
        // Ethanol: C-C-O -> 3 heavy + 6 explicit H = 9 atoms
        let mut eth = ChemMolecule::new();
        let mut c0 = ChemAtom::new(1, 6); c0.x = 0.0; c0.y = 0.0;
        let mut c1 = ChemAtom::new(2, 6); c1.x = 14.4; c1.y = 0.0;
        let mut o2 = ChemAtom::new(3, 8); o2.x = 22.0; o2.y = 12.0;
        eth.add_atom(c0); eth.add_atom(c1); eth.add_atom(o2);
        eth.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        eth.add_bond(11, 1, 2, 1, BondDisplay::Solid);

        let prep = prepare_molecule_for_quantum(&eth, true, true);
        assert_eq!(prep.added_h_count, 6);
        assert_eq!(prep.molecule.atoms.len(), 9);
        assert_eq!(prep.parent_heavy_atom.len(), 9);
        assert_eq!(CcVal::compute(&prep.molecule).implicit_h_counts.iter().sum::<u8>(), 0);

        // Verify all X-H bond lengths are chemically reasonable (0.85..1.25 Å)
        for b in &prep.molecule.bonds {
            let za = prep.molecule.atoms[b.begin_atom].z;
            let zb = prep.molecule.atoms[b.end_atom].z;
            if za == 1 || zb == 1 {
                let d = measure_distance_angstrom(&prep.molecule, b.begin_atom, b.end_atom).unwrap();
                assert!((0.85..=1.25).contains(&d), "X-H distance out of range: {d:.3} Å");
            }
        }
    }
}
