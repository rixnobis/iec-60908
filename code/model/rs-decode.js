'use strict'

/* Errors-and-erasures Reed-Solomon decoder for the CD codes (generator roots
   alpha^0 .. alpha^(2t-1), 2t = 4 parity symbols, so t = 2 correctable errors).

   Convention (verified against rs.js's encoders and read-bits.js's syndrome
   loop): a codeword is an array of n symbols where index 0 is the HIGHEST
   degree term, so the symbol at index i sits at degree (n - 1 - i). Syndromes
   are S_j = R(alpha^j) for j in [0, 2t). An errata at index i has locator value
   X = alpha^(n-1-i), and S_j = sum over errata of e_i * X_i^j.

   The locator search (Berlekamp-Massey + Chien) finds the unknown error
   positions. Magnitudes are then recovered by solving the linear system
   S_j = sum_k e_k * X_k^j directly (a small Vandermonde solve over GF(256)).
   That sidesteps every sign/first-root convention pitfall in the closed-form
   Forney formula: with the locations known, the syndromes ARE a linear system
   in the magnitudes, and the X_k are distinct so it is always solvable. The
   decoder revalidates by recomputing syndromes on the corrected word and
   refuses (returns null) rather than emit a miscorrection. */

const gf = require('./gf-patched')
const ALPHA = 2

function gfpow (base, e) {
  let r = 1
  for (let i = 0; i < e; i++) r = gf.mul(r, base)
  return r
}

// Evaluate a descending-coefficient polynomial (coeffs[0] = highest degree).
function evalDesc (coeffs, x) {
  let acc = 0
  for (let i = 0; i < coeffs.length; i++) acc = gf.add(gf.mul(acc, x), coeffs[i])
  return acc
}

// Evaluate an ascending-coefficient polynomial (p[0] = constant term).
function evalAsc (p, x) {
  let acc = 0
  let xp = 1
  for (let i = 0; i < p.length; i++) {
    acc = gf.add(acc, gf.mul(p[i], xp))
    xp = gf.mul(xp, x)
  }
  return acc
}

// Add two ascending-coefficient polynomials.
function polyAddAsc (a, b) {
  const out = new Array(Math.max(a.length, b.length)).fill(0)
  for (let i = 0; i < a.length; i++) out[i] = gf.add(out[i], a[i])
  for (let i = 0; i < b.length; i++) out[i] = gf.add(out[i], b[i])
  return out
}

// Solve the linear system sum_k m_k * X_k^j = S_j (j = 0..v-1) for the v
// magnitudes m_k, by Gaussian elimination over GF(256). Returns null if the
// system is singular (should not happen for distinct X_k).
function solveMagnitudes (Xs, S) {
  const v = Xs.length
  // Build the v x (v+1) augmented matrix: A[j][k] = X_k^j, rhs = S_j.
  const A = []
  for (let j = 0; j < v; j++) {
    const row = []
    for (let k = 0; k < v; k++) row.push(gfpow(Xs[k], j))
    row.push(S[j])
    A.push(row)
  }
  for (let col = 0; col < v; col++) {
    let piv = -1
    for (let r = col; r < v; r++) if (A[r][col] !== 0) { piv = r; break }
    if (piv < 0) return null
    if (piv !== col) { const t = A[piv]; A[piv] = A[col]; A[col] = t }
    const inv = gf.inv(A[col][col])
    for (let c = col; c <= v; c++) A[col][c] = gf.mul(A[col][c], inv)
    for (let r = 0; r < v; r++) {
      if (r === col || A[r][col] === 0) continue
      const f = A[r][col]
      for (let c = col; c <= v; c++) A[r][c] = gf.add(A[r][c], gf.mul(f, A[col][c]))
    }
  }
  return A.map(row => row[v])
}

/* Decode in place-safe fashion. received: array of n symbols. erasureIdx:
   indices known to be erased (from the EFM layer). Returns
   { corrected, errataIdx } on success, or null if uncorrectable. */
function decode (received, erasureIdx = [], twoT = 4, maxErrors = Infinity) {
  const n = received.length

  const S = []
  let allZero = true
  for (let j = 0; j < twoT; j++) {
    const s = evalDesc(received, gfpow(ALPHA, j))
    S.push(s)
    if (s !== 0) allZero = false
  }
  if (allZero && erasureIdx.length === 0) return { corrected: received.slice(), errataIdx: [] }

  // Erasure locator Gamma(x) = product (1 + X x), ascending coefficients.
  const erasureX = erasureIdx.map(i => gfpow(ALPHA, n - 1 - i))
  let gamma = [1]
  for (const X of erasureX) {
    const ng = new Array(gamma.length + 1).fill(0)
    for (let i = 0; i < gamma.length; i++) {
      ng[i] = gf.add(ng[i], gamma[i])
      ng[i + 1] = gf.add(ng[i + 1], gf.mul(gamma[i], X))
    }
    gamma = ng
  }

  // Forney-modified syndromes T = S * Gamma mod x^(2t).
  const T = new Array(twoT).fill(0)
  for (let i = 0; i < twoT; i++) {
    let acc = 0
    for (let k = 0; k < gamma.length && k <= i; k++) acc = gf.add(acc, gf.mul(S[i - k], gamma[k]))
    T[i] = acc
  }

  // Berlekamp-Massey on the modified syndromes -> error locator (ascending).
  // The error locator is found from the HIGH modified syndromes T[f..2t-1]
  // (offset by the erasure count f); re-index them from 0 as U. With f = 0 this
  // is just T; with no errors BM does not run. The offset only matters in the
  // mixed errors-and-erasures case.
  const f = erasureIdx.length
  const U = T.slice(f)
  const nIter = U.length // == twoT - f
  const maxErr = Math.floor(nIter / 2)
  let Lambda = [1]
  let B = [1]
  let L = 0
  let m = 1
  let b = 1
  for (let i = 0; i < nIter; i++) {
    let delta = U[i]
    for (let j = 1; j <= L; j++) delta = gf.add(delta, gf.mul(Lambda[j] || 0, U[i - j]))
    if (delta === 0) {
      m++
    } else if (2 * L <= i) {
      const Tmp = Lambda.slice()
      const coef = gf.div(delta, b)
      const shifted = new Array(m).fill(0).concat(B.map(v => gf.mul(coef, v)))
      Lambda = polyAddAsc(Lambda, shifted)
      L = i + 1 - L
      B = Tmp
      b = delta
      m = 1
    } else {
      const coef = gf.div(delta, b)
      const shifted = new Array(m).fill(0).concat(B.map(v => gf.mul(coef, v)))
      Lambda = polyAddAsc(Lambda, shifted)
      m++
    }
  }

  // Chien search: error positions are indices i (not already erased) where
  // Lambda(X^-1) == 0, with X = alpha^(n-1-i).
  const errIdx = []
  for (let i = 0; i < n; i++) {
    if (erasureIdx.includes(i)) continue
    const X = gfpow(ALPHA, n - 1 - i)
    if (evalAsc(Lambda, gf.inv(X)) === 0) errIdx.push(i)
  }
  if (errIdx.length > maxErr) return null
  // Strategy cap: a caller (e.g. C1 in a cooperative C1/C2 decoder) may choose
  // to correct only a few errors and flag the rest as erasures for the other
  // line, since low-error decoding rarely miscorrects. Refuse here if asked.
  if (errIdx.length > maxErrors) return null

  const errata = erasureIdx.concat(errIdx)
  if (errata.length === 0) return null
  const Xs = errata.map(i => gfpow(ALPHA, n - 1 - i))
  const mags = solveMagnitudes(Xs, S.slice(0, errata.length))
  if (!mags) return null

  const corrected = received.slice()
  for (let k = 0; k < errata.length; k++) {
    corrected[errata[k]] = gf.add(corrected[errata[k]], mags[k])
  }
  for (let j = 0; j < twoT; j++) {
    if (evalDesc(corrected, gfpow(ALPHA, j)) !== 0) return null
  }
  return { corrected, errataIdx: errata.slice().sort((a, c) => a - c) }
}

module.exports = { decode }
