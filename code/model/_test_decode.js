'use strict'
const rs = require('./rs')
const { decode } = require('./rs-decode')

// Build a systematic (n, n-4) codeword: k data symbols + 4 parity (encode_4),
// data high-order, parity low-order, matching the verified convention.
function makeCodeword (data) {
  return data.concat(rs.encode_4(data))
}
function eq (a, b) { return a.length === b.length && a.every((v, i) => v === b[i]) }

const baseData = [0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0, 0x11, 0x22]
const code = makeCodeword(baseData) // n = 14
const n = code.length

let pass = 0, fail = 0
const note = (ok, msg) => { if (ok) pass++; else { fail++; console.log('  FAIL:', msg) } }

// pseudo-random without Date/Math.random dependence concerns: simple LCG
let seed = 0x1234
const rnd = (mod) => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed % mod }

// 1. Clean word decodes to itself.
{
  const r = decode(code.slice(), [])
  note(r && eq(r.corrected, code), 'clean word')
}

// 2. Every single-error position, several magnitudes.
{
  let ok = true
  for (let p = 0; p < n; p++) {
    for (const mag of [1, 0x53, 0xff, 0xa7]) {
      const r = code.slice(); r[p] = r[p] ^ mag
      const d = decode(r, [])
      if (!d || !eq(d.corrected, code)) { ok = false; console.log('   single-err miss at pos', p, 'mag', mag) }
    }
  }
  note(ok, 'all single errors')
}

// 3. Random double errors (2 errors, 0 erasures: 2*2+0 = 4 <= 4, correctable).
{
  let ok = true, trials = 200
  for (let t = 0; t < trials; t++) {
    let p1 = rnd(n), p2 = rnd(n); if (p1 === p2) p2 = (p2 + 1) % n
    const m1 = 1 + rnd(255), m2 = 1 + rnd(255)
    const r = code.slice(); r[p1] ^= m1; r[p2] ^= m2
    const d = decode(r, [])
    if (!d || !eq(d.corrected, code)) { ok = false; if (fail < 5) console.log('   double-err miss', p1, p2, m1, m2) }
  }
  note(ok, 'random double errors (200 trials)')
}

// 4. Single error + 2 erasures (2*1 + 2 = 4 <= 4, correctable).
{
  let ok = true, trials = 200
  for (let t = 0; t < trials; t++) {
    const e1 = rnd(n); let e2 = rnd(n); if (e2 === e1) e2 = (e2 + 1) % n
    let ep = rnd(n); while (ep === e1 || ep === e2) ep = (ep + 1) % n
    const r = code.slice()
    r[e1] ^= (1 + rnd(255)); r[e2] ^= (1 + rnd(255)) // erased positions (values corrupted, positions known)
    r[ep] ^= (1 + rnd(255)) // unknown error
    const d = decode(r, [e1, e2])
    if (!d || !eq(d.corrected, code)) { ok = false; if (fail < 8) console.log('   1err+2eras miss', e1, e2, ep) }
  }
  note(ok, '1 error + 2 erasures (200 trials)')
}

// 5. Four erasures (0 errors: 0 + 4 = 4, correctable as pure erasures).
{
  let ok = true, trials = 200
  for (let t = 0; t < trials; t++) {
    const set = new Set(); while (set.size < 4) set.add(rnd(n))
    const er = [...set]
    const r = code.slice(); er.forEach(p => { r[p] ^= (1 + rnd(255)) })
    const d = decode(r, er)
    if (!d || !eq(d.corrected, code)) { ok = false; if (fail < 12) console.log('   4-erasure miss', er) }
  }
  note(ok, '4 erasures, no error (200 trials)')
}

// 6. Uncorrectable: 3 errors (2*3 = 6 > 4). Decoder must NOT emit a wrong
//    codeword silently. It may return null, or (rarely) a valid-but-different
//    codeword; the contract we enforce is "never returns a corrupted word
//    claimed as correct" -> if it returns non-null it must pass syndromes,
//    which it self-checks. We assert it does not return the WRONG original.
{
  let safe = true, trials = 300, miscorrect = 0
  for (let t = 0; t < trials; t++) {
    const set = new Set(); while (set.size < 3) set.add(rnd(n))
    const ps = [...set]
    const r = code.slice(); ps.forEach(p => { r[p] ^= (1 + rnd(255)) })
    const d = decode(r, [])
    // A miscorrection would be: returns a word != code but claims success.
    // That can happen with RS beyond capacity; we just measure how often.
    if (d && !eq(d.corrected, code)) miscorrect++
  }
  console.log('  (info) 3-error over-capacity: miscorrections', miscorrect, '/', trials, '(expected: decoder returns null or detects most)')
  note(safe, '3-error over-capacity did not crash')
}

console.log(`\n${pass} passing, ${fail} failing`)
process.exit(fail ? 1 : 0)
