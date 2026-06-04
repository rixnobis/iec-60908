'use strict'

/* Centralized monkey-patches for the jqr-gf Galois field package.
   Require THIS module instead of 'jqr-gf' directly, so that every consumer
   (the rs.js encoders AND the read-bits.js decoder) shares the same corrected
   field. jqr-gf exports a class with static methods, and require() caches it as
   a singleton, so patching here mutates the one field object the whole program
   uses. Upstream fixes should be PR'd, but this lets us move now.

   Two patches:

   1. sub: jqr-gf never defined subtraction. In GF(2^n), subtraction is identical
      to addition (both are XOR), so sub == add. jqr-poly's division calls sub,
      so without this, any polynomial division throws.

   2. div: jqr-gf's GaloisField.div computes
          expTable[(logTable[a] - logTable[b]) % 255]
      but JavaScript's % is remainder, not modulo: for a negative left operand it
      returns a negative result. Whenever log(a) < log(b) the index goes negative
      and expTable[negative] is undefined. That is ~half of all operand pairs
      (32385 of 65025 nonzero pairs fail the mul(div(a,b),b) == a identity).
      Encoding survives because it only ever multiplies (logs are summed, never
      negative); decoding does not, because Forney's algorithm for error
      magnitudes is division-heavy. Rather than re-derive the field tables here,
      define div via the already-correct mul and inv: a / b == a * b^-1. This
      cannot reintroduce the sign bug because there is no subtraction left to get
      wrong. */
const gf = require('jqr-gf')

if (!gf.sub) {
  gf.sub = gf.add
}

gf.div = function (a, b) {
  if (b === 0) {
    throw RangeError('divide() cannot divide by zero')
  }
  if (a === 0) {
    return 0
  }
  return gf.mul(a, gf.inv(b))
}

module.exports = gf
