use crate::arith::{U256, U512};
use crate::fields::{const_fq, FieldElement, Fq};
use bytemuck::{AnyBitPattern, NoUninit};
use core::cmp::Ordering;
use core::ops::{Add, Div, Mul, Neg, Sub};
use rand::Rng;

use super::Sqrt;

cfg_if::cfg_if! {
    if #[cfg(target_os = "zkvm")] {
        use sp1_lib::io::{hint_slice, read_vec};
        use core::convert::TryInto;
        use bytemuck::{cast, cast_ref, cast_mut};
    }
}

#[inline]
fn fq_non_residue() -> Fq {
    // (q - 1) is a quadratic nonresidue in Fq
    // 21888242871839275222246405745257275088696311157297823662689037894645226208582
    const_fq([
        0x3c208c16d87cfd46,
        0x97816a916871ca8d,
        0xb85045b68181585d,
        0x30644e72e131a029,
    ])
}

#[inline]
pub const fn fq2_nonresidue() -> Fq2 {
    Fq2::new(
        Fq::from_raw_unchecked(U256::from_raw_unchecked([9, 0])),
        Fq::from_raw_unchecked(U256::from_raw_unchecked([1, 0])),
    )
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, NoUninit, AnyBitPattern)]
#[repr(C)]
pub struct Fq2 {
    c0: Fq,
    c1: Fq,
}

impl Fq2 {
    pub const fn new(c0: Fq, c1: Fq) -> Self {
        Fq2 { c0, c1 }
    }

    pub fn scale(&self, by: Fq) -> Self {
        Fq2 {
            c0: self.c0 * by,
            c1: self.c1 * by,
        }
    }

    #[inline]
    pub fn mul_by_nonresidue_inp(&mut self) {
        #[cfg(target_os = "zkvm")]
        {
            self.mul_inp(&fq2_nonresidue());
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            *self = *self * fq2_nonresidue();
        }
    }

    #[inline]
    pub fn mul_by_nonresidue(&self) -> Self {
        #[cfg(target_os = "zkvm")]
        {
            let mut res = *self;
            res.mul_inp(&fq2_nonresidue());
            res
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            *self * fq2_nonresidue()
        }
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        if power % 2 == 0 {
            *self
        } else {
            Fq2 {
                c0: self.c0,
                c1: self.c1 * fq_non_residue(),
            }
        }
    }

    pub fn real(&self) -> &Fq {
        &self.c0
    }

    pub fn imaginary(&self) -> &Fq {
        &self.c1
    }

    fn cpu_add(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0.cpu_add(other.c0),
            c1: self.c1.cpu_add(other.c1),
        }
    }

    fn cpu_mul(self, other: Fq2) -> Fq2 {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Karatsuba)

        let aa = self.c0.cpu_mul(other.c0);
        let bb = self.c1.cpu_mul(other.c1);

        Fq2 {
            c0: bb.cpu_mul(fq_non_residue()).cpu_add(aa),
            c1: (self.c0.cpu_add(self.c1))
                .cpu_mul(other.c0.cpu_add(other.c1))
                .cpu_sub(aa)
                .cpu_sub(bb),
        }
    }

    fn cpu_sub(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0.cpu_sub(other.c0),
            c1: self.c1.cpu_sub(other.c1),
        }
    }

    fn cpu_neg(self) -> Fq2 {
        Fq2 {
            c0: self.c0.cpu_neg(),
            c1: self.c1.cpu_neg(),
        }
    }

    #[inline]
    pub(crate) fn add_inp(&mut self, other: &Fq2) {
        #[cfg(target_os = "zkvm")]
        {
            let lhs = cast_mut::<Fq2, [u32; 16]>(self);
            let rhs = cast_ref::<Fq2, [u32; 16]>(&other);
            unsafe {
                sp1_lib::syscall_bn254_fp2_addmod(lhs.as_mut_ptr(), rhs.as_ptr());
            }
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            *self = self.cpu_add(*other);
        }
    }

    #[inline]
    pub(crate) fn sub_inp(&mut self, other: &Fq2) {
        #[cfg(target_os = "zkvm")]
        {
            let lhs = cast_mut::<Fq2, [u32; 16]>(self);
            let rhs = cast_ref::<Fq2, [u32; 16]>(&other);
            unsafe {
                sp1_lib::syscall_bn254_fp2_submod(lhs.as_mut_ptr(), rhs.as_ptr());
            }
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            *self = self.cpu_sub(*other);
        }
    }

    #[inline]
    pub(crate) fn mul_inp(&mut self, other: &Fq2) {
        #[cfg(target_os = "zkvm")]
        {
            let lhs = cast_mut::<Fq2, [u32; 16]>(self);
            let rhs = cast_ref::<Fq2, [u32; 16]>(&other);
            unsafe {
                sp1_lib::syscall_bn254_fp2_mulmod(lhs.as_mut_ptr(), rhs.as_ptr());
            }
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            *self = self.cpu_mul(*other);
        }
    }

    #[inline]
    pub fn square_inp(&mut self) {
        #[cfg(target_os = "zkvm")]
        {
            let lhs = cast_mut::<Fq2, [u32; 16]>(self);
            unsafe {
                sp1_lib::syscall_bn254_fp2_mulmod(lhs.as_mut_ptr(), lhs.as_ptr());
            }
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            *self = self.cpu_mul(*self);
        }
    }

    #[inline]
    pub fn double_inp(&mut self) {
        #[cfg(target_os = "zkvm")]
        {
            let lhs = cast_mut::<Fq2, [u32; 16]>(self);
            unsafe {
                sp1_lib::syscall_bn254_fp2_addmod(lhs.as_mut_ptr(), lhs.as_ptr());
            }
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            *self = self.cpu_add(*self);
        }
    }
}

impl FieldElement for Fq2 {
    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }

    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq2 {
            c0: Fq::random(rng),
            c1: Fq::random(rng),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn squared(&self) -> Self {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Complex squaring)

        let mut out = *self;
        out.square_inp();
        out
    }

    fn inverse(self) -> Option<Self> {
        // "High-Speed Software Implementation of the Optimal Ate Pairing
        // over Barreto–Naehrig Curves"; Algorithm 8

        (self
            .c0
            .cpu_mul(self.c0)
            .cpu_sub((self.c1.cpu_mul(self.c1)).cpu_mul(fq_non_residue())))
        .inverse()
        .map(|t| Fq2 {
            c0: self.c0.cpu_mul(t),
            c1: (self.c1.cpu_mul(t)).cpu_neg(),
        })
    }

    fn inverse_unconstrained(self) -> Option<Self> {
        #[cfg(target_os = "zkvm")]
        {
            // Compute the inverse using the zkvm syscall
            sp1_lib::unconstrained! {
                let mut buf = [0u8; 65];
                self.inverse().map(|inv| {
                    let bytes = cast::<[u128; 4], [u8; 64]>(inv.to_u512().0);
                    buf[0..64].copy_from_slice(&bytes);
                    buf[64] = 1;
                });
                hint_slice(&buf);
            }
            let byte_vec = read_vec();
            let bytes: [u8; 65] = byte_vec.try_into().unwrap();
            match bytes[64] {
                0 => None,
                _ => {
                    let inv = cast::<[u8; 64], Fq2>(bytes[0..64].try_into().unwrap());
                    Some(inv).filter(|inv| !self.is_zero() && self * *inv == Fq2::one())
                }
            }
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            self.inverse()
        }
    }
}

impl Mul for Fq2 {
    type Output = Fq2;

    #[allow(unused_mut)]
    fn mul(mut self, other: Fq2) -> Fq2 {
        #[cfg(target_os = "zkvm")]
        {
            self.mul_inp(&other);
            self
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            // Devegili OhEig Scott Dahab
            //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
            //     Section 3 (Karatsuba)

            let aa = self.c0 * other.c0;
            let bb = self.c1 * other.c1;

            Fq2 {
                c0: bb * fq_non_residue() + aa,
                c1: (self.c0 + self.c1) * (other.c0 + other.c1) - aa - bb,
            }
        }
    }
}

impl Div for Fq2 {
    type Output = Fq2;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs
            .inverse()
            .expect("Failed to compute the inverse of the divisor")
    }
}

impl Sub for Fq2 {
    type Output = Fq2;

    fn sub(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
        }
    }
}

impl Add for Fq2 {
    type Output = Fq2;

    fn add(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

impl Neg for Fq2 {
    type Output = Fq2;

    fn neg(self) -> Fq2 {
        Fq2 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

/// `Fq2` elements are ordered lexicographically.
impl Ord for Fq2 {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl PartialOrd for Fq2 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

lazy_static::lazy_static! {
    static ref FQ: U256 = U256::from([
        0x3c208c16d87cfd47,
        0x97816a916871ca8d,
        0xb85045b68181585d,
        0x30644e72e131a029
    ]);

    static ref FQ_MINUS3_DIV4: Fq =
        Fq::new(3.into()).expect("3 is a valid field element and static; qed").cpu_neg().cpu_mul(
        Fq::new(4.into()).expect("4 is a valid field element and static; qed").inverse()
        .expect("4 has inverse in Fq and is static; qed"));

    static ref FQ_MINUS1_DIV2: Fq =
        Fq::new(1.into()).expect("1 is a valid field element and static; qed").cpu_neg().cpu_mul(
        Fq::new(2.into()).expect("2 is a valid field element and static; qed").inverse()
            .expect("2 has inverse in Fq and is static; qed"));
}

impl Fq2 {
    pub fn i() -> Fq2 {
        Fq2::new(Fq::zero(), Fq::one())
    }

    fn cpu_pow(&self, by: U256) -> Self {
        let mut res = Self::one();

        for i in by.bits() {
            res = res.cpu_mul(res);
            if i {
                res = res.cpu_mul(*self);
            }
        }
        res
    }

    fn cpu_sqrt(&self) -> Option<Self> {
        let a1 = self.cpu_pow((*FQ_MINUS3_DIV4).into());
        let a1a = a1.cpu_mul(*self);
        let alpha = a1.cpu_mul(a1a);
        let a0 = alpha.cpu_pow(*FQ).cpu_mul(alpha);

        if a0 == Fq2::one().cpu_neg() {
            return None;
        }

        if alpha == Fq2::one().cpu_neg() {
            Some(Self::i().cpu_mul(a1a))
        } else {
            let b = (alpha.cpu_add(Fq2::one())).cpu_pow((*FQ_MINUS1_DIV2).into());
            Some(b.cpu_mul(a1a))
        }
    }

    pub fn sqrt(&self) -> Option<Self> {
        #[cfg(target_os = "zkvm")]
        {
            // Compute the square root using the zkvm syscall
            sp1_lib::unconstrained! {
                let mut buf = [0u8; 65];
                self.cpu_sqrt().map(|sqrt| {
                    let bytes = cast::<Fq2, [u8; 64]>(sqrt);
                    buf[0..64].copy_from_slice(&bytes);
                    buf[64] = 1;
                });
                hint_slice(&buf);
            }
            let byte_vec = read_vec();
            let bytes: [u8; 65] = byte_vec.try_into().unwrap();
            match bytes[64] {
                0 => None,
                _ => {
                    let sqrt = cast::<[u8; 64], Fq2>(bytes[0..64].try_into().unwrap());
                    Some(sqrt).filter(|sqrt| *sqrt * *sqrt == *self)
                }
            }
        }
        #[cfg(not(target_os = "zkvm"))]
        {
            self.cpu_sqrt()
        }
    }

    pub fn to_u512(self) -> U512 {
        let c0: U256 = (*self.real()).into();
        let c1: U256 = (*self.imaginary()).into();

        U512::new(&c1, &c0, &FQ)
    }
}

impl Sqrt for Fq2 {
    fn sqrt(&self) -> Option<Self> {
        self.sqrt()
    }
}

#[test]
fn sqrt_fq2() {
    // from zcash test_proof.cpp
    let x1 = Fq2::new(
        Fq::from_str(
            "12844195307879678418043983815760255909500142247603239203345049921980497041944",
        )
        .unwrap(),
        Fq::from_str(
            "7476417578426924565731404322659619974551724117137577781074613937423560117731",
        )
        .unwrap(),
    );

    let x2 = Fq2::new(
        Fq::from_str(
            "3345897230485723946872934576923485762803457692345760237495682347502347589474",
        )
        .unwrap(),
        Fq::from_str(
            "1234912378405347958234756902345768290345762348957605678245967234857634857676",
        )
        .unwrap(),
    );

    assert_eq!(x2.sqrt().unwrap(), x1);

    // i is sqrt(-1)
    assert_eq!(Fq2::one().neg().sqrt().unwrap(), Fq2::i(),);

    // no sqrt for (1 + 2i)
    assert!(
        Fq2::new(Fq::from_str("1").unwrap(), Fq::from_str("2").unwrap())
            .sqrt()
            .is_none()
    );
}
