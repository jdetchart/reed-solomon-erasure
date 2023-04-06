use ark_bls12_381::fr;
use ark_ff::{BigInteger, PrimeField};
use std::ops::{Add, Div, Mul, MulAssign, Sub};

#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct Field;

impl crate::Field for Field {
    //todo: the order is much bigger.
    const ORDER: usize = usize::MAX;
    type Elem = fr::Fr;

    fn add(a: Self::Elem, b: Self::Elem) -> Self::Elem {
        a.add(b)
    }

    fn sub(a: Self::Elem, b: Self::Elem) -> Self::Elem {
        a.sub(b)
    }

    fn mul(a: Self::Elem, b: Self::Elem) -> Self::Elem {
        a.mul(b)
    }

    fn div(a: Self::Elem, b: Self::Elem) -> Self::Elem {
        a.div(b)
    }

    fn exp(a: Self::Elem, n: usize) -> Self::Elem {
        if n == 0 {
            return fr::Fr::from(1);
        }
        if a == fr::Fr::from(0) {
            return fr::Fr::from(0);
        }
        let mut r = a;
        for _ in 1..n {
            r.mul_assign(a);
        }
        r
    }

    fn zero() -> Self::Elem {
        fr::Fr::from(0)
    }

    fn one() -> Self::Elem {
        fr::Fr::from(1)
    }

    fn nth_internal(n: usize) -> Self::Elem {
        fr::Fr::from(n as u64)
    }

    // write each Fr as 32B of data in little endian
    fn serialize(input: &[Self::Elem]) -> Vec<u8> {
        input
            .iter()
            //.flat_map(|&u| u.into_bigint().to_bytes_le())
            .flat_map(|&u| {
                //let nb_bytes = (u.into_bigint().num_bits()+7)/8;
                let mut v = u.into_bigint().to_bytes_le();
                v
            })
            .collect()
    }

    // read data 32B/32B as data is elements of Fr serialized
    fn deserialize(input: Vec<u8>) -> Vec<Self::Elem> {
        let mut output = Vec::new();

        let chunks = input.chunks(((fr::Fr::MODULUS_BIT_SIZE+7) as usize) / 8);
        for chunk in chunks {
            //output.push(fr::Fr::from_random_bytes(chunk).unwrap());
            output.push(fr::Fr::from_le_bytes_mod_order(chunk));
        }

        output
    }

    // read data 31B/31B
    fn from_data(input: Vec<u8>) -> Vec<Self::Elem> {
        let mut output = Vec::new();

        let chunks = input.chunks((fr::Fr::MODULUS_BIT_SIZE as usize) / 8);
        for chunk in chunks {
            //output.push(fr::Fr::from_random_bytes(chunk).unwrap());
            output.push(fr::Fr::from_le_bytes_mod_order(chunk));
        }

        output
    }

    // convert a slice of Fr generated as 31B of data each into elements of 31B ( should be the opposite of from_data)
    fn into_data(input: &[Self::Elem]) -> Vec<u8> {
        input
            .iter()
            //.flat_map(|&u| u.into_bigint().to_bytes_le())
            .flat_map(|&u| {
                //let nb_bytes = (u.into_bigint().num_bits()+7)/8;
                let mut v = u.into_bigint().to_bytes_le();
                v.truncate((fr::Fr::MODULUS_BIT_SIZE as usize) / 8 as usize);
                v
            })
            .collect()
    }
}

pub type ReedSolomon = crate::ReedSolomon<Field>;
pub type ReedSolomonNS = crate::ReedSolomonNonSystematic<Field>;

pub type ShardByShard<'a> = crate::ShardByShard<'a, Field>;
