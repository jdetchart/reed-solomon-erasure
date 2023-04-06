use crate::galois_prime::{Field as PrimeF, ReedSolomon, ReedSolomonNS};
use crate::matrix::Matrix;
use crate::tests::{option_shards_to_shards, shards_to_option_shards};
use crate::Field;
use ark_bls12_381::{fr, Fr};
use ark_ff::{BigInteger, One, PrimeField, Zero};
use rand::prelude::SliceRandom;
use std::ops::{Div, Mul};
use std::str::FromStr;

fn print_shards(shards: &Vec<Vec<Fr>>) {
    for shard in shards {
        print!("shard: ");
        for e in shard {
            print!("{:?}", e.to_string());
        }
        println!(" ");
    }
}

fn display_matrix(mat: &Matrix<PrimeF>) {
    for i in 0..mat.row_count() {
        for j in 0..mat.col_count() {
            print!("0{},", mat.get(i, j));
        }
        println!("");
    }
}

#[test]
fn test_fr() {
    let e = fr::Fr::from(1234);
    let one = fr::Fr::from(1);
    let i = one.div(e);

    let ii = one.div(i);
    println!("{:?}", e);
    println!("{:?}", i.clone());
    println!("{:?}", one);
    println!("{:?}", ii);

    assert_eq!(ii, e);

    let o = e.mul(i);
    assert_eq!(o, one);

    let a = fr::Fr::from(123);
    let b = fr::Fr::from(234);
    let c = fr::Fr::from(345);

    let d = a.mul(b).mul(c); // d= a*b*c
    let e = b.mul(c); // e=b*c

    let f = d.div(e); // f=d/e

    let aaa = f.0;
    println!("{:?}", aaa.to_bytes_le());

    assert_eq!(f, a); // f == a ?
}

#[test]
fn test_vec_fr_small() {
    let mut bytes = Vec::new();

    bytes.push(0);
    bytes.push(0);
    bytes.push(1);
    bytes.push(0);
    let elems = crate::galois_prime::Field::from_data(bytes.clone());
    let mut bytes2 = crate::galois_prime::Field::into_data(elems.as_slice());
    bytes2.truncate(bytes.len());
    assert_eq!(bytes, bytes2);
}

#[test]
fn test_vec_fr_big() {
    let mut bytes = Vec::new();

    bytes.push(1);

    let elems = crate::galois_prime::Field::from_data(bytes.clone());
    let mut bytes2 = crate::galois_prime::Field::into_data(elems.as_slice());
    bytes2.truncate(bytes.len());
    assert_eq!(bytes, bytes2);

    bytes.push(4);
    bytes.push(120);
    bytes.push(45);

    let elems = crate::galois_prime::Field::from_data(bytes.clone());
    let mut bytes2 = crate::galois_prime::Field::into_data(elems.as_slice());
    bytes2.truncate(bytes.len());
    assert_eq!(bytes, bytes2);

    for i in 0..60 {
        bytes.push(i);
    }
    let elems = crate::galois_prime::Field::from_data(bytes.clone());
    println!("{}", elems.len());
    let mut bytes2 = crate::galois_prime::Field::into_data(elems.as_slice());
    bytes2.truncate(bytes.len());
    assert_eq!(bytes, bytes2);
}

#[test]
fn test() {
    let (k, r) = (4, 2);
    let rs = ReedSolomon::new(k, r).expect("cannot create RS_381");

    let mut shards = Vec::with_capacity(k + r);
    for i in 0..k + r {
        let mut s = Vec::with_capacity(1);
        let e = Fr::from(i as u32 + 1);
        s.push(e);
        shards.push(s);
    }

    print_shards(&shards);

    rs.encode(&mut shards).unwrap();
    assert!(rs.verify(&shards).unwrap());

    let master_copy = shards.clone();

    let mut shards = shards_to_option_shards(&shards);

    rs.reconstruct(&mut shards).unwrap();
    {
        let shards = option_shards_to_shards(&shards);
        assert!(rs.verify(&shards).unwrap());
        assert_eq!(&shards, &master_copy);
    }
    shards[0] = None;
    shards[2] = None;
    //shards[4] = None;
    rs.reconstruct(&mut shards).unwrap();
    {
        let shards = option_shards_to_shards(&shards);
        print_shards(&shards);
        assert!(rs.verify(&shards).unwrap());
        assert_eq!(&shards, &master_copy);
    }
}

#[test]
fn test_non_systematic() {
    let (k, n) = (3, 5);
    let rs = ReedSolomonNS::vandermonde(k, n).unwrap();

    let mut shards = Vec::with_capacity(n);
    for i in 0..n {
        let mut s = Vec::with_capacity(1);
        let e = Fr::from(i as u32 + 1);
        s.push(e);
        shards.push(s);
    }

    let master_copy = shards.clone();

    rs.encode(&mut shards).unwrap();

    let mut shards = shards_to_option_shards(&shards);

    shards[1] = None;
    shards[2] = None;
    rs.reconstruct(&mut shards).unwrap();
    let shards = option_shards_to_shards(&shards);
    for i in 0..k {
        assert_eq!(shards.get(i).unwrap(), master_copy.get(i).unwrap());
    }
}

#[test]
fn test_convert() {
    let mut rng = rand::thread_rng();

    let mut s = Vec::with_capacity(100);
    for _ in 0..100 {
        let mut nums: Vec<u8> = (0..31).collect();
        nums.shuffle(&mut rng);
        let e = Fr::from_le_bytes_mod_order(&nums);
        s.push(e);
    }
    let v = crate::galois_prime::Field::serialize(&s);

    let elts = crate::galois_prime::Field::deserialize(v.clone());
    assert_eq!(elts, s);

    let vv = crate::galois_prime::Field::serialize(&elts);
    assert_eq!(vv, v);

    let e = Fr::from(0);
    let f = Fr::zero();
    assert_eq!(e, f);

    let e = Fr::from(1);
    let f = Fr::one();
    assert_eq!(e, f);
}

#[test]
fn test_non_systematic_big() {
    let (k, n) = (100, 400);
    let rs = ReedSolomonNS::vandermonde(k, n).unwrap();

    let mut rng = rand::thread_rng();
    let mut shards = Vec::with_capacity(n);
    for _ in 0..n {
        let mut s = Vec::with_capacity(1);
        let mut nums: Vec<u8> = (0..31).collect();
        nums.shuffle(&mut rng);
        let e = Fr::from_le_bytes_mod_order(&nums);
        println!("{}", e);
        s.push(e);
        shards.push(s);
    }

    let master_copy = shards.clone();

    rs.encode(&mut shards).unwrap();

    let mut shards = shards_to_option_shards(&shards);

    for i in 0..k {
        shards[i] = None;
    }

    rs.reconstruct(&mut shards).unwrap();
    let shards = option_shards_to_shards(&shards);
    for i in 0..k {
        assert_eq!(shards.get(i).unwrap(), master_copy.get(i).unwrap());
    }
}

#[test]
fn test_kzg_interpret() {
    let (k, n) = (9, 16);

    let result = [
        "1",
        "450546001518488004043740862689444221536008393703282834321009581329618042925",
        "10468026038480548076302179832125789037837043786985129940971135420505056612351",
        "19623306155733127681910331647029450587080131057690957914338975907063433086605",
        "5585031919584593081523271106717735964065031941283735935824073952426852060166",
        "19734755435881330559657342592047686213784598574104623986489109888169549046833",
        "42158234230692591021509553361872595711911063218759632923430255714073653963722",
        "1835050255571180876438380995291916676587020443102309563329074505681459593928",
        "13999541685248685924868409598719329183003259454390178360664146922932411377078",
        "18438132830253093969482988563679170302433480715024021486038817893967554183262",
        "18253374773506337245296614778354950310564194920130509571260896177300656180305",
        "44391972966686073740002486674521882893263271329944375685904200539462381966805",
        "43153552156434253219658856461261505522668921394300274764395429315445910276604",
        "11970629992521811888386655719275230342789277541400939550117556550609953372367",
        "50285864025212780532643610729242942504752813429970145631413803854943979290257",
        "15554607430845508155353431713719595382777120897700191346616646132050837540676",
    ];

    let data = vec![
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 255,
    ];

    let d = PrimeF::from_data(data.clone());
    let mut shards = Vec::with_capacity(n);
    for chunk in d.chunks_exact(1) {
        println!(" input {}", chunk[0]);
        shards.push(chunk.to_vec());
    }
    for _ in 0..n - k {
        shards.push(vec![PrimeF::zero(); 1]);
    }
    let rs = ReedSolomonNS::vandermonde(k, n).unwrap();

    rs.encode(&mut shards).unwrap();

    for (i, s) in shards.iter().enumerate() {
        let v1 = s.get(0).unwrap();
        println!("{}", s.get(0).unwrap());
        let v2 = Fr::from_str(result[i]).unwrap();
        assert!(v1.eq(&v2));
    }
}

/*
#[test]
fn test_kzg() {
    let coefficients = vec![1, 2, 3, 1, 1, 17, 32]
        .into_iter()
        .map(Fr::from)
        .collect::<Vec<_>>();

    let degree = 7;

    //setup
    let mut rng = thread_rng();

    let mut secret = [0u8; 32];
    rng.fill_bytes(&mut secret);

    let s = Fr::from_be_bytes_mod_order(&secret);

    let mut points_in_g1 = vec![];

    let g1 = P1::generator();
    for i in 0..=degree {
        let i_as_bigint = BigUint::from_slice(&[i as u32]);
        let s_i_as_bigint = s.pow(&i_as_bigint);

        let mut s_i_bytes = vec![0u8; 32];
        let raw_bytes = s_i_as_bigint.to_bytes_be();
        s_i_bytes[32 - raw_bytes.len()..].copy_from_slice(&raw_bytes);
        let s_i_scalar = Scalar::from_fr_bytes(&s_i_bytes);

        let result = s_i_scalar * g1;
        points_in_g1.push(result);
    }

    // NOTE: `secret` in Fr via prior `assert`.
    let scalar = Scalar::from_fr_bytes(secret);
    let result_in_g2 = scalar * P2::generator();

    let setup = (points_in_g1,result_in_g2);
}*/
