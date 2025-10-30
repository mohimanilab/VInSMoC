use std::{num::ParseIntError, str::FromStr};

use num_traits::{NumCast, PrimInt, Signed};

use nom::{
    branch::alt,
    bytes::complete::{tag, take_while_m_n},
    character::complete::digit1,
    error::VerboseError,
    multi::{many1, separated_list1},
    IResult, Parser,
};

/// Parses a charge given in prefix format into a signed integer.
///
/// # Examples
///
/// ```
/// use parsers::charge::prefix_charge;
///
/// assert_eq!(prefix_charge::<i32>("+1"), Ok(("", 1)));
/// assert_eq!(prefix_charge::<i32>("-1"), Ok(("", -1)));
/// assert_eq!(prefix_charge::<i32>("10"), Ok(("", 10)));
/// assert_eq!(prefix_charge::<i32>("+0"), Ok(("", 0)));
/// assert_eq!(prefix_charge::<i32>("-01"), Ok(("", -1)));
/// assert_eq!(prefix_charge::<i32>("1+"), Ok(("+", 1)));
/// assert!(prefix_charge::<i32>("++1").is_err());
/// assert!(prefix_charge::<i32>("+").is_err());
/// assert!(prefix_charge::<i32>("-").is_err());
/// ```
pub fn prefix_charge<Int: PrimInt + Signed + FromStr<Err = ParseIntError>>(
    input: &str,
) -> IResult<&str, Int, VerboseError<&str>> {
    let (input, sign) = take_while_m_n(0, 1, |c| c == '+' || c == '-')(input)?;
    let sign = match sign {
        "+" | "" => Int::one(),
        "-" => NumCast::from(-1).unwrap(),
        _ => unreachable!(),
    };

    let (input, charge) = digit1(input)?;

    let charge = charge
        .parse::<Int>()
        .expect("failed to parse integer from string")
        .mul(sign);

    Ok((input, charge))
}

/// Parses a charge given in postfix format into a signed integer.
///
/// # Examples
///
/// ```
/// use parsers::charge::postfix_charge;
///
/// assert_eq!(postfix_charge::<i32>("1+"), Ok(("", 1)));
/// assert_eq!(postfix_charge::<i32>("1-"), Ok(("", -1)));
/// assert_eq!(postfix_charge::<i32>("10"), Ok(("", 10)));
/// assert_eq!(postfix_charge::<i32>("0+"), Ok(("", 0)));
/// assert_eq!(postfix_charge::<i32>("01-"), Ok(("", -1)));
/// assert_eq!(postfix_charge::<i32>("1--"), Ok(("-", -1)));
/// assert!(postfix_charge::<i32>("+1").is_err());
/// assert!(postfix_charge::<i32>("++1").is_err());
/// assert!(postfix_charge::<i32>("+").is_err());
/// assert!(postfix_charge::<i32>("-").is_err());
/// ```
pub fn postfix_charge<Int: PrimInt + Signed + FromStr<Err = ParseIntError>>(
    input: &str,
) -> IResult<&str, Int, VerboseError<&str>> {
    let (input, charge) = digit1(input)?;

    let (input, sign) = take_while_m_n(0, 1, |c| c == '+' || c == '-')(input)?;
    let sign = match sign {
        "+" | "" => Int::one(),
        "-" => NumCast::from(-1).unwrap(),
        _ => unreachable!(),
    };

    let charge = charge
        .parse::<Int>()
        .expect("failed to parse integer from string")
        .mul(sign);

    Ok((input, charge))
}

/// Parses a charge given in repetition format (e.g. ++, -, ++++, --).
///
/// # Examples
///
/// ```
/// use parsers::charge::repetition_charge;
///
/// assert_eq!(repetition_charge::<i32>("+"), Ok(("", 1)));
/// assert_eq!(repetition_charge::<i32>("-"), Ok(("", -1)));
/// assert_eq!(repetition_charge::<i32>("+++"), Ok(("", 3)));
/// assert_eq!(repetition_charge::<i32>("---"), Ok(("", -3)));
/// assert_eq!(repetition_charge::<i32>("-+"), Ok(("+", -1)));
/// assert_eq!(repetition_charge::<i32>("+-"), Ok(("-", 1)));
/// assert_eq!(repetition_charge::<i32>("+1"), Ok(("1", 1)));
/// assert!(repetition_charge::<i32>("1+").is_err());
/// ```
pub fn repetition_charge<Int: PrimInt + Signed + FromStr<Err = ParseIntError>>(
    input: &str,
) -> IResult<&str, Int, VerboseError<&str>> {
    let mut parse_repeat = alt((many1(tag("+")), many1(tag("-"))));
    let (input, charges) = parse_repeat(input)?;

    // If many1 succeeded, we have at least one charge
    debug_assert!(!charges.is_empty());

    let sign = match charges[0] {
        "+" => Int::one(),
        "-" => NumCast::from(-1).unwrap(),
        _ => unreachable!(),
    };

    let charge: Int = NumCast::from(charges.len()).unwrap();
    let charge = charge.mul(sign);

    Ok((input, charge))
}

/// Parses a charge given in either prefix or postfix format into a signed integer.
///
/// Attempts to parse using postfix format first since "1+" could be valid in prefix notation
/// if we do not consume the "+" character.
///
/// # Examples
///
/// ```
/// use parsers::charge::anyfix_charge;
///
/// assert_eq!(anyfix_charge::<i32>("+1"), Ok(("", 1)));
/// assert_eq!(anyfix_charge::<i32>("-1"), Ok(("", -1)));
/// assert_eq!(anyfix_charge::<i32>("10"), Ok(("", 10)));
/// assert_eq!(anyfix_charge::<i32>("+0"), Ok(("", 0)));
/// assert_eq!(anyfix_charge::<i32>("-01"), Ok(("", -1)));
/// assert_eq!(anyfix_charge::<i32>("1+"), Ok(("", 1)));
/// assert_eq!(anyfix_charge::<i32>("1-"), Ok(("", -1)));
/// assert_eq!(anyfix_charge::<i32>("10"), Ok(("", 10)));
/// assert_eq!(anyfix_charge::<i32>("0+"), Ok(("", 0)));
/// assert_eq!(anyfix_charge::<i32>("01-"), Ok(("", -1)));
/// assert_eq!(anyfix_charge::<i32>("1--"), Ok(("-", -1)));
/// assert!(anyfix_charge::<i32>("++1").is_err());
/// assert!(anyfix_charge::<i32>("+").is_err());
/// assert!(anyfix_charge::<i32>("-").is_err());
/// ```
pub fn anyfix_charge<Int: PrimInt + Signed + FromStr<Err = ParseIntError>>(
    input: &str,
) -> IResult<&str, Int, VerboseError<&str>> {
    alt((postfix_charge, prefix_charge))(input)
}

/// Parses a charge attempting to use any of the above formats.
///
/// Returns 0 (neutral) and consumes nothing if everything else failed.
///
/// # Examples
///
/// ```
/// use parsers::charge::charge;
///
/// assert_eq!(charge::<i32>("+"), Ok(("", 1)));
/// assert_eq!(charge::<i32>("-"), Ok(("", -1)));
/// assert_eq!(charge::<i32>("+++"), Ok(("", 3)));
/// assert_eq!(charge::<i32>("---"), Ok(("", -3)));
/// assert_eq!(charge::<i32>("-+"), Ok(("+", -1)));
/// assert_eq!(charge::<i32>("+-"), Ok(("-", 1)));
/// assert_eq!(charge::<i32>("+1"), Ok(("", 1)));
/// assert_eq!(charge::<i32>("1+"), Ok(("", 1)));
/// assert_eq!(charge::<i32>("10"), Ok(("", 10)));
/// assert_eq!(charge::<i32>("+0"), Ok(("", 0)));
/// assert_eq!(charge::<i32>("-01"), Ok(("", -1)));
/// assert_eq!(charge::<i32>("1-"), Ok(("", -1)));
/// assert_eq!(charge::<i32>("0+"), Ok(("", 0)));
/// assert_eq!(charge::<i32>("01-"), Ok(("", -1)));
/// assert_eq!(charge::<i32>("1--"), Ok(("-", -1)));
/// assert_eq!(charge::<i32>(""), Ok(("", 0)));
/// assert_eq!(charge::<i32>("cannot parse"), Ok(("cannot parse", 0)));
/// ```
pub fn charge<Int: PrimInt + Signed + FromStr<Err = ParseIntError>>(
    input: &str,
) -> IResult<&str, Int, VerboseError<&str>> {
    let res = alt((anyfix_charge, repetition_charge))(input);

    match res {
        Ok(parsed) => Ok(parsed),
        Err(_) => Ok((input, Int::zero())),
    }
}

/// Constructs a parser to consume a delimited list of charges given a charge parser and a delimiter parser.
///
/// Consumes a list of the form `x1*x2*x3...xn`, where <*> is any string matched by the provided delimiter parser,
/// and where each `xi` can be parsed as a charge. The charges are then returned in a vector of generic type integers.
///
/// # Examples
///
/// ```
/// use parsers::charge::prefix_charge;
/// use parsers::charge::charge_list;
/// use nom::bytes::complete::tag;
///
/// let mut parser = charge_list(prefix_charge, tag(", "));
/// assert_eq!(parser("+1, +2"), Ok(("", vec![1, 2])));
/// assert_eq!(parser("+1, -2"), Ok(("", vec![1, -2])));
/// assert_eq!(parser("+1, +2"), Ok(("", vec![1, 2])));
///
/// // Note that the delimiter must match the provided string exactly to be consumed
/// assert_eq!(parser("+1,+2"), Ok((",+2", vec![1])));
///
/// assert_eq!(parser("+1"), Ok(("", vec![1])));
/// assert_eq!(parser("1, 2"), Ok(("", vec![1, 2])));
/// assert_eq!(parser("1 2"), Ok((" 2", vec![1])));
///
/// let mut parser = charge_list(prefix_charge, tag(" "));
/// assert_eq!(parser("1 2"), Ok(("", vec![1, 2])));
/// assert_eq!(parser("1 2hello"), Ok(("hello", vec![1, 2])));
/// assert_eq!(parser("1 2 hello"), Ok((" hello", vec![1, 2])));
/// ```
pub fn charge_list<'a, Int: PrimInt + Signed + FromStr<Err = ParseIntError>>(
    charge_parser: impl Parser<&'a str, Int, VerboseError<&'a str>>,
    delimiter: impl Parser<&'a str, &'a str, VerboseError<&'a str>>,
) -> impl FnMut(&'a str) -> IResult<&'a str, Vec<Int>, VerboseError<&'a str>> {
    separated_list1(delimiter, charge_parser)
}
