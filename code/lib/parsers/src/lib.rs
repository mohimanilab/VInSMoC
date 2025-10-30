use std::str::FromStr;

use nom::{
    bytes::complete::{take, take_till},
    character::complete::{line_ending, not_line_ending, space0},
    combinator::map_res,
    error::VerboseError,
    sequence::tuple,
    IResult,
};

pub mod charge;
pub mod fasta;

/// A nom parser to consume up to and including whitespace.
///  
/// Consumes and returns everything up to the first whitespace character. If no such character is found,
/// it will still succeed and return everything. __If the terminator is not a newline, it will consume
/// the terminator (and any following spaces) as well__.
///
/// # Examples
///
/// ```
/// use parsers::word;
///
/// assert_eq!(word("hello goodbye"), Ok(("goodbye", "hello")));
/// assert_eq!(word("hello\ngoodbye"), Ok(("\ngoodbye", "hello")));
/// assert_eq!(word("hello   \t  \ngoodbye"), Ok(("\ngoodbye", "hello")));
/// ```
pub fn word(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    let (input, until_white) = take_till(|c: char| c.is_whitespace())(input)?;
    let (input, _) = space0(input)?;
    Ok((input, until_white))
}

/// Consumes and returns everything up to and including the next newline character.
///
/// Will fail if it encounters EOF before the end of a line.
///
/// # Examples
///
/// ```
/// use parsers::take_until_end_of_line;
///
/// assert_eq!(take_until_end_of_line("hello\ngoodbye"), Ok(("goodbye", "hello")));
/// assert_eq!(take_until_end_of_line("hello\n\ngoodbye"), Ok(("\ngoodbye", "hello")));
/// assert_eq!(take_until_end_of_line("hello  \r\ngoodbye"), Ok(("goodbye", "hello  ")));
/// assert!(take_until_end_of_line("hello").is_err());
/// ```
pub fn take_until_end_of_line(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    let (input, (line, _)) = tuple((not_line_ending, line_ending))(input)?;
    Ok((input, line))
}

/// A verbose parser that simply consumes and returns everything as a char Vec for testing purposes.
pub fn identity_parser(input: &str) -> IResult<&str, Vec<char>, VerboseError<&str>> {
    let (input, res) = take_till(|_c| false)(input)?;
    Ok((input, res.chars().collect()))
}

/// Implements the common pattern of taking n characters and then trying to parse them into some type of interest.
///
/// # Examples
/// ```
/// use parsers::take_n_try_parse;
///
/// let mut two_digits = take_n_try_parse::<u8>(2);
///
/// assert_eq!(two_digits("12"), Ok(("", 12)));
/// assert_eq!(two_digits("25hello"), Ok(("hello", 25)));
///
/// // Fails to read two bytes
/// assert!(two_digits("2").is_err());
///
/// // Fails to parse the bytes taken
/// assert!(two_digits("hi").is_err());
/// ```
pub fn take_n_try_parse<T: FromStr>(
    n: usize,
) -> impl FnMut(&str) -> IResult<&str, T, VerboseError<&str>> {
    move |input: &str| map_res(take(n), |s: &str| s.parse::<T>())(input)
}
