use std::{collections::HashMap, iter::FromIterator, path::Path};

use nom::{
    branch::alt,
    bytes::complete::{is_a, is_not, tag},
    character::complete::{alphanumeric1, space0},
    combinator::{all_consuming, eof},
    error::{convert_error, VerboseError},
    multi::{many0, many1},
    sequence::{delimited, pair, terminated, tuple},
    Err, IResult, Parser,
};

use crate::{take_until_end_of_line, word};

#[derive(Debug, thiserror::Error)]
pub enum FastaParseError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("failed parsing:\n{0}")]
    Parsing(String),
}

pub fn parse_fasta<T, U, Seq: FromIterator<U> + IntoIterator<Item = U>>(
    input: &str,
    parse_comment: impl for<'a> Parser<&'a str, T, VerboseError<&'a str>>,
    parse_seq: impl for<'a> Parser<&'a str, Seq, VerboseError<&'a str>>,
) -> Result<Vec<(String, T, Seq)>, FastaParseError> {
    // Convert parse_fasta_record into a parser that can be passed to a nom combinator
    let parser = parse_fasta_record(parse_comment, parse_seq);

    match all_consuming(many1(parser))(input) {
        Err(Err::Error(e)) | Err(Err::Failure(e)) => {
            Err(FastaParseError::Parsing(convert_error(input, e)))
        }
        Ok((_, res)) => {
            let mut out = Vec::<(String, T, Seq)>::new();
            for (a, b, c) in res {
                out.push((a.to_string(), b, c));
            }
            Ok(out)
        }
        _ => unreachable!(),
    }
}

pub fn parse_fasta_path<T, U, Seq: FromIterator<U> + IntoIterator<Item = U>>(
    input: impl AsRef<Path>,
    parse_comment: impl for<'a> Parser<&'a str, T, VerboseError<&'a str>>,
    parse_seq: impl for<'a> Parser<&'a str, Seq, VerboseError<&'a str>>,
) -> Result<Vec<(String, T, Seq)>, FastaParseError> {
    let input = std::fs::read_to_string(input.as_ref())?;
    parse_fasta(input.as_str(), parse_comment, parse_seq)
}

/// Builds a function to parse a fasta string by calling a group of helper parsers in series
fn parse_fasta_record<T, U, Seq: FromIterator<U> + IntoIterator<Item = U>>(
    mut parse_comment: impl for<'a> Parser<&'a str, T, VerboseError<&'a str>>,
    mut parse_seq: impl for<'a> Parser<&'a str, Seq, VerboseError<&'a str>>,
) -> impl for<'a> FnMut(&'a str) -> IResult<&'a str, (&'a str, T, Seq), VerboseError<&'a str>> {
    move |input: &str| {
        let (input, (id, comment, sequence)) =
            tuple((parse_fasta_id, parse_fasta_comment, parse_fasta_body))(input)?;

        // Filter out newlines in the sequence
        let line_parse = terminated(alphanumeric1, alt((is_a("\r\n"), eof)));
        let (_, seq_result_vec) = all_consuming(many1(line_parse))(sequence)?;

        let mut seq_vec = Vec::<Seq>::new();
        for line in seq_result_vec {
            let (_, val) = parse_seq.parse(line)?;
            seq_vec.push(val);
        }
        let seq_result = seq_vec.into_iter().flat_map(|s| s.into_iter()).collect();

        // Apply the provided parsers
        let (_, comment_res) = parse_comment.parse(comment)?;

        Ok((input, (id, comment_res, seq_result)))
    }
}

/// Parses out the id of a FASTA string.
///
/// Consumes and returns the id of a FASTA sring.
fn parse_fasta_id(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    let (input, _) = tag(">")(input)?;
    let (input, id) = word(input)?;

    Ok((input, id))
}

/// Parse the attributes of the FASTA string.
///
/// Useful for parsing FASTA strings that use [attr=value] syntax. Consumes all square bracket attributes
/// and returns them in a hashmap. Assumes we have already parsed the id.
///
/// # Examples
///
/// ```
/// use parsers::fasta::parse_fasta_attributes;
///
/// let result = parse_fasta_attributes("[a=b][c=d] [e=f] rest of comment");
///
/// let (rest, dict) = result.unwrap();
/// assert_eq!(rest, "rest of comment");
/// assert_eq!(dict["a"], "b");
/// assert_eq!(dict["c"], "d");
/// assert_eq!(dict["e"], "f");
/// ```
pub fn parse_fasta_attributes(
    input: &str,
) -> IResult<&str, HashMap<String, String>, VerboseError<&str>> {
    // Consume attributes
    let (input, attr_vector) = many0(square_bracket_attribute)(input)?;

    // Populate output map
    let mut map = HashMap::new();
    for (key, value) in attr_vector {
        map.insert(key.to_string(), value.to_string());
    }

    Ok((input, map))
}

/// Parses the comment of a FASTA string.
///
/// Consumes and returns the comment, if present. Assumes all other header elements have been parsed
/// such that `input` begins with a (potentially empty) string followed by a newline.
fn parse_fasta_comment(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    // Consume comment
    let (input, comment) = take_until_end_of_line(input)?;
    Ok((input, comment))
}

/// Parse the FASTA sequence.
///
/// Parses the sequence of the FASTA entry. Assumes the header has already been parsed.
/// Does not filter linebreaks in the sequence.
fn parse_fasta_body(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    let (input, sequence) = is_not(">")(input)?;
    Ok((input, sequence))
}

/// Parses the contents of square bracket attributes.
///
/// Consumes an attribute of the form: `[attr_i=value_i]` and any optional spaces that follow, and then returns the
/// pair: `(attr_i, value_i)`.
fn square_bracket_attribute(input: &str) -> IResult<&str, (&str, &str), VerboseError<&str>> {
    let (input, key) = delimited(tag("["), is_not("\r\n=]"), tag("="))(input)?;
    let (input, value) = terminated(is_not("\r\n]"), pair(tag("]"), space0))(input)?;
    Ok((input, (key, value)))
}

#[cfg(test)]
mod tests {
    use crate::fasta::{parse_fasta, parse_fasta_record};
    use crate::identity_parser;
    use nom::error::convert_error;
    use nom::multi::many1;
    use std::iter::FromIterator;

    #[test]
    fn parse_basic_test() {
        let fasta = ">foo [a=b] [c=d] bar\nabc\ndefg\n>another";
        let mut parser = parse_fasta_record(identity_parser, identity_parser);

        let res = parser(fasta);
        assert!(res.is_ok());

        let (rest, (id, comment, seq)) = res.unwrap();
        let (comment, seq) = (String::from_iter(comment), String::from_iter(seq));
        assert_eq!(rest, ">another");
        assert_eq!(id, "foo");
        assert_eq!(comment, "[a=b] [c=d] bar");
        assert_eq!(seq, "abcdefg");
    }

    #[test]
    fn parse_test_extra_lines_and_no_comment() {
        let fasta = ">a\nABC\nDEF\n\n>b\nGHI\nJKL";
        let mut parser = parse_fasta_record(identity_parser, identity_parser);

        let res = parser(fasta);
        assert!(res.is_ok());

        let (rest, (id, comment, seq)) = res.unwrap();
        let (comment, seq) = (String::from_iter(comment), String::from_iter(seq));
        assert_eq!(rest, ">b\nGHI\nJKL");
        assert_eq!(id, "a");
        assert_eq!(comment, "");
        assert_eq!(seq, "ABCDEF");

        // Now run the parser again on the second part
        let res = parser(rest);
        assert!(res.is_ok());

        let (rest, (id, comment, seq)) = res.unwrap();
        let (comment, seq) = (String::from_iter(comment), String::from_iter(seq));
        assert_eq!(rest, "");
        assert_eq!(id, "b");
        assert_eq!(comment, "");
        assert_eq!(seq, "GHIJKL");
    }

    #[test]
    fn bad_sequence_test() {
        let fasta = ">header\nA B C D E F";
        let mut parser = parse_fasta_record(identity_parser, identity_parser);

        let res = parser(fasta);
        assert!(res.is_err());

        if let Err(nom::Err::Error(e)) = res {
            // This print is unnecessary, but it illustrates how to display
            // a verbose error
            println!(
                "verbose error from parsing fasta:\n{}",
                convert_error(fasta, e)
            );
        }
    }

    #[test]
    fn no_id() {
        let fasta = ">\nABC\n>\nDEF";
        let parser = parse_fasta_record(identity_parser, identity_parser);

        let res = many1(parser)(fasta);
        assert!(res.is_ok());
        let (_, res) = res.unwrap();

        assert!(res.len() == 2);

        let (id, comment, seq) = res[0].clone();
        let (comment, seq) = (String::from_iter(comment), String::from_iter(seq));
        assert_eq!(id, "");
        assert_eq!(comment, "");
        assert_eq!(seq, "ABC");

        let (id, comment, seq) = res[1].clone();
        let (comment, seq) = (String::from_iter(comment), String::from_iter(seq));
        assert_eq!(id, "");
        assert_eq!(comment, "");
        assert_eq!(seq, "DEF");
    }

    #[test]
    fn no_seq() {
        let fasta = ">a\n";
        let mut parser = parse_fasta_record(identity_parser, identity_parser);

        let res = parser(fasta);
        println!("{:?}", res);
        assert!(res.is_err());
    }

    #[test]
    fn file_with_multiple_entries_test() {
        let result = parse_fasta(
            &std::fs::read_to_string("test_files/fasta/two_entries.fasta").unwrap(),
            identity_parser,
            identity_parser,
        );

        assert!(result.is_ok(), "{}", result.unwrap_err());
        let result = result.unwrap();
        assert!(result.len() == 2);
        let (id, comment, seq) = result[0].clone();
        let (comment, seq) = (String::from_iter(comment), String::from_iter(seq));
        assert_eq!(id, "a");
        assert_eq!(comment, "comment");
        assert_eq!(seq, "ABCD");

        let (id, comment, seq) = result[1].clone();
        let (comment, seq) = (String::from_iter(comment), String::from_iter(seq));
        assert_eq!(id, "b");
        assert_eq!(comment, "");
        assert_eq!(seq, "ABCD");
    }

    #[test]
    fn file_with_extra_garbage() {
        let result = parse_fasta(
            "test_files/fasta/not_all_consuming.fasta",
            identity_parser,
            identity_parser,
        );

        assert!(result.is_err());
    }
}
