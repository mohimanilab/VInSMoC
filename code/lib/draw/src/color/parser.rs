use std::str::FromStr;

use nom::{
    branch::alt,
    bytes::complete::{tag, tag_no_case},
    character::complete::digit1,
    combinator::{all_consuming, map, map_res},
    error::VerboseError,
    number::complete::double,
    sequence::{delimited, preceded, tuple},
    IResult,
};

fn rgb_parser(input: &str) -> IResult<&str, crate::Color, VerboseError<&str>> {
    map(
        preceded(
            tag_no_case("RGB"),
            delimited(
                tag("("),
                tuple((
                    map_res(digit1, u8::from_str),
                    preceded(tag(","), map_res(digit1, u8::from_str)),
                    preceded(tag(","), map_res(digit1, u8::from_str)),
                )),
                tag(")"),
            ),
        ),
        crate::Color::from,
    )(input)
}

fn rgba_parser(input: &str) -> IResult<&str, crate::Color, VerboseError<&str>> {
    map(
        preceded(
            tag_no_case("RGBA"),
            delimited(
                tag("("),
                tuple((
                    map_res(digit1, u8::from_str),
                    preceded(tag(","), map_res(digit1, u8::from_str)),
                    preceded(tag(","), map_res(digit1, u8::from_str)),
                    preceded(tag(","), double),
                )),
                tag(")"),
            ),
        ),
        crate::Color::from,
    )(input)
}

fn hsl_parser(input: &str) -> IResult<&str, crate::Color, VerboseError<&str>> {
    map(
        preceded(
            tag_no_case("HSL"),
            delimited(
                tag("("),
                tuple((
                    double,
                    preceded(tag(","), double),
                    preceded(tag(","), double),
                )),
                tag(")"),
            ),
        ),
        crate::Color::from,
    )(input)
}

pub fn parse(input: &str) -> IResult<&str, crate::Color, VerboseError<&str>> {
    all_consuming(alt((rgb_parser, rgba_parser, hsl_parser)))(input)
}
