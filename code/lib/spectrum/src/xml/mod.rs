use std::{
    convert::TryInto,
    io::{self, Read},
};

use flate2::bufread::ZlibDecoder;
use thiserror::Error;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub(crate) enum Endianness {
    Big,
    Little,
}

/// An error that can occur when decoding mzML or mzXML style spectral data.
///
/// This specifically decodes a vector of given precision, little-endian floating point values from
/// a base64-encoded string, with optional zlib compression. The floating points are always then
/// converted to [`f64`]. This appears to be a common-enough paradigm that we have provided it
/// here.
#[derive(Error, Debug)]
pub enum DecodeDataError {
    /// An error occurred when decoding the original base64 encoded bytes.
    #[error("unable to decode base64 data")]
    Base64(#[from] base64::DecodeError),
    /// An error occurred when decompressing the decoded data via zlib
    #[error("unable to decompress data")]
    Compression(#[from] io::Error),
    /// The provided precision cannot possibly be used to parse the whole input byte slice.
    ///
    /// This will happen if the number of bytes in precision (4 for 32, 8 for 64) does not divide
    /// evenly into the length of the decoded and decompressed byte slice. This is probably due to
    /// an issue with the original input.
    #[error("data with {len} bytes cannot be parsed using precision {precision}")]
    ImpossibleLength { precision: usize, len: usize },
    /// User provided a precision that wasn't 32 or 64 bits
    #[error("invalid floating point precision {0}, must be 32 or 64")]
    FloatPrecision(usize),
}

fn float_parser(precision: usize, endianness: Endianness) -> impl Fn(&[u8]) -> f64 {
    match precision {
        32 => match endianness {
            Endianness::Little => {
                |bytes: &[u8]| f32::from_le_bytes(bytes.try_into().unwrap()) as f64
            }
            Endianness::Big => |bytes: &[u8]| f32::from_be_bytes(bytes.try_into().unwrap()) as f64,
        },
        64 => match endianness {
            Endianness::Little => |bytes: &[u8]| f64::from_le_bytes(bytes.try_into().unwrap()),
            Endianness::Big => |bytes: &[u8]| f64::from_be_bytes(bytes.try_into().unwrap()),
        },
        _ => unreachable!("float precision should have been checked in decode_data_array"),
    }
}

/// Converts a base64-encoded string into a vector of doubles.
///
/// `precision` gives the number of bits expected in each floating-point value and `compressed` is
/// true if data was compressed using zlib. The precision is only used for parsing, the output will
/// always be converted to [`f64`].
///
/// Both mzML and mzXML formats use base64 encoding and optional zlib compression for the actual
/// peak data. For mzML data this will be used to decode either the m/z vector or the intensity
/// vector of the spectrum. For mzXML data this will be used to decode a single vector where each
/// pair of elements is a m/z-intensity pair (the first and third floats are m/z values and the
/// second and fourth floats are the corresponding intensity values).
///
/// # Errors
/// See [`DecodeDataError`] for possible error states.
pub(crate) fn decode_data_array(
    input: &[u8],
    precision: usize,
    compressed: bool,
    byte_order: Endianness,
) -> Result<Vec<f64>, DecodeDataError> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    let mut data = base64::decode(input)?;
    if compressed {
        let mut decompressed = Vec::with_capacity(data.len());
        ZlibDecoder::new(data.as_slice()).read_to_end(&mut decompressed)?;
        data = decompressed;
    }

    if precision != 32 && precision != 64 {
        return Err(DecodeDataError::FloatPrecision(precision));
    }

    let f = float_parser(precision, byte_order);

    let chunk_size = precision / 8;
    if data.len() % chunk_size != 0 {
        return Err(DecodeDataError::ImpossibleLength {
            precision,
            len: data.len(),
        });
    }

    let mut output = Vec::<f64>::with_capacity(data.len() / chunk_size);
    for chunk in data.chunks_exact(chunk_size) {
        output.push(f(chunk));
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_32bit_no_zlib() {
        // Basic 32-bit float no compression test
        let test_array: [f32; 5] = [100.00087, 100.00238, 100.00389, 100.00541, 119.271736];
        let b64_array = b"cgDIQjgByEL+AchCxQLIQiGL7kIjjO5CJY3uQieO7kIpj+5CK5DuQi2R7kIvku5CMpPuQjeU7kI5le5CO5buQj2X7kK6QgJDTUMCQ+FDAkN0RAJDB0UCQ5pFAkMuRgJDwUYCQ1RHAkPoRwJDe0gCQw9JAkOiSQJDj/YGQyr3BkPF9wZDYfgGQ/z4BkOX+QZDM/oGQ876BkNp+wZDBfwGQ6L8BkM9/QZD2P0GQ3T+BkPmBRVDmgYVQ04HFUMCCBVDtggVQ2oJFUMfChVD0woVQ4cLFUM7DBVD8AwVQ6QNFUNYDhVDDA8VQ8EPFUN1EBVDKREVQ/UFFkOrBhZDYQcWQxcIFkPNCBZDgwkWQzkKFkPvChZDpQsWQ1sMFkMRDRZDxw0WQ34OFkM0DxZDw6gZQ4CpGUM9qhlD+aoZQ7arGUNzrBlDL60ZQ+ytGUOprhlDZq8ZQyKwGUPfsBlDnLEZQ0iiGkMGoxpDxaMaQ4OkGkNCpRpDAKYaQ7+mGkN9pxpDPKgaQ/+oGkO9qRpDfKoaQzqrGkO0o0dDy6RHQ+OlR0P6pkdDEahHQympR0NAqkdDWKtHQ2+sR0OHrUdDnq5HQ7avR0PNsEdDaT1RQ5U+UUPAP1FD7EBRQxhCUUNEQ1FDcERRQ5xFUUPHRlFD9EdRQyBJUUNLSlFDd0tRQ3p0bEPidWxDS3dsQ7N4bEMbemxDg3tsQ+t8bENTfmxDvH9sQySBbEONgmxD9YNsQ12FbEPGhmxDDn6EQ+R+hEO5f4RDj4CEQ2WBhEM6goRDEIOEQ+aDhEO7hIRDkYWEQ2eGhEM8h4RDEoiEQ5x6hUN0e4VDTHyFQyR9hUP9fYVD1X6FQ61/hUOFgIVDXYGFQzWChUMNg4VD5YOFQ72EhUOVhYVDlnaGQ3B3hkNLeIZDJXmGQwB6hkPaeoZDtHuGQ498hkNpfYZDRH6GQx5/hkP5f4ZD04CGQ66BhkOIgoZDY4OGQz2EhkMM9oZD6PaGQ8T3hkOf+IZDe/mGQ1f6hkMy+4ZDDvyGQ+r8hkPF/YZDof6GQ33/hkNYAIdDNAGHQxACh0PrAodDG9CIQ/vQiEPb0YhDu9KIQ5vTiEN71IhDXNWIQzzWiEMc14hD/NeIQ9zYiEO82YhDnNqIQ33biENvQYxDV0KMQ0BDjEMpRIxDEUWMQ/pFjEPjRoxDy0eMQ7RIjEOdSYxDhkqMQ25LjENXTIxDT3+MQzmAjEMigYxDC4KMQ/WCjEPeg4xDx4SMQ7CFjEOahoxDg4eMQ2yIjENWiYxDP4qMQyiLjEMSjIxD+4yMQ+WNjEMhAI1DCwGNQ/YBjUPgAo1DywONQ7UEjUOgBY1DigaNQ3UHjUNgCI1DSgmNQzUKjUMfC41DCgyNQ/UMjUP0fY1D4H6NQ8x/jUO3gI1Do4GNQ4+CjUN7g41DZ4SNQ1KFjUM+ho1DKoeNQxaIjUMCiY1D7YmNQ9mKjUO6eY5DqXqOQ5d7jkOFfI5Dc32OQ2F+jkNQf45DPoCOQyyBjkMago5DCYOOQ/eDjkPlhI5D1IWOQ8KGjkOwh45Dn4iOQ3n6jkNo+45DWPyOQ0f9jkM3/o5DJv+OQxYAj0MFAY9D9QGPQ+QCj0PUA49DwwSPQ7IFj0OiBo9DkQePQ4EIj0P7f5VD+4CVQ/uBlUP7gpVD+4OVQ/uElUP7hZVD+4aVQ/uHlUP7iJVD/ImVQ/yKlUP8i5VD/IyVQ/yNlUP8jpVD/I+VQ4AtmkOMLppDmC+aQ6UwmkOxMZpDvTKaQ8kzmkPVNJpD4jWaQ+42mkP6N5pDBjmaQxM6mkMfO5pDbF3DQ+pew0NpYMND52HDQw==";

        let decoded = decode_data_array(b64_array, 32usize, false, Endianness::Little);

        assert!(decoded.is_ok());

        let decoded = decoded.unwrap();
        let decoded = decoded.into_iter().map(|x| x as f32).collect::<Vec<_>>();
        assert_eq!(test_array, decoded[..5]);
    }

    #[test]
    fn test_decode_failures() {
        // Make sure that decoding fails when used improperly in a few simple cases
        let b64_array = b"cgDIQjgByEL+AchCxQLIQiGL7kIjjO5CJY3uQieO7kIpj+5CK5DuQi2R7kIvku5CMpPuQjeU7kI5le5CO5buQj2X7kK6QgJDTUMCQ+FDAkN0RAJDB0UCQ5pFAkMuRgJDwUYCQ1RHAkPoRwJDe0gCQw9JAkOiSQJDj/YGQyr3BkPF9wZDYfgGQ/z4BkOX+QZDM/oGQ876BkNp+wZDBfwGQ6L8BkM9/QZD2P0GQ3T+BkPmBRVDmgYVQ04HFUMCCBVDtggVQ2oJFUMfChVD0woVQ4cLFUM7DBVD8AwVQ6QNFUNYDhVDDA8VQ8EPFUN1EBVDKREVQ/UFFkOrBhZDYQcWQxcIFkPNCBZDgwkWQzkKFkPvChZDpQsWQ1sMFkMRDRZDxw0WQ34OFkM0DxZDw6gZQ4CpGUM9qhlD+aoZQ7arGUNzrBlDL60ZQ+ytGUOprhlDZq8ZQyKwGUPfsBlDnLEZQ0iiGkMGoxpDxaMaQ4OkGkNCpRpDAKYaQ7+mGkN9pxpDPKgaQ/+oGkO9qRpDfKoaQzqrGkO0o0dDy6RHQ+OlR0P6pkdDEahHQympR0NAqkdDWKtHQ2+sR0OHrUdDnq5HQ7avR0PNsEdDaT1RQ5U+UUPAP1FD7EBRQxhCUUNEQ1FDcERRQ5xFUUPHRlFD9EdRQyBJUUNLSlFDd0tRQ3p0bEPidWxDS3dsQ7N4bEMbemxDg3tsQ+t8bENTfmxDvH9sQySBbEONgmxD9YNsQ12FbEPGhmxDDn6EQ+R+hEO5f4RDj4CEQ2WBhEM6goRDEIOEQ+aDhEO7hIRDkYWEQ2eGhEM8h4RDEoiEQ5x6hUN0e4VDTHyFQyR9hUP9fYVD1X6FQ61/hUOFgIVDXYGFQzWChUMNg4VD5YOFQ72EhUOVhYVDlnaGQ3B3hkNLeIZDJXmGQwB6hkPaeoZDtHuGQ498hkNpfYZDRH6GQx5/hkP5f4ZD04CGQ66BhkOIgoZDY4OGQz2EhkMM9oZD6PaGQ8T3hkOf+IZDe/mGQ1f6hkMy+4ZDDvyGQ+r8hkPF/YZDof6GQ33/hkNYAIdDNAGHQxACh0PrAodDG9CIQ/vQiEPb0YhDu9KIQ5vTiEN71IhDXNWIQzzWiEMc14hD/NeIQ9zYiEO82YhDnNqIQ33biENvQYxDV0KMQ0BDjEMpRIxDEUWMQ/pFjEPjRoxDy0eMQ7RIjEOdSYxDhkqMQ25LjENXTIxDT3+MQzmAjEMigYxDC4KMQ/WCjEPeg4xDx4SMQ7CFjEOahoxDg4eMQ2yIjENWiYxDP4qMQyiLjEMSjIxD+4yMQ+WNjEMhAI1DCwGNQ/YBjUPgAo1DywONQ7UEjUOgBY1DigaNQ3UHjUNgCI1DSgmNQzUKjUMfC41DCgyNQ/UMjUP0fY1D4H6NQ8x/jUO3gI1Do4GNQ4+CjUN7g41DZ4SNQ1KFjUM+ho1DKoeNQxaIjUMCiY1D7YmNQ9mKjUO6eY5DqXqOQ5d7jkOFfI5Dc32OQ2F+jkNQf45DPoCOQyyBjkMago5DCYOOQ/eDjkPlhI5D1IWOQ8KGjkOwh45Dn4iOQ3n6jkNo+45DWPyOQ0f9jkM3/o5DJv+OQxYAj0MFAY9D9QGPQ+QCj0PUA49DwwSPQ7IFj0OiBo9DkQePQ4EIj0P7f5VD+4CVQ/uBlUP7gpVD+4OVQ/uElUP7hZVD+4aVQ/uHlUP7iJVD/ImVQ/yKlUP8i5VD/IyVQ/yNlUP8jpVD/I+VQ4AtmkOMLppDmC+aQ6UwmkOxMZpDvTKaQ8kzmkPVNJpD4jWaQ+42mkP6N5pDBjmaQxM6mkMfO5pDbF3DQ+pew0NpYMND52HDQw==";

        // This should fail since the precision is invalid
        let decoded = decode_data_array(b64_array, 10usize, false, Endianness::Little);
        assert!(decoded.is_err());

        // This should fail since the precision is incorrect
        let decoded = decode_data_array(b64_array, 64usize, false, Endianness::Little);
        assert!(decoded.is_err());

        // This should fail since the data were not compressed
        let decoded = decode_data_array(b64_array, 32usize, true, Endianness::Little);
        assert!(decoded.is_err());
    }
}
