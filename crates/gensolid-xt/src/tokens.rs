//! Tokenizer for the text node stream.
//!
//! Classifies whitespace-separated tokens; interpreting them as node fields requires the
//! schema named in the preamble and is done by the strict decoder in `import`.
//! These tokens are only diagnostic counts: XT length-prefixed strings and
//! character fields do not obey ordinary whitespace-separated token framing.

/// Lexical class of a token.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TokenKind {
    /// Optionally signed decimal integer.
    Integer,
    /// Decimal floating-point number with a fraction or exponent.
    Real,
    /// Anything else.
    Word,
}

/// A token and its byte offset in the input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Token<'a> {
    /// Token text.
    pub text: &'a str,
    /// Byte offset of the first character.
    pub offset: usize,
    /// Lexical class.
    pub kind: TokenKind,
}

fn classify(text: &str) -> TokenKind {
    let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
    if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
        TokenKind::Integer
    } else if digits.starts_with(|c: char| c.is_ascii_digit() || c == '.')
        && text.parse::<f64>().is_ok_and(f64::is_finite)
    {
        TokenKind::Real
    } else {
        TokenKind::Word
    }
}

/// Splits `input` into tokens with offsets.
pub fn tokenize(input: &str) -> impl Iterator<Item = Token<'_>> {
    let mut rest = input;
    let mut base = 0;
    std::iter::from_fn(move || {
        let start = rest.find(|c: char| !c.is_whitespace())?;
        let tail = &rest[start..];
        let len = tail.find(char::is_whitespace).unwrap_or(tail.len());
        let text = &tail[..len];
        let token = Token {
            text,
            offset: base + start,
            kind: classify(text),
        };
        base += start + len;
        rest = &tail[len..];
        Some(token)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_tokens_with_offsets() {
        let tokens: Vec<_> = tokenize("  12 -3 0.5 1e3 T51 :\n+7").collect();
        let kinds: Vec<_> = tokens.iter().map(|t| t.kind).collect();
        assert_eq!(
            kinds,
            [
                TokenKind::Integer,
                TokenKind::Integer,
                TokenKind::Real,
                TokenKind::Real,
                TokenKind::Word,
                TokenKind::Word,
                TokenKind::Integer
            ]
        );
        assert_eq!(tokens[0].offset, 2);
        assert_eq!(tokens[6].offset, 22);
        assert_eq!(tokens[6].text, "+7");
    }
}
