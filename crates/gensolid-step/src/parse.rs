//! Independent ISO-10303-21 physical-file and entity parser.
//!
//! Grammar reference: https://www.steptools.com/stds/step/IS_final_p21e3.html .

use std::collections::{HashMap, HashSet};

use gensolid_core::{Error, Result};

// Hide comments and quoted strings while locating section delimiters. Byte offsets
// remain unchanged, including for UTF-8 names and doubled quote escapes.
fn visible_text(text: &str) -> Result<String> {
    let bytes = text.as_bytes();
    let mut mask = bytes.to_vec();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i..].starts_with(b"/*") {
            let Some(end) = text[i + 2..].find("*/") else {
                return Err(Error::parse("unterminated STEP comment"));
            };
            let end = i + end + 4;
            mask[i..end].fill(b' ');
            i = end;
        } else if bytes[i] == b'\'' {
            let start = i;
            i += 1;
            loop {
                if i >= bytes.len() {
                    return Err(Error::parse("unterminated STEP string"));
                }
                if bytes[i] == b'\'' {
                    i += 1;
                    if bytes.get(i) == Some(&b'\'') {
                        i += 1;
                    } else {
                        break;
                    }
                } else {
                    i += 1;
                }
            }
            mask[start..i].fill(b' ');
        } else {
            i += 1;
        }
    }
    String::from_utf8(mask).map_err(|_| Error::parse("invalid STEP text"))
}

#[cfg(test)]
fn data_bounds(text: &str) -> Result<(usize, usize)> {
    let visible = visible_text(text)?;
    let mut offset = 0;
    let mut start = None;
    for statement in visible.split_inclusive(';') {
        let word = statement.trim().trim_end_matches(';').trim();
        if start.is_none() && word.eq_ignore_ascii_case("DATA") {
            start = Some(offset + statement.len());
        } else if let Some(start) = start
            && word.eq_ignore_ascii_case("ENDSEC")
        {
            return Ok((start, offset));
        }
        offset += statement.len();
    }
    Err(Error::parse("STEP file has no closed DATA section"))
}

/// Validates the envelope and header and reads every DATA section in one namespace.
/// External references and signatures require semantics this reader does not implement.
pub(crate) fn parse_document(text: &str) -> Result<HashMap<u32, Instance>> {
    let visible = visible_text(text)?;
    let mut statements = Vec::new();
    let mut offset = 0;
    for statement in visible.split_inclusive(';') {
        let end = offset + statement.len();
        let word = statement.trim();
        if !word.is_empty() {
            if !word.ends_with(';') {
                return Err(Error::parse("STEP statement is missing ';'"));
            }
            statements.push((offset, end, word.trim_end_matches(';').trim()));
        }
        offset = end;
    }
    let mut index = 0;
    let first = statements.get(index).map(|s| s.2);
    if !first.is_some_and(|word| word.eq_ignore_ascii_case("ISO-10303-21")) {
        return Err(Error::parse("STEP physical file header is missing"));
    }
    index += 1;
    if !statements
        .get(index)
        .is_some_and(|s| s.2.eq_ignore_ascii_case("HEADER"))
    {
        return Err(Error::parse("STEP HEADER section is missing"));
    }
    index += 1;
    let mut schemas = Vec::new();
    for required in ["FILE_DESCRIPTION", "FILE_NAME", "FILE_SCHEMA"] {
        let (start, end, _) = *statements
            .get(index)
            .ok_or_else(|| Error::parse(format!("STEP HEADER is missing {required}")))?;
        let (name, args) = record(&text[start..end])?;
        if name != required {
            return Err(Error::parse(format!(
                "STEP HEADER expected {required}, found {name}"
            )));
        }
        match required {
            "FILE_DESCRIPTION" => match args.as_slice() {
                [description, Value::Str(_)] => {
                    string_list(description)?;
                }
                _ => return Err(Error::parse("invalid FILE_DESCRIPTION fields")),
            },
            "FILE_NAME" => match args.as_slice() {
                [
                    Value::Str(_),
                    Value::Str(_),
                    authors,
                    organizations,
                    Value::Str(_),
                    Value::Str(_),
                    Value::Str(_),
                ] => {
                    string_list(authors)?;
                    string_list(organizations)?;
                }
                _ => return Err(Error::parse("invalid FILE_NAME fields")),
            },
            "FILE_SCHEMA" => match args.as_slice() {
                [value] => {
                    schemas = string_list(value)?
                        .into_iter()
                        .map(str::to_ascii_uppercase)
                        .collect()
                }
                _ => return Err(Error::parse("invalid FILE_SCHEMA fields")),
            },
            _ => return Err(Error::internal("unhandled mandatory STEP header record")),
        }
        index += 1;
    }
    loop {
        let (start, end, word) = *statements
            .get(index)
            .ok_or_else(|| Error::parse("STEP HEADER has no ENDSEC"))?;
        index += 1;
        if word.eq_ignore_ascii_case("ENDSEC") {
            break;
        }
        let (name, _) = record(&text[start..end])?;
        if matches!(
            name.as_str(),
            "FILE_DESCRIPTION" | "FILE_NAME" | "FILE_SCHEMA"
        ) {
            return Err(Error::parse(format!("duplicate STEP header record {name}")));
        }
    }

    let mut instances = HashMap::new();
    let mut section_names = HashSet::new();
    let mut section_count = 0;
    let mut unnamed = false;
    loop {
        let (start, end, word) = *statements
            .get(index)
            .ok_or_else(|| Error::parse("STEP physical file trailer is missing"))?;
        index += 1;
        if word.eq_ignore_ascii_case("END-ISO-10303-21") {
            if section_count == 0 {
                return Err(Error::parse("STEP file has no DATA section"));
            }
            if statements
                .get(index)
                .is_some_and(|s| s.2.starts_with("SIGNATURE"))
            {
                return Err(Error::unimplemented(
                    "STEP SIGNATURE section is not implemented",
                ));
            }
            if index != statements.len() {
                return Err(Error::parse("unexpected content after STEP trailer"));
            }
            break;
        }
        if word.eq_ignore_ascii_case("ANCHOR") || word.eq_ignore_ascii_case("REFERENCE") {
            return Err(Error::unimplemented(format!(
                "STEP {word} section is not implemented"
            )));
        }
        if word.eq_ignore_ascii_case("DATA") {
            if section_count != 0 || schemas.len() != 1 {
                return Err(Error::parse(
                    "unnamed STEP DATA needs a single section and schema",
                ));
            }
            unnamed = true;
        } else {
            let (name, args) = record(&text[start..end])?;
            if name != "DATA" {
                return Err(Error::parse(format!(
                    "expected STEP DATA section, found {name}"
                )));
            }
            let [Value::Str(section), value] = args.as_slice() else {
                return Err(Error::parse("invalid STEP DATA section parameters"));
            };
            let section_schemas = string_list(value)?;
            if unnamed
                || section_schemas.len() != 1
                || !schemas.contains(&section_schemas[0].to_ascii_uppercase())
            {
                return Err(Error::parse(
                    "STEP DATA section schema is not declared in FILE_SCHEMA",
                ));
            }
            if !section_names.insert(section.clone()) {
                return Err(Error::parse("duplicate STEP DATA section name"));
            }
        }
        section_count += 1;
        let data_start = end;
        let data_end = loop {
            let (start, _, word) = *statements
                .get(index)
                .ok_or_else(|| Error::parse("STEP DATA section has no ENDSEC"))?;
            index += 1;
            if word.eq_ignore_ascii_case("ENDSEC") {
                break start;
            }
        };
        parse_instances(&text[data_start..data_end], &mut instances)?;
    }
    if instances.is_empty() {
        return Err(Error::parse("STEP DATA sections are empty"));
    }
    let mut ids: Vec<_> = instances.keys().copied().collect();
    ids.sort_unstable();
    for id in ids {
        let entity = instances
            .get(&id)
            .ok_or_else(|| Error::internal("parsed STEP entity is missing"))?;
        validate_references(&entity.args, id, &instances)?;
        let mut components: Vec<_> = entity.components.iter().collect();
        components.sort_by(|a, b| a.0.cmp(b.0));
        for (_, args) in components {
            validate_references(args, id, &instances)?;
        }
    }
    Ok(instances)
}

fn string_list(value: &Value) -> Result<Vec<&str>> {
    let Value::List(values) = value else {
        return Err(Error::parse("expected a STEP list of strings"));
    };
    if values.is_empty() {
        return Err(Error::parse("STEP string list must not be empty"));
    }
    values
        .iter()
        .map(|value| match value {
            Value::Str(text) => Ok(text.as_str()),
            _ => Err(Error::parse("expected a string in STEP list")),
        })
        .collect()
}

fn record(text: &str) -> Result<(String, Vec<Value>)> {
    let mut cursor = Cursor {
        text,
        index: 0,
        depth: 0,
    };
    cursor.skip_ws();
    let name = cursor.type_name()?;
    cursor.skip_ws();
    if !cursor.eat('(') {
        return Err(Error::parse(format!("STEP {name} has no argument list")));
    }
    let args = cursor.arguments()?;
    cursor.skip_ws();
    if !cursor.eat(';') {
        return Err(Error::parse(format!("STEP {name} is missing ';'")));
    }
    cursor.skip_ws();
    if !cursor.done() {
        return Err(Error::parse(format!("unexpected text after STEP {name}")));
    }
    Ok((name, args))
}

fn validate_references(
    values: &[Value],
    owner: u32,
    instances: &HashMap<u32, Instance>,
) -> Result<()> {
    for value in values {
        match value {
            Value::Ref(id) if !instances.contains_key(id) => {
                return Err(Error::parse(format!(
                    "STEP #{owner} refers to missing entity #{id}"
                )));
            }
            Value::List(items) | Value::Typed(_, items) => {
                validate_references(items, owner, instances)?
            }
            _ => (),
        }
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Instance {
    pub type_name: String,
    pub args: Vec<Value>,
    pub components: HashMap<String, Vec<Value>>,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Value {
    Str(String),
    Num(f64),
    Ref(u32),
    Bool(bool),
    Enum(String),
    Typed(String, Vec<Value>),
    Star,
    Unset,
    List(Vec<Value>),
}

#[cfg(test)]
pub(crate) fn parse_data(text: &str) -> Result<HashMap<u32, Instance>> {
    let (start, end) = data_bounds(text)?;
    let mut instances = HashMap::new();
    parse_instances(&text[start..end], &mut instances)?;
    if instances.is_empty() {
        return Err(Error::parse("STEP DATA section is empty"));
    }
    Ok(instances)
}

fn parse_instances(text: &str, instances: &mut HashMap<u32, Instance>) -> Result<()> {
    let mut cursor = Cursor {
        text,
        index: 0,
        depth: 0,
    };
    loop {
        cursor.skip_ws();
        if cursor.done() {
            break;
        }
        if !cursor.eat('#') {
            return Err(Error::parse(format!(
                "expected an entity id near {}",
                cursor.snippet()
            )));
        }
        let id = cursor.integer()?;
        cursor.skip_ws();
        if !cursor.eat('=') {
            return Err(Error::parse(format!("entity #{id} has no '='")));
        }
        cursor.skip_ws();
        let mut components = HashMap::new();
        let (type_name, args) = if cursor.eat('(') {
            loop {
                cursor.skip_ws();
                if cursor.eat(')') {
                    break;
                }
                let name = cursor.type_name()?;
                cursor.skip_ws();
                if !cursor.eat('(') {
                    return Err(Error::parse("complex entity has no argument list"));
                }
                let args = cursor.arguments()?;
                if components.insert(name, args).is_some() {
                    return Err(Error::parse(format!("duplicate component in STEP #{id}")));
                }
            }
            if components.is_empty() {
                return Err(Error::parse("empty complex STEP entity"));
            }
            ("COMPLEX".to_owned(), Vec::new())
        } else {
            let name = cursor.type_name()?;
            cursor.skip_ws();
            if !cursor.eat('(') {
                return Err(Error::parse(format!("{name} #{id} has no argument list")));
            }
            (name, cursor.arguments()?)
        };
        cursor.skip_ws();
        if !cursor.eat(';') {
            return Err(Error::parse(format!("{type_name} #{id} is missing ';'")));
        }
        if id == 0
            || instances
                .insert(
                    id,
                    Instance {
                        type_name,
                        args,
                        components,
                    },
                )
                .is_some()
        {
            return Err(Error::parse(format!(
                "invalid or duplicate STEP entity #{id}"
            )));
        }
    }
    Ok(())
}

struct Cursor<'a> {
    text: &'a str,
    index: usize,
    depth: usize,
}

impl Cursor<'_> {
    fn done(&self) -> bool {
        self.index >= self.text.len()
    }

    fn rest(&self) -> &str {
        &self.text[self.index..]
    }

    fn snippet(&self) -> String {
        let rest = self.rest();
        let end = rest
            .char_indices()
            .nth(24)
            .map(|(i, _)| i)
            .unwrap_or(rest.len());
        rest[..end].replace('\n', " ")
    }

    fn skip_ws(&mut self) {
        while let Some(c) = self.rest().chars().next() {
            if self.rest().starts_with("/*") {
                if let Some(end) = self.rest()[2..].find("*/") {
                    self.index += end + 4;
                    continue;
                }
                // data_bounds already rejects unterminated comments.
                break;
            }
            if c.is_whitespace() {
                self.index += c.len_utf8();
            } else {
                break;
            }
        }
    }

    fn eat(&mut self, expected: char) -> bool {
        if self.rest().starts_with(expected) {
            self.index += expected.len_utf8();
            true
        } else {
            false
        }
    }

    fn integer(&mut self) -> Result<u32> {
        let rest = self.rest();
        let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
        if digits == 0 {
            return Err(Error::parse(format!(
                "expected an integer near {}",
                self.snippet()
            )));
        }
        let value = rest[..digits].parse().map_err(|_| {
            Error::parse(format!(
                "entity id does not fit in u32 near {}",
                self.snippet()
            ))
        })?;
        self.index += digits;
        Ok(value)
    }

    fn type_name(&mut self) -> Result<String> {
        let rest = self.rest();
        if rest.starts_with('!') {
            return Err(Error::unimplemented(
                "STEP user-defined entity keywords are not implemented",
            ));
        }
        if !rest.as_bytes().first().is_some_and(u8::is_ascii_alphabetic) {
            return Err(Error::parse(format!(
                "expected a STEP keyword near {}",
                self.snippet()
            )));
        }
        let len = rest
            .bytes()
            .take_while(|b| b.is_ascii_alphabetic() || b.is_ascii_digit() || *b == b'_')
            .count();
        if len == 0 {
            return Err(Error::parse(format!(
                "expected a type name near {}",
                self.snippet()
            )));
        }
        let name = rest[..len].to_ascii_uppercase();
        self.index += len;
        Ok(name)
    }

    fn arguments(&mut self) -> Result<Vec<Value>> {
        if self.depth >= 128 {
            return Err(Error::parse("STEP argument nesting limit exceeded"));
        }
        self.depth += 1;
        let mut args = Vec::new();
        loop {
            self.skip_ws();
            if self.eat(')') {
                break;
            }
            args.push(self.value()?);
            self.skip_ws();
            if self.eat(',') {
                self.skip_ws();
                if self.rest().starts_with(')') {
                    return Err(Error::parse("trailing STEP argument comma"));
                }
                continue;
            }
            if self.eat(')') {
                break;
            }
            return Err(Error::parse(format!(
                "expected ',' or ')' near {}",
                self.snippet()
            )));
        }
        self.depth -= 1;
        Ok(args)
    }

    fn value(&mut self) -> Result<Value> {
        self.skip_ws();
        let Some(c) = self.rest().chars().next() else {
            return Err(Error::parse("unexpected end of STEP arguments"));
        };
        match c {
            '\'' => Ok(Value::Str(self.string()?)),
            '#' => {
                self.index += 1;
                Ok(Value::Ref(self.integer()?))
            }
            '(' => {
                self.index += 1;
                Ok(Value::List(self.arguments()?))
            }
            '*' => {
                self.index += 1;
                Ok(Value::Star)
            }
            '$' => {
                self.index += 1;
                Ok(Value::Unset)
            }
            '.' if self
                .rest()
                .as_bytes()
                .get(1)
                .is_some_and(u8::is_ascii_digit) =>
            {
                self.number()
            }
            '.' => self.enumeration(),
            '+' | '-' | '0'..='9' => self.number(),
            'A'..='Z' | 'a'..='z' => {
                let name = self.type_name()?;
                self.skip_ws();
                if !self.eat('(') {
                    return Err(Error::parse("typed STEP value has no argument list"));
                }
                let args = self.arguments()?;
                if args.len() != 1 {
                    return Err(Error::parse("typed STEP parameters need exactly one value"));
                }
                Ok(Value::Typed(name, args))
            }
            '"' | '@' | '!' => Err(Error::unimplemented(
                "STEP binary, constant/value references or user-defined parameters are not implemented",
            )),
            _ => Err(Error::parse(format!(
                "unexpected STEP token near {}",
                self.snippet()
            ))),
        }
    }

    fn string(&mut self) -> Result<String> {
        self.index += 1;
        let mut out = String::new();
        loop {
            let Some(c) = self.rest().chars().next() else {
                return Err(Error::parse("unterminated STEP string"));
            };
            self.index += c.len_utf8();
            if c == '\'' {
                if self.rest().starts_with('\'') {
                    self.index += 1;
                    out.push('\'');
                } else {
                    break;
                }
            } else {
                out.push(c);
            }
        }
        Ok(out)
    }

    fn enumeration(&mut self) -> Result<Value> {
        self.index += 1;
        let rest = self.rest();
        let end = rest
            .find('.')
            .ok_or_else(|| Error::parse("unterminated STEP enumeration"))?;
        let name = rest[..end].to_ascii_uppercase();
        self.index += end + 1;
        match name.as_str() {
            "T" => Ok(Value::Bool(true)),
            "F" => Ok(Value::Bool(false)),
            _ if name
                .bytes()
                .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')
                && !name.is_empty() =>
            {
                Ok(Value::Enum(name))
            }
            _ => Err(Error::parse("invalid STEP enumeration")),
        }
    }

    fn number(&mut self) -> Result<Value> {
        let rest = self.rest();
        let mut len = 0;
        let bytes = rest.as_bytes();
        if bytes.first().is_some_and(|b| *b == b'+' || *b == b'-') {
            len += 1;
        }
        let start_digits = len;
        while len < bytes.len() && bytes[len].is_ascii_digit() {
            len += 1;
        }
        if len < bytes.len() && bytes[len] == b'.' {
            len += 1;
            while len < bytes.len() && bytes[len].is_ascii_digit() {
                len += 1;
            }
        }
        if len == start_digits
            || (len == start_digits + 1 && bytes.get(start_digits) == Some(&b'.'))
        {
            return Err(Error::parse(format!(
                "expected a number near {}",
                self.snippet()
            )));
        }
        if len < bytes.len() && (bytes[len] == b'e' || bytes[len] == b'E') {
            let mark = len;
            len += 1;
            if len < bytes.len() && (bytes[len] == b'+' || bytes[len] == b'-') {
                len += 1;
            }
            let exp = len;
            while len < bytes.len() && bytes[len].is_ascii_digit() {
                len += 1;
            }
            if len == exp {
                len = mark;
            }
        }
        let number: f64 = rest[..len]
            .parse()
            .map_err(|_| Error::parse(format!("invalid number near {}", self.snippet())))?;
        self.index += len;
        if !number.is_finite() {
            return Err(Error::parse("STEP number is not finite"));
        }
        Ok(Value::Num(number))
    }
}

pub(crate) fn name_skipped(args: &[Value]) -> Result<&[Value]> {
    match args {
        [Value::Str(_), rest @ ..] => Ok(rest),
        _ => Err(Error::parse("STEP entity is missing its name argument")),
    }
}

pub(crate) fn ref_of(value: &Value) -> Result<u32> {
    match value {
        Value::Ref(id) => Ok(*id),
        _ => Err(Error::parse("expected a STEP entity reference")),
    }
}

pub(crate) fn bool_of(value: &Value) -> Result<bool> {
    match value {
        Value::Bool(flag) => Ok(*flag),
        _ => Err(Error::parse("expected a STEP boolean")),
    }
}

pub(crate) fn nums(value: &Value) -> Result<Vec<f64>> {
    match value {
        Value::List(items) => items
            .iter()
            .map(|item| match item {
                Value::Num(number) => Ok(*number),
                _ => Err(Error::parse("expected a number in a STEP list")),
            })
            .collect(),
        _ => Err(Error::parse("expected a STEP list of numbers")),
    }
}

pub(crate) fn refs(value: &Value) -> Result<Vec<u32>> {
    match value {
        Value::List(items) => items.iter().map(ref_of).collect(),
        Value::Ref(id) => Ok(vec![*id]),
        _ => Err(Error::parse("expected a STEP list of references")),
    }
}
