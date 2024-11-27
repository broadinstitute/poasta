use std::fmt::{self, Display, Write};

/// Represents the DNA strand
#[derive(Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse
}

/// Represents a GFA tag and value
#[derive(Debug)]
pub struct Field {
    pub tag: String,
    pub value: FieldValue
}

impl TryFrom<&str> for Field {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let mut parts = value.trim().splitn(3, ':');
        
        let tag = parts.next().ok_or("No tag")?;
        let value_type = parts.next().ok_or("No type")?;
        let value = parts.next().ok_or("No value")?;
        
        let value = match value_type {
            "A" | "Z" => FieldValue::String(value.to_string()),
            "i" => FieldValue::Integer(value.parse().map_err(|_| "Could not parse integer")?),
            "f" => FieldValue::Float(value.parse().map_err(|_| "Could not parse float")?),
            "J" => FieldValue::Json(value.to_string()),
            "H" => {
                let hex_data: Result<Vec<u8>, _> = value.as_bytes().chunks(2)
                    // TODO: use from_ascii_radix when added to the rust library: 
                    // https://github.com/rust-lang/libs-team/issues/469
                    // Now we need to incur UTF-8 validation overhead...
                    .map(|v| u8::from_str_radix(
                        std::str::from_utf8(v)
                            .map_err(|_| "Could not convert bytes to UTF-8")?,
                        16
                    ).map_err(|_| "Could not parse hex"))
                    .collect();
                
                FieldValue::ByteArray(hex_data?)
            },
            "B" => {
                let mut parts = value.split(',');
                let num_type = parts.next().ok_or("No number type available")?;
                
                match num_type {
                    "c" => FieldValue::NumberListI8(
                        parts.map(|x| x.parse()
                            .map_err(|_| "Could not parse integer"))
                        .collect::<Result<Vec<_>, _>>()?
                    ),
                    "C" => FieldValue::NumberListU8(
                        parts.map(|x| x.parse()
                            .map_err(|_| "Could not parse integer"))
                        .collect::<Result<Vec<_>, _>>()?
                    ),
                    "s" => FieldValue::NumberListI16(
                        parts.map(|x| x.parse()
                            .map_err(|_| "Could not parse integer"))
                        .collect::<Result<Vec<_>, _>>()?
                    ),
                    "S" => FieldValue::NumberListU16(
                        parts.map(|x| x.parse()
                            .map_err(|_| "Could not parse integer"))
                        .collect::<Result<Vec<_>, _>>()?
                    ),
                    "i" => FieldValue::NumberListI32(
                        parts.map(|x| x.parse()
                            .map_err(|_| "Could not parse integer"))
                        .collect::<Result<Vec<_>, _>>()?
                    ),
                    "I" => FieldValue::NumberListU32(
                        parts.map(|x| x.parse()
                            .map_err(|_| "Could not parse integer"))
                        .collect::<Result<Vec<_>, _>>()?
                    ),
                    "f" => FieldValue::NumberListF32(
                        parts.map(|x| x.parse()
                            .map_err(|_| "Could not parse float"))
                        .collect::<Result<Vec<_>, _>>()?
                    ),
                    _ => return Err("Invalid number type")
                }
            },
            _ => return Err("Invalid tag")
        };
        
        Ok(Field { tag: tag.to_string(), value })
    }
}

impl Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let value_type_tag = match self.value {
            FieldValue::String(_) => "Z",
            FieldValue::Integer(_) => "i",
            FieldValue::Float(_) => "f",
            FieldValue::Json(_) => "J",
            FieldValue::ByteArray(_) => "H",
            _ => "B"
        };
        write!(f, "{}:{}:{}", self.tag, value_type_tag, self.value)
    }
}


/// Represents the various types of values a GFA field can have
#[derive(Debug, PartialEq)]
pub enum FieldValue {
    String(String),
    Integer(i64),
    Float(f32),
    Json(String),
    ByteArray(Vec<u8>),
    NumberListI8(Vec<i8>),
    NumberListU8(Vec<u8>),
    NumberListI16(Vec<i16>),
    NumberListU16(Vec<u16>),
    NumberListI32(Vec<i32>),
    NumberListU32(Vec<u32>),
    NumberListF32(Vec<f32>),
}

impl Display for FieldValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FieldValue::String(s) => write!(f, "{}", s),
            FieldValue::Integer(v) => write!(f, "{}", v),
            FieldValue::Float(v) => write!(f, "{}", v),
            FieldValue::Json(v) => write!(f, "{}", v),
            FieldValue::ByteArray(v) => {
                let byte_str = v.iter()
                    .fold(String::new(), |mut output, e| {
                        let _ = write!(output, "{:02X}", e);
                        output
                    });
                
                write!(f, "{}", byte_str)
            },
            FieldValue::NumberListI8(v) => {
                write!(f, "c")?;
                v.iter().try_for_each(|e| write!(f, ",{}", e))
            },
            FieldValue::NumberListU8(v) => {
                write!(f, "C")?;
                v.iter().try_for_each(|e| write!(f, ",{}", e))
            },
            FieldValue::NumberListI16(v) => {
                write!(f, "s")?;
                v.iter().try_for_each(|e| write!(f, ",{}", e))
            },
            FieldValue::NumberListU16(v) => {
                write!(f, "S")?;
                v.iter().try_for_each(|e| write!(f, ",{}", e))
            },
            FieldValue::NumberListI32(v) => {
                write!(f, "i")?;
                v.iter().try_for_each(|e| write!(f, ",{}", e))
            },
            FieldValue::NumberListU32(v) => {
                write!(f, "I")?;
                v.iter().try_for_each(|e| write!(f, ",{}", e))
            },
            FieldValue::NumberListF32(v) => {
                write!(f, "f")?;
                v.iter().try_for_each(|e| write!(f, ",{}", e))
            },
        }
    }
}


/// Represents a GFA header
#[derive(Debug)]
pub struct Header {
    pub version: Option<String>,
    pub fields: Vec<Field>,
}

impl TryFrom<&str> for Header {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let mut parts = value.trim().splitn(3, '\t');
        
        let start = parts.next().ok_or("Empty line?")?;
        if start != "H" {
            return Err("Not a header line");
        }
        
        let version = parts.next()
            .and_then(|v| Field::try_from(v).ok())
            .and_then(|v| {
                if v.tag == "VN" {
                    match v.value {
                        FieldValue::String(v) => Some(v),
                        _ => None
                    }
                } else {
                    None
                }
            });
        
        let fields = if let Some(extra_fields) = parts.next() {
            extra_fields.split('\t')
                .filter_map(|v| Field::try_from(v).ok())
                .collect::<Vec<_>>()
        } else {
            Vec::default()
        };
        
        Ok(Header { version, fields })
    }
}

/// Represents a GFA segment
pub struct Segment {
    pub sid: String,
    pub sequence: Option<String>,
    pub fields: Vec<Field>,
}

impl TryFrom<&str> for Segment {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let mut parts = value.trim().splitn(4, '\t');
        
        let start = parts.next().ok_or("Empty line?")?;
        if start != "S" {
            return Err("Not a segment line");
        }
        
        let sid = parts.next()
            .ok_or("No sid")?;
        
        let sequence = parts.next()
            .and_then(|v| if v != "*" { 
                Some(v.to_ascii_uppercase())
            } else { 
                None 
            });
        
        let fields = if let Some(extra_fields) = parts.next() {
            extra_fields.split('\t')
                .filter_map(|v| Field::try_from(v).ok())
                .collect::<Vec<_>>()
        } else {
            Vec::default()
        };
        
        Ok(Segment { sid: sid.to_string(), sequence, fields })
    }
}


#[derive(Debug)]
pub struct Link {
    pub sid1: String,
    pub sid2: String,
    pub strand1: Strand,
    pub strand2: Strand,
    pub overlap: Option<String>,
    pub fields: Vec<Field>,
}

impl TryFrom<&str> for Link {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let mut parts = value.trim().splitn(6, '\t');

        let start = parts.next().ok_or("Empty line?")?;
        if start != "L" {
            return Err("Not a link line");
        }

        let sid1 = parts.next().ok_or("Missing first segment ID")?;
        let strand1_str = parts.next().ok_or("Missing first strand")?;
        let sid2 = parts.next().ok_or("Missing second segment ID")?;
        let strand2_str = parts.next().ok_or("Missing second strand")?;
        let overlap_str = parts.next().ok_or("Missing overlap")?;

        let strand1 = match strand1_str {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => return Err("Invalid first strand orientation")
        };

        let strand2 = match strand2_str {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => return Err("Invalid second strand orientation")
        };

        let overlap = if overlap_str == "*" {
            None
        } else {
            Some(overlap_str.to_string())
        };

        let fields = if let Some(extra_fields) = parts.next() {
            extra_fields.split('\t')
                .filter_map(|v| Field::try_from(v).ok())
                .collect::<Vec<_>>()
        } else {
            Vec::default()
        };

        Ok(Link {
            sid1: sid1.to_string(),
            sid2: sid2.to_string(),
            strand1,
            strand2,
            overlap,
            fields
        })
    }
}

pub enum GfaLine {
    Header(Header),
    Segment(Segment),
    Link(Link),
    Other(String),
}

impl TryFrom<&str> for GfaLine {
    type Error = &'static str;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let first_char = value.chars().next().ok_or("Empty line")?;

        match first_char {
            'H' => Ok(GfaLine::Header(Header::try_from(value)?)),
            'S' => Ok(GfaLine::Segment(Segment::try_from(value)?)),
            'L' => Ok(GfaLine::Link(Link::try_from(value)?)),
            '#' | 'C' | 'P' | 'J' | 'W' => Ok(GfaLine::Other(value.to_string())),
            _ => Err("Invalid GFA line type")
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_header() {
        let input = "H\tVN:Z:1.0";
        let header = Header::try_from(input).unwrap();
        assert_eq!(header.version, Some("1.0".to_string()));
        assert_eq!(header.fields.len(), 0);
    }

    #[test]
    fn test_parse_segment() {
        let input = "S\tseg1\tACGT\tLN:i:4";
        let segment = Segment::try_from(input).unwrap();
        assert_eq!(segment.sid, "seg1");
        assert_eq!(segment.sequence, Some("ACGT".to_string()));
        assert_eq!(segment.fields.len(), 1);
        
        assert_eq!(segment.fields[0].tag, "LN");
        assert_eq!(segment.fields[0].value, FieldValue::Integer(4));
    }

    #[test]
    fn test_parse_link() {
        let input = "L\tseg1\t+\tseg2\t-\t4M";
        let link = Link::try_from(input).unwrap();
        assert_eq!(link.sid1, "seg1");
        assert_eq!(link.sid2, "seg2");
        assert_eq!(link.strand1, Strand::Forward);
        assert_eq!(link.strand2, Strand::Reverse);
        assert_eq!(link.overlap, Some("4M".to_string()));
        assert_eq!(link.fields.len(), 0);
    }

    #[test]
    fn test_invalid_header() {
        let input = "X\tVN:Z:1.0";
        assert!(Header::try_from(input).is_err());
    }

    #[test]
    fn test_invalid_segment() {
        let input = "X\tseg1\tACGT";
        assert!(Segment::try_from(input).is_err());
    }

    #[test]
    fn test_invalid_link() {
        let input = "X\tseg1\t+\tseg2\t-\t4M";
        assert!(Link::try_from(input).is_err());
    }
}
