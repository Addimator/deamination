use getset::{Getters, Setters};
use std::collections::HashMap;

#[derive(Debug, Getters, Setters, Clone, Hash, PartialEq, Eq)]

pub struct Position {
    #[getset(get = "pub")]
    chrom: String,
    #[getset(get = "pub")]
    pos: u32,
}

impl Position {
    // Benutzerdefinierte `new` Methode, um ein Objekt zu erstellen
    pub fn new(chrom: String, pos: u32) -> Self {
        Position { chrom, pos }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]

pub enum Direction {
    Forward,
    Reverse,
}

impl Direction {
    // Beispielmethode, die die Richtung als String zurückgibt
    pub fn as_str(&self) -> &str {
        match self {
            Direction::Forward => "forward",
            Direction::Reverse => "reverse",
        }
    }
}

pub type NumberOfNucleotides = HashMap<char, usize>;

#[derive(Debug, Setters, Getters, Clone, Eq, PartialEq)]
pub struct MethPos {
    #[getset(get = "pub")]
    position: Position,
    #[getset(get = "pub", set = "pub")]
    methylation: i64,
    #[getset(get = "pub", set = "pub")]
    meth_bases: HashMap<(Direction, u32), NumberOfNucleotides>,
}

impl MethPos {
    // Benutzerdefinierte `new` Methode, um ein Objekt zu erstellen
    pub fn new(position: Position, methylation: i64) -> Self {
        MethPos {
            position,
            methylation,
            meth_bases: HashMap::new(),
        }
    }

    pub fn meth_bases_mut(&mut self) -> &mut HashMap<(Direction, u32), NumberOfNucleotides> {
        &mut self.meth_bases
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct PosType {
    direction: Direction,
    methylation: bool,
    cpg_position: u32,
}

impl PosType {
    // Benutzerdefinierte `new` Methode, um ein Objekt zu erstellen
    pub fn new(direction: Direction, methylation: bool, cpg_position: u32) -> Self {
        PosType {
            direction,
            methylation, // Wenn methylation > 20, dann true, sonst false
            cpg_position,
        }
    }
}
