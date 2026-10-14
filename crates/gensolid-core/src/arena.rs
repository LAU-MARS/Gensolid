//! Generational arenas with typed identifiers.

use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

use crate::error::{Error, ErrorCode, Result};

/// Typed, generational handle to a value stored in an [`Arena`].
///
/// A handle stays valid until its value is removed. After removal the slot may be
/// reused, but with a new generation, so a stale handle never aliases a new value.
pub struct Id<T> {
    index: u32,
    generation: u32,
    marker: PhantomData<fn() -> T>,
}

impl<T> Id<T> {
    const fn new(index: u32, generation: u32) -> Self {
        Self {
            index,
            generation,
            marker: PhantomData,
        }
    }

    /// Slot index, stable for the lifetime of the value.
    pub const fn index(self) -> u32 {
        self.index
    }

    /// Generation of the slot when the value was inserted.
    pub const fn generation(self) -> u32 {
        self.generation
    }
}

impl<T> Clone for Id<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Copy for Id<T> {}

impl<T> PartialEq for Id<T> {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index && self.generation == other.generation
    }
}

impl<T> Eq for Id<T> {}

impl<T> Hash for Id<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
        self.generation.hash(state);
    }
}

impl<T> PartialOrd for Id<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Ord for Id<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.index, self.generation).cmp(&(other.index, other.generation))
    }
}

impl<T> fmt::Debug for Id<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Id({}v{})", self.index, self.generation)
    }
}

#[derive(Clone)]
struct Slot<T> {
    generation: u32,
    value: Option<T>,
}

/// Storage for values addressed by [`Id`].
#[derive(Clone)]
pub struct Arena<T> {
    slots: Vec<Slot<T>>,
    free: Vec<u32>,
    len: usize,
}

impl<T> Default for Arena<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Arena<T> {
    /// Creates an empty arena.
    pub const fn new() -> Self {
        Self {
            slots: Vec::new(),
            free: Vec::new(),
            len: 0,
        }
    }

    /// Number of live values.
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Whether the arena holds no live values.
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Stores `value` and returns its handle.
    pub fn insert(&mut self, value: T) -> Result<Id<T>> {
        if let Some(index) = self.free.pop() {
            let Some(slot) = self.slots.get_mut(index as usize) else {
                return Err(Error::new(
                    ErrorCode::Internal,
                    "arena free list is corrupt",
                ));
            };
            slot.value = Some(value);
            self.len += 1;
            return Ok(Id::new(index, slot.generation));
        }
        let index = u32::try_from(self.slots.len())
            .map_err(|_| Error::new(ErrorCode::CapacityExceeded, "arena is full"))?;
        self.slots.push(Slot {
            generation: 0,
            value: Some(value),
        });
        self.len += 1;
        Ok(Id::new(index, 0))
    }

    /// Returns the value for `id`, or `None` if it was removed or never existed.
    pub fn get(&self, id: Id<T>) -> Option<&T> {
        self.slots
            .get(id.index as usize)
            .filter(|slot| slot.generation == id.generation)
            .and_then(|slot| slot.value.as_ref())
    }

    /// Mutable variant of [`Arena::get`].
    pub fn get_mut(&mut self, id: Id<T>) -> Option<&mut T> {
        self.slots
            .get_mut(id.index as usize)
            .filter(|slot| slot.generation == id.generation)
            .and_then(|slot| slot.value.as_mut())
    }

    /// Whether `id` refers to a live value.
    pub fn contains(&self, id: Id<T>) -> bool {
        self.get(id).is_some()
    }

    /// Removes and returns the value for `id`.
    ///
    /// Generations wrap after 2^32 reuses of one slot.
    pub fn remove(&mut self, id: Id<T>) -> Option<T> {
        let slot = self.slots.get_mut(id.index as usize)?;
        if slot.generation != id.generation {
            return None;
        }
        let value = slot.value.take()?;
        slot.generation = slot.generation.wrapping_add(1);
        self.free.push(id.index);
        self.len -= 1;
        Some(value)
    }

    /// Iterates over live values in slot order.
    pub fn iter(&self) -> impl Iterator<Item = (Id<T>, &T)> {
        self.slots.iter().enumerate().filter_map(|(index, slot)| {
            let index = u32::try_from(index).ok()?;
            slot.value
                .as_ref()
                .map(|value| (Id::new(index, slot.generation), value))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stale_ids_do_not_alias_reused_slots() {
        let mut arena = Arena::new();
        let a = arena.insert("a").unwrap();
        assert_eq!(arena.remove(a), Some("a"));
        let b = arena.insert("b").unwrap();
        assert_eq!(a.index(), b.index());
        assert_ne!(a, b);
        assert_eq!(arena.get(a), None);
        assert_eq!(arena.get(b), Some(&"b"));
        assert_eq!(arena.len(), 1);
    }

    #[test]
    fn iter_skips_removed_values() {
        let mut arena = Arena::new();
        let ids: Vec<_> = (0..4).map(|i| arena.insert(i).unwrap()).collect();
        arena.remove(ids[1]);
        let values: Vec<_> = arena.iter().map(|(_, v)| *v).collect();
        assert_eq!(values, vec![0, 2, 3]);
    }
}
