use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};

use crate::reader::{read_opaque, read_u32};

pub type UnBuffered<'s> = &'s [u8];

pub trait Buffered<'s, 'r, R>: Sized {
    const MIN_BUFFERED_SIZE: usize = 0x500000;
    const BLOCK_SIZE: usize = 0x20000;

    // TODO(buffered): Consider giving the n_bytes from the outside?
    /// Create a new [`Buffer`] reader.
    ///
    /// # Note
    ///
    /// Expects that the first `u32` represents the number of upcoming bytes in the compressed data
    /// stream. If this function is called on a reader that is not at that spot in its stream, the
    /// resulting [`Buffer`] is invalid. This is the same requirement [`read_opaque`] has.
    // We initialize on a Vec<u8> but after preparing this Vec we store the allocation internally
    // as a mutable byte slice, since we do not need to do any Vec-specific operations on it
    // afterwards. When this type is dropped, the ownership of `scratch` is returned since the
    // reference to it dissolves.
    fn new(scratch: &'s mut Vec<u8>, reader: &'r mut R) -> io::Result<Self>;

    /// Get a byte at some index.
    ///
    /// # Panics
    ///
    /// If the `index` exceeds the number of bytes this [`Buffer`] can read, this function panics.
    /// In that case something is seriously wrong anyway.
    fn fetch(&mut self, index: usize) -> u8;

    /// Finish will eat your reader, leaving it at the start of the next frame, and then drops it.
    ///
    /// For an implementation that relies on [`std::io::Seek`] ([`Buffer`] in our case), this
    /// really matters.
    fn finish(self) -> io::Result<()>;
}

/// A specialized buffered reader for the compressed datastream.
pub(crate) struct Buffer<'s, 'r> {
    /// Internal scratch buffer to read into.
    ///
    /// # Warning
    ///
    /// Accessing bytes from this buffer directly is valid iff the index of that byte < `self.idx`.
    scratch: &'s mut [u8],
    /// Points to the next unread/unfilled byte in `scratch`.
    ///
    /// The starting point for reading bytes from `reader` into `scratch`.
    idx: usize,
    reader: &'r mut File,
    // TODO(buffered): Add some notion of a 'rich' heuristic. For instance, if we know there are
    // 1000 atoms, and we only want to read up until the 500th atom, we can pretty safely assume
    // that we can just read (500/1000) * 1.1 * nbytes = 0.55 * nbytes and be fine.
}

impl Buffer<'_, '_> {
    /// Returns the size of this [`Buffer`].
    const fn size(&self) -> usize {
        self.scratch.len()
    }

    // TODO: I have inlined this function wherever possible, so it is basically useless at this
    // point, beyond a debug assert.
    /// Returns a reference to the valid section of this [`Buffer`].
    #[inline(always)]
    fn valid(&self) -> &[u8] {
        &self.scratch[..self.idx]
    }

    /// Returns the number of bytes that are yet to be read by this [`Buffer`].
    const fn left(&self) -> usize {
        self.size() - self.idx
    }

    /// Read enough bytes such that `index` points to a valid byte.
    fn read_to_include(&mut self, index: usize) -> io::Result<()> {
        while self.idx <= index {
            // TODO(buffered): Consider dealing with n_bytes == 0 indicating eof.
            // Read a bunch of bytes limited by the size of the scratch buffer and BLOCK_SIZE.
            // We would rather do a couple more smaller reads (BLOCK_SIZE) than one big one that
            // goes way beyond what we need according to some AtomSelection.
            let until = usize::min(self.size(), index + Self::BLOCK_SIZE);
            self.idx += self.reader.read(&mut self.scratch[self.idx..until])?;
        }
        assert!(
            index < self.idx,
            "index ({index}) must be within than the defined valid range (..{valid})",
            valid = self.idx
        );
        Ok(())
    }
}

impl<'s, 'r> Buffered<'s, 'r, File> for Buffer<'s, 'r> {
    fn new(scratch: &'s mut Vec<u8>, reader: &'r mut File) -> io::Result<Self> {
        let count = read_u32(reader)? as usize;
        let padding = (4 - (count % 4)) % 4; // FIXME: Why, and also, can we do this better?

        // Fill the scratch buffer with a cautionary value.
        scratch.resize(count + padding, 0xff); // FIXME: Is MaybeUninit a good idea here?

        let mut buffer = Self {
            scratch,
            idx: 0,
            reader,
        };

        // In case the buffer size is rather low, it is probably most efficient to just read it all
        // at once, right here.
        if buffer.scratch.len() <= Self::MIN_BUFFERED_SIZE {
            buffer.read_to_include(count.saturating_sub(1))?;
        }

        Ok(buffer)
    }

    #[inline(always)]
    fn fetch(&mut self, index: usize) -> u8 {
        debug_assert!(
            index < self.size(),
            "index ({index}) must be within the defined range of the scratch buffer (..{size})",
            size = self.size()
        );

        // If we're out of bytes, we'll have to read new ones.
        // NOTE: This branch is pretty much singularly responsible for the performance difference
        // between unbuffered and buffered decompression (cf. the impl of this function for
        // `UnBuffered`).
        if index >= self.idx {
            // FIXME(buffered): For now, let's just fuck this up with a terrible unwrap here.
            // Gotta change this to be io::Result at some point? If we can muster the perf hit
            // at least...
            self.read_to_include(index).unwrap();
        }

        self.scratch[index]
    }

    fn finish(self) -> io::Result<()> {
        self.reader.seek(SeekFrom::Current(self.left() as i64))?;
        drop(self);
        Ok(())
    }
}

/// A fallback non-buffered implementation in case [`std::io::Seek`] is not available for `R`.
impl<'s, 'r, R: Read> Buffered<'s, 'r, R> for UnBuffered<'s> {
    fn new(scratch: &'s mut Vec<u8>, reader: &'r mut R) -> io::Result<Self> {
        read_opaque(reader, scratch)?;
        Ok(scratch)
    }

    fn fetch(&mut self, index: usize) -> u8 {
        debug_assert!(
            index < self.len(),
            "index ({index}) must be within the defined range of the scratch buffer (..{size})",
            size = self.len()
        );
        self[index]
    }

    fn finish(self) -> io::Result<()> {
        Ok(()) // Nothing to do, since we already read everything.
    }
}
