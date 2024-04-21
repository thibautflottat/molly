#![allow(non_local_definitions, dead_code)]

use std::io;
use std::num::NonZeroU64;
use std::path::PathBuf;

use molly::selection;
use numpy::ndarray::{Array, Axis, Ix2};
use numpy::{IntoPyArray, Ix1, Ix3, PyArray};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyIterator, PyList, PySlice};

type BoxVec = [[f32; 3]; 3];

#[derive(Debug, Default)]
struct FrameSelection(selection::FrameSelection);
#[derive(Debug, Default)]
struct AtomSelection(selection::AtomSelection);

impl From<FrameSelection> for selection::FrameSelection {
    fn from(val: FrameSelection) -> Self {
        val.0
    }
}

impl From<AtomSelection> for selection::AtomSelection {
    fn from(val: AtomSelection) -> Self {
        val.0
    }
}

impl FromPyObject<'_> for FrameSelection {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        if let Ok(selection) = ob.downcast::<PySlice>() {
            // TODO: This getattr business seems silly, but maybe it's necessary?
            let start = selection.getattr("start")?.extract().ok();
            let end = selection.getattr("stop")?.extract().ok();
            let step = selection.getattr("step")?.extract::<NonZeroU64>().ok();
            let range = selection::Range::new(start, end, step);
            return Ok(FrameSelection(selection::FrameSelection::Range(range)));
        }

        if let Ok(it) = ob.downcast::<PyIterator>() {
            if let Ok(indices) = it
                .iter()?
                .map(|i| i.and_then(PyAny::extract::<usize>))
                .collect::<PyResult<Vec<usize>>>()
            {
                return Ok(FrameSelection(selection::FrameSelection::FrameList(
                    indices,
                )));
            }
        }

        Err(PyTypeError::new_err("Cannot select frames with this type"))
    }
}

impl FromPyObject<'_> for AtomSelection {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        if let Ok(until) = ob.extract::<u32>() {
            return Ok(AtomSelection(selection::AtomSelection::Until(until)));
        }

        if let Ok(list) = ob.downcast::<PyList>() {
            if let Ok(bools) = list
                .iter()
                .map(PyAny::extract::<bool>)
                .collect::<PyResult<Vec<bool>>>()
            {
                return Ok(AtomSelection(selection::AtomSelection::Mask(bools)));
            }
            if let Ok(indices) = list
                .iter()
                .map(PyAny::extract::<u32>)
                .collect::<PyResult<Vec<u32>>>()
            {
                return Ok(AtomSelection(selection::AtomSelection::from_index_list(
                    &indices,
                )));
            }
        }

        Err(PyTypeError::new_err("Cannot select atoms with this type"))
    }
}

#[pyclass]
struct XTCReader {
    inner: molly::XTCReader<std::fs::File>,
    frame: Option<Frame>,
    buffered: bool,
}

#[pymethods]
impl XTCReader {
    /// Open a file as an `XTCReader`.
    #[new]
    #[pyo3(signature = (path, buffered=true))]
    fn open(path: PathBuf, buffered: bool) -> io::Result<Self> {
        let inner = molly::XTCReader::open(path)?;
        Ok(Self {
            inner,
            frame: None,
            buffered,
        })
    }

    #[getter]
    fn get_buffered(&self) -> bool {
        self.buffered
    }

    #[setter]
    fn set_buffered(&mut self, buffered: bool) -> PyResult<()> {
        self.buffered = buffered;
        Ok(())
    }

    #[getter]
    fn get_frame(&self) -> Option<Frame> {
        self.frame.clone() // FIXME: Is there a way around this?
    }

    fn determine_offsets(&mut self, until: Option<usize>) -> io::Result<Vec<u64>> {
        self.inner.determine_offsets(until).map(|l| l.to_vec())
    }

    fn determine_frame_sizes(&mut self, until: Option<usize>) -> io::Result<Vec<u64>> {
        self.inner.determine_frame_sizes(until).map(|l| l.to_vec())
    }

    /// Reset the reading head to the start of the file.
    fn home(&mut self) -> PyResult<()> {
        Ok(self.inner.home()?)
    }

    /// Read a single frame into the `frame` field of the `XTCReader`.
    fn read_frame(&mut self) -> io::Result<()> {
        if self.frame.is_none() {
            self.frame = Some(Frame::default());
        }
        let frame = &mut self.frame.as_mut().unwrap().inner;
        self.inner.read_frame(frame)
    }

    /// Read a single frame and return a copy.
    ///
    /// Calls `read_frame` internally and returns the frame immediately.
    fn pop_frame(&mut self) -> io::Result<Frame> {
        self.read_frame()?;
        Ok(self.frame.clone().unwrap())
    }

    /// Read frames according to the selections and return the frames as a list.
    ///
    /// # Note
    ///
    /// This function can perform the reads in a buffered manner, depending on the value of the
    /// `buffered` attribute.
    fn read_frames(
        &mut self,
        frame_selection: Option<FrameSelection>,
        atom_selection: Option<AtomSelection>,
    ) -> io::Result<Vec<Frame>> {
        let mut frames = Vec::new();
        let frame_selection = frame_selection.unwrap_or_default().into();
        let atom_selection = atom_selection.unwrap_or_default().into();
        match self.buffered {
            true => {
                self.inner
                    .read_frames::<true>(&mut frames, &frame_selection, &atom_selection)?
            }
            false => {
                self.inner
                    .read_frames::<false>(&mut frames, &frame_selection, &atom_selection)?
            }
        };

        Ok(frames.into_iter().map(|frame| frame.into()).collect())
    }

    /// Read all frames into the provided `numpy.ndarray`.
    ///
    /// The `coordinate_array` must have a shape of `(nframes, natoms, 3)`.
    ///
    /// The `boxvec_array` must have a shape of `(nframes, 3, 3)`.
    ///
    /// Returns `True` if the reading operation was successful.
    ///
    /// # Note
    ///
    /// This function can perform the reads in a buffered manner, depending on the value of the
    /// `buffered` attribute.
    fn read_into_array(
        &mut self,
        py: Python<'_>,
        coordinate_array: &PyArray<f32, Ix3>,
        boxvec_array: &PyArray<f32, Ix3>,
        time_array: Option<&PyArray<f32, Ix1>>,
        frame_selection: Option<FrameSelection>,
        atom_selection: Option<AtomSelection>,
    ) -> PyResult<bool> {
        {
            // Verify that the shapes of the arrays are correct.
            let &[nf_coords, na, d] = coordinate_array.shape() else {
                return Err(PyValueError::new_err(
                    ("coordinate array should have 3 dimensions").to_string(),
                ));
            };
            if d != 3 {
                return Err(PyValueError::new_err(format!(
                    "incorrect shape: coordinate array must be of shape (nframes, natoms, 3), found {:?}",
                    (nf_coords, na, d)
                )));
            }
            let &[nf_boxvecs, a, b] = boxvec_array.shape() else {
                return Err(PyValueError::new_err(
                    "boxvec array should have 3 dimensions".to_string(),
                ));
            };
            if d != 3 {
                return Err(PyValueError::new_err(format!(
                    "incorrect shape: boxvec array must be of shape (nframes, 3, 3), found {:?}",
                    (nf_boxvecs, a, b)
                )));
            }
            if nf_coords != nf_boxvecs {
                return Err(PyValueError::new_err(format!(
                    "the number of frames defined in the coordinate and boxvec arrays do not match, {:?} != {:?}",
                    (nf_coords, na, d),
                    (nf_boxvecs, a, b)
                )));
            }
            if let Some(time_array) = time_array {
                let &[nf_times] = time_array.shape() else {
                    return Err(PyValueError::new_err(
                        "time array should have 1 dimension".to_string(),
                    ));
                };
                if nf_times != nf_coords {
                    return Err(PyValueError::new_err(format!(
                        "the number of frames defined in the coordinate and boxvec arrays does not match the number of frames defined in the time array, {:?} != {:?}",
                        (nf_times,),
                        (nf_boxvecs, a, b)
                    )));
                }
            }
        }

        let mut coordinates = coordinate_array.readwrite();
        let mut coordinates = coordinates.as_array_mut();
        let mut boxvecs = boxvec_array.readwrite();
        let mut boxvecs = boxvecs.as_array_mut();
        let mut times = time_array.map(|ts| ts.readwrite());
        let mut times = times.as_mut().map(|ts| ts.as_array_mut());

        let atom_selection: selection::AtomSelection = atom_selection.unwrap_or_default().into();
        let mut frame = molly::Frame::default();
        let until = frame_selection
            .as_ref()
            .and_then(|FrameSelection(selection)| selection.until());
        let offsets = self.inner.determine_offsets(until)?;
        let offsets = offsets.iter().enumerate().filter_map(|(idx, offset)| {
            if let Some(FrameSelection(selection)) = &frame_selection {
                match selection.is_included(idx) {
                    Some(true) => Some(offset),
                    Some(false) => None,
                    None => None,
                }
            } else {
                Some(offset)
            }
        });
        // TODO: Fix up this mess of zips.
        for (i, ((mut array_coordinates, mut array_boxvecs), &offset)) in coordinates
            .axis_iter_mut(Axis(0))
            .zip(boxvecs.axis_iter_mut(Axis(0)))
            .zip(offsets)
            .enumerate()
        {
            py.check_signals()?;
            match self.buffered {
                true => {
                    self.inner
                        .read_frame_at_offset::<true>(&mut frame, offset, &atom_selection)?;
                }
                false => {
                    self.inner.read_frame_at_offset::<false>(
                        &mut frame,
                        offset,
                        &atom_selection,
                    )?;
                }
            };
            // TODO: Check whether the two unwraps here can just be elided somehow.
            array_coordinates
                .rows_mut()
                .into_iter()
                .zip(frame.coords())
                .for_each(|(mut array_coord, frame_coord)| {
                    // Unwrap should be fine here, since we checked the sizes before.
                    frame_coord.write_to_slice(array_coord.as_slice_mut().unwrap())
                });
            array_boxvecs
                .columns_mut()
                .into_iter()
                .zip(frame.boxvec.to_cols_array_2d())
                .for_each(|(array_boxvec, frame_boxvec)| {
                    for (&frame_value, array_value) in frame_boxvec.iter().zip(array_boxvec) {
                        *array_value = frame_value
                    }
                });
            if let Some(ref mut times) = times {
                times[i] = frame.time;
            }
        }

        Ok(true)
    }
}

/// A single trajectory frame.
///
/// All distances are given in nanometers.
#[pyclass]
#[derive(Default, Clone)]
struct Frame {
    inner: molly::Frame,
}

impl From<molly::Frame> for Frame {
    fn from(frame: molly::Frame) -> Self {
        Self { inner: frame }
    }
}

#[pymethods]
impl Frame {
    #[getter]
    fn get_step(&self) -> u32 {
        self.inner.step
    }

    #[getter]
    fn get_time(&self) -> f32 {
        self.inner.time
    }

    /// The box vectors of this frame as an array of columns of a 3×3 matrix.
    #[getter]
    fn get_box(&self) -> BoxVec {
        self.inner.boxvec.to_cols_array_2d()
    }

    #[getter]
    fn get_precision(&self) -> f32 {
        self.inner.precision
    }

    /// Get the positions as an `np.ndarray`.
    ///
    /// Distances in nanometers.
    #[getter]
    fn get_positions<'py>(&self, py: Python<'py>) -> &'py PyArray<f32, Ix2> {
        // Fingers crossed we don't make any position lists of which the length isn't a multiple of
        // 3... This is a guarantee from the implementation, so It's Fine(tm).
        let positions = &self.inner.positions;
        let natoms = positions.len() / 3;
        Array::from_shape_vec((natoms, 3), positions.clone())
            .unwrap()
            .into_pyarray(py)
    }
}

/// Read xtc files, fast.
///
/// Marieke Westendorp, 2024.
#[pymodule]
fn _molly(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(function_name, m)?)?;
    m.add_class::<XTCReader>()?;

    Ok(())
}
