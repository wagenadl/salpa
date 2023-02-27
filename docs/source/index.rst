.. SALPA documentation master file, created by
   sphinx-quickstart on Sun Feb 26 13:42:21 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SALPA
=====

Introduction
-----------------

This document explains how to use the SALPA algorithm from the command line,
or from within Python or Octave/Matlab. The algorithm itself is explained in
Wagenaar and Potter, Real-time multi-channel stimulus artifact suppression by local curve fitting. *J. Neurosci. Methods* **120**, 113–120, available from the
`Wagenaar Lab website <https://www.danielwagenaar.net/pubs.html#02-WP>`_.

Installation
------------

Prerequisites:

- `SCons <https://scons.org>`_
- A C++ compiler
- Python (version 3.8 or later, optional)
- Octave or Matlab (optional)

Download the source from github:

    git clone https://github.com/wagenadl/salpa.git

enter the download folder:

    cd salpa

and build using Scons:

    scons

This creates a single binary "build/salpa", which may be copied to any
convenient location on your $PATH.

Command line usage
------------------

The SALPA program reads and writes raw files comprising 16-bit signed
integers, with arbitrary numbers of channels.

The program understands the following arguments:

- **F** *f*

  The sampling rate is *f*, measured in Hertz. This argument must
  precede any of the arguments that are measured in
  milliseconds. (Default: 30,000 Hz.)

- **c** *n*

  The number of electrode channels in the file is *n*. (Default: equal
  to the total number of channels (see **C**), or 64 if neither is
  explicitly given.)

- **C** *m*

  The total number of channels in the file is *m*. As an example,
  MultiChannel Systems’ binary files contain 60 electrode channels and
  4 auxiliary channels, for a total of 64 channel. Similar
  arrangements are found in some OpenEphys recordings. (Default: equal
  to the number of probe channels (see **c**).)

  Auxiliary channels are copied verbatim from input to output, without
  any processing.

- **t** *d*

  The threshold on the residual error for resuming normal operation
  after detecting an artifact is *d* digital units. (No default; see
  **x**.)

- **x** *x*

  The threshold on the residual error for resuming normal operation
  after detecting an artifact is *x* times the RMS noise of each
  respective channel. (Default: 3.)

  Note that **t** and **x** are mutually exclusive
    
- **l** *τ*

  (That's lowercase ℓ, not upper case 𝐼.) The half-width of the
  template for artifact-shape fitting is *τ*, measured in
  milliseconds. (Default: 3 ms.)

- **a** *t*

  The window over which residual error is calculated is *t*, measured
  in milliseconds. (Can usually be left to its default value, 0.2 ms.)

- **b** *t*

  A segment of *t* milliseconds is blanked before resuming normal
  operation. *t* must be less than the *τ* value given with
  **l**. (Can usually be left to its default value, 0.4 ms.)

- **Z**

  If **Z** is given, “blanking” operations requested by **b** are
  terminated early if the signal crosses zero.

- **r** *a*,*b*

  When the digital value of any channel dips is less than or equal to
  *a*, or greater than or equal to *b*, treat the following segment as
  artifact. (Default: *a* = -32767, *b* = 32767. See also **A**.)
  
- **A** *t*

  Look ahead *t* milliseconds for automatic artifact
  detection. (Default: 0.2 ms. See also **r**.)

- **f** *t*

  Enforce a minimum duration of *t* milliseconds for automatically
  detected artifacts.

- **P** *filename*

  Artifacts are assumed to exist across all channels at the times
  specified in the named file. Each line of the file must contain two
  numbers: the start of a presumed artifact (in samples) and the
  duration of that artifact (also in samples). The lines must appear
  in temporal order.

- **T** *n*

  Use *n* CPU threads for processing. (Default: 8.)

- **S** *n*

  Use a buffer size of *n* scans. (Default: 4096, internally rounded
  down to a power of two.)

- **M** *n*

  Skip first *n* scans from the beginning of the file. (Default: do
  not skip.)

  **Note:** The numbers in the file specified with **P** are measured
  after the skip. That is, if the number *k* appears in the file, the
  artifact is presumed to exist *n* + *k* from the beginning of the
  file.

- **N** *n*

  Only process the first *n* scans of the file. (Not including any
  skipped scans specified with **M**; Default: process entire file.)

- **B**

  If given, subtract baseline before processing. Values specified with
  **r** are relative to subtracted baseline.

- **i** *filename*

  Read input from the named file. (Default: read from *stdin*.)

- **o** *filename*

  Send output to the named file. (Default: write to *stdout*.)

Usage example
^^^^^^^^^^^^^

To process a Neuropixels file recorded by OpenEphys:

  salpa -F 30000 -l 3 -c 384 -r-20000,20000 -i continuous.dat -o clean.dat


Python usage
------------

A Python module named “salpa” is provided to wrap around the binary.

This module defines two functions: :ref:`salpaparams` and :ref:`salparun`. 
 
 
.. toctree::
   :maxdepth: 1
   :caption: TOC

   python
   
