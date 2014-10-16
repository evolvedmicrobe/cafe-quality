cafe-quality
============

Tools to assess and diagnose accuracy issues in PacBio data.

This is a library designed to help analyze PacBio data and assess quality issues.  It is composed of the following items.

Directory Structure
===================

    /lib - header files and dlls needed for analysis.
    /src - source files
        /bio - .NET Bio fork for basic bioinformatics tasks (C#).
        /ConsensusCore - Fork of PacBio consensus core, for compute intense tasks (C++).
        /VariantCaller - Code for loading ZMW reads, calling variants and creating quality metrics (C#).
        /scripts - Scripts to use the library to examine quality (F#)
