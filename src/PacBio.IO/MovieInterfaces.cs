#region Copyright (c) 2010, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// THIS SOFTWARE CONSTITUTES AND EMBODIES PACIFIC BIOSCIENCES’ CONFIDENTIAL
// AND PROPRIETARY INFORMATION.
//
// Disclosure, redistribution and use of this software is subject to the
// terms and conditions of the applicable written agreement(s) between you
// and Pacific Biosciences, where “you” refers to you or your company or
// organization, as applicable.  Any other disclosure, redistribution or
// use is prohibited.
//
// THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#endregion

using System;
using System.Collections.Generic;
using PacBio.HDF;

namespace PacBio.IO
{
    /// <summary>
    /// An object that identifies the Uri of the data that it is presenting data from
    /// </summary>
    public interface ISourceIdentifier
    {
        /// <summary>
        /// The Uri of the data being presented by this object
        /// </summary>
        string FileName   { get; }

        /// <summary>
        /// A string describing the module that generated the data
        /// </summary>
        string SoftwareVersion { get; }

        /// <summary>
        /// Perforce changelist ID of build
        /// </summary>
        string ChangelistID { get; }

        /// <summary>
        /// Date / Time when this data source was written
        /// </summary>
        DateTime DateCreated { get; }
    }
    
    public interface IMovieMetadata
    {
        /// <summary>
        /// The movie identifier string without path or extension information (eg: m081210_101405_Uni_p1_b5)
        /// </summary>
        string MovieName { get; }
        
        /// <summary>
        /// The LIMS run code of the run that contains this movie (eg: 1500011-0004) 
        /// </summary>
        string RunCode { get; }

        /// <summary>
        /// Name of the instrument on which the movie was acquired
        /// </summary>
        string InstrumentName { get; }

        /// <summary>
        /// Id of the instrument on which the movie was acquired
        /// </summary>
        uint InstrumentId { get; }

        /// <summary>
        /// Unique platform ID: 1=Astro, 2=Springfield
        /// </summary>
        uint PlatformId { get; }

        /// <summary>
        /// Platform name, to be replaced with platform product name on release
        /// </summary>
        string PlatformName { get; }

        /// <summary>
        /// Requested camera frame rate (frames per second)
        /// </summary>
        float FrameRate { get; }
        
        /// <summary>
        /// Total number of frames acquired
        /// </summary>
        uint NumFrames { get; }

        /// <summary>
        /// Approximate frame when laser illumination was turned on
        /// </summary>
        int LaserOnFrame { get; }
        
        /// <summary>
        /// Approximate frame when interesting data begins
        /// </summary>
        int HotStartFrame { get; }

        /// <summary>
        /// Dyeset base map in channel order
        /// </summary>
        char[] BaseMap { get; }

        /// <summary>
        /// The ScanData HDF5 group used by this data source
        /// </summary>
        IGroup ScanDataGroup { get; }

        /// <summary>
        /// A string describing the enzyme kit
        /// </summary>
        string BindingKit { get; }

        /// <summary>
        /// A string describing the chemistry kit
        /// </summary>
        string SequencingKit { get; }

        /// <summary>
        /// Version of the basecaller
        /// </summary>
        string BaseCallerChangelistID { get; }

        /// <summary>
        /// Chemistry barcode triple (binding kit, sequencing kit, basecaller changelist id)
        /// </summary>
        PacBio.Data.ChemistryBarcodeTriple ChemistryBarcode { get; }

        /// <summary>
        /// Sequencing chemistry derived from the BindingKit/SequencingKit/BaseCallerChangeListID
        /// </summary>
        string SequencingChemistry { get; }
    }

    /// <summary>
    /// Represents a collection of data objects from a single movie. For example a Trace data source would have a collection
    /// of trace data objects, one for each zmw in the movie. 
    /// IDataSource objects implement IList meaning they can be indexed as if they were an array. This indexes into the list
    /// of data objects available from this source, and will not correspond to zmw hole number.  Use the ByHoleNumber() or ByXY()
    /// methods to get the data for a particular zmw (if it is available from this data source). 
    /// </summary>
    /// <typeparam name="T">The type of object that is supplied by this data source</typeparam>
    public interface IDataSource<T> :  IList<T>, IZmwRangeable<T>, ISourceIdentifier, IDisposable
    {
        /// <summary>
        /// Get the metadata about the movie associated with this data
        /// </summary>
        IMovieMetadata Movie { get; }

        /// <summary>
        /// Access a data object by it's zmw hole number
        /// </summary>
        /// <param name="holeNum"></param>
        /// <returns></returns>
        T ByHoleNumber(int holeNum);

        /// <summary>
        /// Access a data object by it zmw X/Y coordinates
        /// </summary>
        /// <param name="x">One based X position of ZMW</param>
        /// <param name="y">One based Y position of ZMW</param>
        /// <returns></returns>
        T ByXY(int x, int y);

        /// <summary>
        /// Try providing access to the underlying source indexing.
        /// </summary>
        IZmwSource ZmwSource { get; }

        /// <summary>
        /// Validate underlying data source.
        /// </summary>
        void ValidateSource();
    }

    /// <summary>
    /// A ZMW data source that can return data within a specified range.
    /// </summary>
    /// <typeparam name="T">The data object typ.</typeparam>
    public interface IZmwRangeable<T>
    {
        /// <summary>
        /// Access data objects by hole-number range.
        /// </summary>
        /// <param name="range">Range parameter. Count=-1 means everything
        /// after Start. Range is inclusive.</param>
        /// <returns>Enumerable list of data.</returns>
        IEnumerable<T> ByHoleNumberRange(IZmwRange range);
    }
    
    /// <summary>
    /// A data source for accessing basic movie metadata and zmw layout information.  The physical indexing of this
    /// data source may not match that of the other data source, so access using ByHoleNumber() is suggested.
    /// </summary>
    public interface IZmwSource : IDataSource<ISequencingZmw>, IMovieMetadata
    {
        int GetIndexByHoleNumber(int holeNumber);
        
        int GetIndexByHoleXY(int x, int y);

        ZmwIndexer ZmwIndexer { get; }
    }

    /// <summary>
    /// A data source for accessing basecalls
    /// </summary>
    public interface IBaseSource : IDataSource<IZmwBases>
    {
        bool HaveConsensusBases { get; }

        IConsensusBaseSource CreateConsensusBaseSource();

        PacBio.Data.ChemistryBarcodeTriple ChemistryBarcode { get; }

        string SequencingChemistry { get; }
    }

    /// <summary>
    /// A data source for accessing consensus basecalls
    /// </summary>
    public interface IConsensusBaseSource : IDataSource<IZmwConsensusBases>
    {
        
    }
}
