// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: David Alexander

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <limits>
#include <vector>

#include "Matrix/SparseMatrix.hpp"

namespace ConsensusCore {
    // Performance insensitive routines are not inlined

    SparseMatrix::SparseMatrix(int rows, int cols)
        : columns_(cols),  nCols_(cols), nRows_(rows), columnBeingEdited_(-1),
          usedRanges_(cols, Interval(0, 0))
    {
        for (int j = 0; j < nCols_; j++)
        {
            columns_[j] = NULL;
        }
    }

    SparseMatrix::SparseMatrix(const SparseMatrix& other)
        : columns_(other.nCols_),
          nCols_(other.nCols_),
          nRows_(other.nRows_),
          columnBeingEdited_(other.columnBeingEdited_),
          usedRanges_(other.usedRanges_)
    {
        for (int j = 0; j < nCols_; j++)
        {
            if (other.columns_[j] != NULL)
            {
                columns_[j] = new SparseVector(*other.columns_[j]);
            }
            else
            {
                columns_[j] = NULL;
            }
        }
    }

    SparseMatrix::~SparseMatrix()
    {
        for (int j = 0; j < nCols_; j++)
        {
            if (columns_[j] != NULL) delete columns_[j];
        }
    }

    int
    SparseMatrix::UsedEntries() const
    {
        // use column ranges
        int filledEntries = 0;
        for (int col = 0; col < Columns(); ++col)
        {
            int start, end;
            boost::tie(start, end) = UsedRowRange(col);
            filledEntries += (end - start);
        }
        return filledEntries;
    }

    int
    SparseMatrix::AllocatedEntries() const
    {
        int sum = 0;
        for (int j = 0; j < nCols_; j++)
        {
            sum += (columns_[j] != NULL ?
                    columns_[j]->AllocatedEntries() : 0);
        }
        return sum;
    }

    void
    SparseMatrix::ToHostMatrix(double** mat, int* rows, int* cols) const
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        *mat = new double[Rows() * Columns()];
        *rows = Rows();
        *cols = Columns();
        for (int i = 0; i < Rows(); i++) {
            for (int j = 0; j < Columns(); j++) {
                (*mat)[i * Columns() + j] = IsAllocated(i, j) ? Get(i, j) : nan;
            }
        }
    }

    void
    SparseMatrix::CheckInvariants(int column) const
    {
        for (int j = 0; j < nCols_; j++)
         {
             if (columns_[j] != NULL) columns_[j]->CheckInvariants();
         }
    }
}
