
#pragma once

#include <boost/noncopyable.hpp>
#include <string>

// TODO(dalexander): how can we remove this include??
//  We should move all template instantiations out to another
//  header, I presume.
#include "Quiver/SimpleRecursor.hpp"
#include "Types.hpp"
#include "Mutation.hpp"
#include "ContextParameters.hpp"
//#include "Quiver/MultiReadMutationScorer.hpp"

namespace ConsensusCore
{
   
    template<typename R>
    class MutationScorer
    {
        public:
            typedef typename R::MatrixType    MatrixType;
            typedef R                         RecursorType;

        public:
            MutationScorer(const R& recursor)
                throw(AlphaBetaMismatchException);

            MutationScorer(const MutationScorer& other);
            virtual ~MutationScorer();

        public:
            WrappedTemplateParameterPair Template() const;
            void Template(WrappedTemplateParameterPair tpl)
                throw(AlphaBetaMismatchException);

            double Score() const;
            // Must be called from a multiread mutation context.
            double ScoreMutation(const Mutation& m) const;
            /* Every emission move is multiplied by a constant factor to weight moves
             in that direction.  To correct the likelihood, we must add the log of this 
             scaling factor times the number of emitted bases */
        double MatchScalingFactorCorrection() const;
        
        private:
            void DumpMatrix(const MatrixType& mat, const std::string& fname) const;

        public:
            void DumpBetaMatrix() const;
            void DumpAlphaMatrix() const;

        public:
            // Accessors that are handy for debugging.
            const MatrixType* Alpha() const;
            const MatrixType* Beta() const;
            const int NumFlipFlops() const { return numFlipFlops_; }

        private:
             R* recursor_;
            /**
             The Alpha matrix, holding the forward variables
             */
            MatrixType* alpha_;
            MatrixType* beta_;
            /**
             The extension matrix, a clean buffer to work with data where we need it.
             */
            MatrixType* extendBuffer_;
            int numFlipFlops_;
            double cachedMatchScalingCorrectionFactor;
        
            //friend class MultiReadMutationScorer<R>;
    };
    
    typedef MutationScorer<SimpleQvRecursor>       SimpleQvMutationScorer;
    typedef MutationScorer<SimpleQvSumProductRecursor> SimpleQvSumProductMutationScorer;
    typedef MutationScorer<SparseSimpleQvRecursor> SparseSimpleQvMutationScorer;
    typedef MutationScorer<SparseSimpleQvSumProductRecursor> SparseSimpleQvSumProductMutationScorer;
    
}
