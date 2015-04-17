
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

namespace ConsensusCore
{
   
    template<typename R>
    class MutationScorer
    {
    public:
        typedef typename R::MatrixType    MatrixType;
        typedef typename R::EvaluatorType EvaluatorType;
        typedef R                         RecursorType;

    public:
        MutationScorer(const EvaluatorType& evaluator, const R& recursor)
            throw(AlphaBetaMismatchException);

        MutationScorer(const MutationScorer& other);
        virtual ~MutationScorer();

    public:
        TemplateParameterPair Template() const;
        void Template(TemplateParameterPair tpl)
            throw(AlphaBetaMismatchException);

        double Score() const;
        double ScoreMutation(const Mutation& m, const ContextParameters& params) const;

    private:
        void DumpMatrix(const MatrixType& mat, const std::string& fname) const;

    public:
        void DumpBetaMatrix() const;
        void DumpAlphaMatrix() const;

    public:
        // Accessors that are handy for debugging.
        const MatrixType* Alpha() const;
        const MatrixType* Beta() const;
//        const PairwiseAlignment* Alignment() const;
        const EvaluatorType* Evaluator() const;
        const int NumFlipFlops() const { return numFlipFlops_; }

    private:
        EvaluatorType* evaluator_;
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
    };

    typedef MutationScorer<SimpleQvRecursor>       SimpleQvMutationScorer;
    typedef MutationScorer<SimpleQvSumProductRecursor> SimpleQvSumProductMutationScorer;
    typedef MutationScorer<SparseSimpleQvRecursor> SparseSimpleQvMutationScorer;
    typedef MutationScorer<SparseSimpleQvSumProductRecursor> SparseSimpleQvSumProductMutationScorer;
    
}
