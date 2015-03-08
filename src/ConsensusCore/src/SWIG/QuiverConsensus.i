%{
/* Includes the header in the wrapper code */
#include "Sequence.hpp"
#include "Mutation.hpp"
#include "Read.hpp"
#include "ContextParameterProvider.hpp"
#include "ContextParameters.hpp"
#include "TemplateParameterPair.hpp"
#include "TransitionParameters.hpp"
#include "Quiver/QuiverConfig.hpp"
#include "Quiver/QvEvaluator.hpp"
#include "Quiver/MultiReadMutationScorer.hpp"
#include "Quiver/MutationScorer.hpp"
#include "Quiver/SimpleRecursor.hpp"
#include "Quiver/ReadScorer.hpp"
#include "Quiver/Diploid.hpp"
#include "Quiver/QuiverConsensus.hpp"


using namespace ConsensusCore;
using namespace std;
%}

#ifdef SWIGPYTHON

%include "numpy.i"
%numpy_typemaps(float, NPY_FLOAT, int)

%apply (float* IN_ARRAY2, int DIM1, int DIM2)
       { (const float *siteScores, int dim1, int dim2) }

#endif // SWIGPYTHON

 // SWIG now seems to be incorrectly deciding that MultiReadMutationScorer
 // is an abstract class, so we have to tell it otherwise
%feature("notabstract") MultiReadMutationScorer;

#ifdef SWIGCSHARP
%csmethodmodifiers *::ToString() const "public override"
#endif // SWIGCSHARP



%include "Sequence.hpp"
%include "Mutation.hpp"
%include "Read.hpp"
%include "TransitionParameters.hpp" 
%include "ContextParameterProvider.hpp"
%include "ContextParameters.hpp"
%include "TemplateParameterPair.hpp"
%include "Quiver/detail/Combiner.hpp"
%include "Quiver/detail/RecursorBase.hpp"
%include "Quiver/QuiverConfig.hpp"
%include "Quiver/QvEvaluator.hpp"
%include "Quiver/MultiReadMutationScorer.hpp"
%include "Quiver/MutationScorer.hpp"
%include "Quiver/SimpleRecursor.hpp"
%include "Quiver/ReadScorer.hpp"
%include "Quiver/Diploid.hpp"
%include "Quiver/QuiverConsensus.hpp"



namespace ConsensusCore {


    //
    // Dense matrix recursors and such
    //
    %template(QvRecursorBase)           detail::RecursorBase<DenseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SimpleQvRecursor)         SimpleRecursor<DenseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SimpleQvMutationScorer)   MutationScorer<SimpleQvRecursor>;
    
    %template(QvRecursorSumProductBase)           detail::RecursorBase<DenseMatrix, QvEvaluator, detail::SumProductCombiner>;
    %template(SimpleQvSumProductRecursor)         SimpleRecursor<DenseMatrix, QvEvaluator, detail::SumProductCombiner>;
    %template(SimpleQvSumProductMutationScorer)   MutationScorer<SimpleQvSumProductRecursor>;

    //
    // Sparse matrix support
    //
    %template(SparseQvRecursorBase)           detail::RecursorBase<SparseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SparseSimpleQvRecursor)         SimpleRecursor<SparseMatrix, QvEvaluator, detail::ViterbiCombiner>;

    %template(SparseSimpleQvMutationScorer)   MutationScorer<SparseSimpleQvRecursor>;
    %template(SparseQvViterbiMultiReadMutationScorer) MultiReadMutationScorer<SparseSimpleQvRecursor>;
    
    //
    // Sparse matrix sum-product support
    //
    %template(SparseQvSumProductRecursorBase)           detail::RecursorBase<SparseMatrix, QvEvaluator, detail::SumProductCombiner>;
    %template(SparseSimpleQvSumProductRecursor)         SimpleRecursor<SparseMatrix, QvEvaluator, detail::SumProductCombiner>;
    %template(SparseSimpleQvSumProductMutationScorer)   MutationScorer<SparseSimpleQvSumProductRecursor>;


    %template(SparseQvSumProductMultiReadMutationScorer) MultiReadMutationScorer<SparseSimpleQvSumProductRecursor>;


}
