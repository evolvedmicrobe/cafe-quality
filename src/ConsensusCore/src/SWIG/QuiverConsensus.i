%{
/* Includes the header in the wrapper code */
#include "Sequence.hpp"
#include "Mutation.hpp"
#include "Read.hpp"
#include "ContextParameterProvider.hpp"
#include "ContextParameters.hpp"
#include "TransitionParameters.hpp"
#include "TemplateParameterPair.hpp"
#include "Quiver/QuiverConfig.hpp"
#include "Quiver/MultiReadMutationScorer.hpp"
#include "Quiver/MutationScorer.hpp"
#include "Quiver/SimpleRecursor.hpp"
#include "Quiver/ReadScorer.hpp"
#include "Quiver/Diploid.hpp"

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
%include "Quiver/detail/Combiner.hpp"
%include "Quiver/QuiverConfig.hpp"
%include "TemplateParameterPair.hpp"
%include "Quiver/MultiReadMutationScorer.hpp"
%include "Quiver/MutationScorer.hpp"
%include "Quiver/SimpleRecursor.hpp"
%include "Quiver/ReadScorer.hpp"
%include "Quiver/Diploid.hpp"


namespace ConsensusCore {


    //
    // Dense matrix recursors and such
    //
    %template(SimpleQvRecursor)         SimpleRecursor<DenseMatrix, detail::ViterbiCombiner>;
    %template(SimpleQvMutationScorer)   MutationScorer<SimpleQvRecursor>;
    
    %template(SimpleQvSumProductRecursor)         SimpleRecursor<DenseMatrix, detail::SumProductCombiner>;
    %template(SimpleQvSumProductMutationScorer)   MutationScorer<SimpleQvSumProductRecursor>;

    //
    // Sparse matrix support
    //
    %template(SparseSimpleQvRecursor)         SimpleRecursor<SparseMatrix, detail::ViterbiCombiner>;

    %template(SparseSimpleQvMutationScorer)   MutationScorer<SparseSimpleQvRecursor>;
    %template(SparseQvViterbiMultiReadMutationScorer) MultiReadMutationScorer<SparseSimpleQvRecursor>;
    
    //
    // Sparse matrix sum-product support
    //
    %template(SparseSimpleQvSumProductRecursor)         SimpleRecursor<SparseMatrix, detail::SumProductCombiner>;
    %template(SparseSimpleQvSumProductMutationScorer)   MutationScorer<SparseSimpleQvSumProductRecursor>;


    %template(SparseQvSumProductMultiReadMutationScorer) MultiReadMutationScorer<SparseSimpleQvSumProductRecursor>;


}
