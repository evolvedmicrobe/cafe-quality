//
//  TransitionParameters.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/23/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#ifndef __ConsensusCoreXCode__TransitionParameters__
#define __ConsensusCoreXCode__TransitionParameters__

#include <stdio.h>


namespace ConsensusCore {
    
    class TransitionParameters {
    public:
        // NATURAL (NOT LOGGED) SCALE transition parameters
        double Match, Stick, Branch, Deletion;
        TransitionParameters(double match, double stick, double branch, double deletion);
        double CalculateTotal() const;
        void RemoveConstant(double value);
        
        // Define copy and default constructors for SWIG
        TransitionParameters();
        TransitionParameters(const TransitionParameters& other );
        
        TransitionParameters(TransitionParameters&& src) = default;
        TransitionParameters&  operator=(TransitionParameters&& rhs) = default;
        TransitionParameters& operator=(const TransitionParameters& rhs) = default;
        ~TransitionParameters() = default;
        
        
    };
}
#endif /* defined(__ConsensusCoreXCode__TransitionParameters__) */
