//
//  ContextParameters.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/1/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#ifndef __ConsensusCoreXCode__ContextParameters__
#define __ConsensusCoreXCode__ContextParameters__

#include <stdio.h>
#include "ContextParameterProvider.hpp"
#include <unordered_map>
#include <string>
#include <vector>

using namespace std;
namespace ConsensusCore {
    /**
     This class represents a collection of context parameters for a given set of
     SNR values.  It stores the transition parameters for the 
     @param <#parameter#>
     @returns <#retval#>
     @exception <#throws#>
     */
    class ContextParameters {
    public:
        ContextParameters(SNR snrs);
        // Empty constructor for swig
        ContextParameters();
        // Copy constructor for swig
        ContextParameters(const ContextParameters& arg);
        TransitionParameters& GetParametersForContext(char bp1, char bp2);
    private:
        unordered_map<string, TransitionParameters&> param_map;
        vector<string> contexts {"AA", "NA", "CC", "NC", "TT", "NT", "GG", "NG"};
    };
};

#endif /* defined(__ConsensusCoreXCode__ContextParameters__) */
