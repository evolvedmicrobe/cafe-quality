//
//  ContextParameters.cpp
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/1/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#include "ContextParameters.hpp"

using namespace std;
namespace ConsensusCore {
    ContextParameters::ContextParameters(SNR snr)
    {
        for(string ctx : contexts) {
            auto p = ContextParameterProvider::GetTransitionParameters(ctx, snr);
            param_map[ctx] = p;
        }
    }
    
    TransitionParameters&
    ContextParameters::GetParametersForContext(char bp1, char bp2) {
        string s;
        if( bp1 == bp2)
        {
            s +=bp1;
            s += bp2;
        }
        else {
            s += 'N';
            s += bp2;
        }
        return param_map[s];
    }
}