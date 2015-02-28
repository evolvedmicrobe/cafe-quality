//
//  ContextParameterProvider.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/27/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#ifndef __ConsensusCoreXCode__ContextParameterProvider__
#define __ConsensusCoreXCode__ContextParameterProvider__

#include <stdio.h>
#include "TransitionParameters.hpp"
#include <unordered_map>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
using namespace std;
namespace ConsensusCore {
    
    /**
     A data structure to store SNR values.
     */
    class SNR {
    public:
        const double A, C, G, T;
        SNR(double A, double C, double G, double T);
    };
    
    /**
     
     This class is designed to provide the relative transition probabilities for
     a given dinculeotide context at a given SNR value.
     @param <#parameter#>
     @returns <#retval#>
     @exception <#throws#>
     */
    class ContextParameterProvider {
        public:
            static std::shared_ptr<TransitionParameters> GetTransitionParameters(const string& context, const SNR& snrs);
        

       
    };
}

#endif /* defined(__ConsensusCoreXCode__ContextParameterProvider__) */
