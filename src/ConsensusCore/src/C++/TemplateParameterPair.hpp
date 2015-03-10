//
//  TemplateParameterPair.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/3/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#ifndef __ConsensusCoreXCode__TemplateParameterPair__
#define __ConsensusCoreXCode__TemplateParameterPair__

#include <stdio.h>
#include <string>
#include <vector>
#include <memory>

#include "TransitionParameters.hpp"
#include "ContextParameters.hpp"

namespace ConsensusCore {
    
    class TemplateParameterPair {
    public:
        std::string tpl;
        std::vector<TransitionParameters> trans_probs;
        
        TemplateParameterPair(const std::string& tpl_,
                              const std::vector<TransitionParameters>& trans_probs_);
        
        TemplateParameterPair();
        
        // Copy constructor
        TemplateParameterPair(const TemplateParameterPair& other);
        
        TemplateParameterPair(const std::string& tpl_, const ContextParameters& ctx);
        
        /* Get a subsection of this template parameter pair
         TODO: This is a brutal copy operation */
        TemplateParameterPair GetSubSection(int start, int end) const;
        
        // Move assignment operator
        TemplateParameterPair& operator=(TemplateParameterPair&& rhs) = default;
        
        // Destructor
        ~TemplateParameterPair() = default;
        
        // Copy assignment
        TemplateParameterPair& operator=(const TemplateParameterPair& rhs) = default;
        
        //Move Constructor
        TemplateParameterPair(TemplateParameterPair&& src) = default;
        
        
        
    };
}

#endif /* defined(__ConsensusCoreXCode__TemplateParameterPair__) */
