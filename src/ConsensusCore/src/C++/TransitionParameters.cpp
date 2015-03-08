//
//  TransitionParameters.cpp
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/23/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#include "TransitionParameters.hpp"
#include "Quiver/MathUtils.h"

namespace ConsensusCore {
    

    
    
    TransitionParameters::TransitionParameters(double match, double stick,
                                               double branch, double deletion) :
                                               Match(match), Stick(stick),
                                               Branch(branch),Deletion(deletion) {}
    
    
    TransitionParameters::TransitionParameters() :
    Match(0.0) ,Stick(0.0), Branch(0.0), Deletion(0.0) {}
    
    TransitionParameters::TransitionParameters(const ConsensusCore::TransitionParameters &other) :
    Match(other.Match), Stick(other.Stick), Branch(other.Branch), Deletion(other.Deletion) {}
    
    double TransitionParameters::CalculateTotal() const {
        return logsumlog(Match, Stick, Branch, Deletion);
        
    }
    
    void TransitionParameters::RemoveConstant(double value) {
        Match -= value;
        Stick -= value;
        Branch -= value;
        Deletion -= value;
    }
}
