//
//  MatrixTester.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/7/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#pragma once

#include <stdio.h>
#include "Quiver/QuiverConfig.hpp"
#include "Quiver/MutationScorer.hpp"

namespace ConsensusCore {
    class MatrixTester {
    public:
        int TestMutationScorer();
        int TestMultiReadScorer();
    };
}
