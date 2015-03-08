//
//  MatrixTester.cpp
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/7/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#include "MatrixTester.hpp"
#include <fstream>
#include <iostream>

namespace ConsensusCore {
    double  MatrixTester::TestMutationScorer()
    {
        std::ofstream myfile;
        myfile.open("Temp.txt");
        SNR snr(10.0,7.0,5.0,11.0);
        ModelParams mp;
        Read r("tester", "ACGTACGT");
        ContextParameters ctx_params(snr);
        TemplateParameterPair tpp("ACGTCGT", ctx_params);
        QvEvaluator qv(r, tpp, mp);
        SimpleQvSumProductRecursor bo(BandingOptions(4, 12.5));
        myfile << "Made bits";
        myfile.close();
        SimpleQvSumProductMutationScorer t(qv, bo);
        double score = t.Score();
        // C# wants this to be -4.94222030733063
        return score;
    }
}

int main() {
    ConsensusCore::MatrixTester mt;
    std::cout << mt.TestMutationScorer();
    
}