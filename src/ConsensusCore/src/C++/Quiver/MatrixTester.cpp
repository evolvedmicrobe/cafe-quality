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
        ContextParameters ctx_params(snr);
        Read r("tester", "ACGTACGT");
        
        TemplateParameterPair tpp("ACGTCGT", ctx_params);
        QvEvaluator qv(r, tpp, mp);
        SimpleQvSumProductRecursor bo(BandingOptions(4, 12.5));
        myfile << "Made bits";
        myfile.close();
        SimpleQvSumProductMutationScorer t(qv, bo);
        // C# wants this to be -4.94222030733063
        double score = t.Score();
        
        // C# Wants this to be -0.584415070238446
        Mutation m(MutationType::INSERTION, 4, 4, "A");
        auto new_score = t.ScoreMutation(m,ctx_params);
        
        //  Looking for 10.4362503093273
        Mutation m2(MutationType::SUBSTITUTION, 2, 3,"C");
        auto new_score2 = t.ScoreMutation(m2, ctx_params);

        // This should be equal to -9.89216068954291
        Mutation m3(MutationType::DELETION, 4,4);
        auto score3 = t.ScoreMutation(m3, ctx_params);
        
        TemplateParameterPair tpp2("ACCTCGT", ctx_params);
        QvEvaluator qv2(r, tpp2, mp);
        SimpleQvSumProductRecursor bo2(BandingOptions(4, 12.5));
        myfile << "Made bits";
        myfile.close();
        SimpleQvSumProductMutationScorer t2(qv2, bo2);
        // C# wants this to be -10.4362503093273
        double score2 = t2.Score();
        
        
        
        return score2;
    }
}

int main() {
    ConsensusCore::MatrixTester mt;
    std::cout << mt.TestMutationScorer();
    
}