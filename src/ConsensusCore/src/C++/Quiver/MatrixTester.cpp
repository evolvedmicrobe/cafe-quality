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
        // A series of tests, all the correct values are derived from the C# code.
        double eps = .00001; // Maximum percentage difference allowed between C# and C++ (note: log/exp order differs in implementation).
        SNR snr(10.0,7.0,5.0,11.0);
        ModelParams mp;
        ContextParameters ctx_params(snr);
        Read r("tester", "ACGTACGT");
        
        // Create a mutation scorer to test values
        TemplateParameterPair tpp("ACGTCGT", ctx_params);
        QvEvaluator qv(r, tpp, mp);
        SimpleQvSumProductRecursor bo(BandingOptions(4, 12.5));
        SimpleQvSumProductMutationScorer t(qv, bo);
        double score = t.Score();
        assert(fabs(1- score / -4.94222030733063 ) < eps);
        
        // Test an insertion mutation
        Mutation m(MutationType::INSERTION, 4, 'A');
        auto new_score = t.ScoreMutation(m,ctx_params);
        assert(fabs(1- new_score / -0.584415070238446) < eps);
        
        
        // Check the CC Context
        TemplateParameterPair tpp2("ACCTCGT", ctx_params);
        QvEvaluator qv2(r, tpp2, mp);
        SimpleQvSumProductRecursor bo2(BandingOptions(4, 12.5));
        SimpleQvSumProductMutationScorer t2(qv2, bo2);
        double score2 = t2.Score();
        assert(fabs(1- score2 / -10.4362503093273) < eps);
        
        
        // Now get the same value by mutating the original template
        Mutation m2(MutationType::SUBSTITUTION, 2, 'C');
        auto new_score2 = t.ScoreMutation(m2, ctx_params);
        assert(fabs(1- new_score2 / -10.4362503093273) < eps);
        
        // Test deletion near the end (goes through link/alpha beta path).
        Mutation m3(MutationType::DELETION, 4,'-');
        auto score3 = t.ScoreMutation(m3, ctx_params);
        assert(fabs(1- score3 / -9.89216068954291) < eps);
        
        // Test a deletion of the very last base.
        Mutation m4(MutationType::DELETION, 6, '-');
        auto score4 = t.ScoreMutation(m4, ctx_params);
        assert(fabs(1  - score4 / -15.6788158527151) < eps);
        
        // Test an insertion at the very last base
        Mutation m5(MutationType::INSERTION, 7, 'T');
        auto score5 = t.ScoreMutation(m5, ctx_params);
        assert(fabs(1  - score5 / -8.99810225167093) < eps);
        
        // Test a deletion of the first base
        Mutation m6(MutationType::DELETION, 0, '-');
        auto score6 = t.ScoreMutation(m6, ctx_params);
        assert(fabs(1  - score6 / -16.6208180854335) < eps);

       // Test an insertion at the first base
        Mutation m7(MutationType::INSERTION, 0, 'A');
        auto score7 = t.ScoreMutation(m7, ctx_params);
        assert(fabs(1  - score7 / -7.51178602234865) < eps);
        
        

        // Substitution in middle mutations to test link alpha-beta
        Mutation m8(MutationType::SUBSTITUTION, 4, 'A');
        auto score8 = t.ScoreMutation(m8, ctx_params);
        assert(fabs(1  - score8 / -5.23558996122357) < eps);
       
        // Insertion in middle to test link alpha-beta
        Mutation m9(MutationType::INSERTION, 4, 'G');
        auto score9 = t.ScoreMutation(m9, ctx_params);
        assert(fabs(1  - score9 / -6.71553495654471) < eps);
        
        
        return 0.0;
    }
}

int main() {
    ConsensusCore::MatrixTester mt;
    std::cout << mt.TestMutationScorer();
    
}