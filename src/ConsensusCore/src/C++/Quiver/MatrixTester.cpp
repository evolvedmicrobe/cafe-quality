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
#include "MultiReadMutationScorer.hpp"
#include "Read.hpp"

namespace ConsensusCore {
    #define ASSERT(o, e) \
    { \
        if (!(fabs(1 - o/e) < 1e-5)) \
            std::cerr << "failed assert at " __FILE__ ":" << __LINE__ << "! o = " << o << ", e = " << e << std::endl; \
        else \
            std::cerr << "passed test at " __FILE__ ":" << __LINE__ << std::endl; \
    }

    typedef SparseSimpleQvSumProductMutationScorer ScorerType;
    typedef typename ScorerType::RecursorType RecursorType;

    int  MatrixTester::TestMutationScorer()
    {
        // A series of tests, all the correct values are derived from the C# code.
        SNR snr(10.0,7.0,5.0,11.0);
        ModelParams mp;
        ContextParameters ctx_params(snr);
        Read r("tester", "ACGTACGT");
        
        // Create a mutation scorer to test values
        TemplateParameterPair tpp("ACGTCGT", ctx_params);
        QvEvaluator qv(r, tpp, mp);
        RecursorType bo(BandingOptions(4, 12.5));
        ScorerType t(qv, bo);
        double score = t.Score();

        ASSERT(score, -4.94222030733063);
        
        // Test an insertion mutation
        Mutation m(MutationType::INSERTION, 4, 'A');
        auto new_score = t.ScoreMutation(m,ctx_params);
        ASSERT(new_score, -0.584415070238446);
        
        // Check the CC Context
        TemplateParameterPair tpp2("ACCTCGT", ctx_params);
        QvEvaluator qv2(r, tpp2, mp);
        RecursorType bo2(BandingOptions(4, 12.5));
        ScorerType t2(qv2, bo2);
        double score2 = t2.Score();
        ASSERT(score2, -10.4362503093273);
        
        
        // Now get the same value by mutating the original template
        Mutation m2(MutationType::SUBSTITUTION, 2, 'C');
        auto new_score2 = t.ScoreMutation(m2, ctx_params);
        ASSERT(new_score2, -10.4362503093273);
        
        // Test deletion near the end (goes through link/alpha beta path).
        Mutation m3(MutationType::DELETION, 4,'-');
        auto score3 = t.ScoreMutation(m3, ctx_params);
        ASSERT(score3, -9.89216068954291);
        
        // Test a deletion of the very last base.
        Mutation m4(MutationType::DELETION, 6, '-');
        auto score4 = t.ScoreMutation(m4, ctx_params);
        ASSERT(score4, -15.6788158527151);
        
        // Test an insertion at the very last base
        Mutation m5(MutationType::INSERTION, 7, 'T');
        auto score5 = t.ScoreMutation(m5, ctx_params);
        ASSERT(score5, -8.99810225167093);
        
        // Test a deletion of the first base
        Mutation m6(MutationType::DELETION, 0, '-');
        auto score6 = t.ScoreMutation(m6, ctx_params);
        ASSERT(score6, -16.6208180854335);

       // Test an insertion at the first base
        Mutation m7(MutationType::INSERTION, 0, 'A');
        auto score7 = t.ScoreMutation(m7, ctx_params);
        ASSERT(score7, -7.51178602234865);
        
        

        // Substitution in middle mutations to test link alpha-beta
        Mutation m8(MutationType::SUBSTITUTION, 4, 'A');
        auto score8 = t.ScoreMutation(m8, ctx_params);
        ASSERT(score8, -5.23558996122357);
       
        // Insertion in middle to test link alpha-beta
        Mutation m9(MutationType::INSERTION, 4, 'G');
        auto score9 = t.ScoreMutation(m9, ctx_params);
        ASSERT(score9, -6.71553495654471);
        
        return 0;
    }
    
    int MatrixTester::TestMultiReadScorer() {
        BandingOptions bo(3, 18);
        double fast_sçore_threshold = -12.5;
        
        
        
        std::string temp = "AGAGAGATAGCTACTAGTCCTCAGCAAGCTTGATCACACTATATGCGAGCGCGATAGATCGCTCTGCATCGTCACGATGTGTGTATATGACTGAGAGTCATACTATCTCTGCTACGCTCGACGTAGCGCTCATGTCGTCTAGTATGCGTGAGACGACGTAGCAGATACATGAGTGACAGACTCAGCAGTGCGCACAGTCACAGCTGTAGCATCGTACTCTACT";
        SNR snr(15.4944181442261,8.78859329223633,13.521107673645,14.9640893936157);
        ContextParameters ctx_params(snr);
        QuiverConfig qc(ctx_params, bo, fast_sçore_threshold);
        SparseSimpleSumProductMultiReadMutationScorer scorer(qc, temp);
        Read r1("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/19492_19674 RQ=0.848","GAGATACGATGCTACAGGGCTGTGACTGTGCGCACTGGCATGAGATCTGTCACTCCTAATGGTGTATCTTGCTACGCTTCGTACTCTCAGCGCAATACTAGGAACGACAATTGAGCGCGTTACGTCGAGCGTAGCAAGAGGATAGTATGACTCTCAGTCAATAACGACGACAGTTCGGAACG");
        MappedRead mr1(r1, StrandEnum::REVERSE_STRAND, 68, 218, false, false );
        auto result1 = scorer.AddRead(mr1, 1.0);
        Read r2("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/18859_19111 RQ=0.848","GAAGTACGTATTGCTACAGCTAGATGGACTGTGCGCACTGCTGAGTCGTGTCACTCCATGTATCTGCTAACCGCTCGTTCTCACGCATGGACTAGACGACATGAGCCGCTAACGTCGAGCGTCAGCGAGAGATAGTATGACTCTCAATGTCATATGACACACATCGGTGACGAGGTGCATGAGCGATCCTATCGCGCTCAGCATATAGTGGTGGATCAAGCCTTGCCGTGAGGACTAGTAGGCTTCTCTCTC");
        MappedRead mr2(r2, StrandEnum::REVERSE_STRAND, 0, 218, false, true );
        auto result2 = scorer.AddRead(mr2, 1.0);
        Read r3("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/9514_9553 RQ=0.848","TAGAGTAGATGTTAGTCAGCTGTGACTGTGCGCACTGCT");
        MappedRead mr3(r3, StrandEnum::FORWARD_STRAND, 169, 205, false, false );
        auto result3 = scorer.AddRead(mr3, 1.0);
        Read r4("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/16323_16580 RQ=0.848","AGTAGAAGTACTGATGCTACAGCTGCTGACTGTGCGCACTGCTTGAGTCTGTCACTCTATGTATCTGCTACGGTCGTCTCACGCATACTAGAAACGAGCATTGAAAGCGCTACGTGGTCGTAGCGTAGCAGACGAGATAGTATGACTCTACAGTCATAGGTACACACATCGTCGACGATGCCCAGAGCGATCGTAATCGGCGCTCCGAGCATATAGATGTGATCAGCATTGCTGAGGACTAGGTAGCTTCTTCTCTC");
        MappedRead mr4(r4, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result4 = scorer.AddRead(mr4, 1.0);
        Read r5("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/15332_15587 RQ=0.848","CAGAGAGAGGAAGCGTACTAGTTCCTCATGCAAGCTTGATCCAACGATCTATATGCGCAGCGCGTATATGAGTTCGCTCTGCATCGTCCACGATTGTTGGTGTATAATATGACTGAAGAGGTCATACATATCTCCTGCCTACGCTCGACGTAGCGCTCATGTCGTCCTATATGGCGGTGAGACGACGTAGCAGGATACATGAGCTGAAGACTCAGCATGTGCGCAGCAGTCAGCAGCTGTAGCATGCGTACTTAC");
        MappedRead mr5(r5, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result5 = scorer.AddRead(mr5, 1.0);
        Read r6("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/18204_18466 RQ=0.848","AGTAGCAGTACGATGCTACAGGCTGTGACCTGGTCGCACTGCTGTAGATTCTGTCACTGCATGTTCTGGCTACGGTCGTCTCACGCAATGACTAGAGGACTGACATGAGCGCTTACGTCGAGCGTAGCAGAGATAGTATGTTAACTCTCAGTCATATACAACACATTTTCGCTGACGAGTGCAAGAGCGATGCTAATGCGCGGAGCCTCGCAAATATATGTGTTAGATCAACGCTTGCTGAGGATCTAGTAGCTTCTCTCTC");
        MappedRead mr6(r6, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result6 = scorer.AddRead(mr6, 1.0);
        Read r7("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/15692_15951 RQ=0.848","GAGTACGATGGCACAGCTGTGACTGTGCGCAATGTGCGTGTTAGTCCTGGTCAGCTCATGTATCTGCTACGGTCGGTCTCGACGCATACTAGACGACATGAATGGCGATACGTCGCGAGCGTAAAATAGAGATTAGTATGACTCGTCAGTCAATATTACAACACGCACTCGTGGACGATGCACGAGCTGGGATCTATCGCGCTCGGACATTATAGTGTGATCAATAGCTTTGACTAGGACTAGTAGCTTTTCTCATCTC");
        MappedRead mr7(r7, StrandEnum::REVERSE_STRAND, 0, 218, false, true );
        auto result7 = scorer.AddRead(mr7, 1.0);
        Read r8("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/8795_8817 RQ=0.848","AGCTGTGACTGTGCAGCACTGC");
        MappedRead mr8(r8, StrandEnum::FORWARD_STRAND, 186, 204, false, false );
        auto result8 = scorer.AddRead(mr8, 1.0);
        Read r9("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/2226_2455 RQ=0.848","AGAGAGAAAGAAGCTCCTCAGTTCCTCAGCAAGCTTGATCAACTATATGCGAGCGCGATAGATCGGCCTCTGCATCGTCACGAATGTTGTGTATATGAACTGAGAGTCATACTATCTCATGCTACGCCTCGACGTAGCGCTCATGTCGTCTCAAGTTATGGCGTGAGAGCGACTGTTTTAGCAAGATACATGAGTGACAGACTCAGCAATGTGCGGCACCCCGCTCACA");
        MappedRead mr9(r9, StrandEnum::FORWARD_STRAND, 0, 202, true, false );
        auto result9 = scorer.AddRead(mr9, 1.0);
        Read r10("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/16613_16860 RQ=0.848","AGAGAGAAGCTACTGAGTCTCAGCAAGGCTTGTCAACCTACTAATATGCGAGCGCGTGATAGATCGCTTCTGGCATCGTCAACGATGTGTGTATACTGACTGAGAGTCATACTATCTCTGAACTGACGCTCGACGTAGCGCTCATGTCGTCTAGTATGCGTGATGACGACGTAGGCAGACATACAGTGAGTGAATGAGCAGACTCAGCAGTTGCGCACAGTCACAGCGTGTAGCAATCGTACTCTAC");
        MappedRead mr10(r10, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result10 = scorer.AddRead(mr10, 1.0);
        Read r11("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/20400_20644 RQ=0.848","AGAGAAGAGAAAGCTACTAGTCCTCAGCAAGGCTTGATCACACTATATGCGAGAGCCGATAGATCGCTCTGCATCGGTCACGGATGTGTGTATATGACTGAGAGTCATACTATCTCTGGCTACGGCTCGACGTAGCGGCTCATGTCGTCGTAGTATGGCGTGGAAGAAACGACGTAGCAGATACATGATGACAGTACTCAGCAGTGCGACACAAGTCACAGCTGTAGCATCGTAAACTAGCTAC");
        MappedRead mr11(r11, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result11 = scorer.AddRead(mr11, 1.0);
        Read r12("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/6753_6994 RQ=0.848","AGTAGAGTACGAATTGCTACAGCTGTGACTGTGCGCACTGCTGAGCTCTGTCACTCATGTATCTGCTACGTCGTTCTCACGCATACTAGACGACATGAGCGCACGTCGAGCGTAGCGAGAGATAGTATGACTCTCAAGTCCTATACACACATTCGTGACGATGCCCAGAGCGATCTATATCGGGCTCGGCATATAGTGTGGATCAAAGCTTTGCTGAGGAACTAGTAAGCTTTTCTCTCTC");
        MappedRead mr12(r12, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result12 = scorer.AddRead(mr12, 1.0);
        Read r13("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/7969_8009 RQ=0.848","ATGGGCAGTGCGCACAGTTCACAGCCTGTAGCATCGTACT");
        MappedRead mr13(r13, StrandEnum::REVERSE_STRAND, 170, 204, false, false );
        auto result13 = scorer.AddRead(mr13, 1.0);
        Read r14("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/21345_21579 RQ=0.848","AGTAGAAGTTACGATGCTACAGCTGTGACTGTGCGCACTGCTGATCCTGTCACTCAATGTATCTGCTACGTCGTCTCACGCATACTAAGACGACAATGAGCGCTACGTCGAGCGTAGCAGAGAATAGTAATGACCTCTCAGGTCATATACACACATCGTGACGAATGCAGAGCGATCTATCGCGCTCGCATATAGTGTGATCAAGCTTGCTGAGGACTAGTAGCTTCTTCTCTC");
        MappedRead mr14(r14, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result14 = scorer.AddRead(mr14, 1.0);
        Read r15("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/3458_3698 RQ=0.848","AGAGAGAGAAAGCACTAGTCCTCAGCAAGCTTGATAACACTAGTATGCGAGCGCGATAAAAGAGGATTCGCTCTGCATGCGTCACGAAAATGTGTGTATATGACTGAAGAGTCAATACTATCTCTGCTAACGCTGCGACGTAAGCGCTCATTCGTCTAGTTATGCGTGAGACGACGTAGCAGATACATGAGTGCAGACTCAGCTAGTGCGCCAGTCACAGCTTGTAGCATCGTACTCTAC");
        MappedRead mr15(r15, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result15 = scorer.AddRead(mr15, 1.0);
        Read r16("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/6423_6662 RQ=0.848","AGAGGAGAGAAAGCTACTAGTCCCTGCAGCAAAGCTTGATCAGAACTATTGCGAGCGCGATAGATCGCTCGCATCGTCACGATGTGTTTGTATATGACTGAGAAGTCCATACTATGCTCTGCTACGCTCGGACGTAAGCGCTTGCATGTCCGTCTTAGTCATGAGTGAGACGAACGTAGCAAGATACATGAGTGACAGACTCAGCAGTGTGCACAGTCACAGCTGTAGCATGTACTCTA");
        MappedRead mr16(r16, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result16 = scorer.AddRead(mr16, 1.0);
        Read r17("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/17240_17485 RQ=0.848","AGAGAGAGAAAGCTACTTAGTCTCAGCAAGCTTGATCACACTTGATATGCGAGCGCGATAGAGTGCTTCTGGCAATCGTCAACGATGTGTGTATGATGACTGAGAGTCAAACTAATCTCTGCTCAACGCTTCGAGCGTAGCGCTCGCATGTCGTCTAGTATGTCGTGAGACGACGTAGCGAGATACATGAGTGACAGACTCAGCAGTGCGCGCACAGTCACAGGCTGTAGCATCGTACTTTCTAC");
        MappedRead mr17(r17, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result17 = scorer.AddRead(mr17, 1.0);
        Read r18("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/17861_18106 RQ=0.848","AGAGAGAGAAGCTACTAGTCCGGTCAGCAAGCTTGATCACACTGATAAGTCGAGGCGCGATAGATCGCTTGCCTGCGACTCGTTCACGATGTGTGTATATTGAATGAGAAGTCATACTATCTCTGCTACGGCTCGACCGATAGCGCTGCATGTCGTCTAGTATGGTGAGACGACGTAGCAGATACATGAGTGACAAGACTCAGCATGTGCGAACAGTCGCACAGCTGTAGCATCGTATCTCTACC");
        MappedRead mr18(r18, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result18 = scorer.AddRead(mr18, 1.0);
        Read r19("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/17577_17832 RQ=0.848","AAGTAGAGTACGATTGCTACAGCTGTGACTGTGCGCACTGCTGAGTAGTTCTGTCACTTGGCATGTATCTGCTAGTACAGTCGTCTCACGCATACTAGGACGAACATGGAGGCCGCTGAACGTCGAGCGTAGGCCGAGAATAGTTGACTCTCGAGTCATTTACACACATCGTGACGATGCAAAGAGCGATCTAGTCGGCGCTCGCATATATGTGTGAATGGCAGCTTGCTGAAGGATAGTAGCATTCCTCTCTCT");
        MappedRead mr19(r19, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result19 = scorer.AddRead(mr19, 1.0);
        Read r20("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/16968_17208 RQ=0.848","GTAGACGTACGAATGCATACAGTGTGACTGTGCGCCACTGCTGATCTGTCACTCATAGTATCTGCTACGGTCGTCTCACGCTGATACTAAGAGAACGACCATGAGCGCTATGCGTCGAGCGTTAGCAGAGATAGTATAGACTTCAGTCATATACACACATCGTGACGATGCAGGAAGCGATCTATCGCGCTCGCATTATGCGTGTGATCAAGCTTGCTGAGGACTAGTAGCTCGTCTCTC");
        MappedRead mr20(r20, StrandEnum::REVERSE_STRAND, 0, 221, false, true );
        auto result20 = scorer.AddRead(mr20, 1.0);
        Read r21("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/8636_8676 RQ=0.848","ATGCAGTGCGCAGCTCGTCACAGCTGTAGCATCGTACTCT");
        MappedRead mr21(r21, StrandEnum::REVERSE_STRAND, 169, 204, false, false );
        auto result21 = scorer.AddRead(mr21, 1.0);
        Read r22("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/1324_1576 RQ=0.848","AGTAGAGTACCGACTGCTACAGCTGTGCATGTGCGGCCTGTGAGCTCTGTCCACTCCATATGTATCTGCTTAACGTCGTCTCACCGCATACCCTAGACCGACATGAGCGCGCTACGGTCGAGCGTAGCACCGACGATATCCGTATGAACTCTCAGTCATTACCACACCACTCGTGACGATTGCCAGAGCGATCTATCGCGCTCGCAATATAGTGGTGATCAAGCTTGCTGAGACTTAGTAGCTTCTCTCTGC");
        MappedRead mr22(r22, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result22 = scorer.AddRead(mr22, 1.0);
        Read r23("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/6157_6393 RQ=0.848","ATAGGAGTAGGCGATGCTATCAGCTGTGACTTGTGCGCACTGCTGAGTCCTGTCACTCATGTATCTGCTACGTCGTCTGCACGTCATACTAGACGACATGAGCGCTACGTCGAGCCGTAAGCAGAGATAGTATGACTCTCAGTCAATATTACACACATCGGACGAGCAGAGCGAATCTATTCGCGCTCGCATATAGTGGTGATCAAGCTTGCTGGGATAGTAGCTTTTTCTGTCTC");
        MappedRead mr23(r23, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result23 = scorer.AddRead(mr23, 1.0);
        Read r24("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/15983_16224 RQ=0.848","AGAGAGAGAGAGCGTAGCTAGTCCTGCAGCAAAGCTGATCACACTAATAATGCGAGCGCGATAGTCGCTCGATGCACTGCGTCACGATGTGTGTATATGACCTGAAGCAGTCATACTATCTCTGCTACGCTCGACGTAGCGCTCATGGTCGTCTAGTATGCGTGAGAACGACGTAGCAGATACATGAAGTGACAGACTCAGCAGTGCGCACAGTCACAGGCCTGTAGCAATCGTACTCTAC");
        MappedRead mr24(r24, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result24 = scorer.AddRead(mr24, 1.0);
        Read r25("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/1951_2192 RQ=0.848","AGTAGAGTACGAGATGCTACAGCTGTGACTGTGCCGCACTGCTCGATGTCCTGTCAAATCATGTATCCTGCTACATTCGTCGTCCACGCATGACTTAAGATCGACATGCAGCGCTACGTCGAGCGCTAGAAGATAGTATGACCTCTCAGTCATATACAACAACATCGTGACGATGCTACGAGCGATCTACTGCGCTTCGCATCATACGTGTCGATCAAGCTTTGACTAGTAGCTCTCTCTC");
        MappedRead mr25(r25, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result25 = scorer.AddRead(mr25, 1.0);
        Read r26("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/12671_12707 RQ=0.848","ACAGTGCGCATACAGTCACAGCTGTAGCATGCCGTG");
        MappedRead mr26(r26, StrandEnum::REVERSE_STRAND, 169, 204, false, false );
        auto result26 = scorer.AddRead(mr26, 1.0);
        Read r27("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/20743_20984 RQ=0.848","AGTAAGAGTACGATGCATAGGCTTGACTGTGCGCACTGTGCTGAGTCTGTCACATAGTTTCTGCTACGCTGTCTCACGGCATAACTAGTACGAGCATGAGCGCCTACGTCGAGCGTAGCAGAGAAGTAATGACTCTCAGTCATAACTACACGACATCGTGACGATGCCAGCAGCGAATTATCGCGCTTCGGGCATAGTAGTGTTGATCAAAGCTTGCTGAGGACTAGTAGTCTTCTCTCTC");
        MappedRead mr27(r27, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result27 = scorer.AddRead(mr27, 1.0);
        Read r28("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/3802_4037 RQ=0.848","AGAGAGTACGATGCTACAAGCTGTGACTGTGGCGCACTGCTGAGATCTGTCACTCATCGTATCTGCTACGTCGTCTCACGCATAACTAGACACATGAGCGTCTAGCGTCGAGGGCGTAGCAGAGATAGTATGACTCTGCAGTCATATCAGCACATCGTGAGCGATGCAGAGCGCCTCTATCGGCGCTCGCATAATAGTGTGATCAAGCTTGCTGGAGGACTAGTAGCTTCTTCTC");
        MappedRead mr28(r28, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result28 = scorer.AddRead(mr28, 1.0);
        Read r29("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/10553_10784 RQ=0.848","AGAGAGAGAAGCTACTAGTCCTCAGCAGCTTGATCAACACTATATGCGAGCGCGATAGATCGCTCTGCATCGTCAGATGTGTGTTATATGACTGAGAGTCATACTATCTCTGCTACGCCTCGAGTAGGCGCCTCATGTCGTCTAGTAATGCGTGAGACGAACGTAAGCACAGATACATGAGTGACAGACTCAGCAGTGCGCACCAGTCACAGCTGTAGCATCGTACTCTAC");
        MappedRead mr29(r29, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result29 = scorer.AddRead(mr29, 1.0);
        Read r30("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/18502_18748 RQ=0.848","AGAGAGAGCTACTAGTCCTCAAGCCGTCTTGAGTCAACACTCATGCCGATGCGATGCGCGATAAGATCGCATCGCATTTTCGTCACGGATGTGTGTTATATGACTGAGAGTCATACTTCTCGTGGCTAGCTCGACGGTAGCGCTGCATGTCCGTCTAGTATTGCGTGAGAAAGACGTAAGCCAGGATACAGTGAGTGACAGACTCAGGGCAGCTGCGCACAGTCACAGCTGTATGCATCGGTGGAT");
        MappedRead mr30(r30, StrandEnum::FORWARD_STRAND, 0, 218, true, false );
        auto result30 = scorer.AddRead(mr30, 1.0);
        Read r31("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/20131_20370 RQ=0.848","AGTAAACAGATGTACAGCCTTCGCTGTGACTGTGGCACTGTGAGTCTGTCAACTCATGTATCTGCTACGTCGTCTGCACGCATACCTAGAACGACATTGACAGCGGCTACGTCGACGTAGCAGAGAATAGGTATGACTCTCAGGTCATATACCAGCATCGTGACGATGCAGAGCGATCTATACGGCGTGCATATAGTGTGATCAAGCTTGCTGAGAGACTAGCTAGATTCTCTCTAAAT");
        MappedRead mr31(r31, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result31 = scorer.AddRead(mr31, 1.0);
        Read r32("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/4654_4882 RQ=0.848","AGAGAGAGAAGCTTACTAGTCCTCAGCAAAGCCTTTGTCACAACTATATTGCGACGCGATAGATCGCTCTGCATCGTGCACGATGGTGTGTATATGACTGAGAGTCATACTATCTTGCTACGCTCGACGTAGCGCTCATGTCGTCTAGTAGGTGAGACGACGTTAAGCAGATACATGAGTGACAGACTCAGCAGTGCGCACAGTCACAGCTGTACATTCGTACTCTAC");
        MappedRead mr32(r32, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result32 = scorer.AddRead(mr32, 1.0);
        Read r33("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/5573_5804 RQ=0.848","AGTAGAGTACGATGCCTACAGCTGTGACTGTGCGCACTGCTGAGTCTGTCTACTCATGTATCTGCTACGTCGTCTCAACGCATACTAGAGACATGAGGCGCTACGTCGAGCGTAGCAGAGATAGTAATGACTCCTGCAGTCATATACACCAACATCGTTACGATGCCAGAGCGATTATCGGCGCTCGCATATAGTGTGTCAGCTTGCTTGAGGACTAGTAGCTTCTCTCTC");
        MappedRead mr33(r33, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result33 = scorer.AddRead(mr33, 1.0);
        Read r34("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/5243_5464 RQ=0.848","AGAGAGAGAAGCTACTAAGTCTCAGCAAGCTTGATGCACACTATATGCGAGCGCGATAGATCGTCTGCATCGTCAGATGTGTGTATAATGACTGAGAAGTCATACTATCTCTGCTACGCTCGACGTTAGCGCTCATTGTCGTCTAGTATGCGTGAACGACGTAGGCAGATACATGAGTGACAGACTCAGCAGTGCGCAAGTCACAGCTGTAGCATCGTACT");
        MappedRead mr34(r34, StrandEnum::FORWARD_STRAND, 0, 218, true, false );
        auto result34 = scorer.AddRead(mr34, 1.0);
        Read r35("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/4066_4296 RQ=0.848","AGAGAGAAGAAGCTACTAGTCCTCAGCAAGCTTGATCACACTATATGCGAGCGCGATAGATCGCTCTGGCATCGTCACGATGTGTGTATAATGAACTGAGAGTTCATATACTCTGCTAACGCTGCGACGTAGCGCTCATGTCGTCCTCAGTATCCGTGAAGACGACGTGCAGATACATGAGTGACAGACTCAGCATGCGCACAGTCACAGCTGTAGCATCGGTACTCTAC");
        MappedRead mr35(r35, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result35 = scorer.AddRead(mr35, 1.0);
        Read r36("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/10882_11119 RQ=0.848","AAGTAGAGTATGCGATGCTACAGCTGTGACTGTGCGCACTGCTGAGTCTGTCACTCATGTTCTGCTACGTCCGTCTCACGCATACTAAGCGAATGAGCGCTACGTCGACGTAGCAGAGATAGTAATGACTCTCAGTCAATTACACACATCCGTGACGATGGCAGAGCGATCTATCGCGCTGCATATAGTGTGAATCCAAGCTTTGCTGGAGTGACTAGTGATAGCTATATTTCTCTC");
        MappedRead mr36(r36, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result36 = scorer.AddRead(mr36, 1.0);
        Read r37("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/4395_4623 RQ=0.848","AGTAGAGTACGATGCTACAGCTGGTGACTGTGGCGCACTGCTGAGTCCCTGTCACTCATGTATCTGCTACGGTCGTCTCCACGCATAGACGACATGAGCGCTACGGTCGAAGCGTAGCAGAGATAGTATGACTCTCAGTCATATACACACATCGTGACGATGCAAGAGCGATCCTATTCGCGCTCGCATATAGTGTGATCAGCTTGCTGAGGACTATAGCTCTCTCTC");
        MappedRead mr37(r37, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result37 = scorer.AddRead(mr37, 1.0);
        Read r38("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/1608_1846 RQ=0.848","AGAGAGAAGCTACCTAGTCCTCAGCAAGCTTTGATGCAGCAACTATATGCGAGCGCGAATAGATCGCTCTGCATTCGTCACGCCTGTGTGTATTGACTGAGAAAGTCATACTATCTGCTGCTACGGCTCGACGTAGGCGCTCATGTCGTCGCTAGTATGCGTGAGACGACGTAAGAGATACATGAGTGACAGACTCAGGCAGTACGCGCACAGTCACAGCTGTAGCATCGTACTCTAC");
        MappedRead mr38(r38, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result38 = scorer.AddRead(mr38, 1.0);
        Read r39("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/4978_5211 RQ=0.848","AGTAGAAGTACGATGCTTACAGCTGTGACTGGGCGCAACTGCTGAGATCTGTCACTCCATGTATCTGCTACGTCGTTCACGCATACTAACGACATGAGCGCTACGTCGAGCGTGTAGCAGGAGATAGTAAATGCACTCTCAGTCATATACCACATCGTGAACGATGCAGAGCGATCTATCTGGCTCGCATATAGTGTGACAAAGCCTTGCTGAGGACTAGTAGCTTCTCTCTC");
        MappedRead mr39(r39, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result39 = scorer.AddRead(mr39, 1.0);
        Read r40("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/19787_20017 RQ=0.848","AGAGAGAGAGCTACTGTCCTCAGCAAGCTTGATCACACTATATGGCGAGCGCGCGATAGAATCGCTCTGCATCGTCCAGATTGTGTGTAGTTATGACTGAGAGTCAGTACTATCTCTGCTCGCTTCGACGTGCGGCTCATGTCGTCTAGTAAATGGCGGTGAGACTGACGTAAAGCAGGATAACATGAGTTGACAAGACCAGCAGTGCCGCAAGTCAGCAGCTGTAGCAT");
        MappedRead mr40(r40, StrandEnum::FORWARD_STRAND, 0, 212, true, false );
        auto result40 = scorer.AddRead(mr40, 1.0);
        Read r41("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/11156_11400 RQ=0.848","AGAGCGTAGATAGCTACTAGTCTCAGCAAGCTTTGTATGGCACAGCTGATATGCGAAGCGGCGATAGTATATATCTCACGATGTGTGTATAAGTGACGTGAGGAGCTGCTCACATACTATCTCTGCTACGCTCGACGCTAGGCGCTCACTGTCGTCTAGTATGCGTTGAAAGAACGACGTAGGCAGATACATGATGACAGACTCAGCAGTGCGCACAGTCACAGCTGTAGCATCGTTACTCTAC");
        MappedRead mr41(r41, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result41 = scorer.AddRead(mr41, 1.0);
        Read r42("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/5838_6061 RQ=0.848","AGAGAGAAGCTACTAGTCCTTCAGCAAGCTTGATACACTATATGCGAGCGCGATAGATCGCTCTGCATCGTCACGATGTGTGTATATGACTGAGAGTCATACTATCTCTGCTACGCTGACGTAGGCTCAATGTCGTCTAGTATTGCGTGAGACGACGTAGCAGATACATGAGTGACAGACTCAGTCAGTGCGCACAGTCACAGCTGTCGCATCGTACTCTAGC");
        MappedRead mr42(r42, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result42 = scorer.AddRead(mr42, 1.0);
        Read r43("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/21019_21250 RQ=0.848","AGAGAGAGAAGCTACTAGTGCCTCAGCAAGCTTGATCAGCACATTGCGAGCGCGGATAGATGCGCGTCTCATCGTCACGATGGCTGTGTATATCTACTGAAGAGTCCATACTATCCTCTGCTAACGCTCGCGTAGCGCTCAGTCGTCTAGTATGCGTGAGACGAAACGTAGCAGATACATGAGTGACAGACTCAGCTGCGCAACAGTCCACAGCTGTAGCATCGACTCTAC");
        MappedRead mr43(r43, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result43 = scorer.AddRead(mr43, 1.0);
        Read r44("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/2572_2817 RQ=0.848","AGTAGAGTACGATGCTACAGCTGTGACTGTCCGCGCACTGCTGACGTCTGTCACTCAATGTATCTGCCTACGGTCGTCTCATCGCATACTAACGACATGAGCCGCCTACGAGTCGAGCGTACGCTAGAGATAGTATGCTGCTCAGTCATATACACACATCGTCGACGATGCCCAGAGCGATTCTAATCGCGCCATCGCACTATAGTGTGATCAAGCTTGCTGAGGACTAGTAGGCTTTCTCTCTC");
        MappedRead mr44(r44, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result44 = scorer.AddRead(mr44, 1.0);
        Read r45("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/19143_19381 RQ=0.848","AGAGAGAGCTACTAGTATCTCAGAAGCTTGATCACAACTATATGCGAGACGCGATGATCGCTCTGCATTCGTCCACGATGGTGTGTATATGGACGTGAGAGTCATCACTATTCTCTGCTACGGCTCGGACGTAGCGCTCAAATGTCGTCGAGTATGCGTGAGAACGACGTTAGCAGATACATGAGTGACAGAGCTCAGGCAGTGCGCACAGTCACAGCTGTAAGCATCGTACTCTGAC");
        MappedRead mr45(r45, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result45 = scorer.AddRead(mr45, 1.0);
        Read r46("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/11497_11734 RQ=0.848","AGTAGAGTACGATGCTAAGCTGTGACTGTGCGCAACTGCGAGTCTGTCAACTCATGTATCTGCTAACGTCGTCTCACGCATACTAGACGACATGAGCCGCTACGTCGAGCGTAGCAGAGATAGGGTATGACTCGTGCAGTCATATACACATCGATGACGAATTGGCAGAGGCGATCGTACGGTCGGCTCACGCATATATGAAATCATGCTTGCTGGAAGGACTAGTAGCACTTTCTC");
        MappedRead mr46(r46, StrandEnum::REVERSE_STRAND, 0, 222, false, true );
        auto result46 = scorer.AddRead(mr46, 1.0);
        Read r47("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/21610_21833 RQ=0.848","AGAGAGGAAGCCTACTATCCTCAGGCAAGCTTGATCACAACTATATGCGAGCGCGATAGATCGCTTGCATGTCGATGTGTGTATATGACTGAGAGTCATACATCTCTTGCTACCTCGACGTGCGGCTCGATGTCGTCTAGTATGGCGTAGACGACGTAGCAGAATACATGAGTGACAACTGCAGCAGTGCGCACGTGCCACAGCTGTAGCATCGTAACTCTAC");
        MappedRead mr47(r47, StrandEnum::FORWARD_STRAND, 0, 222, true, false );
        auto result47 = scorer.AddRead(mr47, 1.0);
        Read r48("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/2853_3080 RQ=0.848","AGCAGAGAGAAAGCTACCTAGTCCTGCAACAAGCTCTGATCACACTATATCGCCCAGCGCGATAGATCGCTCTGCTATCGTCCACGATGTGTTATACTGACTTGAGAGTCCATACTATCTCTGCTCCGCTCGAGCGTAGCGCTCATGTCGTCTAGTATGCGTAGACGACGTAGGCAGATAACATGAGTGACAGACTCAGCAGTGTCGCCACAGTCCAGCCAGCTGTG");
        MappedRead mr48(r48, StrandEnum::FORWARD_STRAND, 0, 209, true, false );
        auto result48 = scorer.AddRead(mr48, 1.0);
        Read r49("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/3207_3426 RQ=0.848","CACGCTCGCTCGCACTCGCTCGCGCACTGCTGAGTCTGTCACTCATGTAATCTGCTACGTCGTCTCACGCATACTAGACGACATGATGCGCTACGTCGAAGCGTAGCAGAGATAGTATGACTCTCAGTTCATATACACACATCGTGACGATGGCAGAGCGATCTATCGCGCTCGCATATAGTGTGATCAAGCTTGCTGAGGACTAGTAGCTTCTCTCTC");
        MappedRead mr49(r49, StrandEnum::REVERSE_STRAND, 0, 205, false, true );
        auto result49 = scorer.AddRead(mr49, 1.0);
        Read r50("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/15022_15119 RQ=0.848","ATATGATTGATGCGCACATTAGTAGTCTGTCACTCAGTTAGTGATCTGCGCTACGGGTTTTCGTCTCACGCTCGATATATATCGCTGAGAAGTGACG");
        MappedRead mr50(r50, StrandEnum::REVERSE_STRAND, 112, 206, false, false );
        auto result50 = scorer.AddRead(mr50, 1.0);
        Read r51("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/12209_12224 RQ=0.848","GCCACTGACTGGAGA");
        MappedRead mr51(r51, StrandEnum::FORWARD_STRAND, 189, 203, false, false );
        auto result51 = scorer.AddRead(mr51, 1.0);
        Read r52("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/7070_7301 RQ=0.848","TATATTCGAGAAGGTTCATGCGATGCGTGCGATAGATCGCCTTCTGCATCGTCCAAGATGTTGTGTATATGACTGAGAAGTCAATACTATCTTCTGCTACGCTCGACGCTAGCAGCGTCAATGTCGTCCCTAGTATAGTGCGTGGAAGACGACGTAGCACGATACATGAGGTGAACAGACTCGTTAAAAGTGCGCATTAAGGGTCACACAGCTGTAGCCATCGTACTCTAC");
        MappedRead mr52(r52, StrandEnum::FORWARD_STRAND, 11, 222, false, false );
        auto result52 = scorer.AddRead(mr52, 1.0);
        Read r53("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/10405_10514 RQ=0.848","AACGCACCCAGTCGGTGACGAGTGGCAGAGCGATCTATACGCGCTGGTATATTAGTAAGGTGATCAGAGGCTTGCTAGAGGGACTAGTTTAGCCTTTCTTCTCTCGATC");
        MappedRead mr53(r53, StrandEnum::REVERSE_STRAND, 0, 83, false, true );
        auto result53 = scorer.AddRead(mr53, 1.0);
        Read r54("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/5/14719_14894 RQ=0.848","TAGTTGACATGGAGGAGTCACACATATTTGATCTCTTGCTACGCTATATTATAATGCGCGCTCCACTCGTGCGTCCGCTCTAGTAGTGCGTTGTGGACCGCGATCGTAGCAGATACATGAGATGACAGATCTTCAGCAGGTGCGCGCAAGTCACAGCTGTAGCATCGTACTCTAC");
        MappedRead mr54(r54, StrandEnum::FORWARD_STRAND, 83, 222, false, false );
        auto result54 = scorer.AddRead(mr54, 1.0);

        Mutation m(MutationType::INSERTION, 202, 'C');
        auto ress = scorer.Score(m);

        std::cout << ress;
        return 0;
    }
}

int main() {
    ConsensusCore::MatrixTester mt;
    return mt.TestMutationScorer();
    // std::cout << mt.TestMultiReadScorer();
}
