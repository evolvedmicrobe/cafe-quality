// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: David Alexander

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "Poa/PoaConfig.hpp"
#include "Poa/PoaConsensus.hpp"
#include "Utils.hpp"
#include "Mutation.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using boost::erase_all_copy;

using namespace boost::assign;  // NOLINT
using namespace ConsensusCore; // NOLINT

#define MAKE_ALL_PLOTS false

static void plotGraph(const PoaGraph* g, string description, bool REALLY_MAKE_THIS_ONE = false)
{
    if (MAKE_ALL_PLOTS || REALLY_MAKE_THIS_ONE)
    {
        string dotFname = description + ".dot";
        string pngFname = description + ".png";
        string cmd = string("dot -Tpng ") + dotFname + " > " + pngFname;
        g->WriteGraphVizFile(description + ".dot",
                             (PoaGraph::COLOR_NODES | PoaGraph::VERBOSE_NODES));
        cout << cmd << endl;
        system(cmd.c_str());
    }
}

// TEST(PoaGraph, NoReadsTest)
// {
//  // Test that it works with no reads
//  vector<const SequenceFeatures*> reads;
//  const PoaConsensus* pc = PoaConsensus::findConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
//  string dot = pc->getGraph()->toGraphViz();
//  cout << dot << endl;
// }

TEST(PoaGraph, SmallBasicTest)
{
    // Test that it works with a single sequence
    vector<std::string> reads;
    reads += "GGG";
    const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
    string dot = pc->Graph()->ToGraphViz();
    string expectedDot = \
                         "digraph G {"
                         "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                         "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                         "2[shape=Mrecord, label=\"{ G | 1 }\"];"
                         "3[shape=Mrecord, label=\"{ G | 1 }\"];"
                         "4[shape=Mrecord, label=\"{ G | 1 }\"];"
                         "0->2 ;"
                         "2->3 ;"
                         "3->4 ;"
                         "4->1 ;"
                         "}";
    plotGraph(pc->Graph(), "small-basic");
    EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
    EXPECT_EQ("GGG", pc->Sequence());
    delete pc;
}

TEST(PoaGraph, SmallExtraTests)
{
    // Extra at beginning
    {
        vector<std::string> reads;
        reads += "GGG", "TGGG";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        string expectedDot = \
                             "digraph G {"
                             "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                             "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                             "2[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "3[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "4[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "5[shape=Mrecord, label=\"{ T | 1 }\"];"
                             "0->2 ;"
                             "2->3 ;"
                             "3->4 ;"
                             "4->1 ;"
                             "5->2 ;"
                             "0->5 ;"
                             "}";
        plotGraph(pc->Graph(), "extra-at-beginning");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GGG", pc->Sequence());
        delete pc;
    }

    // Extra in middle
    {
        vector<std::string> reads;
        reads += "GGG", "GTGG";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        // this would be easier if we could use the C++0x raw
        // strings feature (in g++ 4.5+)
        string expectedDot = \
                             "digraph G {"
                             "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                             "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                             "2[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "3[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "4[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "5[shape=Mrecord, label=\"{ T | 1 }\"];"
                             "0->2 ;"
                             "2->3 ;"
                             "3->4 ;"
                             "4->1 ;"
                             "5->3 ;"
                             "2->5 ;"
                             "}";
        // plotGraph(pc->getGraph(), "extra-in-middle");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GGG", pc->Sequence());
        delete pc;
    }

    // Extra at end
    {
        vector<std::string> reads;
        reads += "GGG", "GGGT";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        string expectedDot = \
                             "digraph G {"
                             "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                             "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                             "2[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "3[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "4[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "5[shape=Mrecord, label=\"{ T | 1 }\"];"
                             "0->2 ;"
                             "2->3 ;"
                             "3->4 ;"
                             "4->1 ;"
                             "5->1 ;"
                             "4->5 ;"
                             "}";
        plotGraph(pc->Graph(), "extra-at-end");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GGG", pc->Sequence());
        delete pc;
    }
}


TEST(PoaGraph, SmallMismatchTests)
{
    // Mismatch at beginning
    {
        vector<std::string> reads;
        reads += "GGG", "TGG";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        string expectedDot = \
                             "digraph G {"
                             "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                             "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                             "2[shape=Mrecord, label=\"{ G | 1 }\"];"
                             "3[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "4[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "5[shape=Mrecord, label=\"{ T | 1 }\"];"
                             "0->2 ;"
                             "2->3 ;"
                             "3->4 ;"
                             "4->1 ;"
                             "5->3 ;"
                             "0->5 ;"
                             "}";
        plotGraph(pc->Graph(), "mismatch-at-beginning");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GG", pc->Sequence());
        delete pc;
    }

    // Mismatch in middle
    {
        vector<std::string> reads;
        reads += "GGG", "GTG", "GTG";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        string expectedDot = \
                             "digraph G {"
                             "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                             "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                             "2[shape=Mrecord, label=\"{ G | 3 }\"];"
                             "3[shape=Mrecord, label=\"{ G | 1 }\"];"
                             "4[shape=Mrecord, label=\"{ G | 3 }\"];"
                             "5[shape=Mrecord, label=\"{ T | 2 }\"];"
                             "0->2 ;"
                             "2->3 ;"
                             "3->4 ;"
                             "4->1 ;"
                             "5->4 ;"
                             "2->5 ;"
                             "}";
        plotGraph(pc->Graph(), "mismatch-in-middle");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GTG", pc->Sequence());
        delete pc;
    }

    // Mismatch at end
    {
        vector<std::string> reads;
        reads += "GGG", "GGT";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        string expectedDot = \
                             "digraph G {"
                             "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                             "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                             "2[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "3[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "4[shape=Mrecord, label=\"{ G | 1 }\"];"
                             "5[shape=Mrecord, label=\"{ T | 1 }\"];"
                             "0->2 ;"
                             "2->3 ;"
                             "3->4 ;"
                             "4->1 ;"
                             "5->1 ;"
                             "3->5 ;"
                             "}";
        plotGraph(pc->Graph(), "mismatch-at-end");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GG", pc->Sequence());
        delete pc;
    }
}

TEST(PoaGraph, SmallDeletionTests)
{
    // Deletion at beginning
    {
        vector<std::string> reads;
        reads += "GAT", "AT";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        string expectedDot = \
                             "digraph G {"
                             "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                             "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                             "2[shape=Mrecord, label=\"{ G | 1 }\"];"
                             "3[shape=Mrecord, label=\"{ A | 2 }\"];"
                             "4[shape=Mrecord, label=\"{ T | 2 }\"];"
                             "0->2 ;"
                             "2->3 ;"
                             "3->4 ;"
                             "4->1 ;"
                             "0->3 ;"
                             "}";
        plotGraph(pc->Graph(), "deletion-at-beginning");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("AT", pc->Sequence());
        delete pc;
    }

    // Deletion in middle
    {
        vector<std::string> reads;
        reads += "GAT", "GT";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        string expectedDot =
            "digraph G {"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 2 }\"];"
            "3[shape=Mrecord, label=\"{ A | 1 }\"];"
            "4[shape=Mrecord, label=\"{ T | 2 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "2->4 ;"
            "}";
        plotGraph(pc->Graph(), "deletion-in-middle");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        delete pc;
    }

    // Deletion at end
    {
        vector<std::string> reads;
        reads += "GAT", "GA";
        const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
        string dot = pc->Graph()->ToGraphViz();
        string expectedDot = \
                             "digraph G {"
                             "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                             "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                             "2[shape=Mrecord, label=\"{ G | 2 }\"];"
                             "3[shape=Mrecord, label=\"{ A | 2 }\"];"
                             "4[shape=Mrecord, label=\"{ T | 1 }\"];"
                             "0->2 ;"
                             "2->3 ;"
                             "3->4 ;"
                             "4->1 ;"
                             "3->1 ;"
                             "}";
        plotGraph(pc->Graph(), "deletion-at-end");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GA", pc->Sequence());
        delete pc;
    }
}

TEST(PoaConsensus, TestSimple)
{
    vector<std::string> reads;
    reads += \
             "TTTACAGGATAGTCCAGT",
             "ACAGGATACCCCGTCCAGT",
             "ACAGGATAGTCCAGT",
             "TTTACAGGATAGTCCAGTCCCC",
             "TTTACAGGATTAGTCCAGT",
             "TTTACAGGATTAGGTCCCAGT",
             "TTTACAGGATAGTCCAGT";
    const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
    plotGraph(pc->Graph(), "simple");
    EXPECT_EQ("TTTACAGGATAGTCCAGT", pc->Sequence());
    delete pc;
}

TEST(PoaConsensus, TestOverhangSecond)
{
    vector<std::string> reads;
    reads += \
             "TTTACAGGATAGTCCAGT",
             "TTTACAGGATAGTCCAGTAAA",
             "TTTACAGGATAGTCCAGTAAA";
    const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
    EXPECT_EQ("TTTACAGGATAGTCCAGTAAA", pc->Sequence());
    delete pc;
}

// This won't work until local alignment is properly implemented
TEST(PoaConsensus, Tiling)
{
    vector<std::string> reads;
    reads +=  ("CCGTGGGTTGCGATGGTTCAGAGTTTATTCGCACGGCAGT"
               "GCGCTGTTTGGGGCTGCGCACCGCGCCCCAAGGCGGGCCG"
               "TCACATCATGCCAGGCCCGCTGCCTGGGCAATTCATCAGA"
               "AACCTAGGATTACGACGCGGAAATATAGTCCACGGATTGA"
               "TTCCCACCTAGTGAGTGACTGTCACAGGGATGCAACCGCG"
               "GATGTTCCGG"),
              ("CTGCCTGGGCAATTCATCAGAAACCTAGGATTACGACGCG"
               "GAAATATAGTCCACGGATTGATTCCCACCTAGTGAGTGAC"
               "TGTCACAGGGATGCAACCGCGGATGTTCCGGACGGGTAGT"
               "GGGAAATTACCGCAGCGCCTAGGCACCGCGTTACCTACAA"
               "CCGCTCGCTCATGGCCTTCCGCAGGGCTCCGTACTAAAGC"
               "CTAACCCAAGA"),
              ("AGTGGGAAATTACCGCAGCGCCTAGGCACCGCGTTACCTA"
               "CAACCGCTCGCTCATGGCCTTCCGCAGGGCTCCGTACTAA"
               "AGCCTAACCCAAGACCATTGAAGGGCACACGGCATAGCCG"
               "GCTGCAAGGGGGGGGGTGAGCATGCCTGTATGGTGGCCCC"
               "AGAACTACACAAGGCATGTTCC"),
              ("TAAAGCCTAACCCAAGACCATTGAAGGGCACACGGCATAG"
               "CCGGCTGCAAGGGGGGGGGTGAGCATGCCTGTATGGTGGC"
               "CCCAGAACTACACAAGGCATGTTCCGGTTTCTCGAGGCAT"
               "GAGAACCCCGGTGACTAAACCCGAGTTAGAGTTGAGGCTC"
               "TGTCCCCCAAGAAAGCGCGGCTTGGACTTCTTTACCAACC"
               "CTGCGCCTTTCA");
    string expectedResult = \
                            "CCGTGGGTTGCGATGGTTCAGAGTTTATTCGCACGGCAGT"
                            "GCGCTGTTTGGGGCTGCGCACCGCGCCCCAAGGCGGGCCG"
                            "TCACATCATGCCAGGCCCGCTGCCTGGGCAATTCATCAGA"
                            "AACCTAGGATTACGACGCGGAAATATAGTCCACGGATTGA"
                            "TTCCCACCTAGTGAGTGACTGTCACAGGGATGCAACCGCG"
                            "GATGTTCCGGACGGGTAGTGGGAAATTACCGCAGCGCCTA"
                            "GGCACCGCGTTACCTACAACCGCTCGCTCATGGCCTTCCG"
                            "CAGGGCTCCGTACTAAAGCCTAACCCAAGACCATTGAAGG"
                            "GCACACGGCATAGCCGGCTGCAAGGGGGGGGGTGAGCATG"
                            "CCTGTATGGTGGCCCCAGAACTACACAAGGCATGTTCCGG"
                            "TTTCTCGAGGCATGAGAACCCCGGTGACTAAACCCGAGTT"
                            "AGAGTTGAGGCTCTGTCCCCCAAGAAAGCGCGGCTTGGAC"
                            "TTCTTTACCAACCCTGCGCCTTTCA";

    // change this to local alignment:
    const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::LOCAL_ALIGNMENT);
    //    plotGraph(pc->getGraph(), "tiling", true);
    // EXPECT_EQ(expectedResult, pc->getSequence());
    delete pc;
}


// TODO(dalexander): make a test that behaves differently under
// local and global alignment, and test both
TEST(PoaConsensus, SmallLocalTest)
{
    vector<std::string> reads;
    reads +=  "GGTGG", "GGTGG", "T";
    const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::LOCAL_ALIGNMENT);
    plotGraph(pc->Graph(), "small-local");
    string expectedDot = \
                         "digraph G {"
                         "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
                         "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
                         "2[shape=Mrecord, label=\"{ G | 2 }\"];"
                         "3[shape=Mrecord, label=\"{ G | 2 }\"];"
                         "4[shape=Mrecord, label=\"{ T | 3 }\"];"
                         "5[shape=Mrecord, label=\"{ G | 2 }\"];"
                         "6[shape=Mrecord, label=\"{ G | 2 }\"];"
                         "0->2 ;"
                         "2->3 ;"
                         "3->4 ;"
                         "4->5 ;"
                         "5->6 ;"
                         "6->1 ;"
                         "4->1 ;"
                         "0->4 ;"
                         "}";
    string dot = pc->Graph()->ToGraphViz();
    EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
    EXPECT_EQ("GGTGG", pc->Sequence());
    delete pc;
}

TEST(PoaConsensus, SmallTilingTest)
{
    vector<std::string> reads;
    reads += \
            "GGGGAAAA",
            "AAAATTTT",
            "TTTTCCCC",
            "CCCCAGGA";
    const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::LOCAL_ALIGNMENT);
    plotGraph(pc->Graph(), "small-tiling");
    EXPECT_EQ("GGGGAAAATTTTCCCCAGGA", pc->Sequence());
    delete pc;
}


TEST(PoaConsensus, TestVerboseGraphVizOutput)
{
    vector<std::string> reads;
    reads += "GGG", "TGGG";
    const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);
    string dot = pc->Graph()->ToGraphViz(PoaGraph::COLOR_NODES | PoaGraph::VERBOSE_NODES);

    string expectedDot = \
                        "digraph G {"
                        "0[shape=Mrecord, label=\"{ { 0 | ^ } |{ 0 | 1 } |{ 0.00 | 0.00 } }\"];"
                        "1[shape=Mrecord, label=\"{ { 1 | $ } |{ 0 | 0 } |{ 0.00 | 0.00 } }\"];"
                        "2[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" ,"
                                          " label=\"{ { 2 | G } |{ 2 | 2 } |{ 2.00 | 2.00 } }\"];"
                        "3[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" ,"
                                          " label=\"{ { 3 | G } |{ 2 | 2 } |{ 2.00 | 4.00 } }\"];"
                        "4[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" ,"
                                          " label=\"{ { 4 | G } |{ 2 | 0 } |{ 2.00 | 6.00 } }\"];"
                        "5[shape=Mrecord, label=\"{ { 5 | T } |{ 1 | 0 } |{ -0.00 | -0.00 } }\"];"
                        "0->2 ;"
                        "2->3 ;"
                        "3->4 ;"
                        "4->1 ;"
                        "5->2 ;"
                        "0->5 ;}";

    EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
    delete pc;
}


TEST(PoaConsensus, TestMutations)
{
    using ::testing::ElementsAreArray;

    vector<std::string> reads;
    reads += "TGATTACAT",
             "TGATTACAT",
             "TGATTCAT",    // Deletion @ 5
             "TGATTATAT",   // Substitution @ 6
             "TGATTGACAT";  // Insertion @ 5

    const PoaConsensus* pc = PoaConsensus::FindConsensus(reads, PoaConfig::GLOBAL_ALIGNMENT);

    const std::vector<ScoredMutation>* scoredMutations = pc->Mutations();
    std::vector<string> variantDescriptions;
    foreach (const ScoredMutation& scoredMutation, *scoredMutations)
    {
        variantDescriptions.push_back(scoredMutation.ToString());
    }
    sort(variantDescriptions.begin(), variantDescriptions.end());
    const char* expectedDescriptions[] = { "Deletion @5:6 -3.00",
                                           "Insertion (G) @5 -3.00",
                                           "Substitution (T) @6:7 -3.00" };
    ASSERT_THAT(variantDescriptions, ElementsAreArray(expectedDescriptions));
    delete pc;
}
