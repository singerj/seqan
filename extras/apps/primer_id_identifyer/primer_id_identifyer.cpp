// ==========================================================================
//                             primerIdIdentifyer
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@bsse.ethz.ch>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/modifier.h>

#include "ednafull.h"

#include <iostream>
#include <fstream>
#include <algorithm>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // The first (and only) argument of the program is stored here.
    StringSet<CharString> patternFileNames;
    CharString refFileName;
    CharString outputFileName;
    int numThreads;
    unsigned bufferSize;

    AppOptions() :
        numThreads(1),
        bufferSize(100000)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("primerIdIdentifyer");
    // Set short description, version, and date.
    setShortDescription(parser, "Put a Short Description Here");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    ArgParseArgument refArgument(ArgParseArgument::INPUTFILE, "IN");
    setValidValues(refArgument, "fasta fa fna");
    addArgument(parser, refArgument);

    ArgParseArgument patternArg(ArgParseArgument::INPUTFILE, "IN", true);
    setValidValues(patternArg, "fasta fa fna fastq fq fnq");
    addArgument(parser, patternArg);

    addOption(parser, ArgParseOption("o", "outputFileName", "Name of the output file.", ArgParseArgument::OUTPUTFILE, "OUT"));
    setRequired(parser, "o");
    setValidValues(parser, "o", "sam");

    addOption(parser, ArgParseOption("th", "threads", "The number of threads to be used.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("bf", "bufferSize", "The number of hits stored in a buffer before writing them to disk.", ArgParseArgument::INTEGER));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBsuffix_tree_mapping\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getArgumentValue(options.refFileName, parser, 0);
    resize(options.patternFileNames, getArgumentValueCount(parser, 1));
    for (unsigned i =0 ; i < getArgumentValueCount(parser, 1); ++i)
    {
        getArgumentValue(options.patternFileNames[i], parser, 1, i);
    }
    getOptionValue(options.outputFileName, parser, "o");
    getOptionValue(options.numThreads, parser, "th");
    getOptionValue(options.bufferSize, parser, "bf");

    return ArgumentParser::PARSE_OK;
}

template <typename TStream, typename TAlign>
void write(TStream & stream,
           TAlign const & align,
           String<char> const & patternId,
           String<char> const & refId,
           String<char> & cigar,
           bool const reverseComplemented)
{
    getCigarString(cigar,row(align, 0), row(align, 1), 1000);

    stream << patternId << "\t";                                //QNAME
    stream << (unsigned short)0x0002 << "\t";                        //FLAG
    stream << refId << "\t";   //RNAME
    if (reverseComplemented)
        stream << 2 * length(source(row(align,1))) / 5;//  + row(align, 1)._array[0] + 1 << "\t";    //POS
    else
        stream << row(align, 1)._array[0] + 1 << "\t";    //POS
    stream << 0 << "\t";                                             //MAPQ
    stream << cigar << "\t";                                         //CIGAR
    stream << "*" << "\t";                                           //RNEXT
    stream << 0 << "\t";                                             //PNEXT
    stream << 0 << "\t";                                             //TLen
    stream << source(row(align, 1)) << "\t";                                //SEQ
    stream << "*" << "\n";                                           //QUAL
    //std::cerr << stream.rdbuf() << std::endl;
    std::cerr << patternId << "\t";                                //QNAME
    std::cerr << (unsigned short)0x0002 << "\t";                        //FLAG
    std::cerr << refId << "\t";   //RNAME
    if (reverseComplemented)
        std::cerr << 2 * length(source(row(align,0))) / 5  + row(align, 1)._array[0] + 1 << "\t";    //POS
    else
        std::cerr << row(align, 1)._array[0] + 1 << "\t";    //POS
    std::cerr << 0 << "\t";                                             //MAPQ
    std::cerr << cigar << "\t";                                         //CIGAR
    std::cerr << "*" << "\t";                                           //RNEXT
    std::cerr << 0 << "\t";                                             //PNEXT
    std::cerr << 0 << "\t";                                             //TLen
    std::cerr << source(row(align, 1)) << "\t";                                //SEQ
    std::cerr << "*" << "\n";                                           //QUAL
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    typedef String<Dna5> TRefSeq;

    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "primerIdIdentifyer\n"
              << "==================\n\n";

    // Preparation of out stream including SAM header
    std::ofstream stream(toCString(options.outputFileName), std::ios::out | std::ios::app);
    stream << "@HD\tVN:1.5\tSO:coordinate\n";
    // Reference
    SequenceStream seqInRef(toCString(options.refFileName));
    StringSet<String<char> > refIds;
    StringSet<TRefSeq> refSeqsLeft;
    readAll(refIds, refSeqsLeft, seqInRef);

    StringSet<TRefSeq> refSeqsRight;
    resize(refSeqsRight, length(refSeqsLeft));
    for (unsigned i = 0; i < length(refSeqsLeft); ++i)
    {
        stream << "@SQ\tSN:" << refIds[i] << "\tLN:" << length(refSeqsLeft[i]) << "\n";
        refSeqsRight[i] = suffix(refSeqsLeft[i], 2 * length(refSeqsLeft[i]) / 5);
        resize(refSeqsLeft[i], 3 * length(refSeqsLeft[i]) / 5);
    }

    // ReverseComplement References
    StringSet<TRefSeq> refSeqsLeftComp;
    StringSet<TRefSeq> refSeqsRightComp;
    resize(refSeqsLeftComp, length(refSeqsLeft));
    resize(refSeqsRightComp, length(refSeqsRight));
    for (unsigned i = 0; i < length(refSeqsLeft); ++i)
    {
        refSeqsLeftComp[i] = refSeqsLeft[i];
        reverseComplement(refSeqsLeftComp[i]);
        refSeqsRightComp[i] = refSeqsRight[i];
        reverseComplement(refSeqsRightComp[i]);
    }

    // Preparing the align Objects
    // Left side
    String<Align<TRefSeq> > alignObjsLeft;
    resize(alignObjsLeft, length(refSeqsLeft));
    for (unsigned i = 0; i < length(refSeqsLeft); ++i)
    {
        resize(rows(alignObjsLeft[i]), 2);
        assignSource(row(alignObjsLeft[i], 0), refSeqsLeft[i]);
    }
    /*
    String<Align<TRefSeq> > alignObjsLeftComp;
    resize(alignObjsLeftComp, length(refSeqsLeftComp));
    for (unsigned i = 0; i < length(refSeqsLeftComp); ++i)
    {
        resize(rows(alignObjsLeftComp[i]), 2);
        assignSource(row(alignObjsLeftComp[i], 0), refSeqsLeftComp[i]);
    }
    */
    // Right side
    String<Align<TRefSeq> > alignObjsRight;
    resize(alignObjsRight, length(refSeqsRight));
    for (unsigned i = 0; i < length(refSeqsRight); ++i)
    {
        resize(rows(alignObjsRight[i]), 2);
        assignSource(row(alignObjsRight[i], 0), refSeqsRight[i]);
    }
    /*
    String<Align<TRefSeq> > alignObjsRightComp;
    resize(alignObjsRightComp, length(refSeqsRightComp));
    for (unsigned i = 0; i < length(refSeqsRightComp); ++i)
    {
        resize(rows(alignObjsRightComp[i]), 2);
        assignSource(row(alignObjsRightComp[i], 0), refSeqsRightComp[i]);
    }
    */

    // Positioning Segment
    TRefSeq posSegmentLeft = prefix(refSeqsLeft[0], 20);
    TRefSeq posSegmentLeftComp = posSegmentLeft;
    reverseComplement(posSegmentLeftComp);
    //TRefSeq posSegmentLeftComp = suffix(refSeqsLeftComp[0], length(refSeqsLeftComp[0]) - 20);
    TRefSeq posSegmentRight = suffix(refSeqsRight[0], length(refSeqsRight[0]) - 20);
    TRefSeq posSegmentRightComp = posSegmentRight;
    reverseComplement(posSegmentRightComp);
    //TRefSeq posSegmentRightComp = prefix(refSeqsRightComp[0], 20);


    // Pattern
    StringSet<CharString> patternIds;
    StringSet<String<Dna5> > patternSeqs;

    // Cigar String
    CharString cigar;

    // Alignment helpers
    int maxScore, maxId;
    Ednafull scoringScheme(-2, -1);
    AlignConfig<true, false, false, true> alignConfig;

    for (unsigned fileCounter = 0; fileCounter < length(options.patternFileNames); ++fileCounter)
    {
        clear(patternIds);
        clear(patternSeqs);

        // Reading the pattern
        SequenceStream seqInPattern(toCString(options.patternFileNames[fileCounter]));
        SEQAN_OMP_PRAGMA(critical (read_chunk))
        {
            if (!atEnd(seqInPattern))
            {
                if(readBatch(patternIds, patternSeqs, seqInPattern, options.bufferSize) != 0)
                {
                    std::cout << "ERROR: Could not read samples!\n";
                }
            }
        }

       // int lDiag = -30;
       // int uDiag = 30;
        for (unsigned i = 0; i < length(patternSeqs); ++i)
        {

            // Determine if right or left aligned
            String<Dna5> left= prefix(patternSeqs[i], 20);
            String<Dna5> right = suffix(patternSeqs[i], length(patternSeqs[i]) - 20);

            maxId = 0;
            maxScore = globalAlignmentScore(posSegmentLeft, left, MyersBitVector());
            //std::cerr << "posSegmentLeft: " << posSegmentLeft << std::endl;
            //std::cerr << "posSegmentLeft: " << posSegmentLeftComp << std::endl;
            //std::cerr << "posSegmentLeft: " << posSegmentRight << std::endl;
            //std::cerr << "posSegmentLeft: " << posSegmentRightComp << std::endl;
            int tempScore = globalAlignmentScore(posSegmentLeftComp, right, MyersBitVector());
            if (tempScore > maxScore)
            {
                maxScore = tempScore;
                maxId = 1;
            }
            tempScore = globalAlignmentScore(posSegmentRight, right, MyersBitVector());
            if (tempScore > maxScore)
            {
                maxScore = tempScore;
                maxId = 2;
            }
            tempScore = globalAlignmentScore(posSegmentRightComp, left, MyersBitVector());
            if (tempScore > maxScore)
            {
                maxScore = tempScore;
                maxId = 3;
            }


            switch (maxId) 
            {
                case(0) :
                {
                    maxScore = minValue<int>();
                    maxId = -1;
                    for (unsigned j = 0; j < length(alignObjsLeft); ++j)
                    {
                        assignSource(row(alignObjsLeft[j], 1), patternSeqs[i]);
                        int result = globalAlignment(alignObjsLeft[j], scoringScheme, alignConfig);//, lDiag, uDiag);
                        if (maxScore < result)
                        {
                            maxScore = result;
                            maxId = j;
                        }
                    }
                    write(stream, alignObjsLeft[maxId], patternIds[i], refIds[maxId], cigar, false);
                    break;
                }
                case(1) :
                {
                    maxScore = minValue<int>();
                    maxId = -1;
                    for (unsigned j = 0; j < length(alignObjsLeft); ++j)
                    {
                        reverseComplement(patternSeqs[i]);
                        assignSource(row(alignObjsLeft[j], 1), patternSeqs[i]);
                        int result = globalAlignment(alignObjsLeft[j], scoringScheme, alignConfig);//, lDiag, uDiag);
                        if (maxScore < result)
                        {
                            maxScore = result;
                            maxId = j;
                        }
                    }
                    write(stream, alignObjsLeft[maxId], patternIds[i], refIds[maxId], cigar, true);
                    break;
                }
                case(2) :
                {
                    //std::cerr << source(row(alignObjsRight[j], 0)) << std::endl;
                        //char c;
                        //std::cin >> c;
                        //std::cerr << std::endl;
                    maxScore = minValue<int>();
                    maxId = -1;
                    for (unsigned j = 0; j < length(alignObjsRight); ++j)
                    {
                        assignSource(row(alignObjsRight[j], 1), patternSeqs[i]);
                        int result = globalAlignment(alignObjsRight[j], scoringScheme, alignConfig);//, lDiag, uDiag);
                        if (maxScore < result)
                        {
                            maxScore = result;
                            maxId = j;
                        }
                    }
                    write(stream, alignObjsRight[maxId], patternIds[i], refIds[maxId], cigar, false);
                    break;
                }
                case(3) :
                {
                    maxScore = minValue<int>();
                    maxId = -1;
                    for (unsigned j = 0; j < length(alignObjsRight); ++j)
                    {
                        reverseComplement(patternSeqs[i]);
                        assignSource(row(alignObjsRight[j], 1), patternSeqs[i]);
                        int result = globalAlignment(alignObjsRight[j], scoringScheme, alignConfig);//, lDiag, uDiag);
                        if (maxScore < result)
                        {
                            maxScore = result;
                            maxId = j;
                        }
                    }
                    std::cerr << source(row(alignObjsRight[maxId], 0)) << std::endl;
                    std::cerr << patternSeqs[i] << std::endl;
                    std::cerr << maxScore << std::endl << alignObjsRight[maxId] << std::endl;
                    write(stream, alignObjsRight[maxId], patternIds[i], refIds[maxId], cigar, true);
                    char c;
                    //std::cin >> c;
                    std::cerr << std::endl;
                    break;
                }
                default :
                    std::cerr << "A mistake must have happend!" << std::endl;
            }
        }
        stream.close();
    }

    /*
    Dna5String seqH = "CGAAATT";
    Dna5String seqV = "AAATT";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seqH);
    assignSource(row(align, 1), seqV);

    Score<int, Simple> scoringScheme(2, -1, -2);
    AlignConfig<> alignConfig;

    int lDiag = -2;
    int uDiag = 2;

    int result = globalAlignment(align, scoringScheme, alignConfig, lDiag, uDiag);

    std::cout << "Score: " << result << "\n";
    std::cout << "The resulting alignment is\n"
              << align << "\n";

    //WRITING
    record.flag = 0;


    record._qID = MaxValue<__uint32>::VALUE;
    record.rID = BamAlignmentRecord::INVALID_REFID;
    record.beginPos = BamAlignmentRecord::INVALID_POS;
    record.mapQ = 255;
    record.bin = 0;
    clear(record.cigar);
    record.rNextId = BamAlignmentRecord::INVALID_REFID;
    record.pNext = BamAlignmentRecord::INVALID_POS;
    record.tLen = BamAlignmentRecord::INVALID_LEN;
    clear(record.seq);
    clear(record.qual);
    clear(record.tags);

    std::cout << length(row(align, 0)) << std::endl;
    for (int i =0; i < length(row(align, 0)); ++i)
        std::cout << isGap(row(align, 0), i);
    std::cout << std::endl;

    for (int i =0; i < length(row(align, 0)._array); ++i)
        std::cout << row(align, 0)._array[i] << " ";
    std::cout << std::endl;

    for (int i =0; i < length(row(align, 1)._array); ++i)
        std::cout << row(align, 1)._array[i] << " ";
    std::cout << std::endl;
*/


    return 0;
}

