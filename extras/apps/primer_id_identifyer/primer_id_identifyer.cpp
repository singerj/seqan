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
    int minScore;

    AppOptions() :
        numThreads(1),
        bufferSize(100000),
        minScore(minValue<int>())
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
    addOption(parser, ArgParseOption("s", "minScore", "MinScore of alignment in order to be considered.", ArgParseArgument::INTEGER));

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
    getOptionValue(options.minScore, parser, "s");

    return ArgumentParser::PARSE_OK;
}

template <typename TStream, typename TAlign>
void writeLocal(TStream & stream,
           TAlign & align,
           String<char> const & patternId,
           String<char> const & refId,
           String<char> & cigar,
           int score)
{
    int clipBegin = row(align, 1)._array[0];
    int clipEnd = length(row(align, 1)) - row(align, 1)._array[length(row(align, 1)._array) - 1];
    setClippedBeginPosition(row(align, 0), clipBegin);
    setClippedBeginPosition(row(align, 1), clipBegin);
    setClippedEndPosition(row(align, 0), clipEnd);
    setClippedEndPosition(row(align, 1), clipEnd);
    getCigarString(cigar,row(align, 0), row(align, 1), 1000);
    stream << patternId << "\t";                                //QNAME
    stream << (unsigned short)0x0002 << "\t";                        //FLAG
    stream << refId << "\t";   //RNAME
    stream << row(align, 1)._array[0] + 1 << "\t";    //POS
    stream << score << "\t";
    stream << cigar << "\t";                                         //CIGAR
    stream << "*" << "\t";                                           //RNEXT
    stream << 0 << "\t";                                             //PNEXT
    stream << 0 << "\t";                                             //TLen
    stream << source(row(align, 1)) << "\t";                                //SEQ
    stream << "*" << "\n";                                           //QUAL

}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    typedef String<Iupac> TRefSeq;

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
    StringSet<TRefSeq> refSeqs;
    readAll(refIds, refSeqs, seqInRef);

    for (unsigned i = 0; i < length(refSeqs); ++i)
        stream << "@SQ\tSN:" << refIds[i] << "\tLN:" << length(refSeqs[i]) << "\n";
    stream.close();

    // Preparing the align Objects
    String<Align<TRefSeq> > alignObjs, alignObjsRevComp;
    resize(alignObjs, length(refSeqs));
    resize(alignObjsRevComp, length(refSeqs));
    for (unsigned i = 0; i < length(refSeqs); ++i)
    {
        resize(rows(alignObjs[i]), 2);
        resize(rows(alignObjsRevComp[i]), 2);
        assignSource(row(alignObjs[i], 0), refSeqs[i]);
        assignSource(row(alignObjsRevComp[i], 0), refSeqs[i]);
    }

    

    // Alignment helpers
    int maxScore, maxScoreRevComp, maxId, maxIdRevComp;
    Ednafull scoringScheme(-1, -20);
    AlignConfig<true, false, false, true> alignConfig;

    for (unsigned fileCounter = 0; fileCounter < length(options.patternFileNames); ++fileCounter)
    {
        SequenceStream seqInPattern(toCString(options.patternFileNames[fileCounter]));
        // This is the parallelization starting point
        SEQAN_OMP_PRAGMA(parallel num_threads(options.numThreads))
        for (; ;)
        {
            // Pattern
            StringSet<CharString> patternIds;
            StringSet<String<Dna5> > patternSeqs;
            String<Dna5> revComPatternSeq;

            // Cigar String
            CharString cigar;

            std::ostringstream localStream;
            // Reading the pattern
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

            if (length(patternIds) == 0)
            {
                std::cerr << "TEST" << std::endl;
                break;
            }

           // int lDiag = -30;
           // int uDiag = 30;
            for (unsigned i = 0; i < length(patternSeqs); ++i)
            {
                // Determine if right or left aligned
                String<Dna5> left= prefix(patternSeqs[i], 20);
                String<Dna5> right = suffix(patternSeqs[i], length(patternSeqs[i]) - 20);

                maxScore = minValue<int>();
                maxId = -1;
                for (unsigned j = 0; j < length(alignObjs); ++j)
                {
                    assignSource(row(alignObjs[j], 1), patternSeqs[i]);
                    int result = globalAlignment(alignObjs[j], scoringScheme, alignConfig, maxScore);//, lDiag, uDiag);
                    if (maxScore < result)
                    {
                        maxScore = result;
                        maxId = j;
                    }
                }
                revComPatternSeq = patternSeqs[i];
                reverseComplement(revComPatternSeq);
                maxScoreRevComp = minValue<int>();
                maxIdRevComp = -1;
                for (unsigned j = 0; j < length(alignObjsRevComp); ++j)
                {
                    assignSource(row(alignObjsRevComp[j], 1), revComPatternSeq);
                    int result = globalAlignment(alignObjsRevComp[j], scoringScheme, alignConfig, maxScoreRevComp);//, lDiag, uDiag);
                    if (maxScoreRevComp < result)
                    {
                        maxScoreRevComp = result;
                        maxIdRevComp = j;
                    }
                }
                if (maxScore > maxScoreRevComp)
                {
                    int normScore = maxScore / (length(revComPatternSeq) * 5.0) * 255.0;
                    if (normScore > options.minScore)
                        writeLocal(localStream, alignObjs[maxId], patternIds[i], refIds[maxId], cigar, normScore);
                }
                else
                {
                    int normScore = maxScoreRevComp / (length(revComPatternSeq) * 5.0) * 255.0;
                    if(normScore > options.minScore)
                        writeLocal(localStream, alignObjsRevComp[maxIdRevComp], patternIds[i], refIds[maxIdRevComp], cigar, normScore);
                }
            }
            SEQAN_OMP_PRAGMA(critical (write_chunk))
            {
                std::ofstream stream(toCString(options.outputFileName), std::ios::out | std::ios::app);
                stream << localStream.str();
                stream.close();
            }
        }
    }

    return 0;
}

