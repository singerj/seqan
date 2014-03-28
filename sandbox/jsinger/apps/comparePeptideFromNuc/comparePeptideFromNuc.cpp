// ==========================================================================
//                           comparePeptideFromNuc
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/translation.h>

#include <seqan/arg_parse.h>

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

using namespace seqan;

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The first (and only) argument of the program is stored here.
    seqan::CharString goldTextFile;
    
    StringSet<CharString> compTextFiles;

    AppOptions() :
        verbosity(1)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("comparePeptideFromNuc");
    // Set short description, version, and date.
    setShortDescription(parser, "Put a Short Description Here");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT", true));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBcomparePeptideFromNuc\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
    getArgumentValue(options.goldTextFile, parser, 0);

    resize(options.compTextFiles, getArgumentValueCount(parser, 1));
    for (unsigned i =0 ; i < getArgumentValueCount(parser, 1); ++i)
    {
        getArgumentValue(options.compTextFiles[i], parser, 1, i);
        std::cerr << options.compTextFiles[i] << std::endl;
    }

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "EXAMPLE PROGRAM\n"
              << "===============\n\n";

    SequenceStream streamGold(toCString(options.goldTextFile));
    CharString idGold;
    Dna5String seqGold;
 
    int resRead = readRecord(idGold, seqGold, streamGold);
    if (resRead != 0)
    {
        std::cerr << "ERROR: Could not read " << options.goldTextFile << "!\n";
        return 1;
    }

    String<AminoAcid> textGold;
    translate(textGold, seqGold);

    for (unsigned i = 0; i < length(options.compTextFiles); ++i)
    {
        
        SequenceStream streamComp(toCString(options.compTextFiles[i]));
        CharString idComp;
        Dna5String seqComp;
     
        resRead = readRecord(idComp, seqComp, streamComp);
        if (resRead != 0)
        {
            std::cerr << "ERROR: Could not read " << options.compTextFiles[i] << std::endl;
            return 1;
        }

        String<AminoAcid> textComp;
        translate(textComp, seqComp);

        if (textComp != textGold)
        {
            std::cout << "The files: " << options.goldTextFile << " and " << options.compTextFiles[i] << " are NOT equal." << std::endl;
            std::cout << textGold << std::endl;
            std::cout << "vs" << std::endl;
            std::cout << textComp << std::endl;
        }
        else
        {
            std::cout << "The files: " << options.goldTextFile << " and " << options.compTextFiles[i] << " are equal." << std::endl;
        }

    }


    return 0;
}
