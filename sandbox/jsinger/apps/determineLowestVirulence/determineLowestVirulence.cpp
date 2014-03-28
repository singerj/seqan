// ==========================================================================
//                          determineLowestVirulence
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

#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

#include <lemon/list_graph.h>

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
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The first (and only) argument of the program is stored here.
    seqan::CharString text;

    String<char> codonTable;
    String<char> genomeFile;
    String<char> probFile;

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
    seqan::ArgumentParser parser("determineLowestVirulence");
    // Set short description, version, and date.
    setShortDescription(parser, "Put a Short Description Here");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    // We require one argument.
   // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, seqan::ArgParseOption("ct", "codonTable", "Name of the file the codon table is stored in.", ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("i", "genomeTable", "Name of the genome file.", ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("pf", "probFile", "Name of the probability file.", ArgParseArgument::STRING, "TEXT"));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBdetermineLowestVirulence\\fP \\fB-v\\fP \\fItext\\fP",
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
   // seqan::getArgumentValue(options.text, parser, 0);

    getOptionValue(options.codonTable, parser, "ct");
    getOptionValue(options.genomeFile, parser, "i");
    getOptionValue(options.probFile, parser, "pf");

    return seqan::ArgumentParser::PARSE_OK;
}

struct AminoAcidToDna_
{
    String<StringSet<String<char> > > table;
    AminoAcidToDna_(AppOptions const & options);
};

AminoAcidToDna_::AminoAcidToDna_(AppOptions const & options)
{
    String<char, MMap<> > codonString;
    open(codonString, toCString(options.codonTable));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(codonString);

    String<char> codon;
    AminoAcid aminoAcid;
    resize(table, 24);
    while (!atEnd(iter))
    {
        if (value(iter) == '#')
        {
            for (; value(iter) != '\n'; ++iter) ;
            ++iter;
        }
        else
        {
            clear(codon);
            for (; value(iter) != '\t'; ++iter) 
            {
                appendValue(codon, value(iter));
            }
            ++iter;

            aminoAcid = value(iter);
            ++iter;
            ++iter;

            appendValue(table[ordValue(aminoAcid)], codon);
        }
    }

    for (unsigned i = 0; i < length(table); ++i)
    {
        std::cout << AminoAcid(i) << '\t';
        for (unsigned j = 0; j < length(table[i]); ++j)
            std::cout << toCString(table[i][j]) << '\t';
        std::cout << std::endl;
    }
}

template <typename TValue>
unsigned hashCodon(String<TValue> const & codon)
{
    return ordValue(Dna(codon[2])) + 4 * ordValue(Dna(codon[1])) + 16 * ordValue(Dna(codon[0]));
}

bool initTransitionProb(double (&transitionProb)[64][64], AppOptions const & options)
{
    String<char, MMap<> > inString;
    if (!open(inString, toCString(options.probFile)))
        return false;
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(inString);

    String<char> temp;
    String<Dna> tempDna;
    StringSet<String<Dna> > codonOrder;
    while (!atEnd(iter))
    {
        clear(temp);
        clear(tempDna);
        if (value(iter) == '#')
        {
            for (; value(iter) != '\n'; ++iter) ;
            ++iter;
        }
        else if (value(iter) == ';')
        {
            ++iter;
            for (; value(iter) != '\n'; ++iter)
            {
                if(value(iter) == ';')
                {
                    appendValue(codonOrder, tempDna);
                    clear(tempDna);
                }
                else
                    appendValue(tempDna, value(iter));
            }
            appendValue(codonOrder, tempDna);
            ++iter;
        }
        else
        {
            appendValue(tempDna, value(iter));
            ++iter;
            appendValue(tempDna, value(iter));
            ++iter;
            appendValue(tempDna, value(iter));
            ++iter;

            // ;
            ++iter;
            unsigned codonPos = hashCodon(tempDna);
            for (unsigned i = 0; i < 64; ++i)
            {
                clear(temp);
                for (; value(iter) != ';' && value(iter) != '\n'; ++iter)
                {
                    appendValue(temp, value(iter));
                }
                ++iter;
                lexicalCast2(transitionProb[codonPos][hashCodon(codonOrder[i])], temp);
            }
        }
    }
    return true;
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

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
        std::cout << "__OPTIONS____________________________________________________________________\n";

    // reading and creating amino acid table
    AminoAcidToDna_ _table(options);

    typedef double TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef String<String<char> > TProperties;

    // Reading genome sequence.
    TProperties codons;

    SequenceStream seqIO(toCString(options.genomeFile));
    CharString id;
    Dna5String seq;
 
    int resRead = readRecord(id, seq, seqIO);
    if (resRead != 0)
    {
        std::cerr << "ERROR: Could not read record!\n";
        return 1;
    }

    String<AminoAcid> text;
    translate(text, seq);

    // reading codon transition probabilities
    double transitionProb[64][64];
    initTransitionProb(transitionProb, options);

    // creating the graph
    StringSet<String<TVertexDescriptor> > vertices;
    resize(vertices, length(text) + 2);
    resize(vertices[0], 1);
    resize(vertices[length(vertices) - 1], 1);

    TGraph g;

    vertices[0][0] = addVertex(g);
    vertices[length(vertices) - 1][0] = addVertex(g);

    for (unsigned i = 0; i < length(text); ++i)
    {
        resize(vertices[i + 1], length(_table.table[ordValue(text[i])]));
        for (unsigned j = 0; j < length(vertices[i + 1]); ++j)
            vertices[i + 1][j] = addVertex(g);
    }


    resizeVertexMap(g, codons);
    assignProperty(codons, vertices[0][0], "start");
    assignProperty(codons, vertices[length(vertices) - 1][0], "end");
    for (unsigned i = 0; i < length(text); ++i)
    {
        for (unsigned j = 0; j < length(vertices[i + 1]); ++j)
        {
            assignProperty(codons, vertices[i + 1][j], _table.table[ordValue(text[i])][j]);
            for (unsigned k = 0; k < length(vertices[i]); ++k)
            {
                if (i == 0)
                    addEdge(g, vertices[i][k], vertices[i+1][j], 0.0);
                else
                {
                    addEdge(g, vertices[i][k], vertices[i+1][j], transitionProb[hashCodon(getProperty(codons, vertices[i][k]))][hashCodon(getProperty(codons,vertices[i+1][j]))]);
                }
            }
        }
    }

    for (unsigned i = 0; i < length(vertices[length(text)]) ; ++i)
        addEdge(g, vertices[length(text)][i], vertices[length(text) + 1][0], 0.0);

    // Run Dijkstra's algorithm from vertex 0.
    String<unsigned> predMap;
    String<double> distMap;
    InternalMap<TCargo> cargoMap;
    bellmanFordAlgorithm(g, 0, cargoMap, predMap, distMap);

    // Print results to stdout.
    std::cout << "Single-Source Shortest Paths: \n";
    _printPath(g,predMap,(TVertexDescriptor) 0, vertices[length(vertices) - 1][0], codons);
    std::cout << std::endl;

    std::cout << "The score is: " << distMap[1] / (double)(length(text) - 1) << std::endl;

    return 0;
}
