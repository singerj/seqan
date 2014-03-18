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
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

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
    seqan::getArgumentValue(options.text, parser, 0);

    return seqan::ArgumentParser::PARSE_OK;
}

struct AminoAcidToDna_
{
    String<StringSet<String<char> > > table;
    AminoAcidToDna_();
};

AminoAcidToDna_::AminoAcidToDna_()
{
    resize(table, 2);
    resize(table[0], 4);
    table[0][0] = "gct";
    table[0][1] = "gcc"; 
    table[0][2] = "gcg";
    table[0][3] = "gct";
    resize(table[1], 2);
    table[1][0] = "aga";
    table[1][1] = "agg";
}


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    /*seqan::ArgumentParser parser;
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
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "TEXT     \t" << options.text << "\n\n";
    }*/

    AminoAcidToDna_ _table;

    std::cerr << toCString(_table.table[0][0]) << std::endl;
    std::cerr << toCString(_table.table[1][1]) << std::endl;

    typedef double TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef String<char> TCityName;
    typedef String<TCityName> TProperties;

    TProperties codons;

    String<AminoAcid> text = "AR";

    StringSet<String<TVertexDescriptor> > vertices;
    resize(vertices, length(text) + 2);
    resize(vertices[0], 1);
    resize(vertices[length(vertices) - 1], 1);

    TGraph g;

    vertices[0][0] = addVertex(g);
    vertices[length(vertices) - 1] = addVertex(g);

    for (unsigned i = 0; i < length(text); ++i)
    {
        resize(vertices[i + 1], length(_table.table[ordValue(text[i])]));
        for (unsigned j = 0; j < length(vertices[i + 1]); ++j)
        {
            vertices[i + 1][j] = addVertex(g);
            for (unsigned k = 0; k < length(vertices[i]); ++k)
                addEdge(g, vertices[i][k], vertices[i+1][j], 1.0);
        }
    }

    for (unsigned i = 0; i < length(vertices[length(text)]) ; ++i)
        addEdge(g, vertices[length(text)][i], vertices[length(text) + 1][0], 1.0);

    resizeVertexMap(g, codons);
    assignProperty(codons, vertices[0][0], "start");
    assignProperty(codons, vertices[length(vertices) - 1][0], "end");
    for (unsigned i = 0; i < length(text); ++i)
    {
        for (unsigned j = 0; j < length(vertices[i + 1]); ++j)
        {
            //vertices[i + 1][j] = addVertex(g);
            assignProperty(codons, vertices[i + 1][j], _table.table[ordValue(text[i])][j]);
        }
    }

    // Run Dijkstra's algorithm from vertex 0.
    String<unsigned> predMap;
    String<unsigned> distMap;
    InternalMap<TCargo> cargoMap;
    dijkstra(g, 0, cargoMap, predMap, distMap);

    // Print results to stdout.
    std::cout << "Single-Source Shortest Paths: \n";
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Path from 0 to " << getValue(it) << ": ";
        _printPath(g,predMap,(TVertexDescriptor) 0, getValue(it));
        std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")\n";
        goNext(it);
    }
    


    /*TVertexDescriptor vertHamburg = addVertex(g);
    TVertexDescriptor vertHannover = addVertex(g);
    TVertexDescriptor vertMuenchen = addVertex(g);
    addEdge(g, vertBerlin, vertHamburg, 289.0);
    addEdge(g, vertHamburg, vertHannover, 289.0);
    addEdge(g, vertHannover, vertMuenchen, 572.0);
    */
    /*FILE* strmWrite = fopen("graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);
    */
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;

   // TVertexIterator itV(g);
    std::cout << g << std::endl;

    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator itV(g);
    for(;!atEnd(itV);goNext(itV)) {
        std::cout << value(itV) << ':' << getProperty(codons, value(itV)) << std::endl;
    }


    lemon::ListDigraph ga;
    lemon::ListDigraph::Node u = ga.addNode();
    lemon::ListDigraph::Node v = ga.addNode();
    lemon::ListDigraph::Arc  a = ga.addArc(u, v);
    std::cout << "Hello World! This is LEMON library here." << std::endl;
    std::cout << "We have a directed graph with " << countNodes(ga) << " nodes "
       << "and " << countArcs(ga) << " arc." << std::endl;

    return 0;
}
