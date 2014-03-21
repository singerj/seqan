// ==========================================================================
//                          create_blastp_benchmark
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
#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/file.h>
#include <seqan/seq_io.h>
#include <seqan/random.h>
#include <seqan/score.h>
#include <seqan/arg_parse.h>

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

using namespace seqan;

struct AppOptions
{
    // name of database
    String<char> dbFileName;
    // number of sequences in the database
    unsigned numSeqDB;

    // min length of sequences
    unsigned minLengthSeqDB;

    // max length sequences
    unsigned maxLengthSeqDB;

    // name of read file
    String<char> readsFileName;
    // number of sequences in the reads
    unsigned numSeqReads;

    // min length of reads
    unsigned minLengthSeqReads;

    // max length reads
    unsigned maxLengthSeqReads;

    // the error rate
    double errorRate;

    // name of the file to store the alignment qualities
    String<char> alignmentQuality;
    
    // name of the file to store the start position
    String<char> startPosition;

    // name of the file to store the frequency control    
    String<char> frequencyControl;

    AppOptions() :
        dbFileName("db.fasta"),
        numSeqDB(100000),
        minLengthSeqDB(30),
        maxLengthSeqDB(5000),
        readsFileName("reads.fasta"),
        numSeqReads(100000),
        minLengthSeqReads(50),
        maxLengthSeqReads(50),
        errorRate(5.0),
        alignmentQuality("alignmentQualities.txt"),
        startPosition("startPosition.txt"),
        frequencyControl("frequencyControl.txt")
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
    seqan::ArgumentParser parser("create_blastp_benchmark");
    // Set short description, version, and date.
    setShortDescription(parser, "This tool creates a database and reads for a BlastP benchmark.");
    setVersion(parser, "0.1");
    setDate(parser, "Feb 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This app creates a database and reads for a blastp benchmark.");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("nD", "numSeqDB", "Set the number of sequences in the db.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("miD", "minLengthSeqDB", "Set the minimum length of a sequence in the db.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("maD", "maxLengthSeqDB", "Set the maximum length of a sequence in the db.", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption("nR", "numSeqReads", "Set the number of reads.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("miR", "minLengthSeqReads", "Set the minimum length of reads.", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("maR", "maxLengthSeqReads", "Set the maximum length of reads.", ArgParseArgument::INTEGER, "INT"));
    
    addOption(parser, seqan::ArgParseOption("er", "errorRate", "Set the error rate of the reads.", ArgParseArgument::DOUBLE, "ERROR"));
    
    addOption(parser, seqan::ArgParseOption("aq", "alignmentQuality", "Name of the file the alignment qualities should be stored in.", ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("st", "startPosition", "Name of the file the alignment start positions should be stored in.", ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("fc", "frequencyControl", "Name of the file the frequency control should be stored in.", ArgParseArgument::STRING, "TEXT"));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBcreate_blastp_benchmark\\fP \\fB-nd 10000\\fP \\fB-miD 100\\fP \\fB-maD 3000\\fP \\fBreads.fasta\\fP \\fB-nR 1000\\fP \\fB-miR 30\\fP \\fB-maR 100\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.dbFileName, parser, 0); 
    getOptionValue(options.numSeqDB, parser, "numSeqDB");
    getOptionValue(options.minLengthSeqDB, parser, "minLengthSeqDB");
    getOptionValue(options.maxLengthSeqDB, parser, "maxLengthSeqDB");

    getArgumentValue(options.readsFileName, parser, 1); 
    getOptionValue(options.numSeqReads, parser, "numSeqDB");
    getOptionValue(options.minLengthSeqReads, parser, "minLengthSeqReads");
    getOptionValue(options.maxLengthSeqReads, parser, "maxLengthSeqReads");

    getOptionValue(options.errorRate, parser, "errorRate");

    getOptionValue(options.alignmentQuality, parser, "alignmentQuality");
    getOptionValue(options.startPosition, parser, "startPosition");
    getOptionValue(options.frequencyControl, parser, "frequencyControl");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

static String<float> freq;
static StringSet<String<double> > dist;

void initFreq()
{
    // Frequencies according to: http://web.expasy.org/docs/relnotes/relstat.html (20.03.2014)

    appendValue(freq, 0);
    appendValue(freq, 8.25);
    appendValue(freq, 13.78);
    appendValue(freq, 17.84);
    appendValue(freq, 23.29);
    appendValue(freq, 24.66);
    appendValue(freq, 28.59);
    appendValue(freq, 35.34);
    appendValue(freq, 42.41);
    appendValue(freq, 44.68);
    appendValue(freq, 50.63);
    appendValue(freq, 60.29);
    appendValue(freq, 66.13);
    appendValue(freq, 68.54);
    appendValue(freq, 72.4);
    appendValue(freq, 77.11);
    appendValue(freq, 83.68);
    appendValue(freq, 89.02);
    appendValue(freq, 90.1);
    appendValue(freq, 93.02);
    appendValue(freq, 99.88);
}

template <typename TValue>
void initDist(TValue const & scoreMatrix)
{
    resize(dist, 20);
    for (unsigned i = 0; i < 20; ++i)
    {
        resize(dist[i], 21);
        dist[i][0] = 0;
    }

    for (unsigned i = 0; i < 20; ++i)
    {
        int currentScore = score(scoreMatrix, i, i);
        for (unsigned j = 0; j < 20; ++j)
        {
            dist[i][j + 1] = std::abs(score(scoreMatrix, i , j) - currentScore);
            if (dist[i][j + 1] == 0)
                dist[i][j + 1] += dist[i][j];
            else
                dist[i][j + 1] = dist[i][j] + 1 / dist[i][j + 1];
        }
    }

}

AminoAcid getAminoAcidDB(double rng)
{
    for (unsigned i = 1; i < 20; ++i)
        if (freq[i] <= rng)
            continue;
        else
            return i - 1;

    return 19;
}

AminoAcid getAminoAcidRead(AminoAcid as)
{
    Rng<MersenneTwister> rng(42);
    Pdf<Uniform<double> > uniformDouble(0, dist[(unsigned)as][19]);
    double newAmino = pickRandomNumber(rng, uniformDouble);

    for (unsigned i = 1; i < 20; ++i)
        if (dist[i] <= newAmino)
            continue;
        else
            return i - 1;

    return 19;
}


void createDatabase(AppOptions const & options){

    Rng<MersenneTwister> rng(42);

    String<char, MMap<> > db;
    open(db, toCString(options.dbFileName));

    String<unsigned> dbFreq;
    resize(dbFreq, 20, 0);

    AminoAcid temp;
    for (unsigned i = 0; i < options.numSeqDB; ++i)
    {
        appendValue(db, '>');
        append(db, "a\n");

        Pdf<Uniform<int> > uniformInt(options.minLengthSeqDB, options.maxLengthSeqDB - 1);
        unsigned entryLength = pickRandomNumber(rng, uniformInt);

        Pdf<Uniform<double> > uniformDouble(0, 100.09);
        for (unsigned j = 0; j < entryLength; ++j)
        {
            temp = getAminoAcidDB(pickRandomNumber(rng, uniformDouble));
            appendValue(db, temp);
            ++dbFreq[(unsigned)temp];
        }
        append(db, "\n");
    }

    std::ofstream outStream;
    outStream.open(toCString(options.frequencyControl), std::ios::out);
    outStream << "#This file shows for every amino acid its frequency in the db and the expected one." << std::endl;
    for (unsigned i = 0; i < length(dbFreq); ++i)
        outStream << AminoAcid(i) << "\t" << static_cast<double>(dbFreq[i]) / static_cast<double>(lengthSum(db)) << "\t" << freq[i+1] - freq[i] << "\n";
}

void createReads(AppOptions const & options){

    Rng<MersenneTwister> rng(42);

    SequenceStream seqStream(toCString(options.dbFileName));
    StringSet<String<char> > ids;
    StringSet<String<char> > seqs;

    readAll(ids, seqs, seqStream);

    String<char, MMap<> > reads;
    open(reads, toCString(options.readsFileName));

    Pdf<Uniform<int> > uniformSeqIdRng(0, length(seqs) - 1);

    std::ofstream outStream;
    outStream.open(toCString(options.startPosition), std::ios::out);
    outStream << "#Stores the start positions of the reads in the db (id of read/pos within read)." << std::endl;

    std::ofstream outStreamQual;
    outStreamQual.open(toCString(options.alignmentQuality), std::ios::out);
    outStreamQual << "#The alignment score of the mutaten reads and the alignment score of the original read (with itself) are stored here. At the end of the file the averages are stored." << std::endl;
    outStreamQual << "#entryLength\tnumErrors\talignScore\talignOrigScore" << std::endl;

    long globalAlignScore = 0;
    long globalOrigAlignScore = 0;
    for (unsigned i = 0; i < options.numSeqReads; ++i)
    {
        appendValue(reads, '>');
        append(reads, "a\n");

        Pdf<Uniform<int> > uniformInt(options.minLengthSeqReads, options.maxLengthSeqReads);
        unsigned entryLength = pickRandomNumber(rng, uniformInt);

        unsigned int uniformSeqId = pickRandomNumber(rng, uniformSeqIdRng);

        if (length(seqs[uniformSeqId]) - 1 < entryLength)
            --i;
        else
        {
            Pdf<Uniform<int> > uniformSeqStartRng(0, length(seqs[uniformSeqId]) - 1 - entryLength);
            unsigned int uniformSeqStart = pickRandomNumber(rng, uniformSeqStartRng);

            outStream << uniformSeqId << "\t" << uniformSeqStart << std::endl;

            Pdf<Uniform<double> > uniformErrorRate(0, 100);

            unsigned numErrors = 0;
            String<AminoAcid> read;
            for (unsigned i = 0; i < entryLength; ++i)
            {
                if (pickRandomNumber(rng, uniformErrorRate) < options.errorRate)
                {
                    appendValue(read, getAminoAcidRead(seqs[uniformSeqId][uniformSeqStart + i]));
                    ++numErrors;
                }
                else 
                    appendValue(read, seqs[uniformSeqId][uniformSeqStart + i]);
            }
            append(reads, read);
            append(reads, "\n");

            Align<String<AminoAcid> > align;
            resize(rows(align), 2);
            assignSource(row(align, 0), infix(seqs[uniformSeqId], uniformSeqStart, uniformSeqStart + entryLength));
            assignSource(row(align, 1), read);
            Blosum62 scoringScheme;
 
            int alignScore = localAlignment(align, scoringScheme);
            globalAlignScore += alignScore;

            assignSource(row(align, 1), infix(seqs[uniformSeqId], uniformSeqStart, uniformSeqStart + entryLength));
            int alignOrigScore = localAlignment(align, scoringScheme);
            globalOrigAlignScore += alignOrigScore;

            outStreamQual << entryLength << "\t" << numErrors << "\t" << alignScore << "\t" << alignOrigScore << std::endl;
        }
    }
    outStreamQual << globalAlignScore / options.numSeqReads << "\t" << globalOrigAlignScore / options.numSeqReads << std::endl; ;
}

int main(int argc, char const ** argv)
{
    initFreq();
    initDist(Blosum62());

    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "Create benchmark with:\n"
              << "numSeqDB: " << options.numSeqDB << '\n'
              << "minLengthSeqDB: " << options.minLengthSeqDB << '\n'
              << "maxLengthSeqDB: " << options.maxLengthSeqDB << '\n';

    createDatabase(options);
    createReads(options);

    // Print the command line arguments back to the user.

    return 0;
}

