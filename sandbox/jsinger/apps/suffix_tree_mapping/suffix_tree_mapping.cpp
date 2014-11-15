// ==========================================================================
//                            suffix_tree_mapping
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

#define SEQAN_DEBUG_INDEX
#define SEQAN_DEBUG_INDEX
#define DEBUG_INDEX
#define SEQAN_DEBUG

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/gff_io.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/translation.h>

#include <seqan/arg_parse.h>

namespace seqan
{
template <>
struct SAValue<StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > >
{
    typedef Pair<unsigned short, unsigned> Type;
};
}

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.

struct AppOptions
{
    CharString sampleFileName;
    CharString dbFileName;
    CharString outputFileName;
    CharString outputIndexName;
    CharString inputType;
    int numThreads;
    unsigned bufferSize; // #reads read at a time
    bool onlyStoreIndex;
    unsigned geneticCode; // which translation code to use

    AppOptions() :
        numThreads(1),
        bufferSize(100000),
        onlyStoreIndex(false),
        geneticCode(0)
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
    seqan::ArgumentParser parser("suffix_tree_mapping");
    // Set short description, version, and date.
    setShortDescription(parser, "This porgram searches for exact query matches using a suffix tree index.");
    setVersion(parser, "0.1");
    setDate(parser, "July 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-is\\fP \\fISAMPLEFILE\\fP \\fB-ir\\fP \\fIREFERENCE\\fP \\fB-t\\fP \\fITYPE [dna, peptide]\\fP \\fB-o\\fP \\fIOUTFILE\\fP");
    addDescription(parser, "This program determines the locations of the input queries in the reference.");

    addOption(parser, ArgParseOption("is", "inputSample", "Name of the multi-FASTA input.", ArgParseArgument::INPUTFILE, "IN"));
    setRequired(parser, "is");
    setValidValues(parser, "is", "fasta fa fna");
   
    addOption(parser, ArgParseOption("t", "inputType", "The type of the input.", ArgParseArgument::STRING));
    setRequired(parser, "t");
    setValidValues(parser, "t", "dna peptide");

    addOption(parser, ArgParseOption("ir", "inputRefs", "Name of the multi-FASTA input.", ArgParseArgument::INPUTFILE, "IN"));
    setRequired(parser, "ir");
    setValidValues(parser, "ir", "fasta fa fna");

    addOption(parser, ArgParseOption("o", "outputFileName", "Name of the output file.", ArgParseArgument::OUTPUTFILE, "OUT"));
    setRequired(parser, "o");
    setValidValues(parser, "o", "gff");

    addOption(parser, ArgParseOption("oI", "outputIndexName", "Name of the index to store.", ArgParseArgument::OUTPUTFILE, "OUT"));
    addOption(parser, seqan::ArgParseOption("sO", "storeOnly", "Stop after storing the index."));

    addOption(parser, ArgParseOption("th", "threads", "The number of threads to be used.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("bf", "bufferSize", "The number of reads stored in a buffer before writing them to disk.", ArgParseArgument::INTEGER));
    ArgParseOption gcOption("gc", "geneticCode", "The genetic code to be used for translation.", ArgParseArgument::INTEGER);
    setHelpText(gcOption, "There are several different genetic codes available taken from NCBI: 0: Canonical, 1: VertMitochondrial, 2: YeastMitochondrial, 3: MoldMitochondrial, 4: InvertMitochondrial, 5: Ciliate, 6: FlatwormMitochondrial, 7: Euplotid, 8: Prokaryote, 9: AltYeast, 10: AscidianMitochondrial, 11: AltFlatwormMitochondrial, 12: Blepherisma, 13: ChlorophyceanMitochondrial, 14: TrematodeMitochondrial, 15: ScenedesmusMitochondrial, 16: ThraustochytriumMitochondrial, 17: PterobranchiaMitochondrial, 18: Gracilibacteria");
    addOption(parser, gcOption);

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getOptionValue(options.sampleFileName, parser, "is");
    getOptionValue(options.dbFileName, parser, "ir");
    getOptionValue(options.outputFileName, parser, "o");
    options.outputIndexName = options.dbFileName;
    getOptionValue(options.outputIndexName, parser, "oI");
    getOptionValue(options.inputType, parser, "t");
    getOptionValue(options.numThreads, parser, "th");
    getOptionValue(options.bufferSize, parser, "bf");
    getOptionValue(options.geneticCode, parser, "gc");
    options.onlyStoreIndex = isSet(parser, "storeOnly");

    return seqan::ArgumentParser::PARSE_OK;
}

template <typename TNewPos, typename TOldPos>
bool genomicLocationSet(TNewPos & newPos, TOldPos const & pos, unsigned queryLength, unsigned dbHitLength)
{
    unsigned start = pos.i2;

    switch (pos.i1 % 6)
    {
    case (0) :
    {
        newPos.i1 = (start * 3);
        newPos.i2 = newPos.i1 + (queryLength * 3);
        break;
    }
    case (1) :
    {
        newPos.i1 = (start * 3) +1;
        newPos.i2 = newPos.i1 + (queryLength * 3);
        break;
    }
    case (2) :
    {
        newPos.i1 = start * 3 + 2;
        newPos.i2 = newPos.i1 + (queryLength * 3);
        break;
    }
    case (3) :
    {
        int genomicStartR = (start * 3);
        int genomicEndR = genomicStartR + (queryLength * 3);
        newPos.i1 = dbHitLength  - genomicEndR;
        newPos.i2 = dbHitLength  - genomicStartR;
        break;
    }
    case (4) :
    {
        int genomicStartR = (start * 3) + 1;
        int genomicEndR = genomicStartR + (queryLength * 3);
        newPos.i1 = dbHitLength  - genomicEndR;
        newPos.i2 = dbHitLength  - genomicStartR;
        break;
    }
    default :
    {
        int genomicStartR = start * 3 + 2;
        int genomicEndR = genomicStartR + (queryLength * 3);
        newPos.i1 = dbHitLength  - genomicEndR;
        newPos.i2 = dbHitLength  - genomicStartR;
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

    std::cout << "Suffix-Tree Mapper\n"
              << "==================\n\n";

    // DB
    double startTime = sysTime();
    std::cout << "Reading database took ";
    seqan::SequenceStream seqInDB(seqan::toCString(options.dbFileName));
    seqan::StringSet<seqan::CharString> dbIds;
    seqan::StringSet<seqan::String<char>, Owner<seqan::ConcatDirect<> > > dbRaw;

    if (atEnd(seqInDB))
    {
        std::cout << "ERROR: File does not contain any sequences!\n";
        return 1;
    }
    if(seqan::readAll(dbIds, dbRaw, seqInDB) != 0)
    {
        std::cout << "ERROR: Could not read db!\n";
        return 1;
    }
    std::cout << sysTime() - startTime << " seconds" << std::endl;


    seqan::StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > aaSeqs;

    typedef Index<StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > > TIndex;
    typedef Finder<TIndex> TFinder;

    TIndex _index;
    TFinder finder(_index);

    startTime = sysTime();
    std::cout << "Reading/computing the index took ";
    if (!open(_index, toCString(options.outputIndexName)))
    {
        if (options.inputType == "dna")
        {
            switch(options.geneticCode) {
                case (0) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<Canonical>());
                           break;
                case (1) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<VertMitochondrial>());
                           break;
                case (2) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<YeastMitochondrial>());
                           break;
                case (3) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<MoldMitochondrial>());
                           break;
                case (4) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<InvertMitochondrial>());
                           break;
                case (5) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<Ciliate>());
                           break;
                case (6) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<FlatwormMitochondrial>());
                           break;
                case (7) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<Euplotid>());
                           break;
                case (8) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<Prokaryote>());
                           break;
                case (9) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<AltYeast>());
                           break;
                case (10) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<AscidianMitochondrial>());
                           break;
                case (11) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<AltFlatwormMitochondrial>());
                           break;
                case (12) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<Blepherisma>());
                           break;
                case (13) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<ChlorophyceanMitochondrial>());
                           break;
                case (14) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<TrematodeMitochondrial>());
                           break;
                case (15) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<ScenedesmusMitochondrial>());
                           break;
                case (16) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<ThraustochytriumMitochondrial>());
                           break;
                case (17) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<PterobranchiaMitochondrial>());
                           break;
                case (18) : translate(aaSeqs, dbRaw, SIX_FRAME, GeneticCode<Gracilibacteria>());
                           break;
                default : std::cout << "Unknown Gentetic Code!" << std::endl;
                          return 1;
            }
            _index = TIndex(aaSeqs);
        }
        else
        {
            for (unsigned i = 0; i < length(dbRaw); ++i)
                appendValue(aaSeqs, dbRaw[i]);
            _index = TIndex(aaSeqs);
        }
        find(finder, "A"); // dummy to force the index creation.
        save(_index, toCString(options.outputIndexName));

        if (options.onlyStoreIndex)
            return 0;
    }
    std::cout << sysTime() - startTime << " seconds" << std::endl;

    String<unsigned> origLength;
    resize(origLength, length(dbRaw), 0);
    for (unsigned i = 0; i < length(origLength); ++i)
        origLength[i] = length(dbRaw[i]);

    // free some space
    clear(dbRaw);
    shrinkToFit(dbRaw);

    seqan::SequenceStream seqInSample(seqan::toCString(options.sampleFileName));
    std::fstream fout(toCString(options.outputFileName), std::ios::binary | std::ios::out | std::ios::app);

    startTime = sysTime();
    std::cout << "Reading and searching for the queries took ";
    SEQAN_OMP_PRAGMA(parallel num_threads(options.numThreads))
    for (; ;)
    {
        // Sample
        seqan::StringSet<seqan::CharString> sampleIds;
        seqan::StringSet<seqan::String<char> > sampleSeqs;

        SEQAN_OMP_PRAGMA(critical (read_chunk))
        {

            if (!atEnd(seqInSample))
            {
                if(seqan::readBatch(sampleIds, sampleSeqs, seqInSample, options.bufferSize) != 0)
                {

                    std::cout << "ERROR: Could not read samples!\n";
                }
            }
        }

        if (length(sampleIds) == 0)
            break;


        String<GffRecord> records;
        GffRecord record;
        record.type="peptide";
        appendValue(record.tagName, "ID");
        appendValue(record.tagName, "Name");
        resize(record.tagValue, 2);
        Pair<unsigned> correctPos;
        if (options.inputType == "dna")
        {
            for (unsigned i = 0; i < length(sampleSeqs); ++i)
            {
                clear(finder);
                while(find(finder, sampleSeqs[i]))
                {
                    record.ref=dbIds[position(finder).i1 / 6];
                    record.source=options.sampleFileName;
                    genomicLocationSet(correctPos, position(finder), length(sampleSeqs[i]), origLength[position(finder).i1 / 6]);
                    record.beginPos=correctPos.i1;
                    record.endPos=correctPos.i2;
                    (position(finder).i2 % 6 > 2) ? record.strand = '+' : record.strand = '-';
                    record.tagValue[0] = sampleIds[i];
                    record.tagValue[1] = record.source;
                    appendValue(records, record);
                }
                SEQAN_OMP_PRAGMA(critical (write_chunk))
                if (length(records) > options.bufferSize)
                {
                    for (unsigned j = 0; j < length(records); ++j)
                    {
                        writeRecord(fout, records[j], Gff());
                    }
                    clear(records);
                }
            }
        }
        else
        {
            record.strand = '+';
            for (unsigned i = 0; i < length(sampleSeqs); ++i)
            {
                clear(finder);
                while(find(finder, sampleSeqs[i]))
                {
                    record.ref=dbIds[position(finder).i1];
                    record.source=options.sampleFileName;
                    record.beginPos=position(finder).i2;
                    record.endPos=record.beginPos + length(sampleSeqs[i]);
                    record.tagValue[0] = sampleIds[i];
                    record.tagValue[1] = record.source;
                    appendValue(records, record);

                }
                SEQAN_OMP_PRAGMA(critical (write_chunk))
                if (length(records) > options.bufferSize)
                {
                    for (unsigned j = 0; j < length(records); ++j)
                        writeRecord(fout, records[j], Gff());
                    clear(records);
                }
            }
        }
        SEQAN_OMP_PRAGMA(critical (write_chunk))
        for (unsigned j = 0; j < length(records); ++j)
        {
            writeRecord(fout, records[j], Gff());
        }
    }
    std::cout << sysTime() - startTime << " seconds" << std::endl;
    return 0;
}


