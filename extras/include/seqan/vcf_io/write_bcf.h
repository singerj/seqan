// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_WRITE_BCF_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_WRITE_BCF_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Bcf_;
typedef Tag<Bcf_> Bcf;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function write()                                                 [VcfHeader]
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should become writeRecord?

/*!
 * @fn VcfIO#write
 * @headerfile <seqan/vcf_io.h>
 * @brief Write a VcfHeader.
 *
 * @signature int write(target, header, context, Vcf());
 *
 * @param[in,out] target  The StreamConcept to write to.
 * @param[out]    header  The VcfHeader to write.
 * @param[in,out] context VcfIOContext to use.
 *
 * @return int A status code, 0 on success, a different value otherwise.
 */

/**
.Function.VCF I/O#write
..cat:VCF I/O
..summary:Write a @Class.VcfHeader@.
..signature:int write(stream, header, context, Vcf())
..param.stream:The @Concept.StreamConcept@ to write to.
...type:Concept.StreamConcept
..param.header:The @Class.VcfHeader@ to write.
...type:Class.VcfHeader
..param.context:The @Class.VcfIOContext@ to use for writing.
...class:Class.VcfIOContext
..return:$0$ on success, $1$ on failure.
..include:seqan/vcf_io.h
*/

template <typename TTarget, typename TBcfValueType>
void _getBcfDictionaryId(TTarget & target, TBcfValueType const & bcfValue)
{
    typename Iterator<TBcfValueType const>::Type it = begin(bcfValue);
    clear(target);
    for (; !atEnd(it); ++it)
        if (value(it) == '=')
        {
            ++it;
            for (; value(it) != ','; ++it)
                writeValue(target, value(it));
            return;
        }
}

template <typename TTarget>
void
_write(TTarget & target,
      VcfHeader const & header,
      VcfIOContext const & vcfIOContextTag,
      StringSet<String<char> > & contigNameStore,
      NameStoreCache<StringSet<String<char> > > & contigNameStoreCache,
      StringSet<String<char> > & restNameStore,
      NameStoreCache<StringSet<String<char> > > & restNameStoreCache,
      Bcf const & /*tag*/)
{
    typedef Size<VcfHeader>::Type       TSize;

    write(target, "BCF\2\2");

    // create the header including ids for contigs on the one side and 
    // INFO, FORMAT and FILTER on the other side
    String<char> vcfHeader;
    String<char> tempId;
    unsigned int contigId = 0;

    //write(vcfHeader, "##");
    //write(vcfHeader, header.headerRecords[0].key);
    //writeValue(vcfHeader, '=');
    //write(vcfHeader, header.headerRecords[0].value);
    //writeValue(vcfHeader, '\n');

    /*if (header.headerRecords[0].key != "FILTER" || prefix(header.headerRecords[0].value, 9) != "<ID=PASS")
    {
        write(vcfHeader, "##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>\n");
        appendName(restNameStore, "PASS", restNameStoreCache);
        refresh(restNameStoreCache);
    }*/

    for (TSize i = 0; i < length(header.headerRecords); ++i)
    {
        write(vcfHeader, "##");
        write(vcfHeader, header.headerRecords[i].key);
        writeValue(vcfHeader, '=');
        write(vcfHeader, header.headerRecords[i].value);
        if (header.headerRecords[i].key == "contig")
        {
            _getBcfDictionaryId(tempId, header.headerRecords[i].value);
            appendName(contigNameStore, tempId, contigNameStoreCache);
            back(vcfHeader) = ',';
            write(vcfHeader, "IDX=");
            appendNumber(vcfHeader, contigId);
            writeValue(vcfHeader, '>');
            ++contigId;
        }
        else if (header.headerRecords[i].key == "INFO" ||
                 header.headerRecords[i].key == "FORMAT" ||
                 header.headerRecords[i].key == "FILTER")
        {
            unsigned int idx;
            _getBcfDictionaryId(tempId, header.headerRecords[i].value);
            if (!getIdByName(restNameStore, tempId, idx, restNameStoreCache))
            {
                appendName(restNameStore, tempId, restNameStoreCache);
                refresh(restNameStoreCache);
                idx = length(restNameStore) - 1;
                back(vcfHeader) = ',';
                write(vcfHeader, "IDX=");
                appendNumber(vcfHeader, idx);
                write(vcfHeader, '>');
            }
        }
        writeValue(vcfHeader, '\n');
    }

    write(vcfHeader, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (unsigned i = 0; i < length(header.sampleNames); ++i)
    {
        writeValue(vcfHeader, '\t');
        write(vcfHeader, header.sampleNames[i]);
    }
    writeValue(vcfHeader, '\n');
    writeValue(vcfHeader, '\0');
    appendRawPod(target, (__uint32)length(vcfHeader));
    write(target, vcfHeader);
}

template <typename TTarget>
void
write(TTarget & target,
      VcfHeader const & header,
      VcfIOContext const & vcfIOContextTag,
      Bcf const & /*tag*/)
{
    typedef StringSet<String<char> >    TNameStore;
    typedef NameStoreCache<TNameStore>  TNameStoreCache;

    TNameStore contigNameStore;
    TNameStoreCache contigNameStoreCache(contigNameStore);

    TNameStore restNameStore;
    TNameStoreCache restNameStoreCache(restNameStore);

    _write(target, header, vcfIOContextTag, contigNameStore, contigNameStoreCache, restNameStore, restNameStoreCache, Bcf());
}



// ----------------------------------------------------------------------------
// Function writeRecord()                                           [VcfRecord]
// ----------------------------------------------------------------------------


/*!
 * @fn VcfIO#writeRecord
 * @headerfile <seqan/vcf_io.h>
 * @brief Write a VcfRecord.
 *
 * @signature int writeRecord(target, record, context, Vcf());
 *
 * @param[in,out] target  The StreamConcept to write to.
 * @param[out]    record  The VcfRecord to write.
 * @param[in,out] context VcfIOContext to use.
 *
 * @return int A status code, 0 on success, a different value otherwise.
 */

/**
.Function.VCF I/O#writeRecord
..cat:VCF I/O
..summary:Write a @Class.VcfRecord@.
..signature:int writeRecord(stream, record, context, Vcf())
..param.stream:The @Concept.StreamConcept@ to write to.
...type:Concept.StreamConcept
..param.record:The @Class.VcfRecord@ to write.
...type:Class.VcfRecord
..param.context:The @Class.VcfIOContext@ to use for writing.
...class:Class.VcfIOContext
..return:$0$ on success, $1$ on failure.
..include:seqan/vcf_io.h
*/

template <typename TSequence>
unsigned int _countEntriesAlt(TSequence const & seq)
{
    unsigned counter = 0;
    for (unsigned i = 0; i < length(seq); ++i)
        if (seq[i] == ',' || seq[i] == ';' || seq[i] == ':')
            ++counter;

    if (back(seq) == ',' || back(seq) == ';' || back(seq) == ':')
        return counter;
    return counter + 1;
}
template <typename TTarget>
void
_write(TTarget & target,
             VcfRecord const & record,
             VcfIOContext const & vcfIOContext,
             StringSet<String<char> > & contigNameStore,
             NameStoreCache<StringSet<String<char> > > & contigNameStoreCache,
             StringSet<String<char> > & restNameStore,
             NameStoreCache<StringSet<String<char> > > & restNameStoreCache,
             Bcf const & /*tag*/)
{
    String<char> buffer;
    __int32 idx;
    
    // l_shared
    appendRawPod(target, (__uint32)length(buffer));
    
    // l_indiv
    appendRawPod(target, (__uint32)length(buffer));
    
    // CHROM
    getIdByName(contigNameStore, (*vcfIOContext.sequenceNames)[record.rID], idx, contigNameStoreCache);
    std::cerr << idx << " " << toCString((*vcfIOContext.sequenceNames)[record.rID]) << std::endl;
    appendRawPod(target, idx);
    
    // POS
    appendRawPod(target, record.beginPos);
    
    // rlen
    appendRawPod(target, (__int32)length(record.ref));
    
    // QUAL
    appendRawPod(target, record.qual);
    
    // n_allele_info
    idx = (_countEntriesAlt(record.ref) + _countEntriesAlt(record.alt)) << 16 | _countEntriesAlt(record.info);
    appendRawPod(target, idx);
    
    // n_fmt_sample
    idx = _countEntriesAlt(record.format) << 24 | length(record.genotypeInfos);
    appendRawPod(target, idx);
    
    // ID
    int bla = -105; // dont know why, but this is included in the bcftools of samtools
    appendRawPod(target, (char)(bla));
    write(target, record.id);
    
    // allele
    bla = 23; // dont know why, but this is included in the bcftools of samtools
    appendRawPod(target, (char)(bla));
    write(target, record.ref);
    bla = 23; // dont know why, but this is included in the bcftools of samtools
    appendRawPod(target, (char)(bla));
    write(target, record.alt);

    // filter info
    // split the filter info
    StringSet<String<char> > dummy;
    splitString(dummy, record.filter, ';', '"'); 
    for (unsigned i = 0; i < length(dummy); ++i)
    {
        bla = 17; // dont know why, but this is included in the bcftools of samtools
        appendRawPod(target, (char)(bla));
        getIdByName(restNameStore, dummy[i], idx, restNameStoreCache);
        appendRawPod(target, (__int8)idx);
    }

    // filter info
    // split the filter info
    splitString(dummy, record.info, ';', '"'); 
    for (unsigned i = 0; i < length(dummy); ++i)
    {
        bla = 17; // dont know why, but this is included in the bcftools of samtools
        appendRawPod(target, (char)(bla));
        getIdByName(restNameStore, dummy[i], idx, restNameStoreCache);
        appendRawPod(target, (__int8)idx);
    }


    //bla = 17; // dont know why, but this is included in the bcftools of samtools
    //appendRawPod(target, (char)(bla));
    
    char c;
    std::cin>>c;
    std::cerr << std::endl;
    /*
    write(target, (*vcfIOContext.sequenceNames)[record.rID]);
    writeValue(target, '\t');
    appendNumber(target, record.beginPos + 1);
    writeValue(target, '\t');
    if (empty(record.id))
        writeValue(target, '.');
    else
        write(target, record.id);
    writeValue(target, '\t');
    if (empty(record.ref))
        writeValue(target, '.');
    else
        write(target, record.ref);
    writeValue(target, '\t');
    if (empty(record.alt))
        writeValue(target, '.');
    else
        write(target, record.alt);
    writeValue(target, '\t');
    if (record.qual != record.qual)  // only way to test for nan
        writeValue(target, '.');
    else
        appendNumber(target, record.qual);
    writeValue(target, '\t');
    if (empty(record.filter))
        writeValue(target, '.');
    else
        write(target, record.filter);
    writeValue(target, '\t');
    if (empty(record.info))
        writeValue(target, '.');
    else
        write(target, record.info);
    writeValue(target, '\t');
    if (empty(record.format))
        writeValue(target, '.');
    else
        write(target, record.format);
    for (unsigned i = 0; i < length(record.genotypeInfos); ++i)
    {
        writeValue(target, '\t');
        if (empty(record.genotypeInfos[i]))
            writeValue(target, '.');
        else
            write(target, record.genotypeInfos[i]);
    }
    writeValue(target, '\n');
    */
}

// TODO(singer): VcfHeader not const because VcfIOContext can not handle it
template <typename TTarget, typename TRecord>
void
write(TTarget & target,
      VcfHeader & header,
      String<TRecord> const & records,
      Bcf const & /*tag*/)
{
    typedef StringSet<String<char> >    TNameStore;
    typedef NameStoreCache<TNameStore>  TNameStoreCache;

    TNameStore contigNameStore;
    TNameStoreCache contigNameStoreCache(contigNameStore);

    TNameStore restNameStore;
    TNameStoreCache restNameStoreCache(restNameStore);

    seqan::VcfIOContext vcfIOContext(header.sequenceNames, header.sampleNames);

    _write(target, header, vcfIOContext, contigNameStore, contigNameStoreCache, restNameStore, restNameStoreCache, Bcf());
    for (unsigned i = 0; i < length(records); ++i)
        _write(target, records[i], vcfIOContext, contigNameStore, contigNameStoreCache, restNameStore, restNameStoreCache, Bcf());

}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_WRITE_BCF_H_

