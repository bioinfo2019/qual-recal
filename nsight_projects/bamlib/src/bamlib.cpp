/*
 * bamlib.cpp
 *
 *  Created on: Mar 11, 2017
 *      Author: eliot
 */

#include "bamlib.hpp"
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/misc/interval_tree.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/bed_io.h>
#include <seqan/basic.h>
#include <seqan/vcf_io.h>

#include <iterator>
#include <random>

using namespace seqan;



/** Compress a STL string using zlib with given compression level and return
  * the binary data. */
int compress_string(const std::string& str,
                            int compressionlevel = Z_BEST_COMPRESSION)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();           // set the z_stream's input

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw(std::runtime_error(oss.str()));
    }

    return outstring.size();
}


template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


template<typename TK, typename TV>
std::vector<TK> extract_keys(std::map<TK, TV> const& input_map) {
	std::vector<TK> retval;
	for (auto const& element : input_map) {
		retval.push_back(element.first);
	}
	return retval;
}


template<typename TK, typename TV>
std::vector<TV> extract_values(std::map<TK, TV> const& input_map) {
	std::vector<TV> retval;
	for (auto const& element : input_map) {
		retval.push_back(element.second);
	}
	return retval;
}

int xformCoords(int pos, int contigLength, int readLength)
{
    return contigLength - pos - readLength;
}


bool isGoodRead(BamAlignmentRecord record, int refContigLength, std::string refContigName)
{
	// ART read name
	// @SL2.50ch00-5451450|-|18320722|+|3484789/2
	// format is:
	// @read name|read1 strand|read1 start pos|read2 strand|read2 start pos
	//      0          1             2               3            4

	// varsim read name
	// @SL2.50ch00-18999559-,SL2.50ch00-18999587-D,SL2.50ch00-18999560-:SL2.50ch00-18999635--:::::::::::270378:0/1

	int beginPos = record.beginPos;
	int mateBeginPos = record.pNext;


	std::vector<std::string> splits = split(toCString(record.qName), ':');

	// read1 location(s) in splits[0], read2 location(s) in splits[1]. now split locations on comma


	std::vector<std::string> read1_locs = split(splits[0], ',');
	std::vector<std::string> read2_locs = split(splits[1], ',');

	std::string read1Chromosome = "";

	std::vector<int> read1_positions;
	std::vector<std::string> tmp_vec;

	for(int i = 0;i < read1_locs.size();++i)
	{
		tmp_vec = split(read1_locs[i], '-');
		read1Chromosome = tmp_vec[0];
		read1_positions.push_back(std::stoi(tmp_vec[1]));
	}


	std::string read2Chromosome = "";

	std::vector<int> read2_positions;

	for(int i = 0;i < read2_locs.size();++i)
	{
		tmp_vec = split(read2_locs[i], '-');
		read2Chromosome = tmp_vec[0];
		read2_positions.push_back(std::stoi(tmp_vec[1]));
	}

	for(int i = 0;i < read1_positions.size();++i)
	{
		if((fabs(read1_positions[i] - beginPos) < 6) && read1Chromosome == refContigName)
			return true;
	}

	for(int i = 0;i < read2_positions.size();++i)
	{
		if((fabs(read2_positions[i] - beginPos) < 6) && read2Chromosome == refContigName)
			return true;
	}

	return false;

}



std::map<int, double> build_mapq_histogram(std::string bamFileName, int underSamplingFactor) {

	int numMultiMaps = 0;
	int numBadReads = 0;
	int numGoodReads = 0;

	//Open BamFileIn for reading.
	BamFileIn bamFileIn;
	std::map<int, double> mapqHistogram;

	if (!open(bamFileIn, toCString(bamFileName)))
	{
		std::cerr << "ERROR: Could not open " << bamFileName << " for reading.\n";
		return mapqHistogram; // return an empty map
	}

	BamHeader header;
	readHeader(header, bamFileIn);

	typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

	TBamContext const & bamContext = context(bamFileIn);


	// Read records.
	BamAlignmentRecord record;

	while (!atEnd(bamFileIn))
	{
		readRecord(record, bamFileIn);

		if (hasFlagUnmapped(record))
		    continue;

		int mapQ = record.mapQ;

		int tagIdx = 0;
		int tagValInt = 0;


		BamTagsDict tagsDict(record.tags);

		int primaryAlignmentScore;
		int secondaryAlignmentScore = -75;

		if (findTagKey(tagIdx, tagsDict, "XS"))
			extractTagValue(secondaryAlignmentScore, tagsDict, tagIdx);

		if (findTagKey(tagIdx, tagsDict, "AS"))
			extractTagValue(primaryAlignmentScore, tagsDict, tagIdx);

		if (primaryAlignmentScore == secondaryAlignmentScore)
		{

			numMultiMaps++;
			continue;
		}

		std::string refContigName = toCString(contigNames(bamContext)[record.rID]);
		int refContigLength = contigLengths(bamContext)[record.rID];

		if (isGoodRead(record, refContigLength, refContigName))
		{
			mapqHistogram[mapQ] += 1.0;
			numGoodReads += 1;
		}
		else
		{
			numBadReads += 1;
		}

	}

	//at this point, normalize the values in the histogram. Basically divide each value by the
	//number of good reads. Then multiply that number by numBadReads/numGoodReads. When we sample
	//a Bernouli distribution with the above figure as the probability of success for each histogram bin,
	//we will end up sampling approx the same number of good reads as bad reads and with a mapQ
	//distribution that approximates the file we sampled from. And we should have a nice even
	//distribution across all chromosomes.

	std::vector<int> keys = extract_keys(mapqHistogram);
    double total_prob = 0;
	for(int i = 0;i < keys.size();++i) {
		mapqHistogram[keys[i]] = mapqHistogram[keys[i]] / (double)numGoodReads * (double)underSamplingFactor;

		total_prob += mapqHistogram[keys[i]];
		std::cerr << "MapQ " << keys[i] << ": " << mapqHistogram[keys[i]] << "\n";
	}
	std::cerr << "Found a total of " << numGoodReads << " good reads, and " << numBadReads << " bad reads, total prob " << total_prob << "\n";

	return mapqHistogram;

}


std::random_device rd;
long seed = rd();
std::mt19937 mt_rand(seed);

bool random_bool_with_prob( double prob )  // probability between 0.0 and 1.0
{
    std::bernoulli_distribution d(prob);
    return d(mt_rand);
}

float getBadReadFraction(std::string inBamFileName)
{

	int numMultiMaps = 0;
	int numBadReads = 0;
	int numGoodReads = 0;

	//Open BamFileIn for reading.
	BamFileIn bamFileIn;

	if (!open(bamFileIn, toCString(inBamFileName)))
	{
		std::cerr << "ERROR: Could not open " << inBamFileName << " for reading.\n";
		return 1;
	}

	BamHeader header;

	readHeader(header, bamFileIn);

	typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

	TBamContext const & bamContext = context(bamFileIn);

	// Read records.
	BamAlignmentRecord record;

	std::cerr << "Counting reads....\n";

	while (!atEnd(bamFileIn))
	{
		readRecord(record, bamFileIn);

        	if (hasFlagUnmapped(record))
		    continue;

		int mapQ = record.mapQ;

		int tagIdx = 0;

		BamTagsDict tagsDict(record.tags);

		int primaryAlignmentScore = 0;
		int secondaryAlignmentScore = 75;

		if (findTagKey(tagIdx, tagsDict, "XS"))
			extractTagValue(secondaryAlignmentScore, tagsDict, tagIdx);

		if (findTagKey(tagIdx, tagsDict, "AS"))
			extractTagValue(primaryAlignmentScore, tagsDict, tagIdx);


		if (primaryAlignmentScore == secondaryAlignmentScore)
		{
			numMultiMaps++;
			continue;
		}

		std::string refContigName = toCString(contigNames(bamContext)[record.rID]);
		int refContigLength = contigLengths(bamContext)[record.rID];

		if (isGoodRead(record, refContigLength, refContigName))
		{

			numGoodReads += 1;

		}
		else
		{
			numBadReads += 1;

		}

	}

	return ((float)numBadReads)/((float)numGoodReads);

}

int sampleBamFile(std::string inBamFileName, std::string outBamFileName)
{
	float badReadFraction = getBadReadFraction(inBamFileName);

	int numMultiMaps = 0;
	int numBadReads = 0;
	int numGoodReads = 0;

	//Open BamFileIn for reading.
	BamFileIn bamFileIn;

	if (!open(bamFileIn, toCString(inBamFileName)))
	{
		std::cerr << "ERROR: Could not open " << inBamFileName << " for reading.\n";
		return 1;
	}

	BamFileOut bamFileOut(context(bamFileIn), toCString(outBamFileName));

	BamHeader header;

	readHeader(header, bamFileIn);

	writeHeader(bamFileOut, header);

	typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

	TBamContext const & bamContext = context(bamFileIn);

	// Read records.
	BamAlignmentRecord record;

	std::cerr << "Sampling BAM file....\n";

	while (!atEnd(bamFileIn))
	{
		readRecord(record, bamFileIn);

        	if (hasFlagUnmapped(record))
		    continue;

		int mapQ = record.mapQ;

		int tagIdx = 0;

		BamTagsDict tagsDict(record.tags);

		int primaryAlignmentScore = 0;
		int secondaryAlignmentScore = 75;

		if (findTagKey(tagIdx, tagsDict, "XS"))
			extractTagValue(secondaryAlignmentScore, tagsDict, tagIdx);

		if (findTagKey(tagIdx, tagsDict, "AS"))
			extractTagValue(primaryAlignmentScore, tagsDict, tagIdx);


		if (primaryAlignmentScore == secondaryAlignmentScore)
		{
			numMultiMaps++;
			continue;
		}

		std::string refContigName = toCString(contigNames(bamContext)[record.rID]);
		int refContigLength = contigLengths(bamContext)[record.rID];

		if (isGoodRead(record, refContigLength, refContigName))
		{

			if(random_bool_with_prob(badReadFraction))
			{
				writeRecord(bamFileOut, record);
				numGoodReads += 1;
			}

		}
		else
		{
				writeRecord(bamFileOut, record);
				numBadReads += 1;

		}

	}

	std::cerr << "Wrote a total of " << numGoodReads << " good reads, and " << numBadReads << " bad reads to new BAM file\n";
	return 0;

}


int buildIntervalTree(IntervalTree<int, CharString>& tree, std::string bamFileName, int chromosomeIndex, int readLength)
{

	//Open BamFileIn for reading.
	BamFileIn inFile;


	if (!open(inFile, toCString(bamFileName)))
	{
		std::cerr << "ERROR: Could not open " << bamFileName << " for reading.\n";
		return 1;
	}

	std::string baiFileName = bamFileName + ".bai";

	// Read BAI index.
	BamIndex<Bai> baiIndex;
	if (!open(baiIndex, toCString(baiFileName)))
	{
		std::cerr << "ERROR: Could not read BAI index file " << baiFileName << "\n";
		return 1;
	}

	// Read header.
	BamHeader header;
	readHeader(header, inFile);

	typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

	TBamContext const & bamContext = context(inFile);

	// Translate BEGIN and END arguments to number, 1-based to 0-based.
	int beginPos = 1, endPos = contigLengths(bamContext)[chromosomeIndex];

	// 1-based to 0-based.
	beginPos -= 1;
	endPos -= 1;

	// Jump the BGZF stream to this position.
	bool hasAlignments = false;
	if (!jumpToRegion(inFile, hasAlignments, chromosomeIndex, beginPos, endPos, baiIndex))
	{
		std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
		return 1;
	}

	if (!hasAlignments)
		return 0;  // No alignments here.

	// Seek linearly to the selected position.
	BamAlignmentRecord record;

	// fill a string of intervals and keys
	typedef IntervalAndCargo<int, CharString> TInterval;
	String<TInterval> intervals;

	int readCounter = 0;

	while (!atEnd(inFile))
	{

		readRecord(record, inFile);

		// If we are on the next reference or at the end already then we stop.
		if (record.rID == -1 || record.rID > chromosomeIndex || record.beginPos >= endPos)
			break;
		// If we are left of the selected position then we skip this record.
		if (record.beginPos < beginPos)
			continue;

		readCounter += 1;

		appendValue(intervals, TInterval(record.beginPos, record.beginPos + readLength, "read" + std::to_string(readCounter)));


	}

    createIntervalTree(tree, intervals);

    return 0;
}


float calcGcContent(CharString dnaSequence)
{
	float gcCount = 0;

	for (unsigned int i = 0;i < length(dnaSequence);++i)
		if(dnaSequence[i] == 'C' || dnaSequence[i] == 'G')
			gcCount += 1.0;

	return gcCount/((float) length(dnaSequence));
}


int getChromosomeIntervalSequences(std::vector<unsigned char>& readSequenceData, std::string fastaFileName, std::string bamFileName, unsigned chromosomeIdx, int readLength, int offSet)
{

    int intervalStart = 0;
    int intervalCounter = 0;

    // Try to load index and create on the fly if necessary.
    FaiIndex faiIndex;
    int result = 0;


    if (!open(faiIndex, toCString(fastaFileName)) != 0)
    {
        if (build(faiIndex, toCString(fastaFileName)) != 0)
        {
            std::cerr << fastaFileName << " " << "ERROR: Index could not be loaded or built.\n";
            return 1;
        }

        if (save(faiIndex) != 0)  // Name is stored from when reading.
        {
            std::cerr << "ERROR: Index could not be written do disk.\n";
            return 1;
        }
    }

	//Open BamFileIn for reading.
	BamFileIn inFile;


	if (!open(inFile, toCString(bamFileName)))
	{
		std::cerr << "ERROR: Could not open " << bamFileName << " for reading.\n";
		return 1;
	}

	// Read BAI index.
	BamIndex<Bai> baiIndex;
	std::string baiFileName = bamFileName + ".bai";;

	if (!open(baiIndex, toCString(baiFileName)))
	{
		std::cerr << "ERROR: Could not read BAI index file " << baiFileName << "\n";
		return 1;
	}

	// Read header.
	BamHeader header;
	readHeader(header, inFile);

	typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

	TBamContext const & bamContext = context(inFile);


	// Translate BEGIN and END arguments to number, 1-based to 0-based.
	int beginPos = 1, endPos = contigLengths(bamContext)[chromosomeIdx];

	// 1-based to 0-based.
	beginPos -= 1;
	endPos -= 1;

	// Jump the BGZF stream to this position.
	bool hasAlignments = false;
	if (!jumpToRegion(inFile, hasAlignments, chromosomeIdx, beginPos, endPos, baiIndex))
	{
		std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
		return 1;
	}

	if (!hasAlignments)
		return 0;  // No alignments here.

	// Seek linearly to the selected position.
	BamAlignmentRecord record;

	int numMultiMaps = 0;

	std::cerr << "Processing chromosome: "  << toCString(contigNames(bamContext)[chromosomeIdx]) << "\n";

	CharString sequenceInfix;

	while (!atEnd(inFile))
	{
		readRecord(record, inFile);

		// If we are on the next reference or at the end already then we stop.
		if (record.rID == -1 || record.rID > chromosomeIdx || record.beginPos >= endPos)
			break;
		// If we are left of the selected position then we skip this record.
		if (record.beginPos < beginPos)
			continue;

		int leftPos = record.beginPos - offSet;
		int rightPos = record.beginPos + readLength + offSet;

		if (leftPos < 0)
		{
			leftPos = 0;
			rightPos = readLength + 2 * offSet;
		}

		if (rightPos > endPos)
		{
			leftPos = endPos - readLength - 2 * offSet;
			rightPos = endPos;
		}

		readRegion(sequenceInfix, faiIndex, chromosomeIdx, leftPos, rightPos);

		unsigned int k = 0;

		toUpper(sequenceInfix);

		String<char> s = toCString(sequenceInfix);
		String<char> N = toCString("N");

		for (k = 0;k < length(sequenceInfix); ++k)
		{
			readSequenceData.push_back((unsigned char) s[k]);
		}

		intervalCounter += 1;

    }

//	for (int k = 1;k <= readSequenceData.size(); ++k)
//	{
//		std::cout << readSequenceData[k-1];
//		if(k % 500 == 0)
//			std::cout  << std::endl;
//	}

    std::cerr << "Total genomic sequences: " << intervalCounter << std::endl;

    return 0;
}


int getFeatureVectors(std::vector<unsigned char>& readSequenceData, std::vector<unsigned char>& chromosomeIntervalSequences, std::vector<unsigned char>& readQualities, std::vector<std::array<float, 34>>& featureVectors, std::vector<std::string>& readNames, int chromosomeIdx, int readLength, int offSet, std::string bamFileName, std::string sampledBamFileName, std::string fastaFileName, int outputReadClass)
{
	/**
	 * This function will gather features for all reads covering a particular chromosome. The feature vectors are stored in a C++ map, with rean name as key.
	 * Read sequences and qualities are kept in separate vectors for later processing on the GPU. after the GPU has returned the feature vectors derived fro,
	 * invariant calculation and regression of qualities against read position, they will be married up upon return to main.
	 */

    // Try to load index and create on the fly if necessary.
    FaiIndex faiIndex;
    int result = 0;
    CharString sequenceInfix;

    if (!open(faiIndex, toCString(fastaFileName)) != 0)
    {
        if (build(faiIndex, toCString(fastaFileName)) != 0)
        {
            std::cerr << fastaFileName << " " << "ERROR: Index could not be loaded or built.\n";
            return 1;
        }

        if (save(faiIndex) != 0)  // Name is stored from when reading.
        {
            std::cerr << "ERROR: Index could not be written do disk.\n";
            return 1;
        }
    }

	std::array<float, 34> featureVector = {0};

	IntervalTree<int, CharString> readIntervals;

	result = buildIntervalTree(readIntervals, sampledBamFileName, chromosomeIdx, readLength);

	std::string baiFileName = bamFileName + ".bai";
	std::string readSequence = "";

	std::map<std::string, int> chromosomeLengths;

	int readCounter = 0;

	//Open BamFileIn for reading.
	BamFileIn inFile;


	if (!open(inFile, toCString(bamFileName)))
	{
		std::cerr << "ERROR: Could not open " << bamFileName << " for reading.\n";
		return 1;
	}

	// Read BAI index.
	BamIndex<Bai> baiIndex;
	if (!open(baiIndex, toCString(baiFileName)))
	{
		std::cerr << "ERROR: Could not read BAI index file " << baiFileName << "\n";
		return 1;
	}

	// Read header.
	BamHeader header;
	readHeader(header, inFile);

	typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

	TBamContext const & bamContext = context(inFile);


	// Translate BEGIN and END arguments to number, 1-based to 0-based.
	int beginPos = 1, endPos = contigLengths(bamContext)[chromosomeIdx];

	// 1-based to 0-based.
	beginPos -= 1;
	endPos -= 1;

	// Jump the BGZF stream to this position.
	bool hasAlignments = false;
	if (!jumpToRegion(inFile, hasAlignments, chromosomeIdx, beginPos, endPos, baiIndex))
	{
		std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
		return 1;
	}

	if (!hasAlignments)
		return 0;  // No alignments here.

	// Seek linearly to the selected position.
	BamAlignmentRecord record;

	int numMultiMaps = 0;

	std::cerr << "Processing chromosome: "  << toCString(contigNames(bamContext)[chromosomeIdx]) << "\n";

	IntervalTree<int, CharString> fpIntervals;

	while (!atEnd(inFile))
	{

		readRecord(record, inFile);


		int tagIdx = 0;
		int tagValInt = 0;


		BamTagsDict tagsDict(record.tags);

		int primaryAlignmentScore = 0;
		int secondaryAlignmentScore = 75;
		int mateAlignmentScore = 0;

		if (findTagKey(tagIdx, tagsDict, "YS"))
			extractTagValue(mateAlignmentScore, tagsDict, tagIdx);

		featureVector[MATE_ALIGNMENT_SCORE] = mateAlignmentScore;

		if (findTagKey(tagIdx, tagsDict, "XS"))
			extractTagValue(secondaryAlignmentScore, tagsDict, tagIdx);

		if (findTagKey(tagIdx, tagsDict, "AS"))
			extractTagValue(primaryAlignmentScore, tagsDict, tagIdx);

		featureVector[ALIGNMENT_SCORES_DIFF] = secondaryAlignmentScore - primaryAlignmentScore;

		// If we are on the next reference or at the end already then we stop.
		if (record.rID == -1 || record.rID > chromosomeIdx || record.beginPos >= endPos)
			break;
		// If we are left of the selected position then we skip this record.
		if (record.beginPos < beginPos)
			continue;


		if (primaryAlignmentScore == secondaryAlignmentScore)
		{
			numMultiMaps++;
			continue;
		}


		int leftPos = record.beginPos - offSet;
		int rightPos = record.beginPos + readLength + offSet;

		if (leftPos < 0)
		{
			leftPos = 0;
			rightPos = readLength + 2 * offSet;
		}

		if (rightPos > endPos)
		{
			leftPos = endPos - readLength - 2 * offSet;
			rightPos = endPos;
		}

//		readRegion(sequenceInfix, faiIndex, chromosomeIdx, leftPos, rightPos);
		unsigned int k = 0;

//		toUpper(sequenceInfix);

		String<char> s; // = toCString(sequenceInfix);

//		for (k = 0;k < length(sequenceInfix); ++k)
//		{
//			chromosomeIntervalSequences.push_back((unsigned char) s[k]);
//		}

		Byte compr[150];
		uLong comprLen = sizeof(compr);
		uLong Length = (uLong)length(record.seq)+1;

		int ReturnCode;

		ReturnCode = compress(compr, &comprLen, (const Bytef*)toCString(sequenceInfix), Length);

		featureVector[REF_COMPLEXITY_SCORE] = (float) comprLen;


		std::string refContigName = toCString(contigNames(bamContext)[record.rID]);
		int refContigLength = contigLengths(bamContext)[record.rID];

		featureVector[CLASS] = 0.0;

		if(outputReadClass == 1)
		{
			if (isGoodRead(record, refContigLength, refContigName))
				featureVector[CLASS] = 0.0;
			else
				featureVector[CLASS] = 1.0;

			featureVector[SNPREADCLASS] = 0;
		}

		featureVector[CHROMOSOME_IDX] = (float) chromosomeIdx;
		featureVector[READ_BEGIN_POS] = (float) (record.beginPos + 1);
		featureVector[MAPQ] = (float) (record.mapQ);

		if (record.tLen > 0)
			featureVector[INSERT_SIZE] = 1.0; //(float) (record.tLen);
		if (record.tLen < 0)
			featureVector[INSERT_SIZE] = -1.0; //(float) (record.tLen);
		if (record.tLen == 0)
			featureVector[INSERT_SIZE] = 0.0; //(float) (record.tLen);


		int gapExtensions = 0;

		if (findTagKey(tagIdx, tagsDict, "XG"))
			extractTagValue(gapExtensions, tagsDict, tagIdx);

		featureVector[GAP_EXTENSIONS] = (float) gapExtensions;


		int editDistance = 0;

		if (findTagKey(tagIdx, tagsDict, "NM"))
			extractTagValue(editDistance, tagsDict, tagIdx);

		featureVector[EDIT_DISTANCE] = (float) editDistance;


		int misMatches = 0;

		if (findTagKey(tagIdx, tagsDict, "XM"))
			extractTagValue(misMatches, tagsDict, tagIdx);

		featureVector[MISMATCHES] = (float) misMatches;

		int gapOpens = 0;

		if (findTagKey(tagIdx, tagsDict, "XO"))
			extractTagValue(gapOpens, tagsDict, tagIdx);

		featureVector[GAP_OPENS] = (float) gapOpens;


		int alignmentScore = 0;

		if (findTagKey(tagIdx, tagsDict, "AS"))
			extractTagValue(alignmentScore, tagsDict, tagIdx);

		featureVector[ALIGNMENT_SCORE] = (float) alignmentScore;

		featureVector[SECONDARY_ALIGNMENT_SCORE] = secondaryAlignmentScore;

		comprLen = sizeof(compr);
		Length = (uLong)length(record.seq)+1;

		ReturnCode = compress(compr, &comprLen, (const Bytef*)toCString(record.seq), Length);

		featureVector[READ_COMPLEXITY_SCORE] = (float) comprLen;

		// Pair orientation

		if(!hasFlagRC(record) && hasFlagNextRC(record)) // F1R2
			featureVector[PAIR_ORIENTATION] = 0;
		if(!hasFlagRC(record) && !hasFlagNextRC(record)) // F1F2
			featureVector[PAIR_ORIENTATION] = 1.0;
		if(hasFlagRC(record) && hasFlagNextRC(record)) // R1R2
			featureVector[PAIR_ORIENTATION] = 2.0;
		if(hasFlagRC(record) && !hasFlagNextRC(record)) // R1F2
			featureVector[PAIR_ORIENTATION] = 3.0;


		// Pair alignment type

		CharString pairAlignType = "CP";

		if (findTagKey(tagIdx, tagsDict, "YT"))
			extractTagValue(pairAlignType, tagsDict, tagIdx);


		if (pairAlignType == "CP")
			featureVector[PAIR_ALIGNMENT_TYPE] = 0.0;
		if (pairAlignType == "DP")
			featureVector[PAIR_ALIGNMENT_TYPE] = 1.0;
		if (pairAlignType == "UP")
			featureVector[PAIR_ALIGNMENT_TYPE] = 2.0;
		if (pairAlignType == "UU")
			featureVector[PAIR_ALIGNMENT_TYPE] = 3.0;


		// find intervals that overlap the query interval. Gets read depth in original BAM file
		String<CharString> results;
		findIntervals(results, readIntervals, record.beginPos + readLength/2 - 10, record.beginPos + readLength/2 + 10);

		featureVector[DEPTH] = (float) length(results);

		featureVectors.push_back(featureVector);

		std::string readSuffix = "/2";

		if(hasFlagFirst(record))
			readSuffix = "/1";

		readNames.push_back(toCString(record.qName) + readSuffix);

		//s = toCString(record.seq);

		//for (unsigned int k = 0;k < length(record.seq); ++k)
		//{
		//	readSequenceData.push_back((unsigned char) s[k]);
		//}

		s = toCString(record.qual);

		float nLowQualBases = 0;
		float avgLowQualScore = 0;

		for (int k = 0;k < length(record.qual); ++k)
		{
			readQualities.push_back(s[k] - 33);

			if(s[k] - 33 < 20)
			{
				nLowQualBases += 1.0;
				avgLowQualScore += (float)(s[k] - 33);

			}
		}

		if(nLowQualBases == 0)
			avgLowQualScore = 40;
		else
			avgLowQualScore = avgLowQualScore / nLowQualBases;

		featureVector[N_LOW_QUAL_BASES] = nLowQualBases;
		featureVector[AVG_LOW_QUAL_BASE_SCORE] = avgLowQualScore;

	}


	return 0;
}




int getTrainingFeatures(std::vector<unsigned char>& trainingData, int readLength, int nReads, bool noMultiMap, std::string bamFileName)
{
	std::string baiFileName = bamFileName + ".bai";
	std::string readSequence = "";

	std::map<std::string, int> chromosomeLengths;

	int readCounter = 0;

	//Open BamFileIn for reading.
	BamFileIn inFile;


//std::vector<std::string> x = split("SL2.50ch00_maternal-109031|-|18321202", '|');

//std::cout<<x[1]<<" "<<x[2]<<"\n";

//x = split("SL2.50ch00_maternal-109029|+|17175409", '|');

//std::cout<<x[1]<<" "<<x[2]<<"\n";

	if (!open(inFile, toCString(bamFileName)))
	{
		std::cerr << "ERROR: Could not open " << bamFileName << " for reading.\n";
		return 1;
	}

	// Read BAI index.
	BamIndex<Bai> baiIndex;
	if (!open(baiIndex, toCString(baiFileName)))
	{
		std::cerr << "ERROR: Could not read BAI index file " << baiFileName << "\n";
		return 1;
	}

	// Read header.
	BamHeader header;
	readHeader(header, inFile);

	typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

	TBamContext const & bamContext = context(inFile);

	for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
		chromosomeLengths[toCString(contigNames(bamContext)[i])] = contigLengths(bamContext)[i];

	for (unsigned i = 0; i < length(contigNames(bamContext)); ++i) {

		CharString rName = toCString(contigNames(bamContext)[i]);

		// Translate from reference name to rID.
		int rID = 0;
		if (!getIdByName(rID, contigNamesCache(context(inFile)), rName))
		{
			std::cerr << "ERROR: Reference sequence named " << rName << " not known.\n";
			return 1;
		}

		// Translate BEGIN and END arguments to number, 1-based to 0-based.
		int beginPos = 1, endPos = contigLengths(bamContext)[i];

		// 1-based to 0-based.
		beginPos -= 1;
		endPos -= 1;

		// Translate number of elements to print to number.
		int num = 3;

		// Jump the BGZF stream to this position.
		bool hasAlignments = false;
		if (!jumpToRegion(inFile, hasAlignments, rID, beginPos, endPos, baiIndex))
		{
			std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
			return 1;
		}

		if (!hasAlignments)
			return 0;  // No alignments here.

		// Seek linearly to the selected position.
		BamAlignmentRecord record;

		int numMultiMaps = 0;

		std::cerr << "Processing chromosome: "  << toCString(contigNames(bamContext)[i]) << "\n";

		while (!atEnd(inFile))
		{
			readRecord(record, inFile);


			int tagIdx = 0;
			int tagValInt = 0;


			BamTagsDict tagsDict(record.tags);

			int primaryAlignmentScore = 0;
			int secondaryAlignmentScore = -75;

			if (findTagKey(tagIdx, tagsDict, "XS"))
				extractTagValue(secondaryAlignmentScore, tagsDict, tagIdx);

			if (findTagKey(tagIdx, tagsDict, "AS"))
				extractTagValue(primaryAlignmentScore, tagsDict, tagIdx);

			// If we are on the next reference or at the end already then we stop.
			if (record.rID == -1 || record.rID > rID || record.beginPos >= endPos)
				break;
			// If we are left of the selected position then we skip this record.
			if (record.beginPos < beginPos)
				continue;


			if (primaryAlignmentScore == secondaryAlignmentScore)
			{
				numMultiMaps++;
				continue;
			}



			readSequence = "";
			String<char> s = toCString(record.seq);
			if (i > 0){
			readCounter++;
			for (int k = 0;k < length(record.seq); ++k)
			{
				trainingData.push_back((unsigned char) s[k]);
			}
			}
			//trainingData.insert(trainingData.end(), readSequence.begin(), readSequence.end());


			if(readCounter == nReads)
				break;
		}

		if(readCounter == nReads)
			break;
	}

	return 0;
}


int recalibrateBamFile(std::string scoresFile, std::string inBamFileName, std::string outBamFileName)
{

	std::ifstream input(scoresFile);
	std::map<std::string, int> qualScoresMap;
	std::string key;
	int value;

    std::string line;
    while (getline(input, line))
    {
    	std::vector<std::string> splits = split(line.c_str(), '\t');
    	qualScoresMap[splits[1]] = std::stoi(splits[2]);
    }

	//while (input >> key >> value)
		//qualScoresMap[key] = value;

	//std::cerr << qualScoresMap.size() << std::endl;

    // Open input file.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(inBamFileName)))
    {
        std::cerr << "ERROR: Could not open " << inBamFileName << std::endl;
        return 1;
    }


    try
    {
    	BamFileOut bamFileOut(bamFileIn, toCString(outBamFileName));
        // Read header.
        BamHeader header;
        readHeader(header, bamFileIn);

//        if (!open(bamFileOut, toCString(outBamFileName)))
//        {
//            std::cerr << "ERROR: Could not open " << outBamFileName << std::endl;
//            return 1;
//        }
        // Read records.

        writeHeader(bamFileOut, header);

        BamAlignmentRecord record;

        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);

            std::string readName = toCString(record.qName);

            readName += "/1";

            if(!hasFlagFirst(record)) {
            	readName = toCString(record.qName);
            	readName += "/2";
            }

			if (qualScoresMap.count(readName) != 0){

				record.mapQ = qualScoresMap[readName];

			}

            writeRecord(bamFileOut, record);


        }

        close(bamFileOut);  // flushes
    }
    catch (Exception const & e)
    {
        //std::cout << "ERROR: " << e.what() << std::endl;
        std::cout << "ERROR: " << std::endl;
        return 1;
    }




    return 0;


}


















































