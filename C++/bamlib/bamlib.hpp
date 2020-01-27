/*
 * bamlib.hpp
 *
 *  Created on: Mar 11, 2017
 *      Author: eliot
 */


#include <string>
#include <zlib.h>
#include <pthread.h>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <sstream>

#include <stdexcept>
#include <iomanip>
#include <array>

#define INVARIANTS_BUFFER_WIDTH 5

#define CHROMOSOME_IDX 0
#define READ_BEGIN_POS 1
#define MAPQ 2
#define INSERT_SIZE 3
#define GAP_EXTENSIONS 4
#define EDIT_DISTANCE 5
#define MISMATCHES 6
#define GAP_OPENS 7
#define ALIGNMENT_SCORE 8
#define SECONDARY_ALIGNMENT_SCORE 9
#define READ_COMPLEXITY_SCORE 10
#define READ_GC_CONTENT 11
#define PAIR_ORIENTATION 12
#define PAIR_ALIGNMENT_TYPE 13
#define INTERCEPT 14
#define SLOPE 15
#define R_VALUE 16
#define DEPTH 17
#define READ_EIGENVAL_1 18
#define READ_EIGENVAL_2 19
#define READ_TRACE 20
#define READ_DET 21
#define REF_EIGENVAL_1 22
#define REF_EIGENVAL_2 23
#define REF_TRACE 24
#define REF_DET 25
#define REF_COMPLEXITY_SCORE 26
#define REF_GC_CONTENT 27
#define N_LOW_QUAL_BASES 28
#define AVG_LOW_QUAL_BASE_SCORE 29
#define MATE_ALIGNMENT_SCORE 30
#define ALIGNMENT_SCORES_DIFF 31
#define CLASS 32
#define SNPREADCLASS 33

void calc_invariants_cpu(const unsigned int height, const unsigned int width, float* matrix, const unsigned* dna_seqs);
int getTrainingFeatures(std::vector<unsigned char>& trainingData, int readLength, int nReads, bool noMultiMap, std::string bamFileName);
int sampleBamFile(std::string inBamFileName, std::string outBamFileName);
int getFeatureVectors(std::vector<unsigned char>& readSequenceData, std::vector<unsigned char>& chromosomeIntervalSequences, std::vector<unsigned char>& readQualities, std::vector<std::array<float, 34>>& featureVectors, std::vector<std::string>& readNames, int chromosomeIdx, int readLength, int offSet, std::string bamFileName, std::string sampledBamFileName, std::string fastaFileName, int outputReadClass);
int recalibrateBamFile(std::string scoresFile, std::string inBamFileName, std::string outBamFileName);
//int buildSnpIntervalTree(IntervalTree<int, CharString>& tree, std::string vcfFileName, int chromosomeIndex, int readLength);
//int isGoodReadSnp(BamAlignmentRecord record, int readLength, IntervalTree<int, CharString> snpIntervals);
