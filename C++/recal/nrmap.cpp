/*
 * nrmap.cpp
 *
 *  Created on: Mar 11, 2017
 *      Author: eliot
 */


#include <stdlib.h>
#include "bamlib.hpp"
#include "nrmap.h"
#include "OptionPrinter.hpp"

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"

namespace
{
  const size_t ERROR_IN_COMMAND_LINE = 1;
  const size_t SUCCESS = 0;
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace


void print_usage() {
   std::cout << "Usage: nrmap [ap] -l num -n num\n";
}


int main(int argc, char const * argv[])
{

  try
  {
    std::string appName = boost::filesystem::basename(argv[0]);
    int add = 0;
    int like = 0;
    std::vector<std::string> sentence;

	int readLength = 0;
	int nReads = 0;
	std::string inCommand, bamOutFileName, fastaFileName, scoresFileName, depthBamFileName;
	int a = 0;
	int p = 0;
	int option = 0;
        std::string bamFileName;
    /** Define and parse the program options
     */
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Print help messages")
      ("verbose,v", "Print all messages and notifcations")
      //("word,w", po::value<std::vector<std::string> >(&sentence),
      // "words for the sentence, specify multiple times")
      //(",t", "just a temp option that does very little")
      ("command,c", po::value<std::string>(&inCommand)->required(), "Command")
      ("bam_file_in,b", po::value<std::string>(&bamFileName)->required(), "Input BAM file")
      ("bam_file_out,o", po::value<std::string>(&bamOutFileName)->required(), "Output BAM file")
      ("fasta_file_in,f", po::value<std::string>(&fastaFileName)->required(), "FASTA reference")
      ("scores_file_in,s", po::value<std::string>(&scoresFileName)->required(), "TAB delimited input qual scores for recalibration")
      ("depth_bam_in,d", po::value<std::string>(&depthBamFileName)->required(), "BAM file for estimating read depth")
      //("manual,m", po::value<std::string>(), "extract value manually")
      ("read_length,l", po::value<int>(&readLength)->required(), "Length of reads in BAM file (they must all be of equal length!)");

    po::positional_options_description positionalOptions;
    positionalOptions.add("add", 1);
    positionalOptions.add("like", 1);

    po::variables_map vm;

    try
    {
      po::store(po::command_line_parser(argc, argv).options(desc)
                  .positional(positionalOptions).run(),
                vm); // throws on error

      /** --help option
       */
      if ( vm.count("help")  )
      {
        std::cout << "The program extracts training features from a BAM file."
                  << " Some of those features are further post-processed "
                  << " on the GPU." << std::endl << std::endl;
        rad::OptionPrinter::printStandardAppDesc(appName,
                                                 std::cout,
                                                 desc,
                                                 &positionalOptions);
        return SUCCESS;
      }

      po::notify(vm); // throws on error, so do after help in case
                      // there are any problems
    }
    catch(boost::program_options::required_option& e)
    {
      rad::OptionPrinter::formatRequiredOptionError(e);
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      rad::OptionPrinter::printStandardAppDesc(appName,
                                               std::cout,
                                               desc,
                                               &positionalOptions);
      return ERROR_IN_COMMAND_LINE;
    }
    catch(boost::program_options::error& e)
    {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      rad::OptionPrinter::printStandardAppDesc(appName,
                                               std::cout,
                                               desc,
                                               &positionalOptions);
      return ERROR_IN_COMMAND_LINE;
    }



    //examples of handling various command options
    // can do this without fear because it is required to be present
/*
    std::cout << "Necessary = "
              << vm["necessary"].as<std::string>() << std::endl;

    if ( vm.count("verbose") )
    {
      std::cout << "VERBOSE PRINTING" << std::endl;
    }
    if (vm.count("verbose") && vm.count("t"))
    {
      std::cout << "heeeeyahhhhh" << std::endl;
    }

    std::cout << "Required Positional, add: " << add
              << " like: " << like << std::endl;

    if ( sentence.size() > 0 )
    {
      std::cout << "The specified words: ";
      std::string separator = " ";
      if (vm.count("verbose"))
      {
        separator = "__**__";
      }
      for(size_t i=0; i<sentence.size(); ++i)
      {
        std::cout << sentence[i] << separator;
      }
      std::cout << std::endl;

    }

    if ( vm.count("manual") )
    {
      std::cout << "Manually extracted value: "
                << vm["manual"].as<std::string>() << std::endl;
    }
 */



///APP CODE HERE!!!

    if(inCommand == "recal")
    {
    	recalibrateBamFile(scoresFileName, bamFileName, bamOutFileName);
    	return 0;
    }

    if(inCommand == "samplebam")
    {
    	int ret = sampleBamFile(bamFileName, bamOutFileName);
    	return 0;
    }

	std::vector<unsigned char> trainingSet;
	trainingSet.reserve(readLength*nReads);

	std::vector<unsigned char> chromosomeIntervalSequences;


	int offSet = 125;

	std::vector<unsigned char> readQualities;

	std::vector<std::array<float, 34>> featureVectors;
	std::vector<std::string> readNames;
	std::vector<float> floatQuals;

	for(int chromosomeIdx = 0;chromosomeIdx < 13;++chromosomeIdx)
	{
		trainingSet.clear();
		chromosomeIntervalSequences.clear();
		readQualities.clear();
		featureVectors.clear();
		readNames.clear();
		floatQuals.clear();


		int foo = getFeatureVectors(trainingSet, chromosomeIntervalSequences, readQualities, featureVectors, readNames, chromosomeIdx, readLength, offSet, bamFileName, depthBamFileName, fastaFileName, 1);
/*
		for(int i = 0;i<readQualities.size();++i)
			floatQuals.push_back((float)readQualities[i]);

		std::cerr << "Computing least squares.... ";

		float* corrcoeffs = new float[readQualities.size()/readLength];
		float* lsq = new float[floatQuals.size()/readLength*2];

		gpuMatMult(1, readLength, floatQuals.data(), readLength, lsq, 2, floatQuals.size()/readLength, readLength);

		std::cerr << "Done. " << std::endl;

		std::cerr << "Computing correlation coefficients.... ";

		getCorrelationCoefficients(readQualities.size()/readLength, readLength, corrcoeffs, readQualities.data());


		std::cerr << "Done. " << std::endl;


		for (int i=0;i<featureVectors.size();++i)
		{
			featureVectors[i][R_VALUE] = corrcoeffs[i];
			featureVectors[i][SLOPE] = lsq[i*2];
			featureVectors[i][INTERCEPT] = lsq[i*2 + 1];
		}

*/
		for (int i=0;i<featureVectors.size();++i) {
			std::cout << readNames[i] << "\t";
			for (int j=0;j<33;++j)
				std::cout << featureVectors[i][j] << "\t";
			std::cout << std::endl;
		}

	}


  }



  catch(std::exception& e)
  {
    std::cerr << "Unhandled Exception reached the top of main: "
              << e.what() << ", application will now exit" << std::endl;
    return ERROR_UNHANDLED_EXCEPTION;

  }

  return SUCCESS;



}
