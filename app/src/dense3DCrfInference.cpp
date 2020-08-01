/*
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the copyright holder nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURposE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE posSIBILITY OF SUCH DAMAGE.
*/
#include "itkio.h"
#include "miaImageProcessing.h"

#include <cstdio>
#include <cmath>

#include <vector>
#include <string>
#include <sstream> //for ostringstream.
#include <fstream>
#include "densecrf.h"

#include "ioparsing.hpp"
#include "parameters.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <chrono>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace mia;

std::vector<Image> load_multi_channel_image(std::vector<std::string> image_filenames)
{

	std::vector<Image> images;
	images.push_back(itkio::load(image_filenames[0]));

	for (int c = 1; c < image_filenames.size(); c++)
	{
		Image img = itkio::load(image_filenames[c]);
		if (img.size() != images[0].size())
		{
			std::cout << "Image modalities size mismatch! Aborting..." << std::endl;
			exit(0);
		}
		else
			images.push_back(img);
	}
	return images;
}

float cap_prob(float p)
{

	return std::fmax(0.05, std::fmin(p, 0.95));
}

//Loads the nii with the unary potentials and makes a matrix that the rest of the original code needs.
void setupUnaryEnergyMatrix(MatrixXf &unaryCosts, std::vector<Image> probmaps, int numberOfForegroundClasses, bool computeBackgroundProbmap = true)
{

	int numberOfVoxelsInImage = probmaps[0].size();
	int start_class_index = 0;
	int numberOfClasses = numberOfForegroundClasses + 1;
	if (computeBackgroundProbmap)
		start_class_index = 1;

	float eps = 0.001;
	unaryCosts.setZero();
	for (int c = 0; c < numberOfForegroundClasses; ++c)
	{

		int index = 0;
		for (int z = 0; z < probmaps[c].sizeZ(); z++)
		{
			for (int y = 0; y < probmaps[c].sizeY(); y++)
			{
				for (int x = 0; x < probmaps[c].sizeX(); x++)
				{
					float val = probmaps[c](x, y, z);

					val = -log(cap_prob(probmaps[c](x, y, z)));
					unaryCosts(c + start_class_index, index) = val;
					if (computeBackgroundProbmap)
						unaryCosts(0, index) += std::fmax(0, probmaps[c](x, y, z));
					index++;
				}
			}
		}
	}
	if (computeBackgroundProbmap)
	{
		int index = 0;
		for (int z = 0; z < probmaps[0].sizeZ(); z++)
		{
			for (int y = 0; y < probmaps[0].sizeY(); y++)
			{
				for (int x = 0; x < probmaps[0].sizeX(); x++)
				{

					unaryCosts(0, index) = -log(cap_prob(1 - unaryCosts(0, index)));
					index++;
				}
			}
		}
	}
}

void addPairwiseBilateralMultiMod(DenseCRF3D &crf3d, std::vector<float> sXYZ, std::vector<float> sMod, std::vector<Image> &images, LabelCompatibility *function = NULL, KernelType kernel_type = DIAG_KERNEL, NormalizationType normalization_type = NORMALIZE_SYMMETRIC)
{
	int numImages = images.size();
	int sizeX = images[0].sizeX();
	int sizeY = images[0].sizeY();
	int sizeZ = images[0].sizeZ();

	MatrixXf feature(3 + numImages, images[0].size());
	feature.setZero();
	int index = 0;
	for (int z = 0; z < images[0].sizeZ(); z++)
	{
		for (int y = 0; y < images[0].sizeY(); y++)
		{
			for (int x = 0; x < images[0].sizeX(); x++)
			{
				feature(0, index) = x / sXYZ[0];
				feature(1, index) = y / sXYZ[1];
				feature(2, index) = z / sXYZ[2];

				for (int imgIdx = 0; imgIdx < numImages; imgIdx++)
					feature(3 + imgIdx, index) = images[imgIdx](x, y, z) / sMod[imgIdx];

				index++;
			}
		}
	}
	crf3d.addPairwiseEnergy(feature, function, kernel_type, normalization_type);
}

int main(int argc, char *argv[])
{

	Parameters params;

	std::string config_file;

	std::vector<std::string> simages;
	std::vector<std::string> simagelists;
	std::vector<std::string> sprobmaps;
	std::vector<std::string> sprobmaplists;

	std::string output_path = ".";

	int numberOfModalities = 0;
	int numberOfForegroundClasses = 0;
	bool computeBackgroundProbmap = true;
	bool overwrite_output = false;
	bool do_output_probmaps = false;
	bool do_output_unary = false;
	std::vector<std::string> sbilateralXYZStds;
	std::vector<std::string> sposXYZStds;
	std::vector<std::string> sbilateralModStds;

	float minIntensity = -3, maxIntensity = +3;

	//====================================

	try
	{
		// Declare the supported options.
		po::options_description generic("generic options");
		generic.add_options()("help", "produce help message")("config", po::value<std::string>(&config_file), "configuration file")("probs", po::bool_switch(&do_output_probmaps)->default_value(false), "")("unary", po::bool_switch(&do_output_unary)->default_value(false), "")("overwrite", po::bool_switch(&overwrite_output)->default_value(false), "");

		po::options_description config("specific options");
		config.add_options()("output", po::value<std::string>(&output_path), "output path")("image", po::value<std::vector<std::string>>(&simages)->multitoken(), "filename(s) of multi-modality images")("imagelist", po::value<std::vector<std::string>>(&simagelists)->multitoken(), "text file(s) listing images")("probmap", po::value<std::vector<std::string>>(&sprobmaps)->multitoken(), "filename(s) of multi-modality images")("probmaplist", po::value<std::vector<std::string>>(&sprobmaplists)->multitoken(), "text file(s) listing images")("posXYZStds", po::value<std::vector<std::string>>(&sposXYZStds), "posXYZStds")("biModStds", po::value<std::vector<std::string>>(&sbilateralModStds), "bilateralModStds")("biXYZStds", po::value<std::vector<std::string>>(&sbilateralXYZStds), "bilateralXYZStds")("numberOfForegroundClasses", po::value<int>(&numberOfForegroundClasses), "numberOfForegroundClasses")("posW", po::value<float>(&params.posW), "posW")("biW", po::value<float>(&params.bilateralW), "biW")("maxIterations", po::value<int>(&params.maxIterations), "maxIteration");

		po::options_description cmdline_options("options");
		cmdline_options.add(generic).add(config);

		po::options_description config_file_options;
		config_file_options.add(config);

		po::variables_map vm;

		po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
		po::notify(vm);

		if (vm.count("config"))
		{
			std::ifstream ifs(config_file.c_str());
			if (!ifs)
			{
				std::cout << "cannot open config file: " << config_file << std::endl;
				return 0;
			}
			else
			{
				po::store(parse_config_file(ifs, config_file_options), vm);
				po::notify(vm);
			}
		}

		if (vm.count("help"))
		{
			std::cout << cmdline_options << std::endl;
			return 0;
		}
	}
	catch (std::exception &e)
	{
		std::cout << e.what() << std::endl;
		return 1;
	}

	if (!fs::exists(output_path)) fs::create_directories(output_path);

	//================================

	std::vector<std::string> images;
	strings_to_values(simages, images);

	std::vector<std::string> imagelists;
	strings_to_values(simagelists, imagelists);

	strings_to_values(sbilateralXYZStds, params.bilateralXYZStds, true);

	strings_to_values(sposXYZStds, params.posXYZStds, true);

	std::vector<std::string> probmaps;
	strings_to_values(sprobmaps, probmaps);

	std::vector<std::string> probmaplists;
	strings_to_values(sprobmaplists, probmaplists);

	std::vector<std::vector<std::string>> image_filenames = parse_multichannel_filepaths(images, imagelists);
	std::vector<std::vector<std::string>> probmap_filenames = parse_multichannel_filepaths(probmaps, probmaplists);

	int numberOfImages = image_filenames.size();
	int numberOfProbChannels = probmap_filenames.size();

	if (numberOfImages == 0)
	{
		std::cout << "No images provided. Aborting..." << std::endl;
	}
	if (numberOfProbChannels == 0)
	{
		std::cout << "No probmaps provided. Aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (numberOfProbChannels == numberOfForegroundClasses)
		computeBackgroundProbmap = true;
	else if (numberOfForegroundClasses > 0 && numberOfProbChannels == numberOfForegroundClasses + 1)
		computeBackgroundProbmap = false;
	else
	{
		std::cout << "number of probmaps mismatches with provided number of foreground classes" << std::endl;
		exit(EXIT_FAILURE);
	}

	numberOfModalities = image_filenames[0].size();
	int numClasses = numberOfForegroundClasses + 1;

	params.setNumberOfModalities(numberOfModalities);
	strings_to_values(sbilateralModStds, params.bilateralModsStds, true);

	std::vector<std::vector<Image>> intensity_images;

	params.print();

	for (int imgIdx = 0; imgIdx < numberOfImages; imgIdx++)
	{

		fs::path input_path(image_filenames[imgIdx][0]);

		std::string basename = fs::basename(input_path);
		if (fs::extension(basename) != "")
			basename = fs::basename(basename);
		std::stringstream base_output_sfilename;
		std::stringstream labelmap_output_path;

		base_output_sfilename << output_path << "/" << basename;
		labelmap_output_path << base_output_sfilename.str() << "_labelmap.nii.gz";

		if (fs::exists(labelmap_output_path.str()))
		{
			if (!overwrite_output)
			{
				std::cout << labelmap_output_path.str() << " already exists. skipping." << std::endl;
				continue;
			}
			else
				std::cout << labelmap_output_path.str() << " already exists. Overwriting." << std::endl;
		}
		std::cout << labelmap_output_path.str() << std::endl;

		namespace ch = std::chrono;
		auto start = ch::high_resolution_clock::now();

		std::cout << "loading data " << imgIdx + 1 << " of " << image_filenames.size() << " with " << numberOfModalities << " channel(s)" << std::endl;

		std::vector<Image> multiChannelImage = load_multi_channel_image(image_filenames[imgIdx]);
		std::vector<Image> imageProbmaps = load_multi_channel_image(probmap_filenames[imgIdx]);

		int sizeX = multiChannelImage[0].sizeX();
		int sizeY = multiChannelImage[0].sizeY();
		int sizeZ = multiChannelImage[0].sizeZ();

		DenseCRF3D crf3d(sizeX, sizeY, sizeZ, numberOfForegroundClasses + 1);

		std::cout << "********** Setting up Unary. **********" << std::endl;
		MatrixXf unaryCost = MatrixXf(numClasses, multiChannelImage[0].size());
		setupUnaryEnergyMatrix(unaryCost, imageProbmaps, numberOfForegroundClasses, computeBackgroundProbmap);
		crf3d.setUnaryEnergy(unaryCost);

		std::cout << "********** Setting up Pairwise potentials **********" << std::endl;
		crf3d.addPairwiseGaussian(params.posXYZStds[0],
								  params.posXYZStds[1],
								  params.posXYZStds[2],
								  new PottsCompatibility(params.posW));

		addPairwiseBilateralMultiMod(crf3d,
									 params.bilateralXYZStds,
									 params.bilateralModsStds,
									 multiChannelImage,
									 new PottsCompatibility(params.bilateralW));

		std::cout << "++++++++++++++++++++++++++ Performing Inference ++++++++++++++++++++++++++" << std::endl;

		MatrixXf probMapsMatrix = crf3d.inference(params.maxIterations);

		if (do_output_probmaps)
		{
			Image output_probmap = imageProbmaps[0].clone();
			for (int i = 0; i < numClasses; ++i)
			{
				std::stringstream output_filename;
				output_filename << base_output_sfilename.str() << "_probmap" << i << ".nii.gz";
				int index = 0;
				for (int z = 0; z < imageProbmaps[0].sizeZ(); z++)
				{
					for (int y = 0; y < imageProbmaps[0].sizeY(); y++)
					{
						for (int x = 0; x < imageProbmaps[0].sizeX(); x++)
						{
							output_probmap(x, y, z) = probMapsMatrix(i, index++);
						}
					}
				}
				itkio::save(output_probmap, output_filename.str());
			}
		}

		if (do_output_unary)
		{
			Image output_probmap = imageProbmaps[0].clone();
			for (int i = 0; i < numClasses; ++i)
			{
				std::stringstream output_filename;
				output_filename << base_output_sfilename.str() << "_unary" << i << ".nii.gz";
				int index = 0;
				for (int z = 0; z < imageProbmaps[0].sizeZ(); z++)
				{
					for (int y = 0; y < imageProbmaps[0].sizeY(); y++)
					{
						for (int x = 0; x < imageProbmaps[0].sizeX(); x++)
						{
							output_probmap(x, y, z) = unaryCost(i, index++);
						}
					}
				}
				itkio::save(output_probmap, output_filename.str());
			}
		}

		VectorXs segmentationVector = crf3d.currentMap(probMapsMatrix);

		Image output_labelmap = imageProbmaps[0].clone();
		output_labelmap.dataType(mia::USHORT);
		int index = 0;
		for (int z = 0; z < imageProbmaps[0].sizeZ(); z++)
		{
			for (int y = 0; y < imageProbmaps[0].sizeY(); y++)
			{
				for (int x = 0; x < imageProbmaps[0].sizeX(); x++)
				{
					output_labelmap(x, y, z) = segmentationVector[index++];
				}
			}
		}

		itkio::save(output_labelmap, labelmap_output_path.str());

		std::cout << "++++++++++++++++++++++++++ Done. ++++++++++++++++++++++++++" << std::endl;
		auto stop = ch::high_resolution_clock::now();
		std::cout << "done. took " << ch::duration_cast<ch::milliseconds>(stop - start).count() << " ms" << std::endl;
	}
}
