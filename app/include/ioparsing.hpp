#pragma once

#include <vector>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

template <typename T>
void strings_to_values(const std::vector<std::string>& string_seq, std::vector<T>& values, bool overwrite=false)
{
  if(string_seq.size() > 0 && overwrite)
    values.resize(0);
  for (std::vector<std::string>::const_iterator it = string_seq.begin(); it != string_seq.end(); ++it)
  {
    std::cout << *it << std::endl;
    std::stringstream ss(*it);
    std::copy(std::istream_iterator<T>(ss), std::istream_iterator<T>(), back_inserter(values));
  }
}

void params_to_vector(const std::vector<double>& params, Eigen::Vector3d& vector)
{
  if (params.size() > 0)
  {
    vector[0] = (params.size() > 0) ? params[0] : 0;
    vector[1] = (params.size() > 1) ? params[1] : params[0];
    vector[2] = (params.size() > 2) ? params[2] : params[0];
  }
  else
  {
    vector.setZero();
  }
}

template <typename T>
std::vector<std::vector<T>> transpose_2d_vector(std::vector<std::vector<T>>& vector2d) {
	int r = vector2d.size();
	if (r > 0 ){
	int c = vector2d[0].size();

	std::vector<std::vector<T>> tvect(c, std::vector<T>(r));
	for(int i = 0; i < r; i++)
			for(int j = 0; j < c; j++){
			
				tvect[j][i] = vector2d[i][j];
			}
			
	return tvect;
}
else 
return vector2d;
}


std::vector<std::string> parse_image_filepaths(std::string imageListPath, std::string singleImagePath) {

  std::vector<std::string> image_filenames;
  std::string image_filename;
  std::ifstream ifs_image(imageListPath);
  if (imageListPath.size() > 0)
  {
    while (getline(ifs_image, image_filename))
    {
      boost::trim(image_filename);
      if (image_filename != "")
      {
        image_filenames.push_back(image_filename);
      }
    }
  }
  if (singleImagePath != "") image_filenames.push_back(singleImagePath);

  return image_filenames;
}


std::vector<std::vector<std::string>> parse_multichannel_filepaths(std::vector<std::string> images, std::vector<std::string> imagelists) {

  int channels = std::max(imagelists.size(), images.size());

  std::vector<std::vector<std::string>> image_filenames(channels);
  for (int c = 0; c < channels; c++) {
    if (imagelists.size() > 0)
      image_filenames[c] =  parse_image_filepaths(imagelists[c], "");
    if (images.size() == channels) { image_filenames[c].push_back(images[c]); }
  }

  return  transpose_2d_vector(image_filenames);
}


