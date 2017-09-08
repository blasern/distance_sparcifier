#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <queue>
#include <list>
#include <sstream>
#include <unordered_map>
#include <functional>

typedef float value_t;
typedef int64_t index_t;
typedef int16_t coefficient_t;

enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout> class compressed_distance_matrix {
public:
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	void init_rows();

	compressed_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(_distances), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
		assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	template <typename DistanceMatrix>
	compressed_distance_matrix(const DistanceMatrix& mat)
	    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
		init_rows();

		for (index_t i = 1; i < size(); ++i)
			for (index_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
	}

	value_t operator()(const index_t i, const index_t j) const;

	size_t size() const { return rows.size(); }
};

template <> void compressed_distance_matrix<LOWER_TRIANGULAR>::init_rows() {
	value_t* pointer = &distances[0];
	for (index_t i = 1; i < size(); ++i) {
		rows[i] = pointer;
		pointer += i;
	}
}

template <> void compressed_distance_matrix<UPPER_TRIANGULAR>::init_rows() {
	value_t* pointer = &distances[0] - 1;
	for (index_t i = 0; i < size() - 1; ++i) {
		rows[i] = pointer;
		pointer += size() - i - 2;
	}
}

template <> value_t compressed_distance_matrix<UPPER_TRIANGULAR>::operator()(index_t i, index_t j) const {
	if (i > j) std::swap(i, j);
	return i == j ? 0 : rows[i][j];
}

template <> value_t compressed_distance_matrix<LOWER_TRIANGULAR>::operator()(index_t i, index_t j) const {
	if (i > j) std::swap(i, j);
	return i == j ? 0 : rows[j][i];
}

typedef compressed_distance_matrix<LOWER_TRIANGULAR> compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR> compressed_upper_distance_matrix;

class euclidean_distance_matrix {
public:
	std::vector<std::vector<value_t>> points;

	euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points) : points(_points) {}

	value_t operator()(const index_t i, const index_t j) const {
		return std::sqrt(std::inner_product(points[i].begin(), points[i].end(), points[j].begin(), value_t(),
		                                    std::plus<value_t>(),
		                                    [](value_t u, value_t v) { return (u - v) * (u - v); }));
	}

	size_t size() const { return points.size(); }
};

enum file_format { LOWER_DISTANCE_MATRIX, UPPER_DISTANCE_MATRIX, DISTANCE_MATRIX, POINT_CLOUD, DIPHA };

template <typename T> T read(std::istream& s) {
	T result;
	s.read(reinterpret_cast<char*>(&result), sizeof(T));
	return result; // on little endian: boost::endian::little_to_native(result);
}

compressed_lower_distance_matrix read_point_cloud(std::istream& input_stream) {
	std::vector<std::vector<value_t>> points;

	std::string line;
	value_t value;
	while (std::getline(input_stream, line)) {
		std::vector<value_t> point;
		std::istringstream s(line);
		while (s >> value) {
			point.push_back(value);
			s.ignore();
		}
		if (!point.empty()) points.push_back(point);
		assert(point.size() == points.front().size());
	}

	euclidean_distance_matrix eucl_dist(std::move(points));

	index_t n = eucl_dist.size();

	// std::cout << "point cloud with " << n << " points in dimension " << eucl_dist.points.front().size() << std::endl;

	std::vector<value_t> distances;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_upper_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(distances)));
}

compressed_lower_distance_matrix read_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;

	std::string line;
	value_t value;
	for (int i = 0; std::getline(input_stream, line); ++i) {
		std::istringstream s(line);
		for (int j = 0; j < i && s >> value; ++j) {
			distances.push_back(value);
			s.ignore();
		}
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_dipha(std::istream& input_stream) {
	if (read<int64_t>(input_stream) != 8067171840) {
		std::cerr << "input is not a Dipha file (magic number: 8067171840)" << std::endl;
		exit(-1);
	}

	if (read<int64_t>(input_stream) != 7) {
		std::cerr << "input is not a Dipha distance matrix (file type: 7)" << std::endl;
		exit(-1);
	}

	index_t n = read<int64_t>(input_stream);

	std::vector<value_t> distances;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (i > j)
				distances.push_back(read<double>(input_stream));
			else
				read<double>(input_stream);

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_file(std::istream& input_stream, file_format format) {
	switch (format) {
	case LOWER_DISTANCE_MATRIX:
		return read_lower_distance_matrix(input_stream);
	case UPPER_DISTANCE_MATRIX:
		return read_upper_distance_matrix(input_stream);
	case DISTANCE_MATRIX:
		return read_distance_matrix(input_stream);
	case POINT_CLOUD:
		return read_point_cloud(input_stream);
	case DIPHA:
		return read_dipha(input_stream);
	}
}

void print_usage_and_exit(int exit_code) {
	std::cerr << "Usage: "
	          << "distance_sparcifier "
	          << "[options] [filename]" << std::endl
	          << std::endl
	          << "Options:" << std::endl
	          << std::endl
	          << "  --help           print this screen" << std::endl
	          << "  --format         use the specified file format for the input. Options are:" << std::endl
	          << "                     lower-distance (lower triangular distance matrix; default)" << std::endl
	          << "                     upper-distance (upper triangular distance matrix)" << std::endl
	          << "                     distance       (full distance matrix)" << std::endl
	          << "                     point-cloud    (point cloud in Euclidean space)" << std::endl
	          << "                     dipha          (distance matrix in DIPHA file format)" << std::endl
	          << "  --interleaving   interleaving constant >= 1" << std::endl
	          << std::endl;

	exit(exit_code);
}

compressed_lower_distance_matrix sparcify_distance(compressed_lower_distance_matrix dist,
						   float interleaving_const = 1.0,
						   int initial_point_index = 0) {
  /*
  Sparcify distance using farthest point sampling

  Parameters
  ==========
  compressed_lower_distance_matrix dist
    Original distance matrix
  float interleaving_const
    Desired interleaving constant (default: 1)
  int initial_point_index (default: 0)
    Index of point to start farthest point sampling with

  Return value
  ============
  The sparcified compressed_lower_distance_matrix
  */
  if (interleaving_const == 1){
    return dist;
  }
  // constants
  const int n = dist.size();
  const float inf = *max_element(dist.distances.begin(), dist.distances.end());
  const float delta = 1/interleaving_const;

  // initialize indices and insertion radii
  std::vector<int> possible_indices(n-1);
  std::iota(possible_indices.begin(), possible_indices.end(), 1);
  std::vector<int> indices;
  indices.push_back(initial_point_index);
  std::vector<float> insertion_radii;
  std::vector<float> final_radii;
  insertion_radii.push_back(inf);
  final_radii.push_back(2.0 * inf / (1-delta));
  
  // find ordered indices and insertion radii
  while(indices.size() < n){
    std::vector<float> indices_dist;
    for (auto i : possible_indices){
      std::vector<float> min_dist;
      for (auto j : indices){
	min_dist.push_back(dist(i, j));
      }
      float min = *min_element(min_dist.begin(), min_dist.end());
      indices_dist.push_back(min);
    }
    int possible_index = std::distance(indices_dist.begin(),
				       std::max_element(indices_dist.begin(),
							indices_dist.end()));
    int new_index = possible_indices.at(possible_index);
    possible_indices.erase(possible_indices.begin() + possible_index);
    indices.push_back(new_index);
    insertion_radii.push_back(*std::max_element(indices_dist.begin(), indices_dist.end()));
  }

  // update fps distance
  std::vector<float> fps_dist;
  std::vector<int> subindices = indices;
  for (auto i : indices){
    subindices.erase(subindices.begin());
    if (subindices.size() > 0){
      for (auto j : subindices){
	fps_dist.push_back(dist(i, j));
      }
    }
  }

  // update distances
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < i; ++j){
      int ij = j + (i * (i - 1)) / 2;
      float dij = fps_dist.at(ij);
      if (dij > (insertion_radii.at(i) + insertion_radii.at(j))/(1+delta)){
	fps_dist.at(ij) = inf;
      }
    }
  }

  // return sparse distance
  compressed_lower_distance_matrix spdist = compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(fps_dist)));
  return spdist;
}

int main(int argc, char** argv) {

	const char* filename = nullptr;

	file_format format = DISTANCE_MATRIX;

	float interleaving = 1.0;

	for (index_t i = 1; i < argc; ++i) {
		const std::string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--interleaving") {
		  std::string parameter = std::string(argv[++i]);
		  size_t next_pos;
		  interleaving = std::stof(parameter, &next_pos);		 
		  if (next_pos != parameter.size()) print_usage_and_exit(-1);
		  if (interleaving < 1) print_usage_and_exit(-1);
		} else if (arg == "--format") {
			std::string parameter = std::string(argv[++i]);
			if (parameter == "lower-distance")
				format = LOWER_DISTANCE_MATRIX;
			else if (parameter == "upper-distance")
				format = UPPER_DISTANCE_MATRIX;
			else if (parameter == "distance")
				format = DISTANCE_MATRIX;
			else if (parameter == "point-cloud")
				format = POINT_CLOUD;
			else if (parameter == "dipha")
				format = DIPHA;
			else
				print_usage_and_exit(-1);
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}

	compressed_lower_distance_matrix dist = read_file(filename ? file_stream : std::cin, format);

	index_t n = dist.size();

	// std::cout << "distance matrix with " << n << " points" << std::endl;

	auto value_range = std::minmax_element(dist.distances.begin(), dist.distances.end());
	// std::cout << "value range: [" << *value_range.first << "," << *value_range.second << "]" << std::endl;

	// sparcify distance matrix
	compressed_lower_distance_matrix sparse_dist = sparcify_distance(dist,
									 interleaving);

	// output sparse distance matrix
	for (int i = 0; i < n; ++i){
	  for (int j = 0; j < i; ++j){
	    std::cout << sparse_dist(i, j) << ",";
	  }
	  std::cout << std::endl;
	}

}

