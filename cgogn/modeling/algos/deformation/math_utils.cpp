/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/modeling/algos/deformation/math_utils.h>
#include <iostream>

namespace cgogn
{

namespace modeling
{

double getAngleBetweenUnitVectors(const Vec3& a, const Vec3& b)
{
	return 2.0 * asin((a - b).norm() / 2.0);
}

/**
 * Call the Vec3 norm function on the difference of vectors
*/
double distance_vec3(const Vec3& p1, const Vec3& p2)
{
	return (p1 - p2).norm();
}

/**
 * compute main direction from provided set
 * use Principal Component Analysis (PCA)
*/
std::tuple<Vec3, Vec3, Vec3> find_main_directions_from_set(const CMap2& m, 
							const CMap2::Attribute<Vec3>* attribute, 
						cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set,
						const Vec3& center)
{
	Eigen::Matrix3d covariance_matrix;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			covariance_matrix(i, j) = 0.0;

			control_set->foreach_cell([&](CMap2::Vertex v) {
				const Vec3& position = value<Vec3>(m, attribute, v);

					covariance_matrix(i, j) += 
							(center[i] - position[i]) * (center[j] - position[j]);
				});

			covariance_matrix(i, j) /= control_set->size() - 1;
		}
	}

	Eigen::SelfAdjointEigenSolver<
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> 
												eigen_solver(covariance_matrix);
	
	Eigen::Matrix<double, 1, Eigen::Dynamic> eigen_values = 
													eigen_solver.eigenvalues();

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_vectors = 
													eigen_solver.eigenvectors();

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sorted_eigen_vectors = 
								sort_eigen_vectors(eigen_values, eigen_vectors);
	
	Eigen::Vector3d main_eigen_vector = sorted_eigen_vectors.col(0); 
	main_eigen_vector.normalize();

	Eigen::Vector3d second_eigen_vector = sorted_eigen_vectors.col(1); 
	second_eigen_vector.normalize();

	Eigen::Vector3d third_eigen_vector = sorted_eigen_vectors.col(2); 
	third_eigen_vector.normalize();

	const Vec3 cross_product = main_eigen_vector.cross(second_eigen_vector);  
	if (cross_product == third_eigen_vector)
	{
		return std::make_tuple(main_eigen_vector, second_eigen_vector, third_eigen_vector);  
	} 
	else 
	{
		return std::make_tuple(main_eigen_vector, second_eigen_vector, -third_eigen_vector);
	}

}

/**
 * sort eigen vectors by their eigen values
*/
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sort_eigen_vectors(
				const Eigen::Matrix<double, 1, Eigen::Dynamic>& eigen_values,
	const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& eigen_vectors)
{
	std::vector<std::pair<int, int>> eigen_value_index_vector;
	for (int i = 0; i < eigen_values.size(); ++i)
	{
		eigen_value_index_vector.push_back(std::make_pair(eigen_values[i], i));
	}

	std::sort(std::begin(eigen_value_index_vector), 
				std::end(eigen_value_index_vector),
				std::greater<std::pair<int, int>>());

	auto sorted_eigen_values = Eigen::Matrix<double, 1, 
											Eigen::Dynamic>(eigen_values.cols());
	auto sorted_eigen_vectors =
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(
									eigen_vectors.rows(), eigen_vectors.cols());
	for (int i = 0; i < eigen_values.size(); ++i)
	{
		sorted_eigen_values[i] =
			eigen_values[eigen_value_index_vector[i].second]; 
		sorted_eigen_vectors.col(i) = eigen_vectors.col(eigen_value_index_vector[i].second);
	}

	return sorted_eigen_vectors;
}

/**
 * scale the bounding box by extension factor 
 * and relative to the bounding box center
 * @param {Vec3} bb_min 
 * @param {Vec3} bb_max
 * @param {float} extension_factor 
 * @returns a tuple with extended_bounding_box minimum and maximum positions 
 * 			and center
*/
std::tuple<Vec3, Vec3, Vec3> get_extended_bounding_box(const Vec3& bb_min, const Vec3& bb_max,
													   const float& extension_factor)
{

	Vec3 center = (bb_min + bb_max) / Scalar(2);
	Vec3 e_bb_min = ((bb_min - center) * extension_factor) + center;
	Vec3 e_bb_max = ((bb_max - center) * extension_factor) + center;

	return std::make_tuple(e_bb_min, e_bb_max, center);
}

/**
 * loop on the CMap2 vertices to compute the mean of the attribute
 * valid for attribute in Vec3
 * @returns {Vec3} mean value 
*/
Vec3 get_mean_value_attribute_from_set(const CMap2& m, const CMap2::Attribute<Vec3>* attribute,
									   cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set)
{
	Vec3 mean_value = {0.0, 0.0, 0.0};

	control_set->foreach_cell([&](CMap2::Vertex v) {
		const Vec3& pos = value<Vec3>(m, attribute, v);

		mean_value += pos;
	});

	mean_value /= control_set->size();

	return mean_value;
}

/// @brief mean value of std vector of Vec3
/// @param positions 
/// @return mean value 
Vec3 get_mean_value_in_array_Vec3(const std::vector<Vec3> positions)
{
	if (positions.empty())
	{
		return Vec3();
	}

	std::size_t count = positions.size(); 

	Vec3 mean_value = {0.0, 0.0, 0.0}; 
	for (std::size_t p = 0; p < count ; p++)
	{
		mean_value += positions[p]; 
	}

	mean_value /= count; 
	return mean_value; 
}

/**
 * loop on the CMap2 vertices to compute the extrema values of a set
 * valid for attribute in Vec3
 * @returns {pair<Vec3, Vec3>} extrema values
*/
std::pair<Vec3, Vec3> get_border_values_in_set(const CMap2& m, 
								const CMap2::Attribute<Vec3>* attribute, 
						cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set)
{ 
	Vec3 local_min = {1000.0, 1000.0, 1000.0};
	Vec3 local_max = {0.0, 0.0, 0.0};

	control_set->foreach_cell([&](CMap2::Vertex v) {
		const Vec3& position = value<Vec3>(m, attribute, v);

		for (size_t j = 0; j < 3; j++)
		{
			if (position[j] < local_min[j])
			{
				local_min[j] = position[j];
			}

			if (position[j] > local_max[j])
			{
				local_max[j] = position[j];
			}
		}
	});

	return std::make_pair(local_min, local_max);
}

/**
 * get extrema values in provided set, in local frame  
*/
std::pair<Vec3,Vec3> get_border_values_in_set_along_local_frame(const CMap2& m, 
					const CMap2::Attribute<Vec3>* attribute, 
					cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set, 
					const std::tuple<Vec3, Vec3, Vec3>& main_directions)
{
	Vec3 local_min = {1000.0, 1000.0, 1000.0};
	Vec3 local_max = {-1000.0, -1000.0, -1000.0};

	std::vector<Vec3> main_directions_vector(3); 
	main_directions_vector[0] = std::get<0>(main_directions); 
	main_directions_vector[1] = std::get<1>(main_directions); 
	main_directions_vector[2] = std::get<2>(main_directions);

	control_set->foreach_cell([&](CMap2::Vertex v) {
		const Vec3& position = value<Vec3>(m, attribute, v);

		for (size_t j = 0; j < 3; j++)
		{
			const double projection_j = position.dot(main_directions_vector[j]);
			 
			if (projection_j < local_min[j])
			{
				local_min[j] = projection_j;
			}

			if (projection_j > local_max[j])
			{
				local_max[j] = projection_j;
			}
		}
	});

	Vec3 min_position = {0.0, 0.0, 0.0}, 
	max_position = {0.0, 0.0, 0.0};

	for (size_t j = 0; j < 3; j++)
	{
		min_position += local_min[j]*main_directions_vector[j]; 
		max_position += local_max[j]*main_directions_vector[j]; 	
	} 

	return std::make_pair(min_position, max_position);
}

/**
 * compute extrema values from a array of Vec3
 * @returns {pair<Vec3, Vec3>} extrema values
*/
std::pair<Vec3, Vec3> get_border_values_in_array_Vec3(
											const std::vector<Vec3> positions)
{ 
	Vec3 local_min = {1000.0, 1000.0, 1000.0};
	Vec3 local_max = {0.0, 0.0, 0.0};

	for (uint32_t i = 0; i < positions.size(); i++) {
		const Vec3& pos = positions[i];

		for (size_t j = 0; j < 3; j++)
		{
			if (pos[j] < local_min[j])
			{
				local_min[j] = pos[j];
			}

			if (pos[j] > local_max[j])
			{
				local_max[j] = pos[j];
			}
		}
	}

	return std::make_pair(local_min, local_max);
}

/**
 * closest vertex in set from target position 
 * search for closest distance between vertex position and target one
 * useful to position handle at the user's chosen position
 * @returns {CMap2::Vertex} closest vertex
*/
CMap2::Vertex closest_vertex_in_set_from_value(const CMap2& m, 
						const CMap2::Attribute<Vec3>* vertex_position,
						cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set,
											const Vec3& target_position)
{

	CMap2::Vertex closest_vertex;

	double min_dist = 1000000;
	control_set->foreach_cell([&](CMap2::Vertex v) {
		const Vec3& pos = value<Vec3>(m, vertex_position, v);

		double dist = (pos - target_position).squaredNorm();
		if (dist < min_dist)
		{
			min_dist = dist;
			closest_vertex = v;
		}
	});

	return closest_vertex; 
}

/**
 * projection of a point on a segment 
 * @returns {float} distance of A to projection of P on [AB]
*/
float projection_on_segment(const Vec3& A, const Vec3& B, const Vec3& P){

	const Vec3 AB = B - A;
	const Vec3 AP = P - A; 
	return AB.dot(AP) / AB.squaredNorm(); 
}

/// @brief get projection of point on direction 
/// @param point 
/// @param direction 
/// @return dot product of point and direction 
double get_projection_on_direction(const Vec3& point, const Vec3& direction)
{
	return point.dot(direction); 
}

/// @brief check if target projection is inside delimited area
/// @param projection_value 
/// @param min_border 
/// @param max_border 
/// @return true if inside or on the borders of area, false otherwise 
bool check_projection_in_area(const double& projection_value, const double& min_border, const double& max_border)
{
	return (projection_value >= min_border && projection_value <= max_border); 
}

/// @brief get index of virtual cube 
/// cubes already assigned an index 
/// @param projection_values 
/// @param valid_values 
/// @param max_area_values 
/// @return 
size_t get_index_virtual_cube(const std::vector<double> projection_values, const std::vector<bool> valid_values, const std::vector<double> max_area_values)
{
	if (valid_values[0])
		{
			if (valid_values[1])
			{
				if (projection_values[2] > max_area_values[2])
				{
					return 5;
				}
				else
				{
					return 4;
				}
			}
			else if (valid_values[2])
			{
				if (projection_values[1] > max_area_values[1])
				{
					return 3;
				}
				else
				{
					return 2;
				}
			}
			else
			{
				if (projection_values[1] > max_area_values[1])
				{
					if (projection_values[2] > max_area_values[2])
					{
						return 17;
					}
					else
					{
						return 15;
					}
				}
				else
				{
					if (projection_values[2] > max_area_values[2])
					{
						return 16;
					}
					else
					{
						return 14;
					}
				}
			}
		}
		else if (valid_values[1])
		{
			if (valid_values[2])
			{
				if (projection_values[0] > max_area_values[0])
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}
			else
			{
				if (projection_values[0] > max_area_values[0])
				{
					if (projection_values[2] > max_area_values[2])
					{
						return 9;
					}
					else
					{
						return 7;
					}
				}
				else
				{
					if (projection_values[2] > max_area_values[2])
					{
						return 8;
					}
					else
					{
						return 6;
					}
				}
			}
		}
		else if (valid_values[2])
		{
			if (projection_values[1] > max_area_values[1])
			{
				if (projection_values[0] > max_area_values[0])
				{
					return 13;
				}
				else
				{
					return 11;
				}
			}
			else
			{
				if (projection_values[0] > max_area_values[0])
				{
					return 12;
				}
				else
				{
					return 10;
				}
			}
		}
		else
		{
			if (projection_values[0] > max_area_values[0])
			{
				if (projection_values[1] > max_area_values[1])
				{
					if (projection_values[2] > max_area_values[2])
					{
						return 25;
					}
					else
					{
						return 23;
					}
				}
				else
				{
					if (projection_values[2] > max_area_values[2])
					{
						return 21;
					}
					else
					{
						return 19;
					}
				}
			}
			else
			{
				if (projection_values[1] > max_area_values[1])
				{
					if (projection_values[2] > max_area_values[2])
					{
						return 24;
					}
					else
					{
						return 22;
					}
				}
				else
				{
					if (projection_values[2] > max_area_values[2])
					{
						return 20;
					}
					else
					{
						return 18;
					}
				}
			}
		}
}

} // namespace modeling

} // namespace cgogn

