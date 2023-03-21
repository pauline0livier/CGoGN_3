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
 * sort eigen vectors by their eigen values
*/
Eigen::Vector3f sort_eigen_vectors(const Eigen::Matrix<float, 1, Eigen::Dynamic>& eigen_values,
		const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& eigen_vectors)
{
	std::vector<std::pair<int, int>> eigen_value_index_vector;
	for (int i = 0; i < eigen_values.size(); ++i)
	{
		eigen_value_index_vector.push_back(std::make_pair(eigen_values[i], i));
	}

	std::sort(std::begin(eigen_value_index_vector), 
				std::end(eigen_value_index_vector),
			  	std::greater<std::pair<int, int>>());

	auto sorted_eigen_values = Eigen::Matrix<float, 1, 
											Eigen::Dynamic>(eigen_values.cols());
	auto sorted_eigen_vectors =
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(
									eigen_vectors.rows(), eigen_vectors.cols());
	for (int i = 0; i < eigen_values.size(); ++i)
	{
		sorted_eigen_values[i] =
			eigen_values[eigen_value_index_vector[i].second]; 
		sorted_eigen_vectors.col(i) = eigen_vectors.col(eigen_value_index_vector[i].second);
	}

	return sorted_eigen_vectors.col(0);
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
		const Vec3& pos = value<Vec3>(m, attribute, v);

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
	});

	return std::make_pair(local_min, local_max);
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

} // namespace modeling

} // namespace cgogn

