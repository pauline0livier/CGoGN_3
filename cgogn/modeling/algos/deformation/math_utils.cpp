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

double distance_vec3(const Vec3& p1, const Vec3& p2)
{
	// compute Euclidean distance or whatever
	return (p1 - p2).norm();
}

// from https://github.com/diegomazala/pca/blob/master/src/pca.h
Eigen::Vector3f sort_eigen_vectors(const Eigen::Matrix<float, 1, Eigen::Dynamic>& eigen_values,
								   const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& eigen_vectors)
{
	// Stuff below is done to sort eigen values. This can be done in other ways too.
	std::vector<std::pair<int, int>> eigen_value_index_vector;
	for (int i = 0; i < eigen_values.size(); ++i)
	{
		eigen_value_index_vector.push_back(std::make_pair(eigen_values[i], i));
	}

	std::sort(std::begin(eigen_value_index_vector), std::end(eigen_value_index_vector),
			  std::greater<std::pair<int, int>>());

	auto sorted_eigen_values = Eigen::Matrix<float, 1, Eigen::Dynamic>(eigen_values.cols());
	auto sorted_eigen_vectors =
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(eigen_vectors.rows(), eigen_vectors.cols());
	for (int i = 0; i < eigen_values.size(); ++i)
	{
		sorted_eigen_values[i] =
			eigen_values[eigen_value_index_vector[i].second]; // can also be eigen_value_index_vector[i].first
		sorted_eigen_vectors.col(i) = eigen_vectors.col(eigen_value_index_vector[i].second);
	}

	return sorted_eigen_vectors.col(0);
}

std::tuple<Vec3, Vec3, Vec3> get_extended_bounding_box(const Vec3& bb_min, const Vec3& bb_max,
													   const float& extension_factor)
{

	Vec3 center = (bb_min + bb_max) / Scalar(2);
	Vec3 e_bb_min = ((bb_min - center) * extension_factor) + center;
	Vec3 e_bb_max = ((bb_max - center) * extension_factor) + center;

	return std::make_tuple(e_bb_min, e_bb_max, center);
}

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

std::pair<Vec3, Vec3> get_border_values_in_set(const CMap2& m, const CMap2::Attribute<Vec3>* attribute, cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set)
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

CMap2::Vertex closest_vertex_in_set_from_value(const CMap2& m, const CMap2::Attribute<Vec3>* vertex_position,
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

float projection_on_segment(const Vec3& A, const Vec3& B, const Vec3& P){

	const Vec3 AB = B - A;
 	const Vec3 AP = P - A; 
  	return AB.dot(AP) / AB.squaredNorm(); 
	//return AB.dot(AP) / (AB.dot(AB));
}

} // namespace modeling

} // namespace cgogn

/*

Eigen::Matrix3f covariance_matrix;
		// inspired from https://gist.github.com/atandrau/847214/882418ab34737699a6b1394d3a28c66e2cc0856f
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				covariance_matrix(i, j) = 0.0;
				control_set->foreach_cell([&](MeshVertex v) {
					const Vec3& pos = value<Vec3>(m, vertex_position, v);
					covariance_matrix(i, j) += (center[i] - pos[i]) * (center[j] - pos[j]);
				});

				covariance_matrix(i, j) /= control_set->size() - 1;
			}

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> eigen_solver(
			covariance_matrix);
		Eigen::Matrix<float, 1, Eigen::Dynamic> eigen_values = eigen_solver.eigenvalues();
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> eigen_vectors = eigen_solver.eigenvectors();

		Eigen::Vector3f main_eigen_vector = cgogn::modeling::sort_eigen_vectors(eigen_values, eigen_vectors);
		main_eigen_vector.normalize();

		Vec3 main_direction = {main_eigen_vector[0], main_eigen_vector[1], main_eigen_vector[2]};
*/