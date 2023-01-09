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
#include <cgogn/modeling/algos/deformation/deformation_utils.h>
#include <iostream>

namespace cgogn
{

namespace modeling
{

///////////
// CAGE //
///////////

double getAngleBetweenUnitVectors(const Vec3& a, const Vec3& b)
{
	return 2.0 * asin((a - b).norm() / 2.0);
}

float compute_mvc(const Vec3& surface_point, Dart vertex, CMap2& cage, const Vec3& cage_point,
				  CMap2::Attribute<Vec3>* cage_position)
{

	double r = (cage_point - surface_point).norm();

	if (r == 0)
	{
		std::cerr << "Alert" << std::endl;
	}

	double sumU(0.);

	Dart it = vertex;

	do
	{
		Vec3 vi = value<Vec3>(cage, cage_position, CMap2::Vertex(it));
		Vec3 vj = value<Vec3>(cage, cage_position, CMap2::Vertex(phi1(cage, it)));
		Vec3 vk = value<Vec3>(cage, cage_position, CMap2::Vertex(phi_1(cage, it)));

		Vec3 ei = (vi - surface_point).normalized();
		Vec3 ej = (vj - surface_point).normalized();
		Vec3 ek = (vk - surface_point).normalized();

		double Bjk = getAngleBetweenUnitVectors(ej, ek);
		double Bij = getAngleBetweenUnitVectors(ei, ej);
		double Bki = getAngleBetweenUnitVectors(ek, ei);

		Vec3 eiej = ei.cross(ej);
		Vec3 ejek = ej.cross(ek);
		Vec3 ekei = ek.cross(ei);

		Vec3 nij = eiej.normalized();
		Vec3 njk = ejek.normalized();
		Vec3 nki = ekei.normalized();

		double ui = (Bjk + (Bij * (nij.dot(njk))) + (Bki * (nki.dot(njk)))) / (2.f * ei.dot(njk));

		sumU += ui;

		it = phi<2, 1>(cage, it);
	} while (it != vertex);

	return (1.0f / r) * sumU;
}

const double GCTriInt(const Vec3& p, const Vec3& v1, const Vec3& v2, const Vec3& nu)
{
	const Vec3 v2_v1 = v2 - v1;
	const Vec3 p_v1 = p - v1;
	const Vec3 v1_p = v1 - p;
	const Vec3 v2_p = v2 - p;

	const Vec3 p_nu = p - nu;

	const double alpha = std::acos((v2_v1.dot(p_v1)) / (v2_v1.norm() * p_v1.norm()));
	const double beta = std::acos((v1_p.dot(v2_p)) / (v1_p.norm() * v2_p.norm()));

	const double lambda = p_v1.squaredNorm() * std::sin(alpha) * std::sin(alpha);

	const double c = p_nu.squaredNorm();

	const double sqrt_c = sqrt(c);
	const double sqrt_lambda = sqrt(lambda);

	const std::array<double, 2> theta = {M_PI - alpha, M_PI - alpha - beta};
	std::array<double, 2> I;

	for (size_t i = 0; i < 2; ++i)
	{
		const double S = std::sin(theta[i]);
		const double C = std::cos(theta[i]);

		const double SS = S * S;

		const auto sign = S < 0 ? -1.0 : 1.0;

		const auto tan_part = 2 * sqrt_c * std::atan2(sqrt_c * C, sqrt(lambda + SS * c));

		const auto log_part =
			log((2 * sqrt_lambda * SS / (std::pow(1.0 - C, 2)) *
				 (1.0 - ((2 * c * C) / ((c * (1 + C) + lambda + sqrt((lambda * lambda) + (lambda * c * SS))))))));

		I[i] = (-sign * 0.5) * (tan_part + (sqrt_lambda * log_part));
	}

	return (-1.0 / (4.0 * M_PI)) * abs(I[0] - I[1] - sqrt_c * beta);
}

// https://github.com/blaisebundle/green_cage_deformer/blob/master/src/green_cage_deformer.cc
namespace
{
const Vec3 NULL_VECTOR(0.0, 0.0, 0);
constexpr auto TOLERANCE = 1e-6;
} // namespace

namespace MathConstants
{
constexpr auto PI = 3.14159265358979323846;
constexpr auto ONE_OVER_FOUR_PI = 0.07957747154594767;
constexpr auto SQRT8 = 2.828427124746190097603;
}; // namespace MathConstants

const double GCTriInt2(const Vec3& p, const Vec3& v1, const Vec3& v2)
{
	const auto v2_v1 = v2 - v1;
	const auto p_v1 = p - v1;
	const auto v1_p = v1 - p;
	const auto v2_p = v2 - p;

	const auto alpha =
		std::acos(std::min(std::max(((v2_v1).dot((p_v1))) / ((v2_v1).norm() * (p_v1).norm()), -1.0), 1.0));
	if (abs(alpha - MathConstants::PI) < TOLERANCE || abs(alpha) < TOLERANCE)
	{
		return 0.0;
	}

	const auto beta =
		std::acos(std::min(std::max((((v1_p).dot((v2_p)))) / ((v1_p).norm() * (v2_p).norm()), -1.0), 1.0));
	const auto lambda = p_v1.squaredNorm() * std::sin(alpha) * std::sin(alpha);
	const auto c = p.squaredNorm();

	const auto sqrt_c = sqrt(c);
	const auto sqrt_lambda = sqrt(lambda);

	const std::array<double, 2> theta{MathConstants::PI - alpha, MathConstants::PI - alpha - beta};
	std::array<double, 2> I;

	for (size_t i = 0; i < 2; ++i)
	{
		const auto S = std::sin(theta[i]);
		const auto C = std::cos(theta[i]);
		const auto sign = S < 0 ? -1.0 : 1.0;

		const auto SS = S * S;
		const auto half_sign = (-sign * 0.5);
		const auto tan_part = (2 * sqrt_c * std::atan2((sqrt_c * C), sqrt(lambda + (SS * c))));
		const auto log_part =
			log(((2 * sqrt_lambda * SS) / std::pow(1.0 - C, 2)) *
				(1.0 - ((2 * c * C) / ((c * (1 + C) + lambda + sqrt((lambda * lambda) + (lambda * c * SS)))))));

		I[i] = half_sign * (tan_part + (sqrt_lambda * log_part));
	}

	return -MathConstants::ONE_OVER_FOUR_PI * abs(I[0] - I[1] - sqrt_c * beta);
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

Scalar vertex_gradient_divergence(const CMap2& m, CMap2::Vertex v, const CMap2::Attribute<Vec3>* face_gradient,
								  const CMap2::Attribute<Vec3>* vertex_position)
{
	Scalar div = 0.0;
	std::vector<CMap2::Edge> edges = incident_edges(m, v);
	for (uint32 i = 0; i < edges.size(); ++i)
	{
		CMap2::Edge e1 = edges[i];
		CMap2::Edge e2 = edges[(i + 1) % edges.size()];

		CMap2::Face f(e1.dart);

		const Vec3& X = value<Vec3>(m, face_gradient, f);

		const Vec3& p0 = value<Vec3>(m, vertex_position, v);
		const Vec3& p1 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m, e1.dart)));
		const Vec3& p2 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m, e2.dart)));

		Vec3 vecR = p0 - p2;
		Vec3 vecL = p1 - p2;
		Scalar cotValue1 = vecR.dot(vecL) / vecR.cross(vecL).norm();

		vecR = p2 - p1;
		vecL = p0 - p1;
		Scalar cotValue2 = vecR.dot(vecL) / vecR.cross(vecL).norm();

		div += cotValue1 * (p1 - p0).dot(X) + cotValue2 * (p2 - p0).dot(X);
	}
	return div / 2.0;
}

double distance_vec3(const Vec3& p1, const Vec3& p2)
{
    // compute Euclidean distance or whatever
	return (p1 - p2).norm(); 
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