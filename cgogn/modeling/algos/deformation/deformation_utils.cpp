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

		// double Bjk = cgogn::geometry::angle((vj - surface_point), (vk - surface_point));
		// double Bij = cgogn::geometry::angle((vi - surface_point), (vj - surface_point));
		// double Bki = cgogn::geometry::angle((vk - surface_point), (vi - surface_point));

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

} // namespace modeling

} // namespace cgogn
