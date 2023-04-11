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

const double GCTriInt(const Vec3& p, const Vec3& v1, 
						const Vec3& v2, const Vec3& nu)
{
	const Vec3 v2_v1 = v2 - v1;
	const Vec3 p_v1 = p - v1;
	const Vec3 v1_p = v1 - p;
	const Vec3 v2_p = v2 - p;

	const Vec3 p_nu = p - nu;

	const double alpha = 
		std::acos((v2_v1.dot(p_v1)) / (v2_v1.norm() * p_v1.norm()));
	const double beta = 
		std::acos((v1_p.dot(v2_p)) / (v1_p.norm() * v2_p.norm()));

	const double lambda = 
		p_v1.squaredNorm() * std::sin(alpha) * std::sin(alpha);

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

		const auto tan_part = 
			2 * sqrt_c * std::atan2(sqrt_c * C, sqrt(lambda + SS * c));

		const auto log_part =
			log((2 * sqrt_lambda * SS / (std::pow(1.0 - C, 2)) *
			(1.0 - ((2 * c * C) / ((c * (1 + C) + lambda + 
				sqrt((lambda * lambda) + (lambda * c * SS))))))));

		I[i] = (-sign * 0.5) * (tan_part + (sqrt_lambda * log_part));
	}

	return (-1.0 / (4.0 * M_PI)) * abs(I[0] - I[1] - sqrt_c * beta);
}


const double GCTriInt2(const Vec3& p, const Vec3& v1, const Vec3& v2)
{
	const auto v2_v1 = v2 - v1;
	const auto p_v1 = p - v1;
	const auto v1_p = v1 - p;
	const auto v2_p = v2 - p;

	const auto alpha =
		std::acos(std::min(std::max(
				((v2_v1).dot((p_v1))) / ((v2_v1).norm() * (p_v1).norm()), -1.0),
					 1.0));
	if (abs(alpha - M_PI) < TOLERANCE || abs(alpha) < TOLERANCE)
	{
		return 0.0;
	}

	const auto beta =
		std::acos(std::min(std::max(
			(((v1_p).dot((v2_p)))) / ((v1_p).norm() * (v2_p).norm()), -1.0), 1.0));
	const auto lambda = p_v1.squaredNorm() * std::sin(alpha) * std::sin(alpha);
	const auto c = p.squaredNorm();

	const auto sqrt_c = sqrt(c);
	const auto sqrt_lambda = sqrt(lambda);

	const std::array<double, 2> theta{M_PI - alpha, M_PI - alpha - beta};
	std::array<double, 2> I;

	for (size_t i = 0; i < 2; ++i)
	{
		const auto S = std::sin(theta[i]);
		const auto C = std::cos(theta[i]);
		const auto sign = S < 0 ? -1.0 : 1.0;

		const auto SS = S * S;
		const auto half_sign = (-sign * 0.5);
		const auto tan_part = 
			(2 * sqrt_c * std::atan2((sqrt_c * C), sqrt(lambda + (SS * c))));
		const auto log_part =
			log(((2 * sqrt_lambda * SS) / std::pow(1.0 - C, 2)) *
				(1.0 - ((2 * c * C) / ((c * (1 + C) + lambda + 
				sqrt((lambda * lambda) + (lambda * c * SS)))))));

		I[i] = half_sign * (tan_part + (sqrt_lambda * log_part));
	}

	return (-1.0 / (4.0 * M_PI))* abs(I[0] - I[1] - sqrt_c * beta);
}


Scalar vertex_gradient_divergence(const CMap2& m, CMap2::Vertex v, 
			const CMap2::Attribute<Vec3>* face_gradient,
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
		const Vec3& p1 = 
			value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m, e1.dart)));
		const Vec3& p2 = 
			value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m, e2.dart)));

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

/**
 * 
 * @param {Vec3} A
 * @param {Vec3} B
 * @param {Vec3} C
 * @param {Vec3} object_point
*/

/// @brief compute axis-related weights of object_point
/// between 2 bones (A-B and B-C)
/// compute projection of object_point on AB to see 
/// if it is projecting on the segment or not
/// if not the case check for BC 
/// if target is inside a bone (distance of projection to an extremity > 0.2)
/// full deformation of this bone 
/// if on B half contribution of each bone
/// otherwise introduce virtual joint before A and after C
/// to compensate for the local property of the axis 
/// if on virtual bone : part of the deformation = identity 
/// @param A first joint of the partial skeleton
/// @param B middle joint of the partial skeleton
/// @param C end joint of the partial skeleton
/// @param object_point 
/// @return pair of weights for the two bones 
/// 		and boolean to specify if on a virtual bone 
std::pair<Eigen::Vector2d, std::vector<bool>> weight_two_bones(const Vec3& A, 
						const Vec3& B, const Vec3& C, const Vec3& object_point) {

	Eigen::Vector2d weights;  

	std::vector<bool> fixed_point(2); 
	fixed_point[0] = false; 
	fixed_point[1] = false; 
	
	const double resAB = projection_on_segment(A, B, object_point); 
	if (resAB > 0.0 && resAB < 1.0){ 
		if (resAB >= 0.8 || resAB <= 0.2){
			const double delta = std::min(resAB, 1.0 - resAB); 
			const double local_value = (delta*0.5)/0.2; 
			weights[0] = 1.0 - local_value; 
			weights[1] = local_value; 
		} else {
			weights[0] = 1.0; 
			weights[1] = 0.0;  
		}
		
	} else if (resAB <= 0.0){
		const double delta = std::abs(resAB); 
		if (delta > 0.2){
			weights[0] = 0.5; 
			weights[1] = 0.5; 
			fixed_point[0] = true; 
			fixed_point[1] = true; 
		} else {
			const double local_value = (delta*0.5)/0.2; 
			weights[0] = local_value;
			weights[1] = 1 - local_value;
			fixed_point[1] = true; 
		}
	} else if (resAB >= 1.0){
		if (resAB == 1.0){
			weights[0] = 0.5; 
			weights[1] = 0.5; 
		} else {
			const double resBC = projection_on_segment(B, C, object_point); 
			if (resBC > 0.0 && resBC < 1.0){
				if (resBC >= 0.8 || resBC <= 0.2){
					const double delta = std::min(resBC, 1.0 - resBC);
					const double local_value = (delta*0.5)/0.2; 
					weights[0] = local_value; 
					weights[1] = 1.0 - local_value; 
				} else {
					weights[0] = 0.0; 
					weights[1] = 1.0; 
				}
			} else if (resBC >= 1.0){
				weights[0] = 0.0; 
				weights[1] = 1.0;
			} else {
				weights[0] = 0.5; 
				weights[1] = 0.5; 
			}
		}
	} 
	return std::make_pair(weights, fixed_point); 
}

/// @brief compute weights of target point on partial skeleton 
/// Constraints: linear partial skeleton
/// introduction of two virtual bones on the extremities
/// to preserve the continuity 
/// @param axis_positions vector of positions of the joints
/// each successive pair forms a bone, (n-1) bones for a vector of size n
/// @param target_point 
/// @return vector of weights of size n+1 to encode the bones + the virtual ones 
/// 		at the extremities 
Eigen::SparseVector<double> weight_partial_skeleton(const std::vector<Vec3> axis_positions, const Vec3& target_point)
{

	const double threshold = 0.3;
	const double opposite_threshold = 1.0 - threshold; 

	std::size_t number_of_bones = axis_positions.size() - 1; 

	Eigen::SparseVector<double> weights(number_of_bones + 2); 

	auto compute_local_weights = [&](const double& local_distance, 
									const size_t& index) {
		if (local_distance <= opposite_threshold || local_distance >= threshold){
			weights.coeffRef(index) = 1.0; 
		} else if (local_distance >= opposite_threshold){
			const double delta = local_distance - opposite_threshold; 
			const double local_value = (delta*0.5)/threshold; 
			weights.coeffRef(index) = 1.0 - local_value; 
			weights.coeffRef(index+1) = local_value;

		} else if (local_distance <= threshold)
		{
			const double delta = threshold - local_distance; 
			const double local_value = (delta*0.5)/threshold; 
			weights.coeffRef(index) = 1.0 - local_value; 
			weights.coeffRef(index-1) = local_value;
		}
	};

	const double distance_first_joint = projection_on_segment(axis_positions[0], axis_positions[1], target_point);

	if (distance_first_joint <= 0.0){
		const double local_distance = 1.0 - std::abs(distance_first_joint); 
		if (local_distance <= opposite_threshold){
			weights.coeffRef(0) = 1.0;  
		} else {
			const double delta = local_distance - opposite_threshold; 
			const double local_value = (delta*0.5)/threshold; 
			weights.coeffRef(1) = local_value;
			weights.coeffRef(0) = 1 - local_value; 
		}
	} else {
		size_t current_index = 1; 
		if (distance_first_joint > 0.0 && distance_first_joint <= 1.0){ 
			compute_local_weights(distance_first_joint, current_index); 
		
		} else {
			current_index ++; 
			double previous_distance_joint = distance_first_joint; 
			while (true){
				const double current_distance_joint = projection_on_segment(
												axis_positions[current_index -1], 
												axis_positions[current_index], 
												target_point); 

				if (current_distance_joint > 0.0 && current_distance_joint <= 1.0){
					compute_local_weights(current_distance_joint, current_index);
					break; 
				} else if (previous_distance_joint > 1.0 && 
												current_distance_joint <= 0.0){
						weights.coeffRef(current_index) = 0.5; 
						weights.coeffRef(current_index - 1) = 0.5; 
						break; 
				} else {
					if (current_index < number_of_bones){
						previous_distance_joint = current_distance_joint; 
						current_index ++; 
					} else {
						weights.coeffRef(current_index) = 1.0; 
						weights.coeffRef(current_index + 1) = 0.0; 
						break; 
					}
					
				}
			}
		}

	}

	return weights; 
}



} // namespace modeling

} // namespace cgogn
