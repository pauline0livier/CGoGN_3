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

#ifndef CGOGN_MODELING_ALGOS_DEFORMATION_UTILS_H_
#define CGOGN_MODELING_ALGOS_DEFORMATION_UTILS_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/math_utils.h>
#include <Eigen/Sparse>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3; 
using Scalar = geometry::Scalar;

using Graph = cgogn::IncidenceGraph;

namespace
{
extern const Vec3 NULL_VECTOR(0.0, 0.0, 0);
extern constexpr auto TOLERANCE = 1e-6;
} 

/**
 * specific function for Green-based deformation of cage 
 * taken from paper [Green Coordinates, Lipman et al. 2008]
*/
const double GCTriInt(const Vec3& p, const Vec3& v1, 
						const Vec3& v2, const Vec3& nu); 

/**
 * adaptation of the function above without parameter nu 
*/
const double GCTriInt2(const Vec3& p, const Vec3& v1, const Vec3& v2); 
    
/**
 * compute gradient divergence of the vertices
 * useful for the geodesic distance
*/
Scalar vertex_gradient_divergence(const CMap2& m, CMap2::Vertex v, 
						const CMap2::Attribute<Vec3>* face_gradient, 
						const CMap2::Attribute<Vec3>* vertex_position); 


/// @brief Compute weights between 2 bones (A-B and B-C) and a point
/// @param A first joint of the partial skeleton
/// @param B middle joint of the partial skeleton
/// @param C end joint of the partial skeleton
/// @param object_point 
/// @return 
std::pair<Eigen::Vector2d, std::vector<bool>> weight_two_bones(const Vec3& A, 
						const Vec3& B, const Vec3& C, const Vec3& object_point); 


/// @brief compute weights between vector of bones and a target point
/// @param axis_positions array of Vector3
/// @param target_point point on the model to deform 
/// @return weights of the target point relative to the partial skeleton
Eigen::SparseVector<double> weight_partial_skeleton(const std::vector<Vec3> axis_positions, const Vec3& target_point);

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_DEFORMATION_UTILS_H_
