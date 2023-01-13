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

#ifndef CGOGN_MODELING_ALGOS_MATH_UTILS_H_
#define CGOGN_MODELING_ALGOS_MATH_UTILS_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/core/types/cells_set.h>


namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3; 
using Scalar = geometry::Scalar;

double getAngleBetweenUnitVectors(const Vec3& a, const Vec3& b); 

double distance_vec3(const Vec3& p1, const Vec3& p2); 

Eigen::Vector3f sort_eigen_vectors(const Eigen::Matrix<float, 1, Eigen::Dynamic>& eigen_values, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& eigen_vectors);

/*
Input: bounding box data & extension factor
Output: std::tuple containing the bounding box extended of factor extension_factor; res(0) for e_bb_min, res(1) for e_bb_max and res(2) for center
*/ 
std::tuple<Vec3,Vec3,Vec3> get_extended_bounding_box(const Vec3& bb_min, const Vec3& bb_max, const float& extension_factor); 

Vec3 get_mean_value_attribute_from_set(const CMap2& m, const CMap2::Attribute<Vec3>* attribute, cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set); 

std::pair<Vec3,Vec3> get_border_values_in_set(const CMap2& m, const CMap2::Attribute<Vec3>* attribute, cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set); 


CMap2::Vertex closest_vertex_in_set_from_value(const CMap2& m, const CMap2::Attribute<Vec3>* vertex_position, cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set, const Vec3& target_position); 


} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_MATH_UTILS_H_
