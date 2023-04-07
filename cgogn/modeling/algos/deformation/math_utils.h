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

/**
 * angle between two unit vectors
*/
double getAngleBetweenUnitVectors(const Vec3& a, const Vec3& b); 

/**
 * distance between two Vec3
*/
double distance_vec3(const Vec3& p1, const Vec3& p2); 


/**
 * find main direction from user defined set 
*/
std::tuple<Vec3, Vec3, Vec3> find_main_directions_from_set(const CMap2& m, 
                            const CMap2::Attribute<Vec3>* attribute, 
                        cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set,
                        const Vec3& center);

/**
 * sort eigen vectors by minimal eigen values
*/
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sort_eigen_vectors(const Eigen::Matrix<double, 1, 
            Eigen::Dynamic>& eigen_values, 
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& eigen_vectors);

/**
 * scale the bounding box by extension factor
*/
std::tuple<Vec3,Vec3,Vec3> get_extended_bounding_box(const Vec3& bb_min, 
                            const Vec3& bb_max, const float& extension_factor); 

/**
 * get mean value of attribute from set 
*/
Vec3 get_mean_value_attribute_from_set(const CMap2& m, 
                            const CMap2::Attribute<Vec3>* attribute, 
                        cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set); 

/**
 * get mean value of Vec3 array
*/
Vec3 get_mean_value_in_array_Vec3(const std::vector<Vec3> positions); 

/**
 * get extrema values in set 
*/
std::pair<Vec3,Vec3> get_border_values_in_set(const CMap2& m, 
                    const CMap2::Attribute<Vec3>* attribute, 
                    cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set); 

/**
 * get extrema values in provided set, in local frame  
*/
std::pair<Vec3,Vec3> get_border_values_in_set_along_local_frame(const CMap2& m, 
                    const CMap2::Attribute<Vec3>* attribute, 
                    cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set, 
                    const std::tuple<Vec3, Vec3, Vec3>& main_directions); 

/**
 * get extrema values in array of Vec3
*/
std::pair<Vec3, Vec3> get_border_values_in_array_Vec3(
                                            const std::vector<Vec3> positions); 

/**
 * get vertex closest from target position
*/
CMap2::Vertex closest_vertex_in_set_from_value(const CMap2& m, 
                            const CMap2::Attribute<Vec3>* vertex_position, 
                        cgogn::ui::CellsSet<CMap2, CMap2::Vertex>* control_set, 
                                const Vec3& target_position); 

/**
 * return distance of A with projection of P on [AB] 
*/
float projection_on_segment(const Vec3& A, const Vec3& B, const Vec3& P); 


/**
 * projection of point on direction 
*/
double get_projection_on_direction(const Vec3& point, const Vec3& direction); 

/**
 * check if target projection is inside the delimited area
*/
bool check_projection_in_area(const double& projection_value, const double& min_border, const double& max_border); 

/**
 * check if the target point is inside the bounding box 
*/
bool check_triple_projection_in_area(const Vec3& point, const Vec3& min_border, const Vec3& max_border); 

/**
 * get index of virtual cube that contains the target point  
*/
size_t get_index_virtual_cube(const Vec3& projection_values, const std::vector<bool>& valid_values, const Vec3& max_area_values); 

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_MATH_UTILS_H_
