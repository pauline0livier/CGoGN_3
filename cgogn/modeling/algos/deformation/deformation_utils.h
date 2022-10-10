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


namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3; 
using Scalar = geometry::Scalar;


double getAngleBetweenUnitVectors(const Vec3& a, const Vec3& b); 

float compute_mvc(const Vec3& surface_point, Dart vertex, CMap2& cage, const Vec3& cage_point,
					  CMap2::Attribute<Vec3>* cage_position); 

const double GCTriInt(const Vec3& p, const Vec3& v1, const Vec3& v2, const Vec3& nu); 

const double GCTriInt2(const Vec3& p, const Vec3& v1, const Vec3& v2); 
     
Eigen::Vector3f sort_eigen_vectors(const Eigen::Matrix<float, 1, Eigen::Dynamic>& eigen_values, const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& eigen_vectors); 

Scalar vertex_gradient_divergence(const CMap2& m, CMap2::Vertex v, const CMap2::Attribute<Vec3>* face_gradient, const CMap2::Attribute<Vec3>* vertex_position); 

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_DEFORMATION_CAGE_H_
