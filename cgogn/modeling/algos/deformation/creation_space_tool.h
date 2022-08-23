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
#ifndef CGOGN_MODELING_ALGOS_CREATION_SPACE_TOOL_H_
#define CGOGN_MODELING_ALGOS_CREATION_SPACE_TOOL_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>
#include <cgogn/core/functions/mesh_ops/volume.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

using Graph = cgogn::IncidenceGraph;

void create_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max); 

void set_attribute_position_indices(CMap2& cage, CMap2::Attribute<uint32>* position_indices); 

void set_attribute_marked_vertices(CMap2& cage, CMap2::Attribute<bool>* marked_vertices); 

void set_attribute_face_indices(CMap2& cage, CMap2::Attribute<uint32>* face_indices); 

void create_handle(Graph& g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius, const Vec3& center1, const Vec3& center2); 

void create_axis(Graph& g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius, const std::vector<Vec3>& vertices_positions); 

void set_graph_attribute_position_indices(Graph& g, Graph::Attribute<uint32>* position_indices); 

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_CREATION_SPACE_TOOL_H_