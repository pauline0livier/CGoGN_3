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

#include <cgogn/geometry/algos/normal.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

using Graph = cgogn::IncidenceGraph;

/**
 * set CMap2 m as a box of size bb_min and bb_max
 */ 
std::vector<CMap2::Vertex> create_bounding_box(CMap2& m, 
                                        CMap2::Attribute<Vec3>* vertex_position, 
                                        const Vec3& bb_min, const Vec3& bb_max); 

/**
 * update CMap2 m representing a box with new size bb_min bb_max
*/
void update_bounding_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, 
                                    const std::vector<CMap2::Vertex>& vertices, 
                                    const Vec3& bb_min, const Vec3& bb_max); 

/**
 * create box with specified center, normal and dimensions 
*/
void create_cage_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, 
                                    const Vec3& bb_min, const Vec3& bb_max, 
                            const Eigen::Matrix3d& local_frame); 
/**
 * set attribute vertex index for CMap2
*/
void set_attribute_vertex_index(CMap2& cage, 
                                    CMap2::Attribute<uint32>* position_indices);

/**
 * set attribute vertex index for IncidenceGraph
*/
void set_attribute_vertex_index_graph(IncidenceGraph& graph, 
                            IncidenceGraph::Attribute<uint32>* position_indices); 

/**
 * set attribute vertex shared for CMap2
*/
void set_attribute_vertex_shared(CMap2& m, CMap2::Attribute<bool>* vertex_shared); 

/**
 * set attribute face index for CMap2
*/
void set_attribute_face_index(CMap2& cage, 
                                        CMap2::Attribute<uint32>* face_indices); 

/**
 * set attribute face normal for CMap2
*/
void set_attribute_face_normal(CMap2& cage, 
                                        CMap2::Attribute<Vec3>* vertex_position, 
                                        CMap2::Attribute<Vec3>* face_normal);
/**
 * set default graph g as handle with provided properties 
*/
Graph::Vertex create_handle(Graph& g, 
		Graph::Attribute<Vec3>* g_vertex_position, 
		Graph::Attribute<Scalar>* g_vertex_radius, 
		const Vec3& handle_position, const Scalar radius_value); 

/**
 * set default graph g as axis with provided properties
*/
std::vector<Graph::Vertex> create_axis(Graph& g, 
                Graph::Attribute<Vec3>* g_vertex_position, 
                Graph::Attribute<Scalar>* g_vertex_radius, 
                const std::vector<Vec3>& vertices_positions, 
                const Scalar& radius_value); 

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_CREATION_SPACE_TOOL_H_