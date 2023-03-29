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
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>

namespace cgogn
{

namespace modeling
{

/**
 * Set m as a volume of type generalized cube 
 * Assign the positions of the vertices of this volume 
 * from the provided bounding box values
 * @param {CMap2} m 
 * @param {CMap2::Attribute<Vec3>} vertex_position
 * @param {Vec3} bb_min 
 * @param {Vec3} bb_max
*/
std::vector<CMap2::Vertex> create_bounding_box(CMap2& m, 
	CMap2::Attribute<Vec3>* vertex_position, 
	const Vec3& bb_min, const Vec3& bb_max)
{
	CMap2::Volume v = add_prism(m, 4);
	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi_1(m, f2))};

	value<Vec3>(m, vertex_position, vertices[0]) = bb_min;
	value<Vec3>(m, vertex_position, vertices[1]) = 
											{bb_min[0], bb_max[1], bb_min[2]};
	value<Vec3>(m, vertex_position, vertices[2]) = 
											{bb_max[0], bb_max[1], bb_min[2]};
	value<Vec3>(m, vertex_position, vertices[3]) = 
											{bb_max[0], bb_min[1], bb_min[2]};

	value<Vec3>(m, vertex_position, vertices[4]) = 
											{bb_min[0], bb_max[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[5]) = 
											{bb_min[0], bb_min[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[6]) = 
											{bb_max[0], bb_min[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[7]) = 
											{bb_max[0], bb_max[1], bb_max[2]};

	return vertices;
}

/**
 * Update the positions of the vertices of CMap2 
 * of type generalized cube  
 * from the provided bounding box values
 * @param {CMap2} m 
 * @param {CMap2::Attribute<Vec3>} vertex_position
 * @param {Vec3} bb_min 
 * @param {Vec3} bb_max 
*/
void update_bounding_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, 
					const std::vector<CMap2::Vertex>& vertices,
					const Vec3& bb_min, const Vec3& bb_max)
{

	value<Vec3>(m, vertex_position, vertices[0]) = bb_min;
	value<Vec3>(m, vertex_position, vertices[1]) = {bb_min[0], bb_max[1], bb_min[2]};
	value<Vec3>(m, vertex_position, vertices[2]) = {bb_max[0], bb_max[1], bb_min[2]};
	value<Vec3>(m, vertex_position, vertices[3]) = {bb_max[0], bb_min[1], bb_min[2]};

	value<Vec3>(m, vertex_position, vertices[4]) = {bb_min[0], bb_max[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[5]) = {bb_min[0], bb_min[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[6]) = {bb_max[0], bb_min[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[7]) = {bb_max[0], bb_max[1], bb_max[2]};
}

/**
 * Set m as a volume of type generalized cube 
 * Compute local frame of this cube from normal and center
 * Assign the positions of the vertices of this volume 
 * from the provided bounding box values
 * @param {CMap2} m 
 * @param {CMap2::Attribute<Vec3>} vertex_position
 * @param {Vec3} bb_min 
 * @param {Vec3} bb_max
 * @param {Vec3} center
 * @param {Vec3} normal
*/
void create_cage_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, 
					const Vec3& bb_min, const Vec3& bb_max,
					const Vec3& center, const Vec3& normal)
{
	CMap2::Volume v = add_prism(m, 4);
	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), 
		CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), 
		CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi_1(m, f2))};

	float d_min = -(normal.dot(bb_min));

	float alpha_min = -d_min - (normal.dot(center));

	Eigen::Vector3d center_min_plane = center + alpha_min * normal;

	Vec3 local_frame_x = bb_min - center_min_plane;
	local_frame_x.normalize();

	Vec3 local_frame_y = normal.cross(local_frame_x);
	local_frame_y.normalize();

	Eigen::Matrix3d local_frame, frame_inverse;

	local_frame.row(0) = local_frame_x;
	local_frame.row(1) = local_frame_y;
	local_frame.row(2) = normal;
	frame_inverse = local_frame.inverse();

	const Vec3 local_bb_min = local_frame * (bb_min - center);
	const Vec3 local_bb_max = local_frame * (bb_max - center);

	const double radius = local_bb_min.norm();

	double min_n, max_n;
	if (local_bb_min[2] > local_bb_max[2])
	{
		min_n = local_bb_max[2];
		max_n = local_bb_min[2];
	}
	else
	{
		min_n = local_bb_min[2];
		max_n = local_bb_max[2];
	}

	Eigen::Vector3d local_vertex0 = 
		{radius * std::cos(0), radius * std::sin(0), max_n};
	Eigen::Vector3d local_vertex1 =	
		{radius * std::cos(M_PI / 2), radius * std::sin(M_PI / 2), max_n};
	Eigen::Vector3d local_vertex2 = 
		{radius * std::cos(M_PI), radius * std::sin(M_PI), max_n};
	Eigen::Vector3d local_vertex3 = 
		{radius * std::cos(3 * M_PI / 2), radius * std::sin(3 * M_PI / 2), max_n};

	Eigen::Vector3d local_vertex4 = 
		{radius * std::cos(-3 * M_PI / 2), radius * std::sin(-3 * M_PI / 2), min_n};
	Eigen::Vector3d local_vertex5 = 
		{radius * std::cos(0), radius * std::sin(0), min_n};
	Eigen::Vector3d local_vertex6 = 
		{radius * std::cos(-M_PI / 2), radius * std::sin(-M_PI / 2), min_n};
	Eigen::Vector3d local_vertex7 = 
		{radius * std::cos(-M_PI), radius * std::sin(-M_PI), min_n};

	value<Vec3>(m, vertex_position, vertices[0]) = 
					(frame_inverse * local_vertex0) + center;
	value<Vec3>(m, vertex_position, vertices[1]) = 
					(frame_inverse * local_vertex1) + center;
	value<Vec3>(m, vertex_position, vertices[2]) = 
					frame_inverse * local_vertex2 + center;
	value<Vec3>(m, vertex_position, vertices[3]) = 
					frame_inverse * local_vertex3 + center;
	value<Vec3>(m, vertex_position, vertices[4]) = 
					frame_inverse * local_vertex4 + center;
	value<Vec3>(m, vertex_position, vertices[5]) = 
					frame_inverse * local_vertex5 + center;
	value<Vec3>(m, vertex_position, vertices[6]) = 
					frame_inverse * local_vertex6 + center;
	value<Vec3>(m, vertex_position, vertices[7]) = 
					frame_inverse * local_vertex7 + center;
}

/**
 * set cage position_indices by looping on the vertices
 * @param {CMap2} cage
 * @param {CMap2::Attribute<uint32>*} position_indices
*/
void set_attribute_vertex_index(CMap2& cage, 
			CMap2::Attribute<uint32>* position_indices)
{
	uint32 nb_vertices = 0;

	foreach_cell(cage, [&](CMap2::Vertex v) -> bool {
		value<uint32>(cage, position_indices, v) = nb_vertices++;
		return true;
	});
}

/**
 * set graph position_indices by looping on the vertices
 * @param {IncidenceGraph} graph
 * @param {Incidence::Attribute<uint32>*} position_indices
*/
void set_attribute_vertex_index_graph(IncidenceGraph& graph, 
				IncidenceGraph::Attribute<uint32>* position_indices)
{
	uint32 nb_vertices = 0;

	foreach_cell(graph, [&](IncidenceGraph::Vertex v) -> bool {
		value<uint32>(graph, position_indices, v) = nb_vertices++;
		return true;
	});
}

/**
 * set cage position_indices by looping on the vertices
 * and marking all of them as false
 * @param {CMap2} cage
 * @param {CMap2::Attribute<uint32>*} marked_vertices
*/
void set_attribute_marked_vertices(CMap2& cage, 
	CMap2::Attribute<bool>* marked_vertices)
{

	parallel_foreach_cell(cage, [&](CMap2::Vertex v) -> bool {
		value<bool>(cage, marked_vertices, v) = false;
		return true;
	});
}

/**
 * set cage face_indices by looping on the vertices
 * @param {CMap2} cage
 * @param {CMap2::Attribute<uint32>*} face_indices
*/
void set_attribute_face_index(CMap2& cage, CMap2::Attribute<uint32>* face_indices)
{
	uint32 nb_faces = 0;
	foreach_cell(cage, [&](CMap2::Face f) -> bool {
		value<uint32>(cage, face_indices, f) = nb_faces++;
		return true;
	});
}

/**
 * set cage face_normal by looping on the faces
 * and computing the normal of one of the two triangles composing the face
 * @param {CMap2} cage
 * @param {CMap2::Attribute<Vec3>*} vertex_position
 * @param {CMap2::Attribute<Vec3>*} face_normal
*/
void set_attribute_face_normal(CMap2& cage, 
								CMap2::Attribute<Vec3>* vertex_position,
							   CMap2::Attribute<Vec3>* face_normal)
{
	foreach_cell(cage, [&](CMap2::Face fc) -> bool {
		std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(cage, fc);

		// triangle 1

		const std::vector<CMap2::Vertex> triangle1 = 
					{face_vertices_[1], face_vertices_[3], face_vertices_[0]};
		const std::vector<Vec3> t1_values = 
				{value<Vec3>(cage, vertex_position, triangle1[0]),
				value<Vec3>(cage, vertex_position, triangle1[1]),
				value<Vec3>(cage, vertex_position, triangle1[2])};

		value<Vec3>(cage, face_normal, fc) = 
			(cgogn::geometry::normal(t1_values[0], t1_values[1], t1_values[2]))
			.normalized(); 

		return true;
	} );
}

/**
 * set object fixed_position attribute
 * set value as Vec3(0,0,0) by default
 * @param {CMap2} object
 * @param {CMap2::Attribute<Vec3>*} object_fixed_position
*/
void set_attribute_fixed_position(CMap2& object, CMap2::Attribute<Vec3>* object_fixed_position)
{
	foreach_cell(object, [&](CMap2::Vertex v) -> bool {
		value<Vec3>(object, object_fixed_position, v) = {0.0, 0.0, 0.0};
		return true;
	});
}


/// @brief Create handle from default graph 
/// set vertices, their positions and radius
/// @param g default graph 
/// @param g_vertex_position default graph vertices position 
/// @param g_vertex_radius default graph vertices radius
/// @param handle_position position of handle
/// @param radius_value radius of vertex
/// @return handle vertex 
Graph::Vertex create_handle(Graph& g, 
		Graph::Attribute<Vec3>* g_vertex_position, 
		Graph::Attribute<Scalar>* g_vertex_radius, 
		const Vec3& handle_position, const Scalar radius_value)
{
	Graph::Vertex nv = add_vertex(g);

	value<Vec3>(g, g_vertex_position, nv) = handle_position;
	value<Scalar>(g, g_vertex_radius, nv) = radius_value;

	return nv;
}


/// @brief Create axis from default graph 
/// set vertices, their positions and radius
/// @param g default graph 
/// @param g_vertex_position default graph vertices position 
/// @param g_vertex_radius default graph vertices radius 
/// @param vertices_positions positions to set to the graph 
/// @param radius_value radius to set to the graph 
/// @return std::vector<Graph::Vertex> axis_vertices 
/// 	usage of std::vector to store the vertices in the order of the axis
std::vector<Graph::Vertex> create_axis(Graph& g,
	Graph::Attribute<Vec3>* g_vertex_position, 
	Graph::Attribute<Scalar>* g_vertex_radius,
	const std::vector<Vec3>& vertices_positions, const Scalar& radius_value)
{
	std::vector<Graph::Vertex> axis_vertices;
	Graph::Vertex nv = add_vertex(g);

	axis_vertices.push_back(nv);

	value<Vec3>(g, g_vertex_position, nv) = vertices_positions[0];
	value<Scalar>(g, g_vertex_radius, nv) = radius_value; 

	Graph::Vertex lastVertex = nv;

	for (unsigned i = 1; i < vertices_positions.size(); i++)
	{
		Graph::Vertex nv1 = add_vertex(g);
		axis_vertices.push_back(nv1);

		value<Vec3>(g, g_vertex_position, nv1) = vertices_positions[i];
		value<Scalar>(g, g_vertex_radius, nv1) = radius_value; 

		connect_vertices(g, lastVertex, nv1);

		lastVertex = nv1;
	}

	return axis_vertices;
}



} // namespace modeling

} // namespace cgogn