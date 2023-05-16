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

/// @brief Update the positions of the vertices of CMap2 
/// of type generalized cube  
/// from the provided bounding box values
/// @param m 
/// @param vertex_position 
/// @param vertices 
/// @param bb_min 
/// @param bb_max 
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

 
std::vector<Vec3> get_bounding_box_positions(const Vec3& bb_min, const Vec3& bb_max)
{
	std::vector<Vec3> positions(8); 
	positions[0] = bb_min; 
	positions[1] = {bb_min[0], bb_max[1], bb_min[2]}; 
	positions[2] = {bb_max[0], bb_max[1], bb_min[2]};
	positions[3] = {bb_max[0], bb_min[1], bb_min[2]};
	positions[4] = {bb_min[0], bb_max[1], bb_max[2]};
	positions[5] = {bb_min[0], bb_min[1], bb_max[2]};
	positions[6] = {bb_max[0], bb_min[1], bb_max[2]};
	positions[7] = {bb_max[0], bb_max[1], bb_max[2]};

	return positions; 
}

/// @brief Create hexahedron from bounding values and main directions 
/// @param m mesh to set to hexahedron
/// @param m_vertex_position positions of the vertices of the mesh to set
/// @param bb_min Vec3 minimum values in local frame
/// @param bb_max Vec3 maximum values in local frame
/// @param local_frame matrice   
std::vector<CMap2::Vertex> create_cage_box(CMap2& m, CMap2::Attribute<Vec3>* m_vertex_position, 
					const Vec3& bb_min, const Vec3& bb_max,
					const Eigen::Matrix3d& local_frame)
{
	CMap2::Volume v = add_prism(m, 4);
	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), 
		CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), 
		CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi_1(m, f2))};

	Eigen::Matrix3d frame_inverse;

	frame_inverse = local_frame.inverse();

	Vec3 center = (bb_min + bb_max) / Scalar(2); 

	const Vec3 local_bb_min = {0.0, 0.0, 0.0};
	const Vec3 local_bb_max = local_frame * (bb_max - bb_min);

	const Vec3 local_bb_min_center = local_frame * (bb_min - center); 

	const double radius = local_bb_min_center.norm();

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

	std::vector<Eigen::Vector3d> local_positions(8);
	local_positions[0] = {0.0, 0.0, 0.0}; 
	local_positions[1] = {0.0, local_bb_max[1], 0.0}; 
	local_positions[2] = {local_bb_max[0], local_bb_max[1], 0.0}; 
	local_positions[3] = {local_bb_max[0], 0.0, 0.0};
	local_positions[4] = {0.0, local_bb_max[1], local_bb_max[2]};
	local_positions[5] = {0.0, 0.0, local_bb_max[2]};
	local_positions[6] = {local_bb_max[0], 0.0, local_bb_max[2]};
	local_positions[7] = local_bb_max;

	for (size_t p = 0; p < 8; p++){
 
		value<Vec3>(m, m_vertex_position, vertices[p]) = 
					(frame_inverse * local_positions[p]) + bb_min;

	}

	return vertices; 

}

/// @brief create cylinder from provided bounding box dimensions 
/// @param m mesh to turn into cylinder
/// @param vertex_position 
/// @param bb_min 
/// @param bb_max 
void create_cylinder_tool(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, 
					const Vec3& bb_min, const Vec3& bb_max)
{
	CMap2::Volume v = add_prism(m, 6);
	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), 
		CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi<1, 1, 1>(m, f1)), 
		CMap2::Vertex(phi<1, 1, 1,1>(m, f1)),CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), 
		CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi<1, 1, 1>(m, f2)), 
		CMap2::Vertex(phi<1, 1, 1, 1>(m, f2)),CMap2::Vertex(phi_1(m, f2))};

	Vec3 bb_min_ = bb_min;
	Vec3 bb_max_ = bb_max;

	const Vec3 center = (bb_min_ + bb_max)/2.0; 
	Vec3 center_min = {center[0], center[1], bb_min[2]}; 
	Vec3 center_max = {center[0], center[1], bb_max[2]};

	double radius = (bb_min[0] - center_min[0])/std::cos(4*M_PI/3); 

	value<Vec3>(m, vertex_position, vertices[1]) = 
		{bb_min_[0], center_min[1]+radius*std::cos(M_PI), 
			center_min[2]+radius*std::sin(M_PI)};
	value<Vec3>(m, vertex_position, vertices[2]) = 
		{bb_min_[0], center_min[1]+radius*std::cos(2*M_PI/3), 
			center_min[2]+radius*std::sin(2*M_PI/3)};
	value<Vec3>(m, vertex_position, vertices[3]) = 
		{bb_min_[0], center_min[1]+radius*std::cos(M_PI/3), 
			center_min[2]+radius*std::sin(M_PI/3)};
	value<Vec3>(m, vertex_position, vertices[4]) = 
		{bb_min_[0], center_min[1]+radius*std::cos(0), 
			center_min[2]+radius*std::sin(0)};
	value<Vec3>(m, vertex_position, vertices[5]) = 
		{bb_min_[0],center_min[1]+radius*std::cos(5*M_PI/3), 
			center_min[2]+radius*std::sin(5*M_PI/3)};

	value<Vec3>(m, vertex_position, vertices[6]) = 
		{bb_max_[0], center_max[1]+radius*std::cos(M_PI), 
			center_max[2]+radius*std::sin(M_PI)};
	value<Vec3>(m, vertex_position, vertices[7]) = 
		{bb_max_[0],center_max[1]+radius*std::cos(4*M_PI/3), 
			center_max[2]+radius*std::sin(4*M_PI/3)};
	value<Vec3>(m, vertex_position, vertices[8]) = 
		{bb_max_[0],center_max[1]+radius*std::cos(5*M_PI/3), 
			center_max[2]+radius*std::sin(5*M_PI/3)};
	value<Vec3>(m, vertex_position, vertices[9]) = 
		{bb_max_[0],center_max[1]+radius*std::cos(0), 
			center_max[2]+radius*std::sin(0)};
	value<Vec3>(m, vertex_position, vertices[10]) = bb_max; 
	value<Vec3>(m, vertex_position, vertices[11]) = 
		{bb_max_[0],center_max[1]+radius*std::cos(2*M_PI/3), 
			center_max[2]+radius*std::sin(2*M_PI/3)};
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
 * set cage position_indices by running through the vertices array 
 * @param {CMap2} cage
 * @param {CMap2::Attribute<uint32>*} position_indices
 * @param {vector of CMap2::Vertex} vertices
*/
void set_attribute_vertex_index_from_vertices_array(CMap2& cage, 
									CMap2::Attribute<uint32>* position_indices, std::vector<CMap2::Vertex> vertices)
{
	const std::size_t vertices_size = vertices.size(); 
	for (std::size_t i = 0; i < vertices_size; i++)
	{
		value<uint32>(cage, position_indices, vertices[i]) = i;
	}
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

/// @brief set mesh vertex shared properties 
/// @param m mesh 
/// @param vertex_shared attribute 
void set_attribute_vertex_shared(CMap2& m, CMap2::Attribute<bool>* vertex_shared)
{
	foreach_cell(m, [&](CMap2::Vertex v) -> bool {
		value<bool>(m, vertex_shared, v) = false;
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