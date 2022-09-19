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


#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/core/functions/mesh_info.h>


namespace cgogn
{

namespace modeling
{


void create_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max)
{
	CMap2::Volume v = add_prism(m, 4);
	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi_1(m, f2))};

	Vec3 bb_min_ = bb_min;
	Vec3 bb_max_ = bb_max;

	value<Vec3>(m, vertex_position, vertices[0]) = bb_min_;
	value<Vec3>(m, vertex_position, vertices[1]) = {bb_min_[0], bb_max_[1], bb_min_[2]};
	value<Vec3>(m, vertex_position, vertices[2]) = {bb_max_[0], bb_max_[1], bb_min_[2]};
	value<Vec3>(m, vertex_position, vertices[3]) = {bb_max_[0], bb_min_[1], bb_min_[2]};

	value<Vec3>(m, vertex_position, vertices[4]) = {bb_min_[0], bb_max_[1], bb_max_[2]};
	value<Vec3>(m, vertex_position, vertices[5]) = {bb_min_[0], bb_min_[1], bb_max_[2]};
	value<Vec3>(m, vertex_position, vertices[6]) = {bb_max_[0], bb_min_[1], bb_max_[2]};
	value<Vec3>(m, vertex_position, vertices[7]) = {bb_max_[0], bb_max_[1], bb_max_[2]};
}

void create_customed_polyhedron(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max, const Vec3& handle_position1, const Vec3& handle_position2){

	CMap2::Volume v = add_custom_polyhedron(m, 6, true, 4);

	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi<1, 1, 1>(m, f1)), CMap2::Vertex(phi<1, 1, 1,1>(m, f1)),CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi<1, 1, 1>(m, f2)), CMap2::Vertex(phi<1, 1, 1, 1>(m, f2)),CMap2::Vertex(phi_1(m, f2))};

	Vec3 bb_min_ = bb_min;
	Vec3 bb_max_ = bb_max;

	double radius = (bb_min[0] - handle_position1[0])/std::cos(4*M_PI/3); 

	value<Vec3>(m, vertex_position, vertices[0]) = bb_min_;
	value<Vec3>(m, vertex_position, vertices[1]) = {bb_min_[0],handle_position1[1]+radius*std::cos(M_PI), handle_position1[2]+radius*std::sin(M_PI)};
	value<Vec3>(m, vertex_position, vertices[2]) = {bb_min_[0], handle_position1[1]+radius*std::cos(2*M_PI/3), handle_position1[2]+radius*std::sin(2*M_PI/3)};
	value<Vec3>(m, vertex_position, vertices[3]) = {bb_min_[0],handle_position1[1]+radius*std::cos(M_PI/3), handle_position1[2]+radius*std::sin(M_PI/3)};
	value<Vec3>(m, vertex_position, vertices[4]) = {bb_min_[0], handle_position1[1]+radius*std::cos(0), handle_position1[2]+radius*std::sin(0)};
	value<Vec3>(m, vertex_position, vertices[5]) = {bb_min_[0],handle_position1[1]+radius*std::cos(5*M_PI/3), handle_position1[2]+radius*std::sin(5*M_PI/3)};

	value<Vec3>(m, vertex_position, vertices[6]) = {bb_max_[0], handle_position2[1]+radius*std::cos(M_PI), handle_position2[2]+radius*std::sin(M_PI)};
	value<Vec3>(m, vertex_position, vertices[7]) = {bb_max_[0],handle_position2[1]+radius*std::cos(4*M_PI/3), handle_position2[2]+radius*std::sin(4*M_PI/3)};
	value<Vec3>(m, vertex_position, vertices[8]) = {bb_max_[0],handle_position2[1]+radius*std::cos(5*M_PI/3), handle_position2[2]+radius*std::sin(5*M_PI/3)};
	value<Vec3>(m, vertex_position, vertices[9]) = {bb_max_[0],handle_position2[1]+radius*std::cos(0), handle_position2[2]+radius*std::sin(0)};
	value<Vec3>(m, vertex_position, vertices[10]) = bb_max; 
	value<Vec3>(m, vertex_position, vertices[11]) = {bb_max_[0],handle_position2[1]+radius*std::cos(2*M_PI/3), handle_position2[2]+radius*std::sin(2*M_PI/3)};
	
}

void set_attribute_position_indices(CMap2& cage, CMap2::Attribute<uint32>* position_indices)
{
	uint32 nb_vertices = 0;

	foreach_cell(cage, [&](CMap2::Vertex v) -> bool {
		value<uint32>(cage, position_indices, v) = nb_vertices++;
		return true;
	});
}

void set_attribute_marked_vertices(CMap2& cage, CMap2::Attribute<bool>* marked_vertices)
{

	parallel_foreach_cell(cage, [&](CMap2::Vertex v) -> bool {
		value<bool>(cage, marked_vertices, v) = false;
		return true;
	});
}

void set_attribute_face_indices(CMap2& cage, CMap2::Attribute<uint32>* face_indices)
{
	uint32 nb_faces = 0;
	foreach_cell(cage, [&](CMap2::Face f) -> bool {
		value<uint32>(cage, face_indices, f) = nb_faces++;
		return true;
	});
}


void create_handle(Graph& g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius, const Vec3& center1, const Vec3& center2){
    
    Graph::Vertex nv = add_vertex(g);
		
	value<Vec3>(g, vertex_position, nv) = center1;
	value<Scalar>(g, vertex_radius, nv) = Scalar(50);

	Graph::Vertex nv1 = add_vertex(g);
		
	value<Vec3>(g, vertex_position, nv1) = center2;
	value<Scalar>(g, vertex_radius, nv1) = Scalar(50);

	connect_vertices(g, nv, nv1);
}

void create_axis(Graph& g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius, const std::vector<Vec3>& vertices_positions){
    Graph::Vertex nv = add_vertex(g);
		
	value<Vec3>(g, vertex_position, nv) = vertices_positions[0];
	value<Scalar>(g, vertex_radius, nv) = Scalar(50);

	Graph::Vertex lastVertex = nv; 

	for (unsigned i  = 1; i < vertices_positions.size(); i++){
		Graph::Vertex nv1 = add_vertex(g);
		
		value<Vec3>(g, vertex_position, nv1) = vertices_positions[i];
		value<Scalar>(g, vertex_radius, nv1) = Scalar(50);

		connect_vertices(g, lastVertex, nv1);

		lastVertex = nv1;
	}	
}

void set_graph_attribute_position_indices(Graph& g, Graph::Attribute<uint32>* position_indices)
{
	uint32 nb_vertices = 0;

	foreach_cell(g, [&](Graph::Vertex v) -> bool {
		value<uint32>(g, position_indices, v) = nb_vertices++;
		return true;
	});
}

}

}