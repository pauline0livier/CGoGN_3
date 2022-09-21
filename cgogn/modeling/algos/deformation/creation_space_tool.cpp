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

void create_handle_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max, const Vec3& center_min, const Vec3& center_max, Vec3& handle_normal){

	CMap2::Volume v = add_prism(m, 6);

	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi<1, 1, 1>(m, f1)), CMap2::Vertex(phi<1, 1, 1,1>(m, f1)),CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi<1, 1, 1>(m, f2)), CMap2::Vertex(phi<1, 1, 1, 1>(m, f2)),CMap2::Vertex(phi_1(m, f2))};

	Vec3 bb_min_ = bb_min;
	Vec3 bb_max_ = bb_max;

	// plane ax+by+cz+d = 0, with (a,b,c) = handle_normal and (x,y,z) = bb_min (same for bb_max)

	// bbmin belongs to the plane 
	float d_handle_min = -(handle_normal.dot(bb_min));

	// find alpha such that center + alpha*handle_normal belongs to plane
	float alpha_min = -d_handle_min - (handle_normal.dot(center_min));

	Eigen::Vector3d center_min_plane = center_min + alpha_min*handle_normal; 

	Vec3 u_min = bb_min - center_min_plane; 
	u_min.normalize(); 

	Vec3 v_min = handle_normal.cross(u_min); 
	v_min.normalize();  

	Eigen::Matrix3d frame_min;
	frame_min.row(0) = u_min;
	frame_min.row(1) = v_min; 
	frame_min.row(2) = handle_normal; 
	Eigen::Matrix3d frame_min_inv = frame_min.inverse();

	Eigen::Vector3d local_bbmin = frame_min * (bb_min - center_min_plane);

	double radius = local_bbmin.norm(); 

	value<Vec3>(m, vertex_position, vertices[0]) = bb_min;

	Eigen::Vector3d local_vertex1 = {radius*std::cos(-M_PI/3), radius*std::sin(-M_PI/3), local_bbmin[2]}; 
	
	Eigen::Vector3d local_vertex2 = {radius*std::cos(-2*M_PI/3), radius*std::sin(-2*M_PI/3), local_bbmin[2]}; 
	
	Eigen::Vector3d local_vertex3 = {radius*std::cos(-M_PI), radius*std::sin(-M_PI), local_bbmin[2]}; 
	
	Eigen::Vector3d local_vertex4 = {radius*std::cos(-4*M_PI/3), radius*std::sin(-4*M_PI/3), local_bbmin[2]}; 
	
	Eigen::Vector3d local_vertex5 = {radius*std::cos(-5*M_PI/3), radius*std::sin(-5*M_PI/3), local_bbmin[2]}; 
	
	value<Vec3>(m, vertex_position, vertices[1]) = (frame_min_inv*local_vertex1) + center_min_plane;

	value<Vec3>(m, vertex_position, vertices[2]) = frame_min_inv*local_vertex2 + center_min_plane;

	value<Vec3>(m, vertex_position, vertices[3]) = frame_min_inv*local_vertex3 + center_min_plane;

	value<Vec3>(m, vertex_position, vertices[4]) = frame_min_inv*local_vertex4+center_min_plane;

	value<Vec3>(m, vertex_position, vertices[5]) = frame_min_inv*local_vertex5 + center_min_plane;



	// Other side 
	// bbmax belongs to the plane 
	float d_handle_max = -(handle_normal.dot(bb_max));

	// find alpha such that center + alpha*handle_normal belongs to plane
	float alpha_max = -d_handle_max - (handle_normal.dot(center_max));

	Eigen::Vector3d center_max_plane = center_max + alpha_max*handle_normal; 

	Vec3 u_max = bb_max - center_max_plane; 
	u_max.normalize(); 

	Vec3 v_max = handle_normal.cross(u_max); 
	v_max.normalize();  

	Eigen::Matrix3d frame_max;
	frame_max.row(0) = u_max;
	frame_max.row(1) = v_max; 
	frame_max.row(2) = handle_normal; 
	Eigen::Matrix3d frame_max_inv = frame_max.inverse();

	Eigen::Vector3d local_bbmax = frame_max * (bb_max - center_max_plane);

	Eigen::Vector3d local_vertex6 = {radius*std::cos(2*M_PI/3), radius*std::sin(2*M_PI/3), local_bbmax[2]}; 
	

	Eigen::Vector3d local_vertex7 = {radius*std::cos(M_PI), radius*std::sin(M_PI), local_bbmax[2]}; 
	

	Eigen::Vector3d local_vertex8 = {radius*std::cos(4*M_PI/3), radius*std::sin(4*M_PI/3), local_bbmax[2]}; 
	

	Eigen::Vector3d local_vertex9 = {radius*std::cos(5*M_PI/3), radius*std::sin(5*M_PI/3), local_bbmax[2]}; 
	
	
	Eigen::Vector3d local_vertex11 = {radius*std::cos(M_PI/3), radius*std::sin(M_PI/3), local_bbmax[2]};

	value<Vec3>(m, vertex_position, vertices[6]) = frame_max_inv*local_vertex6 + center_max_plane;

	value<Vec3>(m, vertex_position, vertices[7]) = frame_max_inv*local_vertex7 + center_max_plane;

	value<Vec3>(m, vertex_position, vertices[8]) = frame_max_inv*local_vertex8 + center_max_plane;

	value<Vec3>(m, vertex_position, vertices[9]) = frame_max_inv*local_vertex9 + center_max_plane;
	
	value<Vec3>(m, vertex_position, vertices[10]) = bb_max;

	value<Vec3>(m, vertex_position, vertices[11]) = frame_max_inv*local_vertex11 + center_max_plane;
	
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