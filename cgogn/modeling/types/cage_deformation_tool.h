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

#ifndef CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_

#include <cgogn/modeling/types/space_deformation_tool.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>

/**
 * @Class CageDeformationTool
 * Represents the local cage deformation tool
 */
class CageDeformationTool : public SpaceDeformationTool<MESH>
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

	struct Point
	{
		Vec3 position;
		uint32_t control_cage_index;
		bool inside_control_cage;
	};

	struct Triangle
	{
		std::vector<Point> points;
		std::vector<uint32_t> virtual_cage_indices;
		Vec3 normal;
		std::pair<Point, Point> edges;
	};

	struct Virtual_cube
	{
		std::vector<Triangle> triangles; // index of point ?
		std::vector<Point> points;
		std::vector<uint32_t> virtual_cage_indices;
	};

	struct Local_direction_control_planes
	{
		double d_min;
		double d_max;
		double d_gap;

		std::pair<Triangle, Triangle> triangles_d_min;
		std::pair<Triangle, Triangle> triangles_d_max;

		Vec3 direction;
		Vec3 shift_after_d_max;
		Vec3 shift_before_d_min;
	};

public:
	MESH* control_cage_;
	std::shared_ptr<Attribute<Vec3>> control_cage_vertex_position_;
	std::shared_ptr<Attribute<Vec3>> influence_cage_vertex_position_;

	Eigen::VectorXd attenuation_;

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;

	std::vector<Triangle> cage_triangles_;

	MESH* influence_cage_;

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> control_cage_coords_;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> control_cage_normal_coords_;

	CageDeformationTool() : SpaceDeformationTool<MESH>(), m_hFactor(-1.0f), control_cage_vertex_position_(nullptr)
	{
	}

	~CageDeformationTool()
	{
	}

	void create_space_tool(MESH* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max,
						   const Vec3& center, const Vec3& normal)
	{
		control_cage_ = m;

		/*cgogn::modeling::create_cage_box(*m, vertex_position, bb_min, bb_max, center, normal);*/ // not working well
																								   // so far
		cgogn::modeling::create_bounding_box(*m, vertex_position, bb_min, bb_max);

		control_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		control_cage_bb_min_ = bb_min;
		control_cage_bb_max_ = bb_max;

		control_cage_vertex_index_ = cgogn::add_attribute<uint32, Vertex>(*control_cage_, "vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*control_cage_, control_cage_vertex_index_.get());

		init_triangles();

		init_control_cage_plane();

		init_virtual_cubes();
	}

	void set_influence_cage(MESH& object, const CMap2::Attribute<Vec3>* object_vertex_position, MESH* m,
							CMap2::Attribute<Vec3>* vertex_position)
	{
		influence_cage_ = m;

		/*std::pair<Vec3, Vec3> influence_area_borders = cgogn::modeling::get_border_values_in_set(object,
		 * object_vertex_position, this->influence_area_); */

		// first simplification 3x control cage
		std::tuple<Vec3, Vec3, Vec3> res_extended_bounding_box =
			cgogn::modeling::get_extended_bounding_box(control_cage_bb_min_, control_cage_bb_max_, 3.0);

		influence_cage_bb_min_ = std::get<0>(res_extended_bounding_box);
		influence_cage_bb_max_ = std::get<1>(res_extended_bounding_box);

		cgogn::modeling::create_bounding_box(*m, vertex_position, influence_cage_bb_min_, influence_cage_bb_max_);

		influence_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, Vertex>(*influence_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*influence_cage_, marked_vertices.get());

		// update influence area
		foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);

			bool inside_cage = check_mvc_on_point_inside_influence_cage(surface_point);

			if (inside_cage)
			{
				this->influence_area_->select(v);
			}

			return true;
		});
	}

	void bind_mvc(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		control_cage_coords_.resize(nbv_object, nbv_cage);
		control_cage_coords_.setZero();

		this->influence_area_->foreach_cell([&](Vertex v) {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);

			uint32 surface_point_index = value<uint32>(object, object_vertex_index, v);

			const double d_x = -(surface_point.dot(local_x_direction_control_planes_.direction)),

						 d_y = -(surface_point.dot(local_y_direction_control_planes_.direction)),

						 d_z = -(surface_point.dot(local_z_direction_control_planes_.direction));

			const bool valid_x_dir = (d_x <= local_x_direction_control_planes_.d_max &&
									  d_x >= local_x_direction_control_planes_.d_min),

					   valid_y_dir = (d_y <= local_y_direction_control_planes_.d_max &&
									  d_y >= local_y_direction_control_planes_.d_min),

					   valid_z_dir = (d_z <= local_z_direction_control_planes_.d_max &&
									  d_z >= local_z_direction_control_planes_.d_min);

			if (valid_x_dir && valid_y_dir && valid_z_dir)
			{
				compute_mvc_on_point_inside_cage(surface_point, surface_point_index);
			}
			else
			{
				Virtual_cube virtual_cube_target;
				if (valid_x_dir)
				{
					if (valid_y_dir)
					{
						if (d_z > local_z_direction_control_planes_.d_max)
						{
							// faceCube[5]
							// Beware, need to inverse the normal of the triangle1 and triangle2
							virtual_cube_target = face_adjacent_virtual_cube_[5];
						}
						else
						{
							// faceCube[4]
							virtual_cube_target = face_adjacent_virtual_cube_[4];
						}
					}
					else if (valid_z_dir)
					{
						if (d_y > local_y_direction_control_planes_.d_max)
						{
							// faceCube[3]
							virtual_cube_target = face_adjacent_virtual_cube_[3];
						}
						else
						{
							// faceCube[2]
							virtual_cube_target = face_adjacent_virtual_cube_[2];
						}
					}
					else
					{
						if (d_y > local_y_direction_control_planes_.d_max)
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// edgeCube[11]
								virtual_cube_target = edge_adjacent_virtual_cube_[11];
							}
							else
							{
								// edgeCube[9]
								virtual_cube_target = edge_adjacent_virtual_cube_[9];
							}
						}
						else
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// edgeCube[10]
								virtual_cube_target = edge_adjacent_virtual_cube_[10];
							}
							else
							{
								// edgeCube[8]
								virtual_cube_target = edge_adjacent_virtual_cube_[8];
							}
						}
					}
				}
				else if (valid_y_dir)
				{
					if (valid_z_dir)
					{
						if (d_x > local_x_direction_control_planes_.d_max)
						{
							// faceCube[1]
							virtual_cube_target = face_adjacent_virtual_cube_[1];
						}
						else
						{
							// faceCube[0]
							virtual_cube_target = face_adjacent_virtual_cube_[0];
						}
					}
					else
					{
						if (d_x > local_x_direction_control_planes_.d_max)
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// edgeCube[3]
								virtual_cube_target = edge_adjacent_virtual_cube_[3];
							}
							else
							{
								// edgeCube[1]
								virtual_cube_target = edge_adjacent_virtual_cube_[1];
							}
						}
						else
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// edgeCube[2]
								virtual_cube_target = edge_adjacent_virtual_cube_[2];
							}
							else
							{
								// edgeCube[0]
								virtual_cube_target = edge_adjacent_virtual_cube_[0];
							}
						}
					}
				}
				else if (valid_z_dir)
				{
					if (d_y > local_y_direction_control_planes_.d_max)
					{
						if (d_x > local_x_direction_control_planes_.d_max)
						{
							// edgeCube[7]
							virtual_cube_target = edge_adjacent_virtual_cube_[7];
						}
						else
						{
							// edgeCube[5]
							virtual_cube_target = edge_adjacent_virtual_cube_[5];
						}
					}
					else
					{
						if (d_x > local_x_direction_control_planes_.d_max)
						{
							// edgeCube[6]
							virtual_cube_target = edge_adjacent_virtual_cube_[6];
						}
						else
						{
							// edgeCube[4]
							virtual_cube_target = edge_adjacent_virtual_cube_[4];
						}
					}
				}
				else
				{
					if (d_x > local_x_direction_control_planes_.d_max)
					{
						if (d_y > local_y_direction_control_planes_.d_max)
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// vertexCube[7]
								virtual_cube_target = vertex_adjacent_virtual_cube_[7];
							}
							else
							{
								// vertexCube[5]
								virtual_cube_target = vertex_adjacent_virtual_cube_[5];
							}
						}
						else
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// vertexCube[3]
								virtual_cube_target = vertex_adjacent_virtual_cube_[3];
							}
							else
							{
								// vertexCube[1]
								virtual_cube_target = vertex_adjacent_virtual_cube_[1];
							}
						}
					}
					else
					{
						if (d_y > local_y_direction_control_planes_.d_max)
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// vertexCube[6]
								virtual_cube_target = vertex_adjacent_virtual_cube_[6];
							}
							else
							{
								// vertexCube[4]
								virtual_cube_target = vertex_adjacent_virtual_cube_[4];
							}
						}
						else
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// vertexCube[2]
								virtual_cube_target = vertex_adjacent_virtual_cube_[2];
							}
							else
							{
								// vertexCube[0]
								virtual_cube_target = vertex_adjacent_virtual_cube_[0];
							}
						}
					}
				}

				compute_mvc_on_point_outside_cage(surface_point, surface_point_index, virtual_cube_target);
			}
		});
	}

	void set_center_control_cage(Vec3& center)
	{
		control_cage_center_ = center;
	}

	void update_influence_cage_position()
	{
		foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, v);

			value<Vec3>(*(this->influence_cage_), this->influence_cage_vertex_position_, v) =
				((cage_point - control_cage_center_) * 1.5) + control_cage_center_;

			return true;
		});
	}

	void update_deformation_object(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_indices =
			get_attribute<uint32, Vertex>(object, "vertex_index");

		this->influence_area_->foreach_cell([&](Vertex v) -> bool {
			uint32 object_vertex_index = value<uint32>(object, object_vertex_indices, v);

			Vec3 new_pos_ = {0.0, 0.0, 0.0};

			foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, cv);

				uint32 cage_point_index = value<uint32>(*control_cage_, control_cage_vertex_index_, cv);

				new_pos_ += control_cage_coords_(object_vertex_index, cage_point_index) * cage_point;

				return true;
			});

			if ((new_pos_[0] != 0.0) && (new_pos_[1] != 0.0) && (new_pos_[2] != 0.0))
			{
				value<Vec3>(object, object_vertex_position, v) = new_pos_;
			}
			else
			{
				std::cout << "empty " << std::endl;
			}

			return true;
		});
	}

	/*void set_up_attenuation(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, Vertex>(object, "vertex_index");

		uint32 nbv_object = nb_cells<Vertex>(object);

		control_area_validity_.resize(nbv_object);
		control_area_validity_.setZero();

		foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, vertex_position, v);

			//bool inside_cage = local_mvc_pt_control_area(surface_point);

			if (inside_cage)
			{
				uint32 surface_point_index = value<uint32>(object, object_vertex_index, v);
				control_area_validity_(surface_point_index) = 1.0f;
			}

			return true;
		});

		this->attenuation_.resize(nbv_object);
		this->attenuation_.setZero();

		//compute_attenuation_cage(object);

	}*/

private:
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_comp_;

	float m_hFactor;

	Vec3 control_cage_center_;
	Vec3 control_cage_bb_min_;
	Vec3 control_cage_bb_max_;

	Vec3 influence_cage_bb_min_;
	Vec3 influence_cage_bb_max_;

	Local_direction_control_planes local_x_direction_control_planes_;
	Local_direction_control_planes local_y_direction_control_planes_;
	Local_direction_control_planes local_z_direction_control_planes_;

	Eigen::VectorXd control_area_validity_;

	std::vector<VirtualCube> face_adjacent_virtual_cube_;
	std::vector<VirtualCube> edge_adjacent_virtual_cube_;
	std::vector<VirtualCube> vertex_adjacent_virtual_cub_e;

	Eigen::Matrix<Vec2, Eigen::Dynamic, Eigen::Dynamic> normal_weights_;

	std::shared_ptr<Attribute<uint32>> control_cage_vertex_index_;

	std::shared_ptr<Attribute<bool>> control_cage_marked_vertices_;

	void init_triangles()
	{

		foreach_cell(*control_cage_, [&](Face fc) -> bool {
			std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*control_cage_, fc);

			std::vector<Points> points;
			for (std::size_t p = 0; p < 4; p++)
			{
				Point point_i;
				CMap2::Vertex vertex_i = face_vertices_[p];
				point_i.position = value<Vec3>(*control_cage_, control_cage_vertex_position_, vertex_i);
				point_i.control_cage_index = value<uint32_t>(*control_cage_, control_cage_vertex_index_, vertex_i);
				point_i.inside_control_cage = true;
				points.push_back(point_i);
			}

			Triangle triangle1;
			triangle1.points = {points[1], points[3], points[0]};
			triangle1.normal = (cgogn::geometry::normal(triangle1.points[0].position, triangle1.points[1].position,
														triangle1.points[2].position))
								   .normalized();
			cage_triangles_.push_back(triangle1);

			Triangle triangle2;
			triangle2.points = {points[1], points[2], points[3]};
			triangle2.normal = (cgogn::geometry::normal(triangle2.points[0].position, triangle2.points[1].position,
														triangle2.points[2].position))
								   .normalized();
			cage_triangles_.push_back(triangle2);

			return true;
		});
	}

	/**
	 * Delimit the control cage area in terms of planes
	 * Plane of equation ax + by + cz = d
	 * Assign each face (or pair of triangles) to its corresponding direction
	 */
	void init_control_cage_plane()
	{
		// TODO set for local frame
		const Vec3 x_dir = {1.0, 0.0, 0.0}, y_dir = {0.0, 1.0, 0.0}, z_dir = {0.0, 0.0, 1.0};

		double d_x_min = 1000.0, d_x_max = -1000.0, d_y_min = 1000.0, d_y_max = -1000.0, d_z_min = 1000.0,
			   d_z_max = -1000.0;

		std::pair<Triangle, Triangle> face_x_min, face_x_max, face_y_min, face_y_max, face_z_min, face_z_max;

		for (std::size_t i = 0; i < 6; i++)
		{ // for each face
			const Triangle triangle1 = cage_triangles_[2 * i];
			const Triangle triangle2 = cage_triangles_[2 * i + 1];

			const Vec3 normal = triangle1.normal;

			const Vec3 triangle_point = triangle1.points[0].position;

			if (normal.dot(x_dir) != 0.0)
			{
				const double d = triangle_point.dot(x_dir);
				if (d < d_x_min)
				{
					d_x_min = d;
					face_x_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_x_max)
				{
					d_x_max = d;
					face_x_max = std::make_pair(triangle1, triangle2);
				}
			}
			else if (normal.dot(y_dir) != 0.0)
			{
				const double d = triangle_point.dot(y_dir);
				if (d < d_y_min)
				{
					d_y_min = d;
					face_y_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_y_max)
				{
					d_y_max = d;
					face_y_max = std::make_pair(triangle1, triangle2);
				}
			}
			else
			{
				const double d = triangle_point.dot(z_dir);
				if (d < d_z_min)
				{
					d_z_min = d;
					face_z_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_z_max)
				{
					d_z_max = d;
					face_z_max = std::make_pair(triangle1, triangle2);
				}
			}
		}

		local_x_direction_control_planes_.d_min = d_x_min;
		local_x_direction_control_planes_.d_max = d_x_max;
		local_x_direction_control_planes_.d_gap = d_x_max - d_x_min;
		local_x_direction_control_planes_.direction = x_dir;
		local_x_direction_control_planes_.triangles_d_min = face_x_min;
		local_x_direction_control_planes_.triangles_d_max = face_x_max;
		local_x_direction_control_planes_.shift_after_d_max = {2.0 * d_x_max - d_x_min, 0.0, 0.0};	// (d_max + d_gap)
		local_x_direction_control_planes_.shift_before_d_min = {2.0 * d_x_min - d_x_max, 0.0, 0.0}; // (d_min - d_gap)

		local_y_direction_control_planes_.d_min = d_y_min;
		local_y_direction_control_planes_.d_max = d_y_max;
		local_y_direction_control_planes_.d_gap = d_y_max - d_y_min;
		local_y_direction_control_planes_.direction = y_dir;
		local_y_direction_control_planes_.triangles_d_min = face_y_min;
		local_y_direction_control_planes_.triangles_d_max = face_y_max;
		local_y_direction_control_planes_.shift_after_d_max = {0.0, 2.0 * d_y_max - d_y_min, 0.0};	// (d_max + d_gap)
		local_y_direction_control_planes_.shift_before_d_min = {0.0, 2.0 * d_y_min - d_y_max, 0.0}; // (d_min - d_gap)

		local_z_direction_control_planes_.d_min = d_z_min;
		local_z_direction_control_planes_.d_max = d_z_max;
		local_z_direction_control_planes_.d_gap = d_z_max - d_z_min;
		local_z_direction_control_planes_.direction = z_dir;
		local_z_direction_control_planes_.triangles_d_min = face_z_min;
		local_z_direction_control_planes_.triangles_d_max = face_z_max;
		local_z_direction_control_planes_.shift_after_d_max = {0.0, 0.0, 2.0 * d_z_max - d_z_min};	// (d_max + d_gap)
		local_z_direction_control_planes_.shift_before_d_min = {0.0, 0.0, 2.0 * d_z_min - d_z_max}; // (d_min - d_gap)
	}

	void init_virtual_cubes()
	{

		init_face_adjacent_virtual_cubes();

		init_edge_ajacent_virtual_cubes();

		init_vertex_adjacent_virtual_cubes();
	}

	void init_face_adjacent_virtual_cubes()
	{
		// One common face with control cage
		Virtual_cube face_adjacent0 = get_virtual_cube_triangles(local_x_direction_control_planes_.triangles_d_min,
																 local_x_direction_control_planes_.shift_before_d_min);
		face_adjacent_virtual_cube_.push_back(face_adjacent0);

		Virtual_cube face_adjacent1 = get_virtual_cube_triangles(local_x_direction_control_planes_.triangles_d_max,
																 local_x_direction_control_planes_.shift_after_d_max);
		face_adjacent_virtual_cube_.push_back(face_adjacent1);

		Virtual_cube face_adjacent2 = get_virtual_cube_triangles(local_y_direction_control_planes_.triangles_d_min,
																 local_y_direction_control_planes_.shift_before_d_min);
		face_adjacent_virtual_cube_.push_back(face_adjacent2);

		Virtual_cube face_adjacent3 = get_virtual_cube_triangles(local_y_direction_control_planes_.triangles_d_max,
																 local_y_direction_control_planes_.shift_after_d_max);
		face_adjacent_virtual_cube_.push_back(face_adjacent3);

		Virtual_cube face_adjacent4 = get_virtual_cube_triangles(local_z_direction_control_planes_.triangles_d_min,
																 local_z_direction_control_planes_.shift_before_d_min);
		face_adjacent_virtual_cube_.push_back(face_adjacent4);

		Virtual_cube face_adjacent5 = get_virtual_cube_triangles(local_z_direction_control_planes_.triangles_d_max,
																 local_z_direction_control_planes_.shift_after_d_max);
		face_adjacent_virtual_cube_.push_back(face_adjacent5);
	}

	void init_edge_adjacent_virtual_cubes()
	{
		// One common edge with control cage
		std::vector<Vec3> intersect_positions0 = find_intersection_positions_face(
			local_x_direction_control_planes_.triangles_d_min, local_z_direction_control_planes_.triangles_d_min);
		std::pair<Triangle, Triangle> new_face0 =
			get_face_from_intersecting_edge(intersect_positions0, local_x_direction_control_planes_.shift_before_d_min);
		Virtual_cube edge_adjacent0 =
			get_virtual_cube_triangles(new_face0, local_z_direction_control_planes_.shift_before_d_min);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent0);

		std::vector<Vec3> intersect_positions1 = find_intersection_positions_face(
			local_x_direction_control_planes_.triangles_d_max, local_z_direction_control_planes_.triangles_d_min);
		std::pair<Triangle, Triangle> new_face1 =
			get_face_from_intersecting_edge(intersect_positions1, local_x_direction_control_planes_.shift_after_d_max);
		Virtual_cube edge_adjacent1 =
			get_virtual_cube_triangles(new_face1, local_z_direction_control_planes_.shift_before_d_min);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent1);

		std::vector<Vec3> intersect_positions2 = find_intersection_positions_face(
			local_x_direction_control_planes_.triangles_d_min, local_z_direction_control_planes_.triangles_d_max);
		std::pair<Triangle, Triangle> new_face2 =
			get_face_from_intersecting_edge(intersect_positions2, local_x_direction_control_planes_.shift_before_d_min);
		Virtual_cube edge_adjacent2 =
			get_virtual_cube_triangles(new_face2, local_z_direction_control_planes_.shift_after_d_max);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent2);

		std::vector<Vec3> intersect_positions3 = find_intersection_positions_face(
			local_x_direction_control_planes_.triangles_d_max, local_z_direction_control_planes_.triangles_d_max);
		std::pair<Triangle, Triangle> new_face3 =
			get_face_from_intersecting_edge(intersect_positions3, local_x_direction_control_planes_.shift_after_d_max);
		Virtual_cube edge_adjacent3 =
			get_virtual_cube_triangles(new_face3, local_z_direction_control_planes_.shift_after_d_max);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent3);

		std::vector<Vec3> intersect_positions4 = find_intersection_positions_face(
			local_y_direction_control_planes_.triangles_d_min, local_x_direction_control_planes_.triangles_d_min);
		std::pair<Triangle, Triangle> new_face4 =
			get_face_from_intersecting_edge(intersect_positions4, local_y_direction_control_planes_.shift_before_d_min);
		Virtual_cube edge_adjacent4 =
			get_virtual_cube_triangles(new_face4, local_x_direction_control_planes_.shift_before_d_min);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent4);

		std::vector<Vec3> intersect_positions5 = find_intersection_positions_face(
			local_y_direction_control_planes_.triangles_d_max, local_x_direction_control_planes_.triangles_d_min);
		std::pair<Triangle, Triangle> new_face5 =
			get_face_from_intersecting_edge(intersect_positions5, local_y_direction_control_planes_.shift_after_d_max);
		Virtual_cube edge_adjacent5 =
			get_virtual_cube_triangles(new_face5, local_x_direction_control_planes_.shift_before_d_min);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent5);

		std::vector<Vec3> intersect_positions6 = find_intersection_positions_face(
			local_y_direction_control_planes_.triangles_d_min, local_x_direction_control_planes_.triangles_d_max);
		std::pair<Triangle, Triangle> new_face6 =
			get_face_from_intersecting_edge(intersect_positions6, local_y_direction_control_planes_.shift_before_d_min);
		Virtual_cube edge_adjacent6 =
			get_virtual_cube_triangles(new_face6, local_x_direction_control_planes_.shift_after_d_max);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent6);

		std::vector<Vec3> intersect_positions7 = find_intersection_positions_face(
			local_y_direction_control_planes_.triangles_d_max, local_x_direction_control_planes_.triangles_d_max);
		std::pair<Triangle, Triangle> new_face7 =
			get_face_from_intersecting_edge(intersect_positions7, local_y_direction_control_planes_.shift_after_d_max);
		Virtual_cube edge_adjacent7 =
			get_virtual_cube_triangles(new_face7, local_x_direction_control_planes_.shift_after_d_max);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent7);

		std::vector<Vec3> intersect_positions8 = find_intersection_positions_face(
			local_y_direction_control_planes_.triangles_d_min, local_z_direction_control_planes_.triangles_d_min);
		std::pair<Triangle, Triangle> new_face8 =
			get_face_from_intersecting_edge(intersect_positions8, local_y_direction_control_planes_.shift_before_d_min);
		Virtual_cube edge_adjacent8 =
			get_virtual_cube_triangles(new_face8, local_z_direction_control_planes_.shift_before_d_min);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent8);

		std::vector<Vec3> intersect_positions9 = find_intersection_positions_face(
			local_y_direction_control_planes_.triangles_d_max, local_z_direction_control_planes_.triangles_d_min);
		std::pair<Triangle, Triangle> new_face9 =
			get_face_from_intersecting_edge(intersect_positions9, local_y_direction_control_planes_.shift_after_d_max);
		Virtual_cube edge_adjacent9 =
			get_virtual_cube_triangles(new_face9, local_z_direction_control_planes_.shift_before_d_min);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent9);

		std::vector<Vec3> intersect_positions10 = find_intersection_positions_face(
			local_y_direction_control_planes_.triangles_d_min, local_z_direction_control_planes_.triangles_d_max);
		std::pair<Triangle, Triangle> new_face10 = get_face_from_intersecting_edge(
			intersect_positions10, local_y_direction_control_planes_.shift_before_d_min);
		Virtual_cube edge_adjacent10 =
			get_virtual_cube_triangles(new_face10, local_z_direction_control_planes_.shift_after_d_max);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent10);

		std::vector<Vec3> intersect_positions11 = find_intersection_positions_face(
			local_y_direction_control_planes_.triangles_d_max, local_z_direction_control_planes_.triangles_d_max);
		std::pair<Triangle, Triangle> new_face11 =
			get_face_from_intersecting_edge(intersect_positions11, local_y_direction_control_planes_.shift_after_d_max);
		Virtual_cube edge_adjacent11 =
			get_virtual_cube_triangles(new_face11, local_z_direction_control_planes_.shift_after_d_max);

		edge_adjacent_virtual_cube_.push_back(edge_adjacent11);
	}

	void init_vertex_adjacent_virtual_cubes()
	{
		Vec3 intersection_position0 = find_intersection_position(local_x_direction_control_planes_.triangles_d_min,
																 local_y_direction_control_planes_.triangles_d_min,
																 local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face0 = get_face_from_intersecting_vertex(
			intersection_position0, local_x_direction_control_planes_.shift_before_d_min,
			local_y_direction_control_planes_.shift_before_d_min);
		Virtual_cube vertex_adjacent0 =
			get_virtual_cube_triangles(new_face0, local_z_direction_control_planes_.shift_before_d_min);

		vertex_adjacent_virtual_cube_.push_back(vertex_adjacent0);

		Vec3 intersection_position1 = find_intersection_position(local_x_direction_control_planes_.triangles_d_max,
																 local_y_direction_control_planes_.triangles_d_min,
																 local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face1 = get_face_from_intersecting_vertex(
			intersection_position1, local_x_direction_control_planes_.shift_after_d_max,
			local_y_direction_control_planes_.shift_before_d_min);
		Virtual_cube vertex_adjacent1 =
			get_virtual_cube_triangles(new_face1, local_z_direction_control_planes_.shift_before_d_min);

		vertex_adjacent_virtual_cube_.push_back(vertex_adjacent1);

		Vec3 intersection_position2 = find_intersection_position(local_x_direction_control_planes_.triangles_d_min,
																 local_y_direction_control_planes_.triangles_d_min,
																 local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face2 = get_face_from_intersecting_vertex(
			intersection_position2, local_x_direction_control_planes_.shift_before_d_min,
			local_y_direction_control_planes_.shift_before_d_min);
		Virtual_cube vertex_adjacent2 =
			get_virtual_cube_triangles(new_face2, local_z_direction_control_planes_.shift_after_d_max);

		vertex_adjacent_virtual_cube_.push_back(vertex_adjacent2);

		Vec3 intersection_position3 = find_intersection_position(local_x_direction_control_planes_.triangles_d_max,
																 local_y_direction_control_planes_.triangles_d_min,
																 local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face3 = get_face_from_intersecting_vertex(
			intersection_position3, local_x_direction_control_planes_.shift_after_d_max,
			local_y_direction_control_planes_.shift_before_d_min);
		Virtual_cube vertex_adjacent3 =
			get_virtual_cube_triangles(new_face3, local_z_direction_control_planes_.shift_after_d_max);

		vertex_adjacent_virtual_cube_.push_back(vertex_adjacent3);

		Vec3 intersection_position4 = find_intersection_position(local_x_direction_control_planes_.triangles_d_min,
																 local_y_direction_control_planes_.triangles_d_max,
																 local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face4 = get_face_from_intersecting_vertex(
			intersection_position4, local_x_direction_control_planes_.shift_before_d_min,
			local_y_direction_control_planes_.shift_after_d_max);
		Virtual_cube vertex_adjacent4 =
			get_virtual_cube_triangles(new_face4, local_z_direction_control_planes_.shift_before_d_min);

		vertex_adjacent_virtual_cube_.push_back(vertex_adjacent4);

		Vec3 intersection_position5 = find_intersection_position(local_x_direction_control_planes_.triangles_d_max,
																 local_y_direction_control_planes_.triangles_d_max,
																 local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face5 = get_face_from_intersecting_vertex(
			intersection_position5, local_x_direction_control_planes_.shift_after_d_max,
			local_y_direction_control_planes_.shift_after_d_max);
		Virtual_cube vertex_adjacent5 =
			get_virtual_cube_triangles(new_face5, local_z_direction_control_planes_.shift_before_d_min);

		vertex_adjacent_virtual_cube_.push_back(vertex_adjacent5);

		Vec3 intersection_position6 = find_intersection_position(local_x_direction_control_planes_.triangles_d_min,
																 local_y_direction_control_planes_.triangles_d_max,
																 local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face6 = get_face_from_intersecting_vertex(
			intersection_position6, local_x_direction_control_planes_.shift_before_d_min,
			local_y_direction_control_planes_.shift_after_d_max);
		Virtual_cube vertex_adjacent6 =
			get_virtual_cube_triangles(new_face6, local_z_direction_control_planes_.shift_after_d_max);

		vertex_adjacent_virtual_cube_.push_back(vertex_adjacent6);

		Vec3 intersection_position7 = find_intersection_position(local_x_direction_control_planes_.triangles_d_max,
																 local_y_direction_control_planes_.triangles_d_max,
																 local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face7 = get_face_from_intersecting_vertex(
			intersection_position7, local_x_direction_control_planes_.shift_after_d_max,
			local_y_direction_control_planes_.shift_after_d_max);
		Virtual_cube vertex_adjacent7 =
			get_virtual_cube_triangles(new_face7, local_z_direction_control_planes_.shift_after_d_max);

		vertex_adjacent_virtual_cube_.push_back(vertex_adjacent7);
	}

	Virtual_cube get_virtual_cube_triangles(const std::pair<Triangle, Triangle> face, const Vec3& shift_vector)
	{
		Virtual_cube new_virtual_cube;

		std::vector<Triangle> virtual_cube_triangles;

		Triangle triangle1 = face.first;
		triangle1.virtual_cage_indices = {1, 3, 0};
		Triangle triangle2 = face.second;
		triangle2.virtual_cage_indices = {1, 2, 3};
		virtual_cube_triangles.push_back(triangle1);
		virtual_cube_triangles.push_back(triangle2);

		Point point4; 
		point4.position = triangle1.points[2].position + shift_vector; 
		point4.inside_control_cage = false; 

		Point point5; 
		point5.position = triangle1.points[0].position + shift_vector; 
		point5.inside_control_cage = false; 

		Point point6; 
		point.position = triangle2.points[1].position + shift_vector; 
		point6.inside_control_cage = false;

		Point point7; 
		point7.position = triangle1.points[1].position + shift_vector; 
		point7.inside_control_cage = false;

		Triangle triangle3, triangle4; 
		triangle3.points = {point5, point7, point4};
		triangle3.virtual_cage_indices = {5, 7, 4};

		triangle4.points = {point5, point6, point7};
		triangle4.virtual_cage_indices = {5, 6, 7};

		std::vector<Point> virtual_cube_points(8);
		virtual_cube_points[0] = triangle1.points[2]; 
		virtual_cube_points[1] = triangle1.points[0]; 
		virtual_cube_points[2] = triangle2.points[1];
		virtual_cube_points[3] = triangle1.points[1];

		virtual_cube_points[4] = triangle3.points[2]; 
		virtual_cube_points[5] = triangle3.points[0];
		virtual_cube_points[6] = triangle4.points[1];
		virtual_cube_points[7] = triangle3.points[1];

		new_virtual_cube.points = virtual_cube_points;

		Triangle triangle5, triangle6, triangle7, triangle8, triangle9, triangle10, triangle11,
			triangle12;

		triangle5.points = {triangle1.points[0], triangle3.points[2], triangle1.points[2]};
		triangle5.virtual_cage_indices = {triangle1.virtual_cage_indices[0], triangle3.virtual_cage_indices[2],
										  triangle1.virtual_cage_indices[2]};

		triangle6.points = {triangle1.points[0], triangle3.points[0], triangle3.points[2]};
		triangle6.virtual_cage_indices = {triangle1.virtual_cage_indices[0], triangle3.virtual_cage_indices[0],
										  triangle3.virtual_cage_indices[2]};

		triangle7.points = {triangle1.points[2], triangle3.points[1], triangle3.points[2]};
		triangle7.virtual_cage_indices = {triangle1.virtual_cage_indices[2], triangle3.virtual_cage_indices[1],
										  triangle3.virtual_cage_indices[2]};

		triangle8.points = {triangle1.points[2], triangle1.points[1], triangle3.points[1]};
		triangle8.virtual_cage_indices = {triangle1.virtual_cage_indices[2], triangle1.virtual_cage_indices[1],
										  triangle3.virtual_cage_indices[1]};

		triangle9.points = {triangle1.points[0], triangle4.points[1], triangle3.points[0]};
		triangle9.virtual_cage_indices = {triangle1.virtual_cage_indices[0], triangle4.virtual_cage_indices[1],
										  triangle3.virtual_cage_indices[0]};

		triangle10.points = {triangle1.points[0], triangle2.points[1], triangle4.points[1]};
		triangle10.virtual_cage_indices = {triangle1.virtual_cage_indices[0], triangle2.virtual_cage_indices[1],
										   triangle4.virtual_cage_indices[1]};

		triangle11.points = {triangle4.points[1], triangle1.points[1], triangle3.points[1]};
		triangle11.virtual_cage_indices = {triangle4.virtual_cage_indices[1], triangle1.virtual_cage_indices[1],
										   triangle3.virtual_cage_indices[1]};

		triangle12.points = {triangle4.points[1], triangle2.points[1], triangle1.positions[1]};
		triangle12.virtual_cage_indices = {triangle4.virtual_cage_indices[1], triangle2.virtual_cage_indices[1],
										   triangle1.virtual_cage_indices[1]};

		virtual_cube_triangles.push_back(triangle3);
		virtual_cube_triangles.push_back(triangle4);
		virtual_cube_triangles.push_back(triangle5);
		virtual_cube_triangles.push_back(triangle6);
		virtual_cube_triangles.push_back(triangle7);
		virtual_cube_triangles.push_back(triangle8);
		virtual_cube_triangles.push_back(triangle9);
		virtual_cube_triangles.push_back(triangle10);
		virtual_cube_triangles.push_back(triangle11);
		virtual_cube_triangles.push_back(triangle12);

		new_virtual_cube.triangles = virtual_cube_triangles;

		/*
		// TODO: find better solution (concatenate?)
		std::vector<uint32_t> virtual_cube_indices; 
		for (std::size_t t = 0; t < virtual_cube_triangles.size(); t++){
			for (std::size_t i = 0; i < 3; i++){
				virtual_cube_indices.push_back(virtual_cube_triangles[t].virtual_cage_indices[i]); 
			}
		}

		new_virtual_cube.virtual_cage_indices = virtual_cage_indices; */

		return virtual_cube;
	}

	bool compute_mvc_on_point_inside_cage(const Vec3& surface_point, const uint32& surface_point_index)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*control_cage_, "vertex_index");

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_control_cage_coords_;

		w_control_cage_coords_.resize(nbv_cage);
		w_control_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, v);

			uint32 cage_point_index = value<uint32>(*control_cage_, control_cage_vertex_index_, v);

			d[cage_point_index] = (surface_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				control_cage_coords_(surface_point_index, cage_point_index) = 1.0;
				return true;
			}

			u[cage_point_index] = (cage_point - surface_point) / d[cage_point_index];

			return true;
		});

		double l[3];
		double theta[3];
		double w[3];
		double c[3];
		double s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{

			std::vector<uint32> triangle_index = 
				{cage_triangles_[t].points[0].control_cage_index, 
				cage_triangles_[t].points[1].control_cage_index, 
				cage_triangles_[t].points[2].control_cage_index}

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = (u[triangle_index[(i + 1) % 3]] - u[triangle_index[(i + 2) % 3]]).norm();
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				theta[i] = 2.0 * asin(l[i] / 2.0);
			}

			double h = (theta[0] + theta[1] + theta[2]) / 2.0;
			if (M_PI - h < epsilon)
			{
				for (std::size_t i = 0; i < 3; i++)
				{
					w[i] = sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_control_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_control_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_control_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / (sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) - 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || fabs(s[2]) < epsilon)
			{
				continue; // eta is on the same plane, outside t  ->  ignore triangle t :
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = (theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] - c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					   (2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_control_cage_coords_[triangle_index[0]] += w[0];
			w_control_cage_coords_[triangle_index[1]] += w[1];
			w_control_cage_coords_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_index = value<uint32>(*control_cage_, control_cage_vertex_index_, v);

			control_cage_coords_(surface_point_index, cage_point_index) =
				w_control_cage_coords_[cage_point_index] / sumWeights;

			return true;
		});

		return false;
	}

	bool compute_mvc_on_point_outside_cage(const Vec3& surface_point, const uint32& surface_point_index,
										   const Virtual_cube virtual_cube_target)
	{
		uint32 nbv_cage = 8; // nb_cells<Vertex>(*control_cage_);

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_control_cage_coords_;

		w_control_cage_coords_.resize(nbv_cage);
		w_control_cage_coords_.setZero();

		Eigen::VectorXd virtual_cage_coords_;
		virtual_cage_coords_.resize(nbv_cage);
		virtual_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);


		for (std::size_t p = 0; p < virtual_cube_target.points.size(); p++)
		{
			const Vec3& cage_point = virtual_cube_target[p].position;

			uint32 cage_point_index = p;

			d[cage_point_index] = (surface_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				virtual_cage_coords_(surface_point_index, cage_point_index) = 1.0;
				return true;
			}

			u[cage_point_index] = (cage_point - surface_point) / d[cage_point_index];
		}

		double l[3];
		double theta[3];
		double w[3];
		double c[3];
		double s[3];

		for (std::size_t t = 0; t < virtual_cube_triangles.size(); t++)
		{

			std::vector<uint32> triangle_index = virtual_cube_triangles[t].virtual_cage_indices;

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = (u[triangle_index[(i + 1) % 3]] - u[triangle_index[(i + 2) % 3]]).norm();
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				theta[i] = 2.0 * asin(l[i] / 2.0);
			}

			double h = (theta[0] + theta[1] + theta[2]) / 2.0;
			if (M_PI - h < epsilon)
			{
				for (std::size_t i = 0; i < 3; i++)
				{
					w[i] = sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_control_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_control_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_control_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / (sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) - 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || fabs(s[2]) < epsilon)
			{
				continue; // eta is on the same plane, outside t  ->  ignore triangle t :
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = (theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] - c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					   (2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_control_cage_coords_[triangle_index[0]] += w[0];
			w_control_cage_coords_[triangle_index[1]] += w[1];
			w_control_cage_coords_[triangle_index[2]] += w[2];
		}

		for (std::size_t p = 0; p < 8; p++)
		{
			uint32 cage_point_index = p;

			virtual_cage_coords_(surface_point_index, cage_point_index) =
				w_control_cage_coords_[cage_point_index] / sumWeights;
		}

		

		return false;
	}

	void bind_mvc_customed_vertices(const Vec3 surface_point, uint32 surface_point_index,
									std::vector<Vec3> position_vertices)
	{
	}

	bool check_mvc_on_point_inside_influence_cage(Vec3 point)
	{
		std::shared_ptr<Attribute<bool>> cage_vertex_marked =
			get_attribute<bool, Vertex>(*influence_cage_, "marked_vertices");

		DartMarker dm(*influence_cage_);

		bool checked = true;
		for (Dart d = influence_cage_->begin(), end = influence_cage_->end(); d != end; d = influence_cage_->next(d))
		{
			Vertex cage_vertex = CMap2::Vertex(d);
			bool vc_marked = value<bool>(*influence_cage_, cage_vertex_marked, cage_vertex);

			if (!dm.is_marked(d) && !vc_marked)
			{

				const Vec3& cage_point = value<Vec3>(*influence_cage_, influence_cage_vertex_position_, cage_vertex);

				float mvc_value =
					compute_mvc(point, d, *influence_cage_, cage_point, influence_cage_vertex_position_.get());

				dm.mark(d);

				value<bool>(*influence_cage_, cage_vertex_marked, cage_vertex) = true;

				if (mvc_value < 0)
				{
					checked = false;
					break;
				}
			}
		}

		parallel_foreach_cell(*influence_cage_, [&](Vertex vc) -> bool {
			value<bool>(*influence_cage_, cage_vertex_marked, vc) = false;
			return true;
		});

		return checked;
	}

	std::vector<Vec3> find_intersection_positions_face(const std::pair<Triangle, Triangle>& face1,
													   const std::pair<Triangle, Triangle>& face2)
	{

		std::vector<Vec3> intersect_positions;
		double epsilon = 0.00000001;

		std::vector<Vec3> positions_face_1;
		positions_face_1.push_back(face1.first.positions[2]);
		positions_face_1.push_back(face1.first.positions[0]);
		positions_face_1.push_back(face1.second.positions[1]);
		positions_face_1.push_back(face1.first.positions[1]);

		std::vector<Vec3> positions_face_2;
		positions_face_2.push_back(face2.first.positions[2]);
		positions_face_2.push_back(face2.first.positions[0]);
		positions_face_2.push_back(face2.second.positions[1]);
		positions_face_2.push_back(face2.first.positions[1]);

		// TODO: find alternative to double loops
		for (std::size_t i = 0; i < 4; i++)
		{
			Vec3 target_position_face1 = positions_face_1[i];
			for (std::size_t j = 0; j < 4; j++)
			{
				Vec3 target_position_face2 = positions_face_2[j];

				const double distance_p1_p2 = (target_position_face1 - target_position_face2).norm();
				if (distance_p1_p2 < epsilon)
				{
					intersect_positions.push_back(target_position_face1);
					continue;
				}
			}
		}

		if (!intersect_positions.size() == 2)
		{
			std::cout << "Warning, size not 2: " << intersect_positions.size() << std::endl;
		}

		return intersect_positions;
	}

	Vec3 find_intersection_position(const std::pair<Triangle, Triangle>& face1,
									const std::pair<Triangle, Triangle>& face2,
									const std::pair<Triangle, Triangle>& face3)
	{
		Vec3 intersection_position;

		double epsilon = 0.00000001;

		std::vector<Vec3> position_face_1;
		position_face_1.push_back(face1.first.positions[2]);
		position_face_1.push_back(face1.first.positions[0]);
		position_face_1.push_back(face1.second.positions[1]);
		position_face_1.push_back(face1.first.positions[1]);

		std::vector<Vec3> position_face_2;
		position_face_2.push_back(face2.first.positions[2]);
		position_face_2.push_back(face2.first.positions[0]);
		position_face_2.push_back(face2.second.positions[1]);
		position_face_2.push_back(face2.first.positions[1]);

		std::vector<Vec3> position_face_3;
		position_face_3.push_back(face3.first.positions[2]);
		position_face_3.push_back(face3.first.positions[0]);
		position_face_3.push_back(face3.second.positions[1]);
		position_face_3.push_back(face3.first.positions[1]);

		// TODO: find alternative to double loops
		for (std::size_t i = 0; i < 4; i++)
		{
			Vec3 target_position_face1 = position_face_1[i];
			for (std::size_t j = 0; j < 4; j++)
			{
				Vec3 target_position_face2 = position_face_2[j];
				const double distance_p1_p2 = (target_position_face1 - target_position_face2).norm();
				if (distance_p1_p2 < epsilon)
				{
					for (std::size_t k = 0; k < 4; k++)
					{
						Vec3 target_position_face3 = position_face_3[k];
						const double distance_p1_p3 = (target_position_face1 - target_position_face3).norm();
						if (distance_p1_p3 < epsilon)
						{
							return target_position_face1;
						}
					}
				}
			}
		}

		return intersection_position;
	}

	std::pair<Triangle, Triangle> get_face_from_intersecting_edge(const std::vector<Vec3>& intersect_positions,
																  const Vec3& shiftVector)
	{
		// create new face composed of the intersecting edge and this edge shifted by y
		const Vec3 face_position0 = intersect_positions[0];
		const Vec3 face_position1 = intersect_positions[1];

		const Vec3 face_position2 = face_position1 + shiftVector;
		const Vec3 face_position3 = face_position0 + shiftVector;

		Triangle local_triangle1, local_triangle2;
		local_triangle1.positions = {face_position1, face_position3, face_position0};
		local_triangle2.positions = {face_position1, face_position2, face_position3};

		return std::make_pair(local_triangle1, local_triangle2);
	}

	std::pair<Triangle, Triangle> get_face_from_intersecting_vertex(const Vec3& intersection_position,
																	const Vec3& shiftVector1, const Vec3& shiftVector2)
	{
		// create new face composed of the intersection vertex and shifted in two directions
		const Vec3 face_position0 = intersection_position;

		const Vec3 face_position1 = face_position0 + shiftVector1;

		const Vec3 face_position2 = (face_position1 + shiftVector1) + shiftVector2;
		const Vec3 face_position3 = face_position0 + shiftVector2;

		Triangle local_triangle1, local_triangle2;
		local_triangle1.positions = {face_position1, face_position3, face_position0};
		local_triangle2.positions = {face_position1, face_position2, face_position3};

		return std::make_pair(local_triangle1, local_triangle2);
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_

/*void compute_attenuation_cage(MESH& object)
{

std::shared_ptr<Attribute<uint32>> cage_face_indices =
	add_attribute<uint32, Face>(*(this->influence_cage_), "face_indices");

cgogn::modeling::set_attribute_face_index(*(this->influence_cage_),
											cage_face_indices.get());

std::shared_ptr<Attribute<uint32>> object_vertex_index =
	get_attribute<uint32, Vertex>(object, "vertex_index");

std::shared_ptr<Attribute<uint32>> cage_vertex_index =
	get_attribute<uint32, Vertex>(*control_cage_, "vertex_index");

std::shared_ptr<Attribute<uint32>> i_cage_vertex_index =
	get_attribute<uint32, Vertex>(*(this->influence_cage_), "vertex_index");

uint32 nbf_cage = 2 * nb_cells<Face>(*(this->influence_cage_));
uint32 nbv_cage = nb_cells<Vertex>(*(this->influence_cage_));

// first loop to find h
//
float h = 0.0f;
float max_dist = 0.0f;
std::vector<Vec2> attenuation_points;
this->influence_area_->foreach_cell([&](Vertex v) {
	uint32 surface_point_index = value<uint32>(object, object_vertex_index, v);

	float i_dist = this->cage_influence_distance(surface_point_index, nbf_cage, nbv_cage);

	this->attenuation_(surface_point_index) = (float)sin(0.5*M_PI * (i_dist ));
	/*if (control_area_validity_(surface_point_index) == 1.0f)
	{

		if (i_dist > h)
		{
			h = i_dist;
		}

		this->attenuation_(surface_const double point_index) = 1.0f;
	}
	else
	{

		if (i_dist > max_dist){
			max_dist = i_dist;
		}
		attenuation_points.push_back({surface_point_index, i_dist});
	}*/
//});

/*for (unsigned int i = 0; i < attenuation_points.size(); i++)
{
	//this->attenuation_(attenuation_points[i][0]) = 0.5f * ((float)sin(M_PI * ((attenuation_points[i][1] / h) -
0.5f))) + 0.5f; this->attenuation_(attenuation_points[i][0]) = (float)sin(0.5*M_PI * (attenuation_points[i][1] /
h));
}*/
//}

/*bool local_mvc_pt_control_area(Vec3 pt)
{
	std::shared_ptr<Attribute<uint32>> i_vertex_index =
		get_attribute<uint32, Vertex>(*control_cage_, "vertex_index");

	std::shared_ptr<Attribute<bool>> cage_vertex_marked =
		get_attribute<bool, Vertex>(*control_cage_, "marked_vertices");

	DartMarker dm(*control_cage_);

	bool checked = true;
	for (Dart d = control_cage_->begin(), end = control_cage_->end(); d != end; d = control_cage_->next(d))
	{
		Vertex cage_vertex = CMap2::Vertex(d);
		bool vc_marked = value<bool>(*control_cage_, cage_vertex_marked, cage_vertex);

		if (!dm.is_marked(d) && !vc_marked)
		{

			const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, cage_vertex);

			float mvc_value = cgogn::modeling::compute_mvc(pt, d, *control_cage_, cage_point,
control_cage_vertex_position_.get());

			dm.mark(d);

			value<bool>(*control_cage_, cage_vertex_marked, cage_vertex) = true;

			if (mvc_value < 0)
			{
				checked = false;
				break;
			}
		}
	}

	parallel_foreach_cell(*control_cage_, [&](Vertex vc) -> bool {
		value<bool>(*control_cage_, cage_vertex_marked, vc) = false;
		return true;
	});

	return checked;
}*/
