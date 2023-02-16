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

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/modeling/types/space_deformation_tool.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class CageDeformationTool : public SpaceDeformationTool<MESH>
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

public:
	MESH* control_cage_;
	std::shared_ptr<Attribute<Vec3>> control_cage_vertex_position_;
	std::shared_ptr<Attribute<Vec3>> influence_cage_vertex_position_;

	Eigen::VectorXd attenuation_;

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;

	std::vector<std::vector<CMap2::Vertex>> cage_triangles_;
	std::vector<Vec3> cage_triangles_normal_;
	std::vector<std::pair<Vec3, Vec3>> cage_triangles_edge_;

	MESH* influence_cage_;

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

		foreach_cell(*control_cage_, [&](Face fc) -> bool {
			std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*control_cage_, fc);

			// triangle 1
			std::vector<CMap2::Vertex> triangle1_vertex(3);
			triangle1_vertex[0] = face_vertices_[1];
			triangle1_vertex[1] = face_vertices_[3];
			triangle1_vertex[2] = face_vertices_[0];

			cage_triangles_.push_back(triangle1_vertex);

			std::vector<CMap2::Vertex> triangle2_vertex(3);
			triangle2_vertex[0] = face_vertices_[1];
			triangle2_vertex[1] = face_vertices_[2];
			triangle2_vertex[2] = face_vertices_[3];

			cage_triangles_.push_back(triangle2_vertex);

			return true;
		});

		
		// TODO set for local frame
		const Vec3 x_dir = {1.0, 0.0, 0.0}, y_dir = {0.0, 1.0, 0.0}, z_dir = {0.0, 0.0, 1.0};

		double d_x_min = 1000.0, d_x_max = -1000.0, d_y_min = 1000.0, d_y_max = -1000.0, d_z_min = 1000.0,
			   d_z_max = -1000.0;

		/*foreach_cell(*control_cage_, [&](Face f) -> bool {
			const Vec3 normal = value<Vec3>(*control_cage_, cage_face_normal_, f);

			std::vector<Vertex> vertices = incident_vertices(*control_cage_, f);
			const Vec3 surface_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, vertices[0]);

			if (normal.dot(x_dir) != 0.0)
			{

				double d = -(surface_point.dot(x_dir));
				if (d < d_x_min)
				{
					d_x_min = d;
					control_cage_face_d_x_min_ = f;
				}

				if (d > d_x_max)
				{
					d_x_max = d;
					control_cage_face_d_x_max_ = f;
				}
			}
			else if (normal.dot(y_dir) != 0.0)
			{
				double d = -(surface_point.dot(y_dir));
				if (d < d_y_min)
				{
					d_y_min = d;
					control_cage_face_d_y_min_ = f;
				}

				if (d > d_y_max)
				{
					d_y_max = d;
					control_cage_face_d_y_min_ = f;
				}
			}
			else
			{
				double d = -(surface_point.dot(z_dir));
				if (d < d_z_min)
				{
					d_z_min = d;
					control_cage_face_d_z_min_ = f;
				}

				if (d > d_z_max)
				{
					d_z_max = d;
					control_cage_face_d_z_min_ = f;
				}
			}

			return true;
		});*/

		//local_x_direction_control_planes_ = std::make_tuple(x_dir, d_x_min, d_x_max);

		//local_y_direction_control_planes_ = std::make_tuple(y_dir, d_y_min, d_y_max);

		//local_z_direction_control_planes_ = std::make_tuple(z_dir, d_z_min, d_z_max);
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

			bool inside_cage = local_mvc_pt_surface(surface_point);

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

			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			const double d_x = -(surface_point.dot(std::get<0>(local_x_direction_control_planes_))),

						 d_y = -(surface_point.dot(std::get<0>(local_y_direction_control_planes_))),

						 d_z = -(surface_point.dot(std::get<0>(local_z_direction_control_planes_)));

			const bool valid_x_dir = (d_x <= std::get<2>(local_x_direction_control_planes_) &&
									  d_x >= std::get<1>(local_x_direction_control_planes_)),

					   valid_y_dir = (d_y <= std::get<2>(local_y_direction_control_planes_) &&
									  d_y >= std::get<1>(local_y_direction_control_planes_)),

					   valid_z_dir = (d_z <= std::get<2>(local_z_direction_control_planes_) &&
									  d_z >= std::get<1>(local_z_direction_control_planes_));

			if (valid_x_dir)
			{
				if (valid_y_dir && valid_z_dir)
				{
					// inside cube
					bind_mvc_inside_cube(surface_point, surface_point_idx);
				}
				else if (valid_y_dir)
				{
					if (d_z > std::get<2>(local_z_direction_control_planes_))
					{
						std::vector<Vertex> closest_face_vertices =
							incident_vertices(*control_cage_, control_cage_face_d_z_max_);

						// find plane of "fake" cube
						const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
												   std::get<1>(local_z_direction_control_planes_);

						const double new_plane_d = std::get<2>(local_z_direction_control_planes_) + gap_z_plane;

						const Vec3 shift = {0.0, 0.0, -new_plane_d};

						/*std::vector<Vec3> position_vertices;
						// add the outside vertices
						for (uint32 i = 0; i < closest_face_vertices.size(); i++)
						{
							const Vec3 position =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, closest_face_vertices[i]);
							position_vertices.push_back(position);

							const Vec3 shifted_position = position - shift;
							position_vertices.push_back(shifted_position);
						}*/

						
						//bind_mvc_customed_vertices(surface_point, surface_point_idx, triangle_set);
					}
					else
					{
					}
				}
				else if (valid_z_dir)
				{
					if (d_y > std::get<2>(local_y_direction_control_planes_))
					{
					}
					else
					{
					}
				}
			}
			else if (valid_y_dir)
			{
				if (valid_z_dir)
				{
					if (d_x > std::get<2>(local_x_direction_control_planes_))
					{
					}
					else
					{
					}
				}
				else
				{
					if (d_x > std::get<2>(local_x_direction_control_planes_))
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
						}
						else
						{
						}
					}
					else
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
						}
						else
						{
						}
					}
				}
			}
			else if (valid_z_dir)
			{
				if (d_y > std::get<2>(local_y_direction_control_planes_))
				{
					if (d_x > std::get<2>(local_x_direction_control_planes_))
					{
					}
					else
					{
					}
				}
				else
				{
					if (d_x > std::get<2>(local_x_direction_control_planes_))
					{
					}
					else
					{
					}
				}
			}
			else
			{
				if (d_x > std::get<2>(local_x_direction_control_planes_))
				{
					if (d_y > std::get<2>(local_y_direction_control_planes_))
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
						}
						else
						{
						}
					}
					else
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
						}
						else
						{
						}
					}
				}
				else
				{
					if (d_y > std::get<2>(local_y_direction_control_planes_))
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
						}
						else
						{
						}
					}
					else
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
						}
						else
						{
						}
					}
				}
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
		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		this->influence_area_->foreach_cell([&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(object, object_vertex_index, v);

			Vec3 new_pos_ = {0.0, 0.0, 0.0};

			foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, cv);

				uint32 cage_point_idx = value<uint32>(*control_cage_, control_cage_vertex_index_, cv);

				new_pos_ += control_cage_coords_(vidx, cage_point_idx) * cage_point;

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
				uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);
				control_area_validity_(surface_point_idx) = 1.0f;
			}

			return true;
		});

		this->attenuation_.resize(nbv_object);
		this->attenuation_.setZero();

		//compute_attenuation_cage(object);

	}*/

private:
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> control_cage_coords_;

	Eigen::VectorXd mixing_factor_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_comp_;

	float m_hFactor;

	Vec3 control_cage_center_;
	Vec3 control_cage_bb_min_;
	Vec3 control_cage_bb_max_;

	Vec3 influence_cage_bb_min_;
	Vec3 influence_cage_bb_max_;

	Face control_cage_face_d_x_min_;
	Face control_cage_face_d_x_max_;

	Face control_cage_face_d_y_min_;
	Face control_cage_face_d_y_max_;

	Face control_cage_face_d_z_min_;
	Face control_cage_face_d_z_max_;

	// local frame direction, mininum and maximum positions of the planes from control cage (d from ax+by+cz+d = 0)
	std::tuple<Vec3, double, double> local_x_direction_control_planes_;
	std::tuple<Vec3, double, double> local_y_direction_control_planes_;
	std::tuple<Vec3, double, double> local_z_direction_control_planes_;

	Eigen::VectorXd control_area_validity_;

	Eigen::Matrix<Vec2, Eigen::Dynamic, Eigen::Dynamic> normal_weights_;

	std::shared_ptr<Attribute<uint32>> control_cage_vertex_index_;

	std::shared_ptr<Attribute<bool>> control_cage_marked_vertices_;

	bool local_mvc_pt_surface(Vec3 pt)
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

				float mvc_value = cgogn::modeling::compute_mvc(pt, d, *influence_cage_, cage_point,
															   influence_cage_vertex_position_.get());

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

	void bind_mvc_inside_cube(const Vec3 surface_point, uint32 surface_point_idx)
	{
		DartMarker dm(*control_cage_);
		float sumMVC = 0.0;

		for (Dart d = control_cage_->begin(), end = control_cage_->end(); d != end; d = control_cage_->next(d))
		{
			Vertex cage_vertex = CMap2::Vertex(d);
			bool vc_marked = value<bool>(*control_cage_, control_cage_marked_vertices_, cage_vertex);

			if (!dm.is_marked(d) && !vc_marked)
			{
				const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, cage_vertex);
				uint32 cage_point_idx = value<uint32>(*control_cage_, control_cage_vertex_index_, cage_vertex);

				float mvc_value =
					compute_mvc(surface_point, d, *control_cage_, cage_point, control_cage_vertex_position_.get());

				control_cage_coords_(surface_point_idx, cage_point_idx) = mvc_value;

				dm.mark(d);

				value<bool>(*control_cage_, control_cage_marked_vertices_, cage_vertex) = true;

				sumMVC += mvc_value;
			}
		}

		parallel_foreach_cell(*control_cage_, [&](Vertex vc) -> bool {
			uint32 cage_point_idx2 = value<uint32>(*control_cage_, control_cage_vertex_index_, vc);

			control_cage_coords_(surface_point_idx, cage_point_idx2) =
				control_cage_coords_(surface_point_idx, cage_point_idx2) / sumMVC;

			value<bool>(*control_cage_, control_cage_marked_vertices_, vc) = false;

			return true;
		});
	}

	void bind_mvc_customed_vertices(const Vec3 surface_point, uint32 surface_point_idx,
									std::vector<Vec3> position_vertices)
	{
	}

	// https://stackoverflow.com/questions/21037241/how-to-determine-a-point-is-inside-or-outside-a-cube
	bool point_inside_cage(const Vec3& point)
	{
		bool resX = (point[0] >= control_cage_bb_min_[0]) && (point[0] <= control_cage_bb_max_[0]);
		bool resY = (point[1] >= control_cage_bb_min_[1]) && (point[1] <= control_cage_bb_max_[1]);
		bool resZ = (point[2] >= control_cage_bb_min_[2]) && (point[2] <= control_cage_bb_max_[2]);

		// simple case - axis aligned cube
		return (resX && resY && resZ);
	}

	// TODO : check property visible == dist_ci < dist_center
	// https://www.gamedev.net/forums/topic/105286-how-to-calculate-which-side-of-a-cube-is-closest-to-a-point/
	std::vector<Vec3> find_visible_cage_points(const Vec3 point)
	{
		std::shared_ptr<Attribute<Vec3>> cage_face_normal =
			cgogn::get_attribute<Vec3, Face>(*control_cage_, "face_normal");

		std::vector<Vec3> visible_points;

		const Vec3 p_to_center = control_cage_center_ - point;

		foreach_cell(*control_cage_, [&](Face fc) -> bool {
			const Vec3 face_normal = value<Vec3>(*control_cage_, cage_face_normal, fc);

			double res = face_normal.dot(p_to_center);

			if (res < 0.0)
			{
				std::vector<Vertex> face_vertices = incident_vertices(*control_cage_, fc);

				for (uint32_t i = 0; i < face_vertices.size(); i++)
				{
					const Vec3 cage_point =
						value<Vec3>(*control_cage_, control_cage_vertex_position_, face_vertices[i]);

					if (!(std::find(visible_points.begin(), visible_points.end(), cage_point) != visible_points.end()))
					{
						visible_points.push_back(cage_point);
					}
				}
			}

			return true;
		});

		// std::cout << "size visible points " << visible_points.size() << std::endl;

		return visible_points;
	}

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
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			float i_dist = this->cage_influence_distance(surface_point_idx, nbf_cage, nbv_cage);

			this->attenuation_(surface_point_idx) = (float)sin(0.5*M_PI * (i_dist ));
			/*if (control_area_validity_(surface_point_idx) == 1.0f)
			{

				if (i_dist > h)
				{
					h = i_dist;
				}

				this->attenuation_(surface_const double point_idx) = 1.0f;
			}
			else
			{

				if (i_dist > max_dist){
					max_dist = i_dist;
				}
				attenuation_points.push_back({surface_point_idx, i_dist});
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
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_

/*

void bind_mvc_old(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*control_cage_, "vertex_index");

		std::shared_ptr<Attribute<bool>> cage_marked_vertices =
			get_attribute<bool, Vertex>(*control_cage_, "marked_vertices");

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		control_cage_coords_.resize(nbv_object, nbv_cage);
		control_cage_coords_.setZero();

		this->influence_area_->foreach_cell([&](Vertex v) {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);

			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			if (point_inside_cage(surface_point))
			{
				DartMarker dm(*control_cage_);
				float sumMVC = 0.0;

				for (Dart d = control_cage_->begin(), end = control_cage_->end(); d != end; d = control_cage_->next(d))
				{
					Vertex cage_vertex = CMap2::Vertex(d);
					bool vc_marked = value<bool>(*control_cage_, cage_marked_vertices, cage_vertex);

					if (!dm.is_marked(d) && !vc_marked)
					{
						const Vec3& cage_point =
							value<Vec3>(*control_cage_, control_cage_vertex_position_, cage_vertex);
						uint32 cage_point_idx = value<uint32>(*control_cage_, cage_vertex_index, cage_vertex);

						float mvc_value = compute_mvc(surface_point, d, *control_cage_, cage_point,
													  control_cage_vertex_position_.get());

						control_cage_coords_(surface_point_idx, cage_point_idx) = mvc_value;

						dm.mark(d);

						value<bool>(*control_cage_, cage_marked_vertices, cage_vertex) = true;

						sumMVC += mvc_value;
					}
				}

				parallel_foreach_cell(*control_cage_, [&](Vertex vc) -> bool {
					uint32 cage_point_idx2 = value<uint32>(*control_cage_, cage_vertex_index, vc);

					control_cage_coords_(surface_point_idx, cage_point_idx2) = control_cage_coords_(surface_point_idx,
cage_point_idx2) / sumMVC;

					value<bool>(*control_cage_, cage_marked_vertices, vc) = false;

					return true;
				});
			}
			else
			{
				const std::vector<Vec3> control_cage_points = find_visible_cage_points(surface_point);
			}
		});
	}

	/*void bind_mvc(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, Vertex>(object, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*influence_cage_, "vertex_index");

		std::shared_ptr<Attribute<bool>> cage_marked_vertices =
			get_attribute<bool, Vertex>(*influence_cage_, "marked_vertices");

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*influence_cage_);

		coords_.resize(nbv_object, nbv_cage);
		coords_.setZero();

		influence_area_->foreach_cell([&](Vertex v) {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			DartMarker dm(*influence_cage_);
			float sumMVC = 0.0;

			for (Dart d = influence_cage_->begin(), end = influence_cage_->end(); d != end;
				 d = influence_cage_->next(d))
			{
				Vertex cage_vertex = CMap2::Vertex(d);
				bool vc_marked = value<bool>(*influence_cage_, cage_marked_vertices, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{
					const Vec3& cage_point =
						value<Vec3>(*influence_cage_, influence_cage_vertex_position_, cage_vertex);
					uint32 cage_point_idx = value<uint32>(*influence_cage_, cage_vertex_index, cage_vertex);

					float mvc_value = compute_mvc(surface_point, d, *influence_cage_, cage_point,
												  influence_cage_vertex_position_.get());

					coords_(surface_point_idx, cage_point_idx) = mvc_value;

					dm.mark(d);

					value<bool>(*influence_cage_, cage_marked_vertices, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			parallel_foreach_cell(*influence_cage_, [&](Vertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(*influence_cage_, cage_vertex_index, vc);

				coords_(surface_point_idx, cage_point_idx2) = coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				value<bool>(*influence_cage_, cage_marked_vertices, vc) = false;

				return true;
			});
		});

	}*/
