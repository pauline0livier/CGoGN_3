/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_CAGE_DEFORMATION_H_
#define CGOGN_MODULE_CAGE_DEFORMATION_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/graph_selection.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/geometry/ui_modules/surface_selectionPO.h>
#include <cgogn/modeling/ui_modules/graph_deformation.h>
#include <cgogn/modeling/ui_modules/surface_deformation.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>

#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/angle.h>

#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

#include <cgogn/modeling/types/axis_deformation_tool.h>
#include <cgogn/modeling/types/cage_deformation_tool.h>
#include <cgogn/modeling/types/handle_deformation_tool.h>

#include <boost/synapse/connect.hpp>

#include <iostream>

#include <string>

namespace cgogn
{

namespace ui
{

using Vec2 = geometry::Vec2;
using Vec3 = geometry::Vec3;
using Mat3 = geometry::Mat3;
using Scalar = geometry::Scalar;

template <typename MESH, typename GRAPH>

class SpaceDeformation : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "SpaceDeformation can only be used with meshes of dimension 2");

	template <typename T>
	using GraphAttribute = typename mesh_traits<GRAPH>::template Attribute<T>;

	template <typename T>
	using MeshAttribute = typename mesh_traits<MESH>::template Attribute<T>;

	using GraphVertex = typename mesh_traits<GRAPH>::Vertex;
	using GraphEdge = typename mesh_traits<GRAPH>::Edge;

	using MeshVertex = typename mesh_traits<MESH>::Vertex;
	using MeshEdge = typename mesh_traits<MESH>::Edge;
	using MeshFace = typename mesh_traits<MESH>::Face;

	enum SelectionTool
	{
		Handle = 0,
		Axis,
		Cage
	};

	struct Parameters
	{
		Parameters() : vertex_position_(nullptr), new_cage_(false), nb_cage(0), nb_tool_(0)
		{
		}

		~Parameters()
		{
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<MeshAttribute<Vec3>> vertex_position_;

		std::vector<MESH*> list_cage_;

		std::vector<Vec3> init_position_;

		std::shared_ptr<MeshAttribute<Vec3>> gammaColor;

		std::vector<Scalar> current_pos_factor_; 
		std::vector<std::vector<Scalar>> weight_matrix_;
		std::vector<Vec3> deformation_vector_; 
		std::vector<Scalar> init_pos_factor_;  

		/*Eigen::VectorXd current_pos_factor_; 
		Eigen::MatrixXd weight_matrix_; 
		Eigen::MatrixXd deformation_matrix_; 
		Eigen::VectorXd init_pos_factor_; 
		Eigen::MatrixXd init_position_; */

		bool new_cage_;

		std::shared_ptr<boost::synapse::connection> cells_set_connection_;

		int nb_cage;

		int nb_tool_; 
	};

public:
	std::unordered_map<std::string, std::shared_ptr<cgogn::modeling::AxisDeformationTool<MESH>>> axis_container_;
	std::unordered_map<std::string, std::shared_ptr<cgogn::modeling::HandleDeformationTool<MESH>>> handle_container_;
	std::unordered_map<std::string, std::shared_ptr<cgogn::modeling::CageDeformationTool<MESH>>> cage_container_;

	SpaceDeformation(const App& app)
		: Module(app, "SpaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selected_cage_(nullptr), selected_handle_(nullptr), influence_cage_mode_(true)
	{
	}
	~SpaceDeformation()
	{
	}

	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.nb_tool_ = 100; 
	}

	void set_selected_mesh(MESH& m)
	{
		selected_mesh_ = &m;
	}

	void set_vertex_position(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		p.vertex_position_ = vertex_position;

		uint32 nbv_m = nb_cells<MeshVertex>(m);
		p.init_position_.resize(nbv_m);

		std::shared_ptr<MeshAttribute<uint32>> vertex_index = get_attribute<uint32, MeshVertex>(m, "vertex_index");

		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(m, vertex_position, v);
			uint32 surface_point_idx = value<uint32>(m, vertex_index, v);

			p.init_position_[surface_point_idx] = surface_point;
			return true;
		});

		p.current_pos_factor_.resize(nbv_m);
		std::fill(p.current_pos_factor_.begin(), p.current_pos_factor_.end(), 0.0);
		
		p.init_pos_factor_.resize(nbv_m);
		std::fill(p.init_pos_factor_.begin(), p.init_pos_factor_.end(), 1.0);

		p.deformation_vector_.resize(p.nb_tool_);
		for (uint32_t i = 0; i < p.deformation_vector_.size(); i++){
			p.deformation_vector_[i] = {0.0, 0.0, 0.0}; 
		}
		
		p.weight_matrix_ = std::vector<std::vector<Scalar>>(nbv_m, std::vector<Scalar>(p.nb_tool_));
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 surface_point_idx = value<uint32>(m, vertex_index, v);

			std::fill(std::begin(p.weight_matrix_[surface_point_idx]), std::end(p.weight_matrix_[surface_point_idx]), 0);;
			return true;
		});


		//std::fill(vect.begin(), vect.end(), value);

	}

	template <typename FUNC>
	void foreach_cage(const FUNC& f)
	{
		static_assert(is_ith_func_parameter_same<FUNC, 0, MESH&>::value, "Wrong function parameter type");
		static_assert(is_ith_func_parameter_same<FUNC, 1, const std::string&>::value, "Wrong function parameter type");
		for (auto& [name, cdt] : cage_container_)
			f(*(cdt->control_cage_), name);
	}

	template <typename FUNC>
	void foreach_handle(const FUNC& f)
	{
		// static_assert(is_ith_func_parameter_same<FUNC, 0, GRAPH&>::value, "Wrong function parameter type");
		// static_assert(is_ith_func_parameter_same<FUNC, 1, const std::string&>::value, "Wrong function parameter
		// type");
		for (auto& [name, hdt] : handle_container_)
			f(*(hdt->control_handle_), name);
	}

	/////////////
	// SIGNALS //
	/////////////
	using cage_added = struct cage_added_ (*)(cgogn::modeling::CageDeformationTool<MESH>* c);
	using handle_added = struct handle_added_ (*)(cgogn::modeling::HandleDeformationTool<MESH>* h);
	using axis_added = struct axis_added_ (*)(cgogn::modeling::AxisDeformationTool<MESH>* a);

private:
	// creation deformation tools
	MESH* generate_global_cage(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position, int i)
	{

		Parameters& p = parameters_[&m];

		MESH* cage = mesh_provider_->add_mesh("cage" + std::to_string(i)); // (m_name + "_cage");
		std::shared_ptr<MeshAttribute<Vec3>> cage_vertex_position =
			cgogn::add_attribute<Vec3, MeshVertex>(*cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*cage, cage_vertex_position);

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		cgogn::modeling::create_bounding_box(*cage, cage_vertex_position.get(), md.bb_min_, md.bb_max_);

		mesh_provider_->emit_connectivity_changed(*cage);
		mesh_provider_->emit_attribute_changed(*cage, cage_vertex_position.get());

		std::shared_ptr<MeshAttribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, MeshVertex>(*cage, "vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*cage, vertex_index.get());

		p.list_cage_.push_back(cage);
		// p.list_weights_.push_back(Weights());
		// cd.cage_vertex_position_ = cage_vertex_position;

		ui::View* v1 = app_.current_view();

		surface_render_->set_vertex_position(*v1, *cage, cage_vertex_position);
		surface_render_->set_render_faces(*v1, *cage, false);

		return cage;
	}

	MESH* generate_local_cage(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position,
							  CellsSet<MESH, MeshVertex>* control_set)
	{
		int cage_number = cage_container_.size();
		std::string cage_name = "local_cage" + std::to_string(cage_number);
		MESH* l_cage = mesh_provider_->add_mesh(cage_name);

		std::shared_ptr<MeshAttribute<Vec3>> l_cage_vertex_position =
			cgogn::add_attribute<Vec3, MeshVertex>(*l_cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*l_cage, l_cage_vertex_position);

		auto mesh_vertex_normal = get_attribute<Vec3, MeshVertex>(m, "normal");

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 local_min = {1000.0, 1000.0, 1000.0};
		Vec3 local_max = {0.0, 0.0, 0.0};

		Vec3 center = {0.0, 0.0, 0.0};
		Vec3 normal = {0.0, 0.0, 0.0};

		control_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			const Vec3& norm = value<Vec3>(m, mesh_vertex_normal, v);

			center += pos;
			normal += norm;

			for (size_t j = 0; j < 3; j++)
			{
				if (pos[j] < local_min[j])
				{
					local_min[j] = pos[j];
				}

				if (pos[j] > local_max[j])
				{
					local_max[j] = pos[j];
				}
			}
		});

		center /= control_set->size();
		normal /= control_set->size();

		Vec3 ray = normal;
		ray.normalize();

		const auto [it, inserted] =
			cage_container_.emplace(cage_name, std::make_shared<cgogn::modeling::CageDeformationTool<MESH>>());
		cgogn::modeling::CageDeformationTool<MESH>* cdt = it->second.get();

		if (inserted)
		{
			cdt->create_space_tool(l_cage, l_cage_vertex_position.get(), local_min, local_max, center, ray);

			mesh_provider_->emit_connectivity_changed(*l_cage);
			mesh_provider_->emit_attribute_changed(*l_cage, l_cage_vertex_position.get());

			View* v1 = app_.current_view();

			surface_render_->set_vertex_position(*v1, *l_cage, l_cage_vertex_position);
			surface_render_->set_render_faces(*v1, *l_cage, false);

			surface_selection_->set_vertex_position(*l_cage, l_cage_vertex_position);

			std::shared_ptr<MeshAttribute<Vec3>> l_cage_vertex_normal =
				cgogn::add_attribute<Vec3, MeshVertex>(*l_cage, "normal");

			surface_diff_pptes_->compute_normal(*l_cage, l_cage_vertex_position.get(), l_cage_vertex_normal.get());

			MESH* i_cage = mesh_provider_->add_mesh("i_cage" + cage_number);

			std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
				cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "position");

			mesh_provider_->set_mesh_bb_vertex_position(*i_cage, i_cage_vertex_position);

			Vec3 i_center = (local_min + local_max) / Scalar(2);
			Vec3 bb_min_ = ((local_min - i_center) * 2.f) + i_center;
			Vec3 bb_max_ = ((local_max - i_center) * 2.f) + i_center;

			cdt->set_center_control_cage(i_center);
			cdt->set_influence_cage(i_cage, i_cage_vertex_position.get(), bb_min_, bb_max_);

			mesh_provider_->emit_connectivity_changed(*i_cage);
			mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

			surface_render_->set_vertex_position(*v1, *i_cage, i_cage_vertex_position);
			surface_render_->set_render_faces(*v1, *i_cage, false);

			CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();
			cdt->influence_area_ = &i_set;

			//  cd.smoothing_factor_ = 0.3;

			// computeCagesAttenuationFactor(*l_cage);
			cdt->update_influence_area(m, vertex_position.get());

			mesh_provider_->emit_cells_set_changed(m, cdt->influence_area_);

			boost::synapse::emit<cage_added>(this, cdt);
		}

		return l_cage;
	}

	std::shared_ptr<cgogn::modeling::HandleDeformationTool<MESH>> generate_handle_tool(MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position,
						   CellsSet<MESH, MeshVertex>* handle_set)
	{
		int handle_number = handle_container_.size();
		std::string handle_name = "local_handle" + std::to_string(handle_number);
		GRAPH* handle = graph_provider_->add_mesh(handle_name);

		auto handle_vertex_position = add_attribute<Vec3, GraphVertex>(*handle, "position");
		auto handle_vertex_radius = add_attribute<Scalar, GraphVertex>(*handle, "radius");

		auto mesh_vertex_normal = get_attribute<Vec3, MeshVertex>(m, "normal");

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 local_min = {1000.f, 1000.f, 1000.f};
		Vec3 local_max = {0.f, 0.f, 0.f};

		Vec3 center = {0.0, 0.0, 0.0};
		Vec3 normal = {0.0, 0.0, 0.0};

		handle_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			const Vec3& norm = value<Vec3>(m, mesh_vertex_normal, v);

			center += pos;
			normal += norm;

			for (size_t j = 0; j < 3; j++)
			{
				if (pos[j] < local_min[j])
				{
					local_min[j] = pos[j];
				}

				if (pos[j] > local_max[j])
				{
					local_max[j] = pos[j];
				}
			}
		});

		center /= handle_set->size();
		normal /= handle_set->size();

		Vec3 ray = normal;
		ray.normalize();

		// const Vec3 handle_position = center;
		const Vec3 handle_position = center; //{center[0] + 0.1f * ray[0], center[1] + 0.1f * ray[1], center[2] + 0.1f * ray[2]};

		const Vec3 inner_handle_position = {center[0] - 2.f * ray[0], center[1] - 2.f * ray[1],
											center[2] - 2.f * ray[2]};

		 
		MeshVertex closest_vertex; 
		double min_dist = 1000000;
		handle_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			
			double dist = (pos - handle_position).squaredNorm(); 
			if (dist < min_dist){
				min_dist = dist; 
				closest_vertex = v; 
			}

		});


		const auto [it, inserted] =
			handle_container_.emplace(handle_name, std::make_shared<cgogn::modeling::HandleDeformationTool<MESH>>());
		cgogn::modeling::HandleDeformationTool<MESH>* hdt = it->second.get();

		if (inserted)
		{

			hdt->create_space_tool(handle, handle_vertex_position.get(), handle_vertex_radius.get(), handle_position,
								   inner_handle_position, ray);

			hdt -> set_handle_mesh_vertex(closest_vertex); 

			graph_provider_->emit_connectivity_changed(*handle);
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_position.get());
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_radius.get());

			View* v1 = app_.current_view();
			graph_render_->set_vertex_position(*v1, *handle, handle_vertex_position);

			graph_selection_->set_vertex_position(*handle, handle_vertex_position);

			graph_deformation_->set_vertex_position(*handle, handle_vertex_position);

			graph_deformation_->set_displacement_normal(*handle, ray);

			auto object_geodesic = get_or_add_attribute<Scalar, MeshVertex>(*selected_mesh_, "geodesic_distance");
			
			hdt->set_geodesic_distance(m, vertex_position);
			mesh_provider_->emit_attribute_changed(m, object_geodesic.get());
			

			boost::synapse::emit<handle_added>(this, hdt);
		}

		return handle_container_[handle_name];
	}

	/*GRAPH* generate_handle(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position,
						   CellsSet<MESH, MeshVertex>* control_set)
	{
		int handle_number = handle_container_.size();
		std::string handle_name = "local_handle" + std::to_string(handle_number);
		GRAPH* handle = graph_provider_->add_mesh(handle_name);

		auto handle_vertex_position = add_attribute<Vec3, GraphVertex>(*handle, "position");
		auto handle_vertex_radius = add_attribute<Scalar, GraphVertex>(*handle, "radius");

		auto mesh_vertex_normal = get_attribute<Vec3, MeshVertex>(m, "normal");

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 local_min = {1000.f, 1000.f, 1000.f};
		Vec3 local_max = {0.f, 0.f, 0.f};

		Vec3 center = {0.0, 0.0, 0.0};
		Vec3 normal = {0.0, 0.0, 0.0};

		control_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			const Vec3& norm = value<Vec3>(m, mesh_vertex_normal, v);

			center += pos;
			normal += norm;

			for (size_t j = 0; j < 3; j++)
			{
				if (pos[j] < local_min[j])
				{
					local_min[j] = pos[j];
				}

				if (pos[j] > local_max[j])
				{
					local_max[j] = pos[j];
				}
			}
		});

		center /= control_set->size();
		normal /= control_set->size();

		Vec3 ray = normal;
		ray.normalize();

		// const Vec3 handle_position = center;
		const Vec3 handle_position = {center[0] + 0.1f * ray[0], center[1] + 0.1f * ray[1], center[2] + 0.1f * ray[2]};

		const Vec3 inner_handle_position = {center[0] - 2.f * ray[0], center[1] - 2.f * ray[1],
											center[2] - 2.f * ray[2]};

		const auto [it, inserted] =
			handle_container_.emplace(handle_name, std::make_shared<cgogn::modeling::HandleDeformationTool<MESH>>());
		cgogn::modeling::HandleDeformationTool<MESH>* hdt = it->second.get();

		if (inserted)
		{

			hdt->create_space_tool(handle, handle_vertex_position.get(), handle_vertex_radius.get(), handle_position,
								   inner_handle_position);

			graph_provider_->emit_connectivity_changed(*handle);
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_position.get());
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_radius.get());

			View* v1 = app_.current_view();
			graph_render_->set_vertex_position(*v1, *handle, handle_vertex_position);

			graph_selection_->set_vertex_position(*handle, handle_vertex_position);

			MESH* i_cage = mesh_provider_->add_mesh("i_cage" + handle_number);

			std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
				cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "position");

			std::shared_ptr<MeshAttribute<Vec3>> i_cage_local_vertex_position =
				cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "local_position");

			mesh_provider_->set_mesh_bb_vertex_position(*i_cage, i_cage_vertex_position);

			Vec3 i_center = (local_min + local_max) / Scalar(2);
			Vec3 bb_min_ = ((local_min - i_center) * 1.5f) + i_center;
			Vec3 bb_max_ = ((local_max - i_center) * 1.5f) + i_center;

			//Vec3 bb_min_ = ((local_min - i_center) * 2.5f) + i_center;
			//Vec3 bb_max_ = ((local_max - i_center) * 2.5f) + i_center;

			hdt->set_influence_cage_handle(i_cage, i_cage_vertex_position.get(), i_cage_local_vertex_position.get(),
										   bb_min_, bb_max_, handle_position, ray);

			mesh_provider_->emit_connectivity_changed(*i_cage);
			mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

			surface_render_->set_vertex_position(*v1, *i_cage, i_cage_vertex_position);
			surface_render_->set_render_faces(*v1, *i_cage, false);

			CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();
			hdt->influence_area_ = &i_set;

			hdt->update_influence_area(m, vertex_position.get());

			mesh_provider_->emit_cells_set_changed(m, hdt->influence_area_);

			boost::synapse::emit<handle_added>(this, hdt);
		}

		return handle;
	}*/

	GRAPH* generate_axis(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position,
						 CellsSet<MESH, MeshVertex>* control_set)
	{
		int axis_number = axis_container_.size();
		std::string axis_name = "axis" + std::to_string(axis_number);
		GRAPH* axis = graph_provider_->add_mesh(axis_name);

		auto axis_vertex_position = add_attribute<Vec3, GraphVertex>(*axis, "position");
		auto axis_vertex_radius = add_attribute<Scalar, GraphVertex>(*axis, "radius");

		auto mesh_vertex_normal = get_attribute<Vec3, MeshVertex>(m, "normal");

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 local_min = {1000.f, 1000.f, 1000.f};
		Vec3 local_max = {0.f, 0.f, 0.f};

		Vec3 center = {0.0, 0.0, 0.0};
		Vec3 normal = {0.0, 0.0, 0.0};

		control_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			const Vec3& norm = value<Vec3>(m, mesh_vertex_normal, v);

			center += pos;
			normal += norm;

			for (size_t j = 0; j < 3; j++)
			{
				if (pos[j] < local_min[j])
				{
					local_min[j] = pos[j];
				}

				if (pos[j] > local_max[j])
				{
					local_max[j] = pos[j];
				}
			}
		});

		center /= control_set->size();
		normal /= control_set->size();

		Vec3 ray = normal;
		ray.normalize();

		Eigen::Matrix3f covariance_matrix; 
		// inspired from https://gist.github.com/atandrau/847214/882418ab34737699a6b1394d3a28c66e2cc0856f
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++) {
				covariance_matrix(i,j) = 0.0;
				control_set->foreach_cell([&](MeshVertex v) {
					const Vec3& pos = value<Vec3>(m, vertex_position, v);
					covariance_matrix(i,j) += (center[i] - pos[i]) * 
												(center[j] - pos[j]);
				}); 
				
			covariance_matrix(i,j) /= control_set->size() -1;
		}	

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> eigen_solver(covariance_matrix);
		Eigen::Matrix<float, 1, Eigen::Dynamic> eigen_values = eigen_solver.eigenvalues();
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> eigen_vectors = eigen_solver.eigenvectors();

		Eigen::Vector3f main_eigen_vector = cgogn::modeling::sort_eigen_vectors(eigen_values, eigen_vectors); 
		main_eigen_vector.normalize();

		Vec3 main_direction = {main_eigen_vector[0], main_eigen_vector[1], main_eigen_vector[2]}; 

		// loop on control set to find extremum points and height 
		double min_x = 1000, max_x = 0.0, max_normal = center.dot(ray); 
		control_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);

			double x_value = pos.dot(main_direction); 
			double n_value = pos.dot(ray); 

			if (x_value > max_x){
				max_x = x_value; 
			} 

			if (x_value < min_x){
				min_x = x_value; 
			}

			if (n_value > max_normal){
				max_normal = n_value; 
			}
		}); 

		const Vec3 axis_center = center + 2.0*(max_normal - center.dot(ray))*ray;
		const Vec3 extrem_min = axis_center + (min_x - center.dot(main_direction))*main_direction; 
		const Vec3 extrem_max = axis_center + (max_x - center.dot(main_direction))*main_direction;

		const double width = (extrem_max - extrem_min).norm()/3.0; 

		std::vector<Vec3> axis_vertices;
		axis_vertices.push_back(extrem_min); 
		axis_vertices.push_back(axis_center);
		axis_vertices.push_back(extrem_max); 

		const auto [it, inserted] =
			axis_container_.emplace(axis_name, std::make_shared<cgogn::modeling::AxisDeformationTool<MESH>>());
		cgogn::modeling::AxisDeformationTool<MESH>* adt = it->second.get();
		if (inserted)
		{

			adt->create_space_tool(axis, axis_vertex_position.get(), axis_vertex_radius.get(), axis_vertices);

			graph_provider_->emit_connectivity_changed(*axis);
			graph_provider_->emit_attribute_changed(*axis, axis_vertex_position.get());
			graph_provider_->emit_attribute_changed(*axis, axis_vertex_radius.get());

			View* v1 = app_.current_view();

			graph_render_->set_vertex_position(*v1, *axis, axis_vertex_position);

			graph_selection_->set_vertex_position(*axis, axis_vertex_position);

			MESH* i_cage = mesh_provider_->add_mesh("i_cage" + axis_number);

			std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
				cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "position");

			std::shared_ptr<MeshAttribute<Vec3>> i_cage_local_vertex_position =
				cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "local_position");

			std::shared_ptr<MeshAttribute<uint32>> i_cage_local_skeleton =
				cgogn::add_attribute<uint32, MeshVertex>(*i_cage, "local_skeleton");

			mesh_provider_->set_mesh_bb_vertex_position(*i_cage, i_cage_vertex_position);

			Vec3 i_center = (local_min + local_max) / Scalar(2);
			Vec3 bb_min_ = ((local_min - i_center) * 1.5f) + i_center;
			Vec3 bb_max_ = ((local_max - i_center) * 1.5f) + i_center;

			adt->set_influence_cage_axis(i_cage, i_cage_vertex_position.get(), i_cage_local_vertex_position.get(), i_cage_local_skeleton.get(), bb_min_, bb_max_, main_direction, ray, width);

			mesh_provider_->emit_connectivity_changed(*i_cage);
			mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

			surface_render_->set_vertex_position(*v1, *i_cage, i_cage_vertex_position);
			surface_render_->set_render_faces(*v1, *i_cage, false);

			CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();
			adt->influence_area_ = &i_set;

			adt->update_influence_area(m, vertex_position.get());

			mesh_provider_->emit_cells_set_changed(m, adt->influence_area_);

			boost::synapse::emit<axis_added>(this, adt);
		}

		return axis;
	}


/// BIND 

	void bind_influence_cage_mvc(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
								 MESH& control_cage, const std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position)

	{
		selected_cdt_->bind_mvc(object, object_vertex_position);
		selected_cdt_->set_up_attenuation(object, object_vertex_position);

		displayGammaColor(object);

		Parameters& p = parameters_[&object];

		selected_cdt_->cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				&control_cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cage_vertex_position.get() == attribute)
					{

						selected_cdt_->update_influence_cage_position();

						MESH* i_cage = selected_cdt_->influence_cage_;

						std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
							get_attribute<Vec3, MeshVertex>(*i_cage, "position");

						mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

						selected_cdt_->update_deformation_object_mvc(object, object_vertex_position, p.init_position_);

						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					}
				});
	}

	void bind_influence_cage_mvc_handle(MESH& object,
										const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
										GRAPH& control_handle,
										const std::shared_ptr<GraphAttribute<Vec3>>& handle_vertex_position)

	{
		selected_hdt_->bind_mvc(object, object_vertex_position);
		selected_hdt_->set_up_attenuation(object, object_vertex_position);

		displayGammaColor(object);

		Parameters& p = parameters_[&object];

		selected_hdt_->handle_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<GRAPH>::template attribute_changed_t<Vec3>>(
				&control_handle, [&](GraphAttribute<Vec3>* attribute) {
					if (handle_vertex_position.get() == attribute)
					{
							selected_hdt_->update_influence_cage_position();

							MESH* i_cage = selected_hdt_->influence_cage_;

							std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
								get_attribute<Vec3, MeshVertex>(*i_cage, "position");

							mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

							selected_hdt_->update_deformation_object_mvc(object, object_vertex_position,
							p.init_position_);

							mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					} else {
						std::cout << " no moving " << std::endl; 
					}
				});
	}

	void bind_handle_influence_area(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position, CellsSet<MESH, MeshVertex>* influence_set, GRAPH& control_handle,
	const std::shared_ptr<GraphAttribute<Vec3>>& handle_vertex_position){

		MeshData<MESH>& md = mesh_provider_->mesh_data(object);
		CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();
		selected_hdt_->influence_area_ = &i_set;

		selected_hdt_->set_influence_area(object, object_vertex_position, influence_set);

		mesh_provider_->emit_cells_set_changed(object, selected_hdt_->influence_area_);

		displayGammaColor(*selected_mesh_);

		Parameters& p = parameters_[&object];

		selected_hdt_->handle_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<GRAPH>::template attribute_changed_t<Vec3>>(
				&control_handle, [&](GraphAttribute<Vec3>* attribute) {
					if (handle_vertex_position.get() == attribute)
					{

						selected_hdt_->update_deformation_object(object, object_vertex_position, p.init_position_);

						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());

					} else {
						std::cout << " no moving " << std::endl; 
					}
				});


	}

	void updateDeformationObject(MESH& object)
	{
		
	}

	void displayGammaColor(MESH& object)
	{

		std::shared_ptr<MeshAttribute<Vec3>> gamma_color =
			cgogn::get_or_add_attribute<Vec3, MeshVertex>(object, "color_gamma");
		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		Parameters& p = parameters_[&object];

		foreach_cell(object, [&](MeshVertex v) -> bool {
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			double gamma_value = 0.0;

			/*for (auto& [name, cdt] : cage_container_)
			{
				gamma_value += cdt->attenuation_(surface_point_idx);
			}*/

			for (auto& [name, hdt] : handle_container_)
			{
				gamma_value += hdt->attenuation_(surface_point_idx);
			
			}

			/*for (auto& [name, adt] : axis_container_)
			{
				gamma_value += adt->attenuation_(surface_point_idx);
			}*/

			Vec3 color;
			if (gamma_value <= 0.0)
				color = {0.0, 0.0, 0.0};
			else
			{
				gamma_value = (gamma_value > 1.0) ? 1.0 : gamma_value;
				color[0] = (1.0 + std::sin((gamma_value - 0.5) * M_PI)) / 2.0;
				color[1] = std::sin((gamma_value + 2.0) * M_PI);
				color[2] = (1.0 + std::sin((gamma_value + 0.5) * M_PI)) / 2.0;
			}

			value<Vec3>(object, gamma_color, v) = color;

			return true;
		});

		View* v1 = app_.current_view();
		surface_render_->set_vertex_color(*v1, object, gamma_color);
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &SpaceDeformation<MESH, GRAPH>::init_mesh));

		graph_provider_ = static_cast<ui::MeshProvider<GRAPH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<GRAPH>::name} + ")"));

		surface_render_ = static_cast<ui::SurfaceRender<MESH>*>(
			app_.module("SurfaceRender (" + std::string{mesh_traits<MESH>::name} + ")"));

		graph_render_ = static_cast<ui::SurfaceRender<GRAPH>*>(
			app_.module("SurfaceRender (" + std::string{mesh_traits<GRAPH>::name} + ")"));

		surface_diff_pptes_ = static_cast<ui::SurfaceDifferentialProperties<MESH>*>(
			app_.module("SurfaceDifferentialProperties (" + std::string{mesh_traits<MESH>::name} + ")"));

		surface_selection_ = static_cast<ui::SurfaceSelectionPO<MESH>*>(
			app_.module("SurfaceSelectionPO (" + std::string{mesh_traits<MESH>::name} + ")"));

		graph_selection_ = static_cast<ui::GraphSelection<GRAPH>*>(
			app_.module("GraphSelection (" + std::string{mesh_traits<GRAPH>::name} + ")"));

		surface_deformation_ = static_cast<ui::SurfaceDeformation<MESH>*>(
			app_.module("SurfaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));

		graph_deformation_ = static_cast<ui::GraphDeformation<GRAPH>*>(
			app_.module("GraphDeformation (" + std::string{mesh_traits<GRAPH>::name} + ")"));
	}

	void left_panel() override
	{
		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Object", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);
			Parameters& p = parameters_[selected_mesh_];

			imgui_combo_attribute<MeshVertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
													[&](const std::shared_ptr<MeshAttribute<Vec3>>& attribute) {
														set_vertex_position(*selected_mesh_, attribute);
													});
			if (p.vertex_position_)
			{
				ImGui::Separator();
				ImGui::Text("Global");
				if (ImGui::Button("Generate global cage"))
				{
					generate_global_cage(*selected_mesh_, p.vertex_position_, p.nb_cage);

					p.nb_cage++;
				}

				ImGui::Separator();
				ImGui::Text("Local");

				ImGui::RadioButton("Handle", reinterpret_cast<int*>(&selected_tool_), Handle);
				ImGui::SameLine();
				ImGui::RadioButton("Axis", reinterpret_cast<int*>(&selected_tool_), Axis);
				ImGui::SameLine();
				ImGui::RadioButton("Cage", reinterpret_cast<int*>(&selected_tool_), Cage);

				if (selected_tool_ == Cage)
				{
					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					imgui_combo_cells_set(md, control_set, "Control ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });

					bool newCage = false;

					if (control_set && control_set->size() > 0)
					{
						MESH* l_cage = generate_local_cage(*selected_mesh_, p.vertex_position_, control_set);
					}

					ImGui::Separator();
					ImGui::Text("Binding");

					if (cage_container_.size() > 0)
					{

						SpaceDeformation<MESH, GRAPH>* space_deformation = static_cast<SpaceDeformation<MESH, GRAPH>*>(
							app_.module("SpaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));
						if (imgui_cage_selector(space_deformation, selected_cage_, "Cage",
												 [&](MESH& m) { selected_cage_ = &m; }))
						{
							
							const std::string cage_name = mesh_provider_->mesh_name(*selected_cage_);

							std::string prefix = "local_cage";

							if (cage_name.size() > 0 && std::mismatch(prefix.begin(), prefix.end(), cage_name.begin()).first == prefix.end())
							{

								// inspired from https://github.com/ocornut/imgui/issues/1658
								const char* items[] = {"MVC", "Green"};
								std::string current_item = "MVC";
								ImGuiComboFlags flags = ImGuiComboFlags_NoArrowButton;

								ImGuiStyle& style = ImGui::GetStyle();
								float w = ImGui::CalcItemWidth();
								float spacing = style.ItemInnerSpacing.x;
								float button_sz = ImGui::GetFrameHeight();
								ImGui::PushItemWidth(w - spacing * 2.0f - button_sz * 2.0f);
								if (ImGui::BeginCombo("##custom combo", current_item.c_str(),
													  ImGuiComboFlags_NoArrowButton))
								{
									for (int n = 0; n < IM_ARRAYSIZE(items); n++)
									{
										bool is_selected = (current_item == items[n]);
										if (ImGui::Selectable(items[n], is_selected))
											current_item = items[n];
										if (is_selected)
											ImGui::SetItemDefaultFocus();
									}
									ImGui::EndCombo();
								}

								std::shared_ptr<cgogn::modeling::CageDeformationTool<MESH>> cdt =
								 cage_container_[cage_name];
								selected_cdt_ = cage_container_[cage_name];
								if (ImGui::Button("Bind object"))
								{

									if (current_item == "MVC")
									{
										bind_influence_cage_mvc(*selected_mesh_, p.vertex_position_, *selected_cage_,selected_cdt_->control_cage_vertex_position_); 
									
									}

									else if (current_item == "Green")
									{
										/*if (!cd.local_def)
										{
											bind_object_green(*selected_mesh_, p.vertex_position_, *selected_cage_,
															  cd.cage_vertex_position_);
										}
										else
										{
											bind_local_green(*selected_mesh_, p.vertex_position_, *selected_cage_,
															 cd.cage_vertex_position_);
										}*/

										std::cout << "green" << std::endl;
									}
									 else
									{
										std::cout << "not available yet" << std::endl;
									}
								}

								/*ImGui::Separator();
									float new_smoothing_factor = selected_cdt_->smoothing_factor_;
									ImGui::SliderFloat("Attenuation factor", &new_smoothing_factor, 0.0f, 1.0f);

									if (new_smoothing_factor != selected_cdt_->smoothing_factor_){
										selected_cdt_->smoothing_factor_ = new_smoothing_factor;
										selected_cdt_->update_attenuation(*selected_mesh_);
									}*/
							}
						}
					}
				}

				if (selected_tool_ == Handle)
				{

					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					imgui_combo_cells_set(md, control_set, "Handle ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });

					if (control_set && control_set -> size() > 0 ){
						selected_hdt_ = generate_handle_tool(*selected_mesh_, p.vertex_position_, control_set);

						control_set = nullptr; 
					}			  

					ImGui::Separator();
					ImGui::Text("Binding");

					CellsSet<MESH, MeshVertex>* influence_set = nullptr;
					
					imgui_combo_cells_set(md, influence_set, "Influence ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { influence_set = cs; });

					if (selected_hdt_ && influence_set)
					{
						bind_handle_influence_area(*selected_mesh_, p.vertex_position_, influence_set, *selected_hdt_->control_handle_,
							selected_hdt_->control_handle_vertex_position_); 
						
						influence_set = nullptr; 

						surface_selection_->clear_selected_vertices_set(*selected_mesh_);
					}

					ImGui::Separator();
					/*ImGui::Text("Binding");


					if (handle_container_.size() > 0)
					{
						SpaceDeformation<MESH, GRAPH>* space_deformation = static_cast<SpaceDeformation<MESH, GRAPH>*>(
							app_.module("SpaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));
						if (!imgui_handle_selector(space_deformation, selected_handle_, "Handle",
												   [&](GRAPH& g) { selected_handle_ = &g; }))
						{
							//std::cout << " there is a bug" << std::endl;
						}
						else
						{
	 
							const std::string handle_name = graph_provider_->mesh_name(*selected_handle_);

							std::string prefix = "local_handle";
					
							if (handle_name.size() > 0 &&
								 std::mismatch(prefix.begin(), prefix.end(), handle_name.begin()).first == prefix.end())
							{ 
								const char* items[] = {"MVC", "Green"};
								std::string current_item = "MVC";
								ImGuiComboFlags flags = ImGuiComboFlags_NoArrowButton;

								ImGuiStyle& style = ImGui::GetStyle();
								float w = ImGui::CalcItemWidth();
								float spacing = style.ItemInnerSpacing.x;
								float button_sz = ImGui::GetFrameHeight();
								ImGui::PushItemWidth(w - spacing * 2.0f - button_sz * 2.0f);
								if (ImGui::BeginCombo("##custom combo", current_item.c_str(),
													  ImGuiComboFlags_NoArrowButton))
								{
									for (int n = 0; n < IM_ARRAYSIZE(items); n++)
									{
										bool is_selected = (current_item == items[n]);
										if (ImGui::Selectable(items[n], is_selected))
											current_item = items[n];
										if (is_selected)
											ImGui::SetItemDefaultFocus();
									}
									ImGui::EndCombo();
								}

								selected_hdt_ = handle_container_[handle_name];
								if (ImGui::Button("Bind object"))
								{

									if (current_item == "MVC")
									{
										bind_influence_cage_mvc_handle(*selected_mesh_, p.vertex_position_,
																	   *selected_handle_,
																	   selected_hdt_->control_handle_vertex_position_);
									}

									else
									{
										std::cout << "not available yet" << std::endl;
									}
								}
						}

						
						}
					}*/
				}

				if (selected_tool_ == Axis)
				{
					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					imgui_combo_cells_set(md, control_set, "Control ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });

					if (control_set && control_set->size() > 0)
					{

						generate_axis(*selected_mesh_, p.vertex_position_, control_set);
					}
				}
			}
		}
	}

private:
	MESH* selected_mesh_;
	MESH* selected_cage_;
	GRAPH* selected_handle_;

	std::shared_ptr<cgogn::modeling::CageDeformationTool<MESH>> selected_cdt_;
	std::shared_ptr<cgogn::modeling::HandleDeformationTool<MESH>> selected_hdt_;
	std::unordered_map<const MESH*, Parameters> parameters_;

	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	MeshProvider<MESH>* mesh_provider_;
	MeshProvider<GRAPH>* graph_provider_;

	SurfaceRender<MESH>* surface_render_;
	SurfaceRender<GRAPH>* graph_render_;

	SurfaceDifferentialProperties<MESH>* surface_diff_pptes_;

	SurfaceSelectionPO<MESH>* surface_selection_;
	GraphSelection<GRAPH>* graph_selection_;

	SurfaceDeformation<MESH>* surface_deformation_;
	GraphDeformation<GRAPH>* graph_deformation_;

	SelectionTool selected_tool_;

	bool influence_cage_mode_; 
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_CAGE_DEFORMATION_H_
