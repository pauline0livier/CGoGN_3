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

#ifndef CGOGN_MODULE_MULTI_TOOLS_DEFORMATION_H_
#define CGOGN_MODULE_MULTI_TOOLS_DEFORMATION_H_

#include <cgogn/core/ui_modules/graph_provider.h>
#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/rendering/types.h>
#include <cgogn/rendering/ui_modules/graph_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/geometry/algos/normal.h>

#include <GLFW/glfw3.h>

#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

#include <cgogn/modeling/types/axis_deformation_tool.h>
#include <cgogn/modeling/types/cage_deformation_tool.h>
#include <cgogn/modeling/types/global_cage_deformation_tool.h>
#include <cgogn/modeling/types/handle_deformation_tool.h>

#include <cgogn/modeling/types/handle_parameters.h>
#include <cgogn/modeling/types/axis_parameters.h>
#include <cgogn/modeling/types/mesh_parameters.h>

#include <boost/synapse/connect.hpp>

#include <iostream>

#include <algorithm>
#include <string>

namespace cgogn
{

namespace ui
{

enum SelectionTool
{
	Handle = 0,
	Axis,
	Cage, 
};

using Vec2 = geometry::Vec2;
using Vec3 = geometry::Vec3;
using Mat3 = geometry::Mat3;
using Scalar = geometry::Scalar;

template <typename MESH, typename GRAPH>

/**
 * @class multiToolsDeformation
 * insertion of tools of type: Handle, Axis, Cage
 * deformation of these tools to deform the mesh
*/
class MultiToolsDeformation : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension == 2,
		"MultiToolsDeformation can only be used with meshes of dimension 2");

	template <typename T>
	using GraphAttribute = typename mesh_traits<GRAPH>::
														template Attribute<T>;

	template <typename T>
	using MeshAttribute = typename mesh_traits<MESH>::template Attribute<T>;

	using GraphVertex = typename mesh_traits<GRAPH>::Vertex;
	using GraphEdge = typename mesh_traits<GRAPH>::Edge;

	using MeshVertex = typename mesh_traits<MESH>::Vertex;
	using MeshEdge = typename mesh_traits<MESH>::Edge;
	using MeshFace = typename mesh_traits<MESH>::Face;

public:
	std::unordered_map<std::string, 
		std::shared_ptr<modeling::AxisDeformationTool<MESH>>> axis_container_;
	std::unordered_map<std::string, 
		std::shared_ptr<modeling::HandleDeformationTool<MESH>>> handle_container_;
	std::unordered_map<std::string, 
		std::shared_ptr<modeling::CageDeformationTool<MESH>>> cage_container_;
	std::unordered_map<std::string, 
	std::shared_ptr<modeling::GlobalCageDeformationTool<MESH>>> global_cage_container_;

	MultiToolsDeformation(const App& app)
		: ViewModule(app, 
		"MultiToolsDeformation (" + std::string{mesh_traits<MESH>::name} + ")"), 
			model_(nullptr), selected_mesh_(nullptr), 
			selected_handle_(nullptr), selected_axis_(nullptr), 
			selected_cage_(nullptr)
	{
	}


	~MultiToolsDeformation()
	{
	}

	
	/// @brief initialize mesh 
	/// set mesh_ for mesh_parameter
	/// create synapse connection when attribute or cell sets changed 
	/// create a cell container to insert the selected elements  
	/// @param m mesh
	void init_mesh(MESH* m)
	{
		parameters_[m] = new modeling::Parameters<MESH>();
		modeling::Parameters<MESH>& p = *parameters_[m];
		p.mesh_ = m;

		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::
									template attribute_changed_t<Vec3>>(
				m, [this, m](MeshAttribute<Vec3>* attribute) {
					modeling::Parameters<MESH>& p = *parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = 
							float32(geometry::mean_edge_length(*m, 
											p.vertex_position_.get()) / 6);
						p.update_selected_vertices_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));

		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::
								template cells_set_changed<MeshVertex>>(
				m, [this, m](CellsSet<MESH, MeshVertex>* set) {
					modeling::Parameters<MESH>& p = *parameters_[m];
					if (p.selected_vertices_set_ == set && 
														p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));

		MeshData<MESH>& md = mesh_provider_->mesh_data(*m);
		md.template add_cells_set<MeshVertex>();
	}

	/// @brief initialize graph
	/// call corresponding function 
	/// @param g graph 
	void init_graph(GRAPH* g)
	{
		GraphData<GRAPH>& gd = graph_provider_->graph_data(*g);

		if (gd.graph_name_.find('h') != std::string::npos)
		{
			init_handle_graph(g); 
		} 
		else {
			init_axis_graph(g); 
		}
	}

	

	/**
	 * set mesh to deform as model for this module 
	 * @param {MESH} m mesh to deform
	 * @param {shared_ptr<MeshAttribute<Vec3>>} vertex_position 
	 * @param {shared_ptr<MeshAttribute<Vec3>>} vertex_normal
	*/

	/// @brief set model to deform 
	/// @param m model to deform 
	/// @param vertex_position position of the vertices of the mesh
	/// @param vertex_normal normal position of the vertices of the mesh
	void set_model(MESH& m,
				const std::shared_ptr<MeshAttribute<Vec3>>& m_vertex_position,
				const std::shared_ptr<MeshAttribute<Vec3>>& m_vertex_normal)
	{
		model_ = &m;

		geometry::compute_normal<MeshVertex>(m, m_vertex_position.get(), 
												m_vertex_normal.get());

		mesh_provider_->emit_attribute_changed(m, m_vertex_normal.get());

		set_vertex_position(m, m_vertex_position);
	}


	/// @brief  loop on the cage of cage_container and global_cage_container
	/// @tparam FUNC function type
	/// @param f function to apply on each cage
	template <typename FUNC>
	void foreach_cage(const FUNC& f)
	{
		static_assert(is_ith_func_parameter_same<FUNC, 0, MESH&>::
										value, "Wrong function parameter type");
		static_assert(is_ith_func_parameter_same<FUNC, 1, 
			const std::string&>::value, "Wrong function parameter type");

		for (auto& [name, gcdt] : global_cage_container_)
			f(*(gcdt->global_cage_), name);

		for (auto& [name, cdt] : cage_container_)
			f(*(cdt->control_cage_), name);
	}


	/// @brief loop on the local cage of cage_container
	/// @tparam FUNC function type
	/// @param f function to apply on each local cage
	template <typename FUNC>
	void foreach_local_cage(const FUNC& f)
	{
		static_assert(is_ith_func_parameter_same<FUNC, 0, MESH&>::
								value, "Wrong function parameter type");
		static_assert(is_ith_func_parameter_same<FUNC, 1, 
			const std::string&>::value, "Wrong function parameter type");

		for (auto& [name, cdt] : cage_container_)
			f(*(cdt->control_cage_), name);
	}


	/// @brief loop on the handle of handle_container
	/// @tparam FUNC function type
	/// @param f function to apply on each handle
	template <typename FUNC>
	void foreach_handle(const FUNC& f)
	{
		for (auto& [name, hdt] : handle_container_)
			f(*(hdt->control_handle_), name);
	}


	/// @brief loop on the axis of axis_container
	/// @tparam FUNC function type
	/// @param f function to apply on each axis
	template <typename FUNC>
	void foreach_axis(const FUNC& f)
	{
		for (auto& [name, adt] : axis_container_)
			f(*(adt->control_axis_), name);
	}

	/////////////
	// SIGNALS //
	/////////////
	using handle_added = struct handle_added_ (*)
						(modeling::HandleDeformationTool<MESH>* h);
	using global_cage_added = struct global_cage_added_ (*)
						(modeling::GlobalCageDeformationTool<MESH>* gcdt);
	using cage_added = struct cage_added_ (*)
						(modeling::CageDeformationTool<MESH>* cdt);
	using axis_added = struct axis_added_ (*)
						(modeling::AxisDeformationTool<MESH>* a);

private:

	/**
	 * set mesh parameters vertex position from the mesh vertex position 
	 * set vertex base size relatively to edge length  
	 * create 
	 * @param {MESH} m mesh to deform
	 * @param {shared_ptr<MeshAttribute<Vec3>>} vertex_position 
	*/

	/// @brief set mesh_parameter vertex position from the mesh ones
	/// set vertex base size
	/// @param m mesh
	/// @param m_vertex_position position of the vertices of the mesh
	void set_vertex_position(MESH& m, 
				const std::shared_ptr<MeshAttribute<Vec3>>& m_vertex_position)
	{
		modeling::Parameters<MESH>& p = *parameters_[&m];
		p.vertex_position_ = m_vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = 
				float32(geometry::mean_edge_length(m, 
										p.vertex_position_.get()) / 6); 

			p.update_selected_vertices_vbo();
		}

		for (View* v : linked_views_)
			v->request_update();
	}


	/// @brief set handle_parameter vertex position from the graph ones
	/// set vertex base size
	/// @param g graph  
	/// @param g_vertex_position position of the vertices of the graph
	void set_handle_vertex_position(const GRAPH& g, 
			const std::shared_ptr<GraphAttribute<Vec3>>& g_vertex_position)
	{
		modeling::HandleParameters<GRAPH>& p = *handle_parameters_[&g];
		p.vertex_position_ = g_vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = 5.0;
			p.update_selected_vertices_vbo();
		}

		for (View* v : linked_views_)
			v->request_update();
	}

	/// @brief set axis_parameter vertex position from the graph ones
	/// set vertex base size
	/// @param g graph  
	/// @param g_vertex_position position of the vertices of the graph
	void set_axis_vertex_position(const GRAPH& g, 
			const std::shared_ptr<GraphAttribute<Vec3>>& g_vertex_position)
	{
		modeling::AxisParameters<GRAPH>& p = *axis_parameters_[&g];
		p.vertex_position_ = g_vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = 5.0;
			p.update_selected_vertices_vbo();
		}

		for (View* v : linked_views_)
			v->request_update();
	}



	/// @brief create global cage around model 
	/// Create GlobalCageDeformationTool
	/// @param object model to deform 
	void create_global_cage_tool(const MESH& object)
	{
		std::string cage_name = "global_cage";

		const auto [it, inserted] =
			global_cage_container_.emplace(cage_name, 
			std::make_shared<modeling::GlobalCageDeformationTool<MESH>>());

		modeling::GlobalCageDeformationTool<MESH>* gcdt = it->second.get();

		if (inserted)
		{
			MESH* cage = mesh_provider_->add_mesh(cage_name);
			std::shared_ptr<MeshAttribute<Vec3>> cage_vertex_position =
				add_attribute<Vec3, MeshVertex>(*cage, "position");

			mesh_provider_->set_mesh_bb_vertex_position(*cage, 
													cage_vertex_position);

			set_vertex_position(*cage, cage_vertex_position);
			MeshData<MESH>& md = mesh_provider_->mesh_data(object);

			std::tuple<Vec3, Vec3, Vec3> extended_bounding_box =
				modeling::get_extended_bounding_box(md.bb_min_, 
													md.bb_max_, 1.2);

			gcdt->create_global_cage(cage, cage_vertex_position.get(), 
									std::get<0>(extended_bounding_box),
									std::get<1>(extended_bounding_box));

			mesh_provider_->emit_connectivity_changed(*cage);
			mesh_provider_->emit_attribute_changed(*cage, 
												cage_vertex_position.get());

			MeshData<MESH>& c_md = mesh_provider_->mesh_data(*cage);

			ui::View* v1 = app_.current_view();

			surface_render_->set_vertex_position(*v1, *cage, 
												cage_vertex_position);
			surface_render_->set_render_faces(*v1, *cage, false);
		}
	}
	
	/// @brief create local handle on model 
	/// Create HandleDeformationTool
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param handle_set user selected set 
	void create_handle_tool(MESH& object, 
				const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
				CellsSet<MESH, MeshVertex>* handle_set)
	{
		int handle_number = handle_container_.size();
		const std::string handle_name = "local_handle" + 
										std::to_string(handle_number);

		const auto [it, inserted] =
			handle_container_.emplace(handle_name, 
				std::make_shared<modeling::HandleDeformationTool<MESH>>());
		modeling::HandleDeformationTool<MESH>* hdt = it->second.get();

		if (inserted)
		{
			GRAPH* handle = graph_provider_->add_graph(handle_name);

			auto handle_vertex_position = 
					add_attribute<Vec3, GraphVertex>(*handle, "position");
			auto handle_vertex_radius = 
					add_attribute<Scalar, GraphVertex>(*handle, "radius");

			set_handle_vertex_position(*handle, handle_vertex_position);

			auto mesh_vertex_normal = 
							get_attribute<Vec3, MeshVertex>(object, "normal");

			const Vec3 handle_position =
				modeling::get_mean_value_attribute_from_set(object, 
													object_vertex_position.get(), 
													handle_set);

			Vec3 normal = modeling::get_mean_value_attribute_from_set(object, 
									mesh_vertex_normal.get(), handle_set);
			normal.normalize();

			modeling::HandleParameters<GRAPH>& handle_p = 
												*handle_parameters_[handle];

			handle_p.name_ = handle_name; 

			hdt->create_space_tool(handle, handle_vertex_position.get(), 
								handle_vertex_radius.get(), handle_position,
								normal);

			CMap2::Vertex closest_vertex =
				modeling::closest_vertex_in_set_from_value(object, 
													object_vertex_position.get(), 
													handle_set, 
													handle_position);

			hdt->set_handle_mesh_vertex(closest_vertex);

			handle_p.normal_ = normal;

			graph_provider_->emit_connectivity_changed(*handle);
			graph_provider_->emit_attribute_changed(*handle, 
											handle_vertex_position.get());
			graph_provider_->emit_attribute_changed(*handle, 
											handle_vertex_radius.get());

			View* v1 = app_.current_view();
			graph_render_->set_vertex_position(*v1, *handle, 
											handle_vertex_position);

			graph_render_->set_vertex_radius(*v1, *handle, 
											handle_vertex_radius);

			hdt->set_geodesic_distance(object, object_vertex_position);

			modeling::Parameters<MESH>& model_p = *parameters_[model_];

			model_p.handle_mesh_vertex_ = closest_vertex; 
			model_p.selection_for_handle_ = true; 

			boost::synapse::emit<handle_added>(this, hdt);
		}
	}

	/// @brief create local axis on model 
	/// Create AxisDeformationTool 
	/// @param object model to deform  
	/// @param object_vertex_position position of the vertices of the model
	void create_axis_tool(const MESH& object, 
			const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position)
	{
		int axis_number = axis_container_.size();
		std::string axis_name = "local_axis" + std::to_string(axis_number);

		const auto [it, inserted] =
			axis_container_.emplace(axis_name, 
				std::make_shared<modeling::AxisDeformationTool<MESH>>());
		modeling::AxisDeformationTool<MESH>* adt = it->second.get();
		if (inserted)
		{
			GRAPH* axis = graph_provider_->add_graph(axis_name);

			auto axis_vertex_position = 
					add_attribute<Vec3, GraphVertex>(*axis, "position");
			auto axis_vertex_radius = 
					add_attribute<Scalar, GraphVertex>(*axis, "radius");

			set_axis_vertex_position(*axis, axis_vertex_position);

			auto mesh_vertex_normal = 
					get_attribute<Vec3, MeshVertex>(object, "normal");

			modeling::Parameters<MESH>& model_p = *parameters_[model_];

			std::vector<Vec3> inside_axis_position;
			std::vector<Vec3> axis_vertices_position;
			std::vector<Vec3> axis_normals;

			for (std::size_t i = 0; 
						i < model_p.selected_depth_vertices_.size(); i++)
			{
				std::pair<MeshVertex, MeshVertex> vertices_set = 
										model_p.selected_depth_vertices_[i];

				const Vec3 front_position = value<Vec3>(object, 
									object_vertex_position, vertices_set.first);
				axis_vertices_position.push_back(front_position);

				const Vec3 back_position = value<Vec3>(object, 
									object_vertex_position, vertices_set.second);

				inside_axis_position.push_back(
									(front_position + back_position) / 2.0);

				const Vec3& normal = 
					value<Vec3>(object, mesh_vertex_normal, vertices_set.first);
				axis_normals.push_back(normal);
			}

			adt->create_space_tool(axis, axis_vertex_position.get(), 
						axis_vertex_radius.get(), axis_vertices_position,
						axis_normals, inside_axis_position);

			modeling::AxisParameters<GRAPH>& axis_p = 
												*axis_parameters_[axis];
			axis_p.normal_ = {0.0, 0.0, 1.0}; 
			axis_p.set_number_of_handles(axis_vertices_position.size()); 

			axis_p.name_ = axis_name; 

			graph_provider_->emit_connectivity_changed(*axis);
			graph_provider_->emit_attribute_changed(*axis, 
												axis_vertex_position.get());
			graph_provider_->emit_attribute_changed(*axis, 
													axis_vertex_radius.get());

			View* v1 = app_.current_view();

			graph_render_->set_vertex_position(*v1, *axis, 
													axis_vertex_position);

			graph_render_->set_vertex_radius(*v1, *axis, 
													axis_vertex_radius);

			boost::synapse::emit<axis_added>(this, adt);
		}
	}


	/// @brief create local cage on model 
	/// create CageDeformationTool
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param control_set user selected vertices 
	void create_local_cage_tool(const MESH& object, 
			const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
			CellsSet<MESH, MeshVertex>* control_set)
	{
		int cage_number = cage_container_.size();
		std::string cage_name = "local_cage" + std::to_string(cage_number);

		const auto [it, inserted] =
			cage_container_.emplace(cage_name, 
				std::make_shared<modeling::CageDeformationTool<MESH>>());
		modeling::CageDeformationTool<MESH>* cdt = it->second.get();

		if (inserted)
		{
			MESH* l_cage = mesh_provider_->add_mesh(cage_name);

			std::shared_ptr<MeshAttribute<Vec3>> l_cage_vertex_position =
				add_attribute<Vec3, MeshVertex>(*l_cage, "position");

			mesh_provider_->set_mesh_bb_vertex_position(*l_cage, 
												l_cage_vertex_position);

			set_vertex_position(*l_cage, l_cage_vertex_position);

			std::shared_ptr<MeshAttribute<Vec3>> mesh_vertex_normal = 
							get_attribute<Vec3, MeshVertex>(object, "normal");

			Vec3 center = 
				modeling::get_mean_value_attribute_from_set(object, 
									object_vertex_position.get(), control_set); 

			std::tuple<Vec3, Vec3, Vec3> main_directions = 
				modeling::find_main_directions_from_set(object, 
										object_vertex_position.get(), 
										control_set, center); 

			std::pair<Vec3, Vec3> local_boundaries = 
				modeling::get_border_values_in_set_along_local_frame(object, 
										object_vertex_position.get(), 
										control_set, main_directions); 

			std::tuple<Vec3, Vec3, Vec3> extended_boundaries =
				modeling::get_extended_bounding_box(local_boundaries.first, 
											local_boundaries.second, 1.1f);

			cdt->create_space_tool(l_cage, l_cage_vertex_position.get(), 
									std::get<0>(extended_boundaries),
									std::get<1>(extended_boundaries), 
									main_directions);

			mesh_provider_->emit_connectivity_changed(*l_cage);
			mesh_provider_->emit_attribute_changed(*l_cage, 
											l_cage_vertex_position.get());

			View* v1 = app_.current_view();

			surface_render_->set_vertex_position(*v1, *l_cage, 
												l_cage_vertex_position);

			surface_render_->set_render_faces(*v1, *l_cage, false);

			std::shared_ptr<MeshAttribute<Vec3>> l_cage_vertex_normal =
				add_attribute<Vec3, MeshVertex>(*l_cage, "normal");

			geometry::compute_normal<MeshVertex>(*l_cage, 
											l_cage_vertex_position.get(), 
											l_cage_vertex_normal.get());

			cdt->set_center_control_cage(std::get<2>(extended_boundaries));

			modeling::Parameters<MESH>& cage_p = 
												*parameters_[l_cage];
			cage_p.name_ = cage_name; 

			boost::synapse::emit<cage_added>(this, cdt);
		}
	}


	/// @brief bind global cage to model and existing tools 
	/// initialize connection when global cage vertices are updated 
	/// 	changes propagating on the model and existing tools
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param global_cage 
	/// @param cage_vertex_position position of the vertices of the global cage
	/// @param binding_type chosen binding type (MVC or Green)
	void bind_global_cage(MESH& object, 
		const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
		MESH& global_cage, 
		std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position,
		const std::string& binding_type)
	{

		std::shared_ptr<modeling::GlobalCageDeformationTool<MESH>> gcdt = 
									global_cage_container_["global_cage"];

		gcdt->set_deformation_type(binding_type);

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		gcdt->init_bind_object(object, object_vertex_position.get(), 
												object_vertex_index.get());

		if (handle_container_.size() > 0)
		{
			for (auto& [name, hdt] : handle_container_)
			{
				if (gcdt->local_handle_weights_.count(name) == 0){
					Vec3 handle_position = hdt->get_handle_position(); 

					gcdt->init_bind_handle(name, handle_position);
				} 
			}
		}

		if (axis_container_.size())
		{
			for (auto& [name, adt] : axis_container_)
			{
				modeling::AxisParameters<GRAPH>& p_axis = 
								*axis_parameters_[adt->control_axis_];

				std::vector<GraphVertex> axis_vertices = 
												adt->get_axis_skeleton();

				gcdt->init_bind_axis(*(adt->control_axis_), name, 
									p_axis.vertex_position_, axis_vertices);
			}
		}

		gcdt->cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::
										template attribute_changed_t<Vec3>>(

				&global_cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cage_vertex_position.get() == attribute)
					{
						if (deformed_tool_ == Cage)
						{
							std::shared_ptr<modeling::
							GlobalCageDeformationTool<MESH>> current_gcdt =
									global_cage_container_["global_cage"];

							std::shared_ptr<MeshAttribute<uint32>> 
							object_vertex_index = 
								get_attribute<uint32, MeshVertex>(object, 
																"vertex_index");

							current_gcdt->deform_object(object, 
												object_vertex_position.get(),
												object_vertex_index.get());

							mesh_provider_->emit_attribute_changed(object, 
												object_vertex_position.get());

							if (current_gcdt->local_handle_weights_.size() > 0)
							{
								for (auto& [name, vw] : 
										current_gcdt->local_handle_weights_)
								{
									std::shared_ptr<modeling::
										HandleDeformationTool<MESH>> local_hdt =
											handle_container_[name];

									modeling::HandleParameters<GRAPH>& p_handle =
									*handle_parameters_[local_hdt->control_handle_];

									current_gcdt->deform_handle(
										*(local_hdt->control_handle_), name,
										p_handle.vertex_position_,
										local_hdt->get_handle_vertex());

									local_hdt->update_handle_position_variable(); 

									graph_provider_->emit_attribute_changed(
										*(local_hdt->control_handle_),
									local_hdt->control_handle_vertex_position_.get());
								}
							}

							if (current_gcdt->local_axis_weights_.size() > 0)
							{
								for (auto& [name, mw] : 
										current_gcdt->local_axis_weights_)
								{
									std::shared_ptr<modeling::
									AxisDeformationTool<MESH>> local_adt =
										axis_container_[name];

									modeling::AxisParameters<GRAPH>& p_axis =
									*axis_parameters_[local_adt->control_axis_];

									current_gcdt->deform_axis(
										*(local_adt->control_axis_), name,
										p_axis.vertex_position_, 
										local_adt->get_axis_skeleton());

									graph_provider_->emit_attribute_changed(
										*(local_adt->control_axis_), 
										local_adt->control_axis_vertex_position_.get());
								}
							}
						}
					}
				});
	}

	/// @brief update deformation type for global cage
	/// @param object 
	/// @param object_vertex_position 
	/// @param global_cage 
	/// @param cage_vertex_position 
	/// @param binding_type 
	void update_bind_global_cage(MESH& object, 
		const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
		MESH& global_cage, 
		std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position,
		const std::string& binding_type)
	{
		std::shared_ptr<modeling::GlobalCageDeformationTool<MESH>> gcdt = 
									global_cage_container_["global_cage"];

		gcdt->set_deformation_type(binding_type);

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		gcdt->bind_object(object, object_vertex_position.get(), 
												object_vertex_index.get());

		if (handle_container_.size() > 0)
		{
			for (auto& [name, hdt] : handle_container_)
			{
				Vec3 handle_position = hdt->get_handle_position();

				gcdt->bind_handle(name, handle_position);
			}
		}

		if (axis_container_.size())
		{
			for (auto& [name, adt] : axis_container_)
			{
				modeling::AxisParameters<GRAPH>& p_axis = 
								*axis_parameters_[adt->control_axis_];

				std::vector<GraphVertex> axis_vertices = 
												adt->get_axis_skeleton();

				gcdt->bind_axis(*(adt->control_axis_), name, 
									p_axis.vertex_position_, axis_vertices);
			}
		}
	}

	/// @brief bind local cage to model and surrounding tools 
	/// initialize connection when local cage vertices are updated 
	/// 	changes propagating on the model and surrounding tools
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param local_cage 
	/// @param cage_vertex_position position of the vertices of the local cage
	/// @param binding_type chosen binding type (MVC or Green)
	void bind_local_cage(MESH& object, 
		const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
		MESH& local_cage, 
		std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position,
		std::string binding_type)
	{

		modeling::Parameters<MESH>& p_cage = *parameters_[&local_cage];

		std::shared_ptr<modeling::
			CageDeformationTool<MESH>> cdt = cage_container_[p_cage.name_];

		cdt->set_deformation_type(binding_type);

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		cdt->init_bind_object(object, object_vertex_position, 
								object_vertex_index);

		if (handle_container_.size() > 0)
		{
			for (auto& [name, hdt] : handle_container_)
			{
				if (cdt->local_handle_data_.count(name) == 0){
					Vec3 handle_position = hdt->get_handle_position(); 

					cdt->init_bind_handle(name, handle_position);
				} 
			}
		}

		cdt->cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::
							template attribute_changed_t<Vec3>>(
				&local_cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cage_vertex_position.get() == attribute)
					{
						std::shared_ptr<modeling::
						CageDeformationTool<MESH>> current_cdt =
							cage_container_[p_cage.name_];

						std::shared_ptr<MeshAttribute<uint32>> 
						object_vertex_index = 
						get_attribute<uint32, MeshVertex>(object, 
														"vertex_index");

						current_cdt->update_virtual_cages(); 
						
						current_cdt->deform_object(object, 
											object_vertex_position.get(), 
											object_vertex_index.get());

						mesh_provider_->emit_attribute_changed(object, 
											object_vertex_position.get());

						if (current_cdt->local_handle_data_.size() > 0)
							{
								for (auto& [name, handle_data] : 
										current_cdt->local_handle_data_)
								{
									std::shared_ptr<modeling::
										HandleDeformationTool<MESH>> local_hdt =
											handle_container_[name];

									modeling::HandleParameters<GRAPH>& p_handle =
									*handle_parameters_[local_hdt->control_handle_];

									current_cdt->deform_handle(
										*(local_hdt->control_handle_), name,
										p_handle.vertex_position_,
										local_hdt->get_handle_vertex());

									local_hdt->update_handle_position_variable(); 

									// see update handle weights on object 

									graph_provider_->emit_attribute_changed(
										*(local_hdt->control_handle_),
									local_hdt->control_handle_vertex_position_.get());
								}
							}
					}
				});
	}

	/// @brief update bind local cage to model and surrounding tools 
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param local_cage 
	/// @param cage_vertex_position position of the vertices of the local cage
	/// @param binding_type chosen binding type (MVC or Green)
	void update_bind_local_cage(MESH& object, 
		const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
		MESH& local_cage, 
		std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position,
		std::string binding_type)
	{

		modeling::Parameters<MESH>& p_cage = *parameters_[&local_cage];

		std::shared_ptr<modeling::
			CageDeformationTool<MESH>> cdt = cage_container_[p_cage.name_];

		cdt->set_deformation_type(binding_type);

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		cdt->bind_object(object, object_vertex_position, 
								object_vertex_index);

		if (handle_container_.size() > 0)
		{
			for (auto& [name, hdt] : handle_container_)
			{
				Vec3 handle_position = hdt->get_handle_position();

				cdt->bind_handle(name, handle_position);
			}
		}
	}

	/// @brief bind local handle to model and surrounding tools 
	/// initialize connection when the handle vertex is updated 
	/// 	changes propagating on the model and surrounding tools
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param control_handle local handle
	/// @param handle_vertex_position position of the vertex of the local handle
	/// @param binding_type chosen binding type (round or spike)
	void bind_local_handle(MESH& object, 
		const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
		GRAPH& control_handle, 
		const std::shared_ptr<GraphAttribute<Vec3>>& handle_vertex_position,
		std::string binding_type)
	{

		modeling::HandleParameters<GRAPH>& p_handle = 
									*handle_parameters_[&control_handle];

		std::shared_ptr<modeling::HandleDeformationTool<MESH>> hdt = 
										handle_container_[p_handle.name_];

		hdt->set_deformation_type(binding_type);

		hdt->init_bind_object(object); 

		if (global_cage_container_.size() > 0)
		{
			std::shared_ptr<modeling::GlobalCageDeformationTool<MESH>> gcdt =
									global_cage_container_["global_cage"];

			
			if (gcdt->local_handle_weights_.count(p_handle.name_) == 0){
				Vec3 handle_position = hdt->get_handle_position(); 

				gcdt->init_bind_handle(p_handle.name_, handle_position);
			} 
		}

		if (cage_container_.size() > 0)
		{
			for (auto& [name, cdt] : cage_container_)
			{
			
				if (cdt->local_handle_data_.count(p_handle.name_) == 0)
				{
					Vec3 handle_position = hdt->get_handle_position(); 

					cdt->init_bind_handle(p_handle.name_, handle_position);
				} 
			}
			
		}

		MeshData<MESH>& md = mesh_provider_->mesh_data(object);

		hdt->handle_attribute_update_connection_ =
			boost::synapse::connect<typename GraphProvider<GRAPH>::
				template attribute_changed_t<Vec3>>(
				&control_handle, [&](GraphAttribute<Vec3>* attribute) {
					if (handle_vertex_position.get() == attribute)
					{
						if (deformed_tool_ == Handle)
						{
							std::shared_ptr<modeling::
							HandleDeformationTool<MESH>> current_hdt =
								handle_container_[p_handle.name_];

							std::shared_ptr<MeshAttribute<uint32>> 
								object_vertex_index =
								get_attribute<uint32, MeshVertex>(object, 
															"vertex_index");

							current_hdt->deform_object(object, 
											object_vertex_position.get(), 
											object_vertex_index.get());

							check_deformation_handle_shared_vertex(p_handle.name_); 

							mesh_provider_->emit_attribute_changed(object, 
											object_vertex_position.get()); 

							if (cage_container_.size() > 0)
							{
								for (auto& [name, cdt] : cage_container_)
								{

									int cage_index = cdt->local_handle_data_[p_handle.name_].cage_index_; 

									if (cage_index == 27)
									{
										// check if outside the cage 
										const Vec3 handle_position = current_hdt->get_handle_position(); 

										const Vec3 local_position = cdt->local_frame_*handle_position;  

										if (!(modeling::check_triple_projection_in_area(local_position, cdt->cage_local_bb_min_, cdt->cage_local_bb_max_)))
										{
											// need resize cage 
											//cdt->update_virtual_cages(); 
										}

										
									} else {
										// check if outside virtual cage 
									}
								}

							}
							
							md.update_bb();
							std::tuple<Vec3, Vec3, Vec3> 
							extended_bounding_box =
								modeling::get_extended_bounding_box(
									md.bb_min_, md.bb_max_, 1.2);

							Vec3 e_bb_min = 
								std::get<0>(extended_bounding_box);
							Vec3 e_bb_max = 
								std::get<1>(extended_bounding_box);

							if (global_cage_container_.size() > 0)
							{
								std::shared_ptr<modeling::
							GlobalCageDeformationTool<MESH>> gcdt =
									global_cage_container_["global_cage"];

									MESH* global_cage = gcdt->global_cage_;
									MeshData<MESH>& cmd = 
									mesh_provider_->mesh_data(*global_cage);

									if (!(e_bb_min == cmd.bb_min_) || 
											!(e_bb_max == cmd.bb_max_))
									{
										gcdt->update_global_cage(e_bb_min,
																 e_bb_max);

										mesh_provider_->
											emit_attribute_changed(
											*global_cage, 
											gcdt->global_cage_vertex_position_.get());

										// check if init bind object has been called first 
										gcdt->bind_object(object, 
												object_vertex_position.get(),
												object_vertex_index.get());

										if (gcdt->local_handle_weights_
																.size() > 0)
										{
											for (auto& [name_handle, vw] : 
												gcdt->local_handle_weights_)
											{
												std::shared_ptr<modeling::
													HandleDeformationTool<MESH>> local_hdt =
														handle_container_[name_handle];
												
												modeling::HandleParameters<GRAPH>& 
													p_handle =
													*handle_parameters_
														[local_hdt->control_handle_];

												Vec3 handle_position = 
													local_hdt->get_handle_position(); 

												gcdt->bind_handle(
													p_handle.name_, 
													handle_position);

												graph_provider_->emit_attribute_changed(
													*(local_hdt->control_handle_),
													local_hdt->control_handle_vertex_position_.get());
											}
										}

										/*if (gcdt->local_axis_weights_
																.size() > 0)
										{
											for (auto& [name_axis, mw] : 
												gcdt->local_axis_weights_)
											{
												std::shared_ptr<modeling::
													AxisDeformationTool<MESH>> local_adt =
													axis_container_[name_axis];
												
												modeling::AxisParameters<GRAPH>& p_axis =
													*axis_parameters_[local_adt->control_axis_];

												std::vector<GraphVertex> 
													axis_vertices = 
													local_adt->get_axis_skeleton();

												gcdt->bind_axis(
													*(local_adt->control_axis_), 
													p_axis.name_,
													p_axis.vertex_position_, 
													axis_vertices);

												graph_provider_->emit_attribute_changed(
													*(local_adt->control_axis_),
													local_adt->control_axis_vertex_position_.get());
											}
										}*/

									}
								}
							}
						}
					});
	}

	/// @brief update bind local handle to model and surrounding tools 
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param control_handle local handle
	/// @param handle_vertex_position position of the vertex of the local handle
	/// @param binding_type chosen binding type (round or spike)
	void update_bind_local_handle(MESH& object, 
		const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
		GRAPH& control_handle, 
		const std::shared_ptr<GraphAttribute<Vec3>>& handle_vertex_position,
		std::string binding_type)
	{

		modeling::HandleParameters<GRAPH>& p_handle = 
									*handle_parameters_[&control_handle];

		std::shared_ptr<modeling::HandleDeformationTool<MESH>> hdt = 
										handle_container_[p_handle.name_];

		hdt->set_deformation_type(binding_type);

		hdt->bind_object(); 
	}


	/// @brief bind local axis to model and surrounding tools 
	/// initialize connection when local axis vertices are updated 
	/// 	changes propagating on the model and surrounding tools
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param control_axis local axis  
	/// @param axis_vertex_position position of the vertices of the local axis
	/// @param binding_type chosen binding type (LBS or DQS)
	void bind_local_axis(MESH& object, 
		const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
		GRAPH& control_axis, 
		const std::shared_ptr<GraphAttribute<Vec3>>& axis_vertex_position,
		std::string binding_type)
	{

		modeling::AxisParameters<GRAPH>& p_axis = 
										*axis_parameters_[&control_axis];

		std::shared_ptr<modeling::AxisDeformationTool<MESH>> adt = 
											axis_container_[p_axis.name_];

		adt->set_deformation_type(binding_type);

		MeshData<MESH>& md = mesh_provider_->mesh_data(object);

		adt->init_bind_object(object, object_vertex_position);

		adt->axis_attribute_update_connection_ =
			boost::synapse::connect<typename GraphProvider<GRAPH>::
								template attribute_changed_t<Vec3>>(
				&control_axis, [&](GraphAttribute<Vec3>* attribute) {
					if (axis_vertex_position.get() == attribute)
					{
						if (deformed_tool_ == Axis)
						{
							std::shared_ptr<modeling::
								AxisDeformationTool<MESH>> current_adt =
								axis_container_[p_axis.name_];

							std::shared_ptr<MeshAttribute<uint32>> 
							object_vertex_index =
								get_attribute<uint32, MeshVertex>(object, 
															"vertex_index");

							if (current_adt->deformation_type_ == "LBS"){
								current_adt->set_axis_transformation(
													p_axis.transformations_);
							}

							if (current_adt ->deformation_type_ == "DQS")
							{
								current_adt->set_dual_quaternions_transformation(
													p_axis.dual_quaternion_transformations_);
							}
							

							current_adt->deform_object(object, 
											object_vertex_position.get(), 
											object_vertex_index.get());

							mesh_provider_->emit_attribute_changed(object, 
											object_vertex_position.get());
						}
					}
				});
	}


	/// @brief update bind local axis to model and surrounding tools 
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param control_axis local axis  
	/// @param axis_vertex_position position of the vertices of the local axis
	/// @param binding_type chosen binding type (LBS or DQS)
	void update_bind_local_axis(MESH& object, 
		const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
		GRAPH& control_axis, 
		const std::shared_ptr<GraphAttribute<Vec3>>& axis_vertex_position,
		std::string binding_type)
	{

		modeling::AxisParameters<GRAPH>& p_axis = 
										*axis_parameters_[&control_axis];

		std::shared_ptr<modeling::AxisDeformationTool<MESH>> adt = 
											axis_container_[p_axis.name_];

		adt->set_deformation_type(binding_type);

		adt->bind_object(object, object_vertex_position);
	}


protected:

	/// @brief initialize the module 
	/// load the necessary exterior modules 
	/// mesh_provider: storage of meshes and their mesh_data
	/// graph_provider: storage of graphs and their graph_data 
	/// surface_render: render module for mesh 
	/// graph_render: render module for graph 
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{
											mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) 
													{ init_mesh(&m); });
		connections_.push_back(boost::synapse::
				connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, 
				&MultiToolsDeformation<MESH, GRAPH>::init_mesh));

		graph_provider_ = static_cast<ui::GraphProvider<GRAPH>*>(
			app_.module("GraphProvider (" + std::string{
										mesh_traits<GRAPH>::name} + ")"));
		graph_provider_->foreach_graph([this](GRAPH& g, const std::string&)
													{ init_graph(&g); });
		connections_.push_back(boost::synapse::
			connect<typename GraphProvider<GRAPH>::graph_added>(
			graph_provider_, this,
				&MultiToolsDeformation<MESH, GRAPH>::init_graph));

		surface_render_ = static_cast<ui::SurfaceRender<MESH>*>(
			app_.module("SurfaceRender (" + 
				std::string{mesh_traits<MESH>::name} + ")"));

		graph_render_ = static_cast<ui::GraphRender<GRAPH>*>(
			app_.module("GraphRender (" + 
				std::string{mesh_traits<GRAPH>::name} + ")"));
	}


	/// @brief mouse pressed event
	/// call the corresponding events in the parameters
	/// pressed Shift to select vertices on selected mesh (model or cage)
	/// pressed Ctrl to selected vertices on selected graph (handle or axis)
	/// @param view current view
	/// @param button left (select) or right (unselect)
	/// @param x retrieved from mouse position
	/// @param y retrieved from mouse position
	void mouse_press_event(View* view, 
						int32 button, int32 x, int32 y) override
	{
		if (selected_mesh_ && view->shift_pressed())
		{

			modeling::Parameters<MESH>& p = *parameters_[selected_mesh_];

			p.select_vertices_set(view, button, x, y);
			mesh_provider_->emit_cells_set_changed(*selected_mesh_, 
												p.selected_vertices_set_);
		}

		if (selected_handle_ && view->control_pressed())
		{

			modeling::HandleParameters<GRAPH>& p = 
										*handle_parameters_[selected_handle_];

			if (p.vertex_position_)
			{
				p.select_vertices_set(view, button, x, y);
				graph_provider_->emit_cells_set_changed(*selected_handle_, 
												p.selected_vertices_set_);
			}
		} 
		
		if (selected_axis_ && view->control_pressed()) {
			modeling::AxisParameters<GRAPH>& p = 
										*axis_parameters_[selected_axis_];

			if (p.vertex_position_)
			{
				p.select_vertices_set(view, button, x, y);
				graph_provider_->emit_cells_set_changed(*selected_axis_, 
												p.selected_vertices_set_);
			}
			
		}
	}


	/// @brief key pressed event 
	/// call the corresponding events in the parameters 
	/// pressed D to deform selected mesh (model and cages)
	/// pressed H to deform selected handle 
	/// pressed A or Q to deform selected axis 
	/// @param view current view 
	/// @param key_code 
	void key_press_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_D)
		{
			if (selected_mesh_)
			{
				modeling::Parameters<MESH>& p = 
											*parameters_[selected_mesh_];
				p.key_pressed_D_event(view);
			}
		}

		if (key_code == GLFW_KEY_H)
		{ 
			if (selected_handle_)
			{
				modeling::HandleParameters<GRAPH>& p = 
										*handle_parameters_[selected_handle_];
				p.key_pressed_H_event(view);
			}
		}

		if (key_code == GLFW_KEY_Q || key_code == GLFW_KEY_A)
		{
			if (selected_axis_)
			{
				modeling::AxisParameters<GRAPH>& p = 
										*axis_parameters_[selected_axis_];
				p.key_pressed_A_event(view);
			}
		}
	}

	/// @brief key released event
	/// call the corresponding events in the parameters 
	/// released D on selected mesh (model and cage)
	/// released H on selected handle
	/// released A or Q on selected axis 
	/// @param view current view
	/// @param key_code 
	void key_release_event(View* view, int32 key_code) override
	{
		unused_parameters(view);
		if (key_code == GLFW_KEY_D)
		{
			modeling::Parameters<MESH>& p = *parameters_[selected_mesh_];
			p.key_release_D_event(view);
		}
		else if (key_code == GLFW_KEY_H)
		{
			modeling::HandleParameters<GRAPH>& p = 
										*handle_parameters_[selected_handle_];
			p.key_release_handle_event(view);
		}
		else if (key_code == GLFW_KEY_Q || key_code == GLFW_KEY_A)
		{
			modeling::AxisParameters<GRAPH>& p = 
										*axis_parameters_[selected_axis_];
			p.key_release_axis_event(view);
		}
	}


	/// @brief mouse move event
	/// call the corresponding events in the parameters 
	/// pressed shift to select mesh vertices (model and cages)
	/// pressed ctrl to select graph vertices
	/// @param view current view 
	/// @param x retrieved from mouse position 
	/// @param y retrieved from mouse position
	void mouse_move_event(View* view, int32 x, int32 y) override
	{

		if (selected_mesh_)
		{
			modeling::Parameters<MESH>& p = 
										*parameters_[selected_mesh_];
			p.mouse_displacement(view, x, y);

			mesh_provider_->emit_attribute_changed(*selected_mesh_, 
												p.vertex_position_.get());
		}

		if (selected_handle_)
		{
			modeling::HandleParameters<GRAPH>& p = 
										*handle_parameters_[selected_handle_];
			p.mouse_displacement(view, x, y);
			graph_provider_->emit_attribute_changed(*selected_handle_, 
												p.vertex_position_.get());
		}

		if (selected_axis_)
		{
			modeling::AxisParameters<GRAPH>& p = 
										*axis_parameters_[selected_axis_];
			p.mouse_displacement(view, x, y);
			graph_provider_->emit_attribute_changed(*selected_axis_, 
												p.vertex_position_.get());
		}
	}

	/// @brief update rendering of selected vertices 
	/// @param view current view 
	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			p->local_draw(view);
		}

		for (auto& [m, p] : handle_parameters_)
		{
			p->local_draw(view);
		}

		for (auto& [m, p] : axis_parameters_)
		{
			p->local_draw(view);
		}

	}

	/// @brief left panel used for user interface
	/// refreshed constantly 
	/// call the corresponding functions relative to user interaction
	void left_panel() override
	{
		imgui_mesh_selector(mesh_provider_, model_, "Object", 
							[&](MESH& m) {
			model_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		bool need_update = false;

		selected_mesh_ = model_;

		if (model_)
		{
			MeshData<MESH>& model_md = mesh_provider_->mesh_data(*model_);

			modeling::Parameters<MESH>& model_p = *parameters_[model_];

			auto object_vertex_position = 
					get_attribute<Vec3, MeshVertex>(*model_, "position");

			set_vertex_position(*model_, object_vertex_position);

			if (model_p.vertex_position_)
			{
				ImGui::Separator();
				ImGui::Text("Create spatial tool");
				ImGui::Separator();

				ImGui::Text("Global");
				if (ImGui::Button("Create global cage"))
				{
					create_global_cage_tool(*model_);
				}

				ImGui::Separator();
				ImGui::Text("Local");
				ImGui::RadioButton("Axis", 
						reinterpret_cast<int*>(&selected_tool_), Axis);
				ImGui::SameLine();
				ImGui::RadioButton("Cage", 
						reinterpret_cast<int*>(&selected_tool_), Cage);
				ImGui::SameLine();
				ImGui::RadioButton("Handle", 
						reinterpret_cast<int*>(&selected_tool_), Handle);

				if (selected_tool_ == Handle)
				{

					selected_mesh_ = model_;
					ImGui::Separator();

					imgui_combo_cells_set(model_md, 
									model_p.selected_vertices_set_, 
									"Set_for_handle",
									[&](CellsSet<MESH, MeshVertex>* cs) {
										model_p.selected_vertices_set_ = cs;
										model_p.update_selected_vertices_vbo();
										need_update = true;
									});

					ImGui::RadioButton("Set_handle", 
						reinterpret_cast<int*>(&model_p.selection_method_),
						(int)modeling::SelectionMethod::SingleCell);

					ImGui::SameLine();
					ImGui::RadioButton("Set_influence_area", 
						reinterpret_cast<int*>(&model_p.selection_method_),
						(int)modeling::SelectionMethod::WithinSphere);

					if (model_p.selection_method_ == 
								modeling::SelectionMethod::SingleCell)
					{
						if (model_p.selected_vertices_set_)
						{

							if (model_p.selected_vertices_set_->size() > 0)
							{
								CellsSet<MESH, MeshVertex>* control_set = 
											model_p.selected_vertices_set_;

								create_handle_tool(*model_, 
									model_p.vertex_position_, control_set);

								model_p.selected_vertices_set_->clear();
								mesh_provider_->emit_cells_set_changed(
									*model_, model_p.selected_vertices_set_);
							}
						}
					}

					if (model_p.selection_method_ == 
									modeling::SelectionMethod::WithinSphere)
					{

						if (model_p.selected_vertices_set_)
						{
							ImGui::Text("(nb elements: %d)", 
									model_p.selected_vertices_set_->size());
							if (ImGui::Button("Clear handle area ##vertices_set"))
							{
								model_p.selected_vertices_set_->clear();
								mesh_provider_->emit_cells_set_changed(
									*model_, model_p.selected_vertices_set_);
							}
						}
						ImGui::TextUnformatted("Drawing parameters");
						need_update |=
							ImGui::ColorEdit3("color_handle##vertices", 
								model_p.param_point_sprite_->color_.data(),
								ImGuiColorEditFlags_NoInputs);

						if (need_update || model_p.object_update_)
						{
							for (View* v : linked_views_)
								v->request_update();
						}

						if (ImGui::Button("Accept handle influence area##vertices_set"))
						{

							accept_handle_influence_area(); 

						}

						model_p.object_update_ = false;
					}
				}

				if (selected_tool_ == Axis)
				{
					selected_mesh_ = model_;
					ImGui::Separator();

					model_p.back_selection_ = true;

					imgui_combo_cells_set(model_md, 
							model_p.selected_vertices_set_, "Set_for_axis",
							[&](CellsSet<MESH, MeshVertex>* cs) {
								model_p.selected_vertices_set_ = cs;
								model_p.update_selected_vertices_vbo();
								need_update = true;
							});

					ImGui::RadioButton("Set_axis", 
						reinterpret_cast<int*>(&model_p.selection_method_),
						(int)modeling::SelectionMethod::SingleCell);

					ImGui::SameLine();
					ImGui::RadioButton("Set_influence_area", 
						reinterpret_cast<int*>(&model_p.selection_method_),
						(int)modeling::SelectionMethod::WithinSphere);

					if (model_p.selection_method_ == 
						modeling::SelectionMethod::SingleCell)
					{
						if (model_p.selected_vertices_set_)
						{

							if (ImGui::Button("Accept axis ##vertices_set"))
							{

								create_axis_tool(*model_, 
												model_p.vertex_position_);

								model_p.selected_vertices_set_->clear();
								mesh_provider_->emit_cells_set_changed(
									*model_, model_p.selected_vertices_set_);

								model_p.selected_depth_vertices_.clear(); 
							}
						}
					}

					if (model_p.selection_method_ == 
								modeling::SelectionMethod::WithinSphere)
					{
						ImGui::SliderFloat("Sphere_radius", 
							&(model_p.sphere_scale_factor_), 1.0f, 50.0f);

						if (model_p.selected_vertices_set_)
						{
							ImGui::Text("(nb elements: %d)", 
									model_p.selected_vertices_set_->size());
							if (ImGui::Button("Clear handle area ##vertices_set"))
							{
								model_p.selected_vertices_set_->clear();
								mesh_provider_->emit_cells_set_changed(
									*model_, model_p.selected_vertices_set_);
							}
						}
						ImGui::TextUnformatted("Drawing parameters");
						need_update |=
							ImGui::ColorEdit3("color_handle##vertices", 
								model_p.param_point_sprite_->color_.data(),
											ImGuiColorEditFlags_NoInputs);

						if (need_update || model_p.object_update_)
						{
							for (View* v : linked_views_)
								v->request_update();
						}

						if (ImGui::Button("Accept axis influence area##vertices_set"))
						{
							model_p.selected_vertices_set_->foreach_cell(
								[&](MeshVertex v) -> bool {
									influence_set_.push_back(v);
									return true;
								});

							std::string last_axis_name = "local_axis" + 
										std::to_string(axis_container_.size() - 1);

							std::shared_ptr<modeling::
							AxisDeformationTool<MESH>> current_adt =
									axis_container_[last_axis_name];
							
							current_adt->object_influence_area_ = influence_set_;

							influence_set_.clear();

							model_p.selected_vertices_set_->clear();
							mesh_provider_->emit_cells_set_changed(
									*model_, model_p.selected_vertices_set_);
						}

						model_p.object_update_ = false;
					}
				}

				if (selected_tool_ == Cage)
				{
					selected_mesh_ = model_;
					ImGui::Separator();

					imgui_combo_cells_set(model_md, 
							model_p.selected_vertices_set_, "Set_for_cage",
							[&](CellsSet<MESH, MeshVertex>* cs) {
								model_p.selected_vertices_set_ = cs;
								model_p.update_selected_vertices_vbo();
								need_update = true;
							});

					model_p.selection_method_ = 
						modeling::SelectionMethod::WithinSphere;
					model_p.back_selection_ = true;

					ImGui::SliderFloat("Sphere radius", 
							&(model_p.sphere_scale_factor_), 10.0f, 100.0f);

					if (model_p.selected_vertices_set_)
					{
						ImGui::Text("(nb elements: %d)", 
							model_p.selected_vertices_set_->size());
						if (ImGui::Button("Clear##vertices_set"))
						{
							model_p.selected_vertices_set_->clear();
							mesh_provider_->emit_cells_set_changed(
									*model_, model_p.selected_vertices_set_);
						}
					}
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("color##vertices", 
								model_p.param_point_sprite_->color_.data(),
									ImGuiColorEditFlags_NoInputs);

					if (need_update || model_p.object_update_)
					{
						for (View* v : linked_views_)
							v->request_update();
					}

					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					if (ImGui::Button("Accept control set##vertices_set"))
					{
						control_set = model_p.selected_vertices_set_;

						create_local_cage_tool(*model_, 
									model_p.vertex_position_, control_set);

						model_p.selected_vertices_set_->clear();
						
						mesh_provider_->emit_cells_set_changed(*model_, 
											model_p.selected_vertices_set_);
					
					}

					model_p.object_update_ = false;

					ImGui::Separator();
					ImGui::Text("Binding");

								
				}

				ImGui::Separator();
				ImGui::Separator();
				ImGui::Text("Deform one tool");
				ImGui::RadioButton("Deform Axis", 
					reinterpret_cast<int*>(&deformed_tool_), Axis);
				ImGui::SameLine();
				ImGui::RadioButton("Deform Cage", 
					reinterpret_cast<int*>(&deformed_tool_), Cage);
				ImGui::SameLine();
				ImGui::RadioButton("Deform Handle", 
					reinterpret_cast<int*>(&deformed_tool_), Handle);

				if (deformed_tool_ == Cage)
				{

					if ((cage_container_.size() + 
								global_cage_container_.size()) > 0)
					{
						MultiToolsDeformation<MESH, GRAPH>* 
							multi_tools_deformation =
							static_cast<MultiToolsDeformation<MESH, GRAPH>*>(
								app_.module("MultiToolsDeformation (" + 
									std::string{mesh_traits<MESH>::name} + ")"));

						imgui_cage_selector(multi_tools_deformation, 
							selected_cage_, "Cage", [&](MESH& m) {
								if (selected_cage_)
								{
									modeling::Parameters<MESH>& old_p = 
												*parameters_[selected_cage_];
									old_p.selected_vertices_set_ = nullptr;
								}
								selected_cage_ = &m;
								mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
							});

						if (selected_cage_)
						{
							selected_handle_ = nullptr;
							selected_axis_ = nullptr;

							modeling::Parameters<MESH>& old_p = 
												*parameters_[selected_mesh_];
							old_p.selected_vertices_set_ = nullptr;

							selected_mesh_ = selected_cage_;
							modeling::Parameters<MESH>& cage_p = 
												*parameters_[selected_mesh_];

							const std::string cage_name = 
								mesh_provider_->mesh_name(*selected_cage_);

							std::string prefix = "local_cage";
							if (cage_name.size() > 0)
							{
								if (std::mismatch(prefix.begin(), 
										prefix.end(), 
								cage_name.begin()).first == prefix.end())
								{
									selected_cdt_ = cage_container_[cage_name];

									if (ImGui::Button("Reset deformation"))
									{
										selected_cdt_->reset_deformation(); 
										mesh_provider_->
											emit_attribute_changed(
											*(selected_cdt_->control_cage_), 
											selected_cdt_->control_cage_vertex_position_.get());
									} 

									ImGui::Separator();
									if (!selected_cdt_->deformation_type_.empty()){
										const char* extension = selected_cdt_->deformation_type_.c_str();

										ImGui::Text("Current binding:");
										ImGui::SameLine();
										ImGui::Text(extension);
									}

								ImGui::Separator();
									static ImGuiComboFlags flags = 0;

								const char* items[] = { "MVC", "Green" };
								static const char* item_current = items[0];       

								if (ImGui::BeginCombo("combo local cage", 
															item_current, flags)) 
								{
									for (int n = 0; n < IM_ARRAYSIZE(items); n++)
									{
										bool is_selected = (item_current == items[n]);
										if (ImGui::Selectable(items[n], is_selected))
											item_current = items[n];
										if (is_selected)
											ImGui::SetItemDefaultFocus();   
										}
									ImGui::EndCombo();
								}
								

								if (ImGui::Button("Bind local cage"))
								{

									if (selected_cdt_->deformation_type_.empty()){
										bind_local_cage(*model_, 
													model_p.vertex_position_, 
													*(selected_cdt_->control_cage_),
											selected_cdt_->control_cage_vertex_position_,
													item_current);
									}
									else if (selected_gcdt_->deformation_type_ != item_current)
									{
										update_bind_local_cage(*model_, 
													model_p.vertex_position_, 
													*(selected_cdt_->control_cage_),
											selected_cdt_->control_cage_vertex_position_,
													item_current);
									}
								}

								}
								else
								{
									selected_gcdt_ = global_cage_container_[cage_name];

									if (ImGui::Button("Reset deformation"))
									{
										selected_gcdt_->reset_deformation(); 
										mesh_provider_->
											emit_attribute_changed(
											*(selected_gcdt_->global_cage_), 
											selected_gcdt_->global_cage_vertex_position_.get());
									} 

									ImGui::Separator();
									if (!selected_gcdt_->deformation_type_.empty()){
										const char* extension = selected_gcdt_->deformation_type_.c_str();

										ImGui::Text("Current binding:");
										ImGui::SameLine();
										ImGui::Text(extension);
									}

								ImGui::Separator();

									static ImGuiComboFlags flags = 0;
									const char* items[] = {"MVC", "Green"};

									static const char* current_item = "MVC";
				
									if (ImGui::BeginCombo("##custom combo", current_item, flags))
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

									if (ImGui::Button("Bind global cage"))
									{

										if (selected_gcdt_->deformation_type_.empty()){
											bind_global_cage(*model_, model_p.vertex_position_, *(selected_gcdt_->global_cage_),
										selected_gcdt_->global_cage_vertex_position_, current_item);
										
										} else if (selected_gcdt_->deformation_type_ != current_item)
										{
											update_bind_global_cage(*model_, model_p.vertex_position_, *(selected_gcdt_->global_cage_),
										selected_gcdt_->global_cage_vertex_position_, current_item); 
	
										}
									}
								}

								MeshData<MESH>& cage_md = 
								mesh_provider_->mesh_data(*selected_cage_);

								imgui_combo_cells_set(cage_md, 
									cage_p.selected_vertices_set_, 
									"Cage D sets",
									[&](CellsSet<MESH, MeshVertex>* cs) {
										cage_p.selected_vertices_set_ = cs;
										cage_p.update_selected_vertices_vbo();
										need_update = true;
									});

								if (cage_p.selected_vertices_set_)
								{
									ImGui::Text("(nb elements: %d)", 
									cage_p.selected_vertices_set_->size());

									if (ImGui::Button("Clear cage selection##vertices_set"))
									{
										cage_p
											.selected_vertices_set_->clear();
										mesh_provider_->emit_cells_set_changed(
											*selected_cage_,
											cage_p.selected_vertices_set_);

										selected_cage_ = nullptr; 
									}
								}
								ImGui::TextUnformatted("Drawing parameters");
								need_update |=
									ImGui::ColorEdit3("color##vertices", 
									cage_p.param_point_sprite_->color_.data(),
									ImGuiColorEditFlags_NoInputs);

								if (need_update || cage_p.object_update_)
								{
									for (View* v : linked_views_)
										v->request_update();
								}
							}
						}
					}
				}

				if (deformed_tool_ == Handle)
				{

					if (handle_container_.size() > 0)
					{
						MultiToolsDeformation<MESH, GRAPH>* 
						multi_tools_deformation =
							static_cast<MultiToolsDeformation<MESH, GRAPH>*>(
								app_.module("MultiToolsDeformation (" + 
								std::string{mesh_traits<MESH>::name} + ")"));

						imgui_handle_selector(multi_tools_deformation, 
									selected_handle_, "Handle",
								[&](GRAPH& g) 
								{ 
									if (selected_handle_ )
									{
										modeling::HandleParameters<GRAPH>& old_p = 
										*handle_parameters_[selected_handle_];

										modeling::HandleParameters<GRAPH>& new_p = 
										*handle_parameters_[&g];

										if (old_p.name_ != new_p.name_){
											old_p.selected_vertices_set_ = nullptr;
										}
									
									}
									
									selected_handle_ = &g; 
								});

						

						if (selected_handle_)
						{
							selected_mesh_ = nullptr;
						
							modeling::HandleParameters<GRAPH>& handle_p = 
										*handle_parameters_[selected_handle_];

							const std::string handle_name = 
							graph_provider_->graph_name(*selected_handle_);

							if (handle_name.size() > 0)
							{
								selected_hdt_ = handle_container_
															[handle_name];

								if (ImGui::Button("Reset deformation"))
								{
									selected_hdt_->reset_deformation(); 

									graph_provider_->emit_attribute_changed(
										*(selected_hdt_->control_handle_),
									selected_hdt_->control_handle_vertex_position_.get());
								} 

								ImGui::Separator();
								if (!selected_hdt_->deformation_type_.empty()){
									const char* extension = selected_hdt_->deformation_type_.c_str();

									ImGui::Text("Current binding:");
									ImGui::SameLine();
									ImGui::Text(extension);
								}

								ImGui::Separator();
								static ImGuiComboFlags flags = 0;

								const char* items[] = { "Spike", "Round" };

								static const char* item_current = items[0];       

								if (ImGui::BeginCombo("combo handle", 
															item_current, flags)) 
								{
									for (int n = 0; n < IM_ARRAYSIZE(items); n++)
									{
										bool is_selected = (item_current == items[n]);
										if (ImGui::Selectable(items[n], is_selected))
											item_current = items[n];
										if (is_selected)
											ImGui::SetItemDefaultFocus();   
										}
									ImGui::EndCombo();
								}

								if (ImGui::Button("Bind handle"))
								{

									if (selected_hdt_->deformation_type_.empty()){
										bind_local_handle(*model_, 
										model_p.vertex_position_, 
										*(selected_hdt_->control_handle_),
										selected_hdt_->control_handle_vertex_position_, 
										item_current);

									} else if (selected_hdt_->deformation_type_ != item_current)
									{
										update_bind_local_handle(*model_, 
										model_p.vertex_position_, 
										*(selected_hdt_->control_handle_),
										selected_hdt_->control_handle_vertex_position_, 
										item_current);
									}
								}

								GraphData<GRAPH>& handle_gd = 
								graph_provider_->graph_data(*selected_handle_);

								imgui_combo_cells_set_graph(
									handle_gd, 
									handle_p.selected_vertices_set_, 
									"Handle_sets",
									[&](CellsSet<GRAPH, GraphVertex>* cs) {
										handle_p.selected_vertices_set_ = cs;
										handle_p.update_selected_vertices_vbo();
										need_update = true;
									});

								if (handle_p.selected_vertices_set_)
								{
									ImGui::Text("(nb elements: %d)", 
									handle_p.selected_vertices_set_->size());

									if (ImGui::Button("Clear handle selection##vertices_set"))
									{
										handle_p.selected_vertices_set_->clear();
										graph_provider_->emit_cells_set_changed(
											*selected_handle_,
											handle_p.selected_vertices_set_);

										selected_handle_ = nullptr; 
									}
								}
								ImGui::TextUnformatted("Drawing parameters");
								need_update |=
									ImGui::ColorEdit3("color##vertices", 
									handle_p.param_point_sprite_->color_.data(),
									ImGuiColorEditFlags_NoInputs);

								if (need_update || handle_p.object_update_)
								{
									for (View* v : linked_views_)
										v->request_update();
								}
							}
						}
					}
				}

				if (deformed_tool_ == Axis)
				{

					if (axis_container_.size() > 0)
					{
						MultiToolsDeformation<MESH, GRAPH>* 
						multi_tools_deformation =
							static_cast<MultiToolsDeformation<MESH, GRAPH>*>(
								app_.module("MultiToolsDeformation (" + 
								std::string{mesh_traits<MESH>::name} + ")"));

						imgui_axis_selector(multi_tools_deformation, 
								selected_axis_, "Axis", 
								[&](GRAPH& g) {
									if (selected_axis_)
									{
										modeling::AxisParameters<GRAPH>& old_p = 
										*axis_parameters_[selected_axis_];

										modeling::AxisParameters<GRAPH>& new_p = 
										*axis_parameters_[&g];

										if (old_p.name_ != new_p.name_){
											old_p.selected_vertices_set_ = nullptr;
										}

									}
									selected_axis_ = &g;
								});

						if (selected_axis_)
						{
							selected_mesh_ = nullptr;

							modeling::AxisParameters<GRAPH>& axis_p = 
										*axis_parameters_[selected_axis_];

							const std::string axis_name = 
								graph_provider_->graph_name(*selected_axis_);

							std::string prefix = "local_axis";
							if (axis_name.size() > 0)
							{
								selected_adt_ = axis_container_[axis_name];

								if (ImGui::Button("Reset deformation"))
								{

									selected_adt_->reset_deformation(); 

									graph_provider_->emit_attribute_changed(
										*(selected_adt_->control_axis_), 
										selected_adt_->control_axis_vertex_position_.get());
								} 

								ImGui::Separator();
								if (!selected_adt_->deformation_type_.empty()){
									const char* extension = selected_adt_->deformation_type_.c_str();

									ImGui::Text("Current binding:");
									ImGui::SameLine();
									ImGui::Text(extension);
								}

								ImGui::Separator();

								static ImGuiComboFlags flags = 0;

								const char* items[] = { "LBS", "DQS" };
								static const char* item_current = items[0];       

								if (ImGui::BeginCombo("combo local axis", 
															item_current, flags)) 
								{
									for (int n = 0; n < IM_ARRAYSIZE(items); n++)
									{
										bool is_selected = (item_current == items[n]);
										if (ImGui::Selectable(items[n], is_selected))
											item_current = items[n];
										if (is_selected)
											ImGui::SetItemDefaultFocus();   
										}
									ImGui::EndCombo();
								}

								if (ImGui::Button("Bind axis"))
								{

									if (selected_adt_->deformation_type_.empty()){
										bind_local_axis(*model_, 
												model_p.vertex_position_, 
												*(selected_adt_->control_axis_),
										selected_adt_->control_axis_vertex_position_, 
												item_current);

									} else if (selected_adt_->deformation_type_ != item_current){

										update_bind_local_axis(*model_, 
												model_p.vertex_position_, 
												*(selected_adt_->control_axis_),
										selected_adt_->control_axis_vertex_position_, 
												item_current);
									}
								}

								GraphData<GRAPH>& axis_gd = 
								graph_provider_->graph_data(*selected_axis_);

								imgui_combo_cells_set_graph(axis_gd, axis_p.selected_vertices_set_, "Axis_sets",
									[&](CellsSet<GRAPH, GraphVertex>* cs) {
										axis_p.selected_vertices_set_ = cs;
										axis_p.update_selected_vertices_vbo();
										need_update = true;
									});

								if (axis_p.selected_vertices_set_)
								{
									ImGui::Text("(nb elements: %d)", 
									axis_p.selected_vertices_set_->size());
									if (ImGui::Button("Clear axis selection##vertices_set"))
									{
										axis_p.selected_vertices_set_->clear();
										graph_provider_->emit_cells_set_changed(
											*selected_axis_,
											axis_p.selected_vertices_set_);

										selected_axis_ =  nullptr; 
									}
								}
								ImGui::TextUnformatted("Drawing parameters");
								need_update |=
									ImGui::ColorEdit3("color##vertices", 
									axis_p.param_point_sprite_->color_.data(),
										ImGuiColorEditFlags_NoInputs);

								if (need_update || axis_p.object_update_)
								{
									for (View* v : linked_views_)
										v->request_update();

								}
							}
						}
					}
				}
			}
		}
	}

private:
	MESH* model_;

	MESH* selected_mesh_;
	MESH* selected_cage_;

	GRAPH* selected_handle_;
	GRAPH* selected_axis_;

	std::shared_ptr<modeling::CageDeformationTool<MESH>> selected_cdt_;
	std::shared_ptr<modeling::GlobalCageDeformationTool<MESH>> 
															selected_gcdt_;
	std::shared_ptr<modeling::HandleDeformationTool<MESH>> selected_hdt_;
	std::shared_ptr<modeling::AxisDeformationTool<MESH>> selected_adt_;

	std::unordered_map<const MESH*, modeling::Parameters<MESH>*> parameters_;

	std::unordered_map<const GRAPH*, modeling::HandleParameters<GRAPH>*> 
														handle_parameters_;

	std::unordered_map<const GRAPH*, modeling::AxisParameters<GRAPH>*> 
														axis_parameters_;

	std::unordered_map<uint32,modeling::SharedVertexData> activation_map_; 

	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;

	std::unordered_map<const MESH*, 
	std::vector<std::shared_ptr<boost::synapse::connection>>> 
														mesh_connections_;

	std::unordered_map<const GRAPH*, 
	std::vector<std::shared_ptr<boost::synapse::connection>>> 
														graph_connections_;

	MeshProvider<MESH>* mesh_provider_;
	GraphProvider<GRAPH>* graph_provider_;

	SurfaceRender<MESH>* surface_render_;
	GraphRender<GRAPH>* graph_render_;

	SelectionTool selected_tool_;
	SelectionTool deformed_tool_;

	std::vector<MeshVertex> influence_set_;

	/// @brief initialize graph of type handle
	/// set graph_ for handle_parameter
	/// create synapse connection when attribute or cell sets changed 
	/// create a cell container to insert the selected elements  
	/// @param g graph 
	void init_handle_graph(GRAPH* g)
	{
		handle_parameters_[g] = new modeling::HandleParameters<GRAPH>();
		modeling::HandleParameters<GRAPH>& p = *handle_parameters_[g];
		p.graph_ = g;

		graph_connections_[g].push_back(
			boost::synapse::connect<typename GraphProvider<GRAPH>::
									template attribute_changed_t<Vec3>>(
				g, [this, g](GraphAttribute<Vec3>* attribute) {

					modeling::HandleParameters<GRAPH>& p = 
													*handle_parameters_[g];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = 3.5;
						p.update_selected_vertices_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));

		graph_connections_[g].push_back(
			boost::synapse::connect<typename GraphProvider<GRAPH>::
								template cells_set_changed<GraphVertex>>(
				g, [this, g](CellsSet<GRAPH, GraphVertex>* set) {
					modeling::HandleParameters<GRAPH>& p = 
													*handle_parameters_[g];
					if (p.selected_vertices_set_ == set && 
														p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));

		GraphData<GRAPH>& gd = graph_provider_->graph_data(*g);
		gd.template add_cells_set<GraphVertex>();
	}

	/// @brief initialize graph of type handle
	/// set graph_ for axis_parameter
	/// create synapse connection when attribute or cell sets changed 
	/// create a cell container to insert the selected elements  
	/// @param g graph 
	void init_axis_graph(GRAPH* g)
	{
		axis_parameters_[g] = new modeling::AxisParameters<GRAPH>();
		modeling::AxisParameters<GRAPH>& p = *axis_parameters_[g];
		p.graph_ = g;

		graph_connections_[g].push_back(
			boost::synapse::connect<typename GraphProvider<GRAPH>::
									template attribute_changed_t<Vec3>>(
				g, [this, g](GraphAttribute<Vec3>* attribute) {

					modeling::AxisParameters<GRAPH>& p = 
													*axis_parameters_[g];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = 3.5;
						p.update_selected_vertices_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));

		graph_connections_[g].push_back(
			boost::synapse::connect<typename GraphProvider<GRAPH>::
								template cells_set_changed<GraphVertex>>(
				g, [this, g](CellsSet<GRAPH, GraphVertex>* set) {
					modeling::AxisParameters<GRAPH>& p = 
													*axis_parameters_[g];
					if (p.selected_vertices_set_ == set && 
														p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));

		GraphData<GRAPH>& gd = graph_provider_->graph_data(*g);
		gd.template add_cells_set<GraphVertex>();
	}

	void accept_handle_influence_area()
	{
		modeling::Parameters<MESH>& model_p = *parameters_[model_];

		model_p.selected_vertices_set_->foreach_cell([&](MeshVertex v) -> bool {
			influence_set_.push_back(v);
			return true;
		});
 
		model_p.selection_for_handle_ = false; 
							
		std::string last_handle_name = "local_handle" + 
							std::to_string(handle_container_.size() - 1);

		std::shared_ptr<modeling::HandleDeformationTool<MESH>> current_hdt =
									handle_container_[last_handle_name];

		std::shared_ptr<MeshAttribute<uint32>> model_vertex_index =
								get_attribute<uint32, MeshVertex>(*model_, 
															"vertex_index");

		current_hdt->set_object_influence_area(*model_, model_vertex_index.get(),
												influence_set_); 
							
		influence_set_.clear();

		model_p.selected_vertices_set_->clear();
		mesh_provider_->emit_cells_set_changed(*model_, 
												model_p.selected_vertices_set_);

		if (handle_container_.size() < 1)
		{
			for ( const auto &myPair : current_hdt->object_influence_area_ ) {
				uint32 vertex_index = myPair.first;
				MeshVertex v = myPair.second.vertex; 

				if (activation_map_.find(vertex_index) != activation_map_.end()){
					const size_t number_of_handles_ = 
								activation_map_[vertex_index].tool_names.size(); 

					if (number_of_handles_ == 1)
					{
						for (std::string h_name : activation_map_[vertex_index].tool_names) 
						{
							std::shared_ptr<modeling::HandleDeformationTool<MESH>> other_hdt =
											handle_container_[h_name];
							other_hdt->object_influence_area_[vertex_index].shared = true; 
							other_hdt->shared_vertex[vertex_index] = v;
						}
					} 

					activation_map_[vertex_index].tool_names.insert(last_handle_name);
					current_hdt->object_influence_area_[vertex_index].shared = true; 
					current_hdt->shared_vertex[vertex_index] = v; 
				} else {
					modeling::SharedVertexData new_data; 
					new_data.tool_names = {}; 
					new_data.current_max_translation = 0.0; 

					const auto [it, inserted] =
						activation_map_.emplace(vertex_index, new_data);

					if (inserted)
					{
						activation_map_[vertex_index].tool_names.insert(last_handle_name); 
						
					}
				}
			}
		}
	}

	void check_deformation_handle_shared_vertex(const std::string& handle_name)
	{
		std::shared_ptr<MeshAttribute<Vec3>> model_vertex_position =
								get_attribute<Vec3, MeshVertex>(*model_, 
															"position");
	
		std::shared_ptr<modeling::HandleDeformationTool<MESH>> current_hdt =
									handle_container_[handle_name];

		for ( const auto &myPair : current_hdt->shared_vertex ) {
				uint32 vertex_index = myPair.first;
				MeshVertex v = myPair.second; 

				if (current_hdt->object_influence_area_[vertex_index].max_local_translation > 
						activation_map_[vertex_index].current_max_translation)
				{

					if (activation_map_[vertex_index].current_max_translation == 0.0)
					{
						value<Vec3>(*model_, model_vertex_position.get(), v) += 
						current_hdt->object_influence_area_[vertex_index].local_translation; 		

					} 
					else 
					{
						std::string old_max_handle_name = activation_map_[vertex_index].name_max_handle; 
						
						std::shared_ptr<modeling::HandleDeformationTool<MESH>> old_hdt =
									handle_container_[old_max_handle_name];

						Vec3 current_position = value<Vec3>(*model_, model_vertex_position.get(), v); 

						Vec3 reset_position =  current_position - 
						old_hdt->get_deformation_from_norm(activation_map_[vertex_index].current_max_translation); 


						value<Vec3>(*model_, model_vertex_position.get(), v) = 
							reset_position + current_hdt->get_deformation_from_norm(current_hdt->object_influence_area_[vertex_index].max_local_translation); 

						std::shared_ptr<MeshAttribute<uint32>> model_vertex_index =
								get_attribute<uint32, MeshVertex>(*model_, 
															"vertex_index");

						const uint32 handle_mesh_vertex_index =
							old_hdt->get_handle_mesh_vertex_index(*model_, model_vertex_index.get());

						if (vertex_index == handle_mesh_vertex_index)
						{
							std::cout << "same vertex as handle" << std::endl; 
							const Vec3 new_position = value<Vec3>(*model_, model_vertex_position.get(), v); 
							old_hdt->update_handle_position(new_position) ; 
							old_hdt->bind_object(); 

						} else {
							old_hdt->bind_isolated_vertex(vertex_index); 
						}
					}	

					activation_map_[vertex_index].name_max_handle = handle_name; 
					activation_map_[vertex_index].current_max_translation = 
								current_hdt->object_influence_area_[vertex_index].max_local_translation;
					
				}

		}

		mesh_provider_->emit_attribute_changed(*model_, 
												model_vertex_position.get());


	}



};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MULTI_TOOLS_DEFORMATION_H_
