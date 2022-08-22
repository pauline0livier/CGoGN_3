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

#ifndef CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cells_set.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>


namespace cgogn
{

namespace modeling
{

template <typename MESH>
class SpaceDeformationTool 
{

    template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;

    using Vec2 = geometry::Vec2;
    using Vec3 = geometry::Vec3;

public:

    MESH* influence_cage_; 
    cgogn::ui::CellsSet<MESH, Vertex>* influence_area_; 

    Eigen::VectorXd attenuation_;

    SpaceDeformationTool(): influence_cage_(nullptr), influence_area_(nullptr), influence_cage_vertex_position_(nullptr)
	{
	
	}

    virtual ~SpaceDeformationTool()
    {

    }

    virtual void create_space_tool()
    {
        
    }

    /*void set_influence_area(CellsSet<MESH, MeshVertex>& influence_area){
        influence_area_ = influence_area; 
    }*/

protected:
    std::shared_ptr<Attribute<Vec3>> influence_cage_vertex_position_; 
    
    

}; 

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_
