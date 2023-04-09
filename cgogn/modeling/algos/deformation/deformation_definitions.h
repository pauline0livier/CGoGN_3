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

#ifndef CGOGN_MODELING_ALGOS_DEFORMATION_DEFINITIONS_H_
#define CGOGN_MODELING_ALGOS_DEFORMATION_DEFINITIONS_H_

namespace cgogn
{

namespace modeling
{

enum class SelectionMethod : int
{
	SingleCell = 0,
	WithinSphere,
	ConnectedComponent
};

enum class SelectingCell : int
{
	VertexSelect = 0,
	EdgeSelect,
	FaceSelect
};

struct VectorWeights
{
	Eigen::VectorXd position_;
	Eigen::VectorXd normal_;
};

struct MatrixWeights
{
	Eigen::MatrixXd position_;
	Eigen::MatrixXd normal_;
};

struct SharedVertexData
{
	std::set<std::string> tool_names; 
	double current_max_translation; 
}; 

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_DEFORMATION_DEFINITIONS_H_
