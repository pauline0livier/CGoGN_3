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

#include <cgogn/ui/app_PO.h>
#include <cgogn/ui/module_PO.h>

namespace cgogn
{

namespace ui
{

/*****************************************************************************/
// Module
/*****************************************************************************/

Module_PO::Module_PO(const App_PO& app, const std::string& name) : app_(app), name_(name)
{
	app.modules_.push_back(this);
}

Module_PO::~Module_PO()
{
}

void Module_PO::init()
{
}

void Module_PO::main_menu()
{
}

void Module_PO::left_panel()
{
}

void Module_PO::popups()
{
}

void Module_PO::close_event()
{
}

/*****************************************************************************/
// ViewModule
/*****************************************************************************/

ViewModule_PO::ViewModule_PO(const App_PO& app, const std::string& name) : Module_PO(app, name)
{
}

ViewModule_PO::~ViewModule_PO()
{
}

void ViewModule_PO::mouse_press_event(View_PO*, int32, int32, int32)
{
}
void ViewModule_PO::mouse_release_event(View_PO*, int32, int32, int32)
{
}
void ViewModule_PO::mouse_dbl_click_event(View_PO*, int32, int32, int32)
{
}
void ViewModule_PO::mouse_move_event(View_PO*, int32, int32)
{
}
void ViewModule_PO::mouse_wheel_event(View_PO*, int32, int32)
{
}
void ViewModule_PO::key_press_event(View_PO*, int32)
{
}
void ViewModule_PO::key_release_event(View_PO*, int32)
{
}

void ViewModule_PO::draw(View_PO*)
{
}

/*****************************************************************************/
// ProviderModule
/*****************************************************************************/

ProviderModule_PO::ProviderModule_PO(const App_PO& app, const std::string& name) : Module_PO(app, name)
{
}

ProviderModule_PO::~ProviderModule_PO()
{
}

} // namespace ui

} // namespace cgogn
