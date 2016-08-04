//
//// Copyright (c) 2016 CNRS
//// Authors: Anna Seppala
////
//// This file is part of hpp-intersect
//// hpp-intersect is free software: you can redistribute it
//// and/or modify it under the terms of the GNU Lesser General Public
//// License as published by the Free Software Foundation, either version
//// 3 of the License, or (at your option) any later version.
////
//// hpp-intersect is distributed in the hope that it will be
//// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//// General Lesser Public License for more details.  You should have
//// received a copy of the GNU Lesser General Public License along with
//// hpp-intersect  If not, see
//// <http://www.gnu.org/licenses/>.
//
//
#ifndef HPP_INTERSECT_INTERSECT_HH
#define HPP_INTERSECT_INTERSECT_HH

#include <hpp/intersect/fwd.hh>

namespace hpp {
    namespace intersect {
    // free functions ?i
    
    /// \addtogroup intersect
    /// \{
        std::vector<fcl::Vec3f> getIntersectionPoints (const fcl::CollisionObjectPtr_t& rom,
               const fcl::CollisionObjectPtr_t& affordance);
    /// \}
    
    } // namespace intersect
} // namespace hpp
    
#endif // HPP_INTERSECT_INTERSECT_HH

