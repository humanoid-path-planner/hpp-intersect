//
//// Copyright (c) 2016 CNRS
//// Authors: Anna Seppala
////
//// This file is part of hpp-affordance
//// hpp-affordance is free software: you can redistribute it
//// and/or modify it under the terms of the GNU Lesser General Public
//// License as published by the Free Software Foundation, either version
//// 3 of the License, or (at your option) any later version.
////
//// hpp-affordance is distributed in the hope that it will be
//// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//// General Lesser Public License for more details.  You should have
//// received a copy of the GNU Lesser General Public License along with
//// hpp-affordance  If not, see
//// <http://www.gnu.org/licenses/>.


#ifndef HPP_INTERSECT_FWD_HH
#define HPP_INTERSECT_FWD_HH

#include <vector>
#include <map>
#include <boost/smart_ptr.hpp>
#include <hpp/fcl/fwd.hh>
#include <hpp/fcl/BVH/BVH_model.h>
#include <hpp/fcl/math/vec_3f.h>

namespace hpp {
      namespace intersect {
          typedef fcl::BVHModel<fcl::OBBRSS> BVHModelOB;
          typedef boost::shared_ptr<BVHModelOB> BVHModelOB_Ptr_t;
          typedef boost::shared_ptr<const BVHModelOB> BVHModelOBConst_Ptr_t;
          typedef std::pair <fcl::CollisionObjectPtr_t, fcl::CollisionObjectPtr_t>
              CollisionPair_t;

      } // namespace intersect
} // namespace hpp

#endif // HPP_INTERSECT_FWD_HH
          
