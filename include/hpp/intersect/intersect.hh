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
    
    /// \addtogroup intersect
    /// \{
       
        /// helper class for stacked inequalities.
        struct Inequality
        {
          Inequality(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                   const Eigen::MatrixXd& N, const Eigen::MatrixXd& V):
                   A_(A), b_(b), N_(N), V_(V) {}
          Eigen::MatrixXd A_;
          Eigen::VectorXd b_;
          Eigen::MatrixXd N_;
          Eigen::MatrixXd V_;
        };

        /// Compute radius and rotation of an elliptic or circular shape
        /// from given vector of parameters of the conic function.
        /// Rotation \param tau is given for an ellipse as the angle of its
        /// major axis from the positive X-axis. Rotation for circles is by default zero.
        /// Also the centroid of the shape is computed.
        /// \param params vector of parameters of the conic function defining a shape
        /// \param centroid the 2d-centroid of the elliptic or circular shape
        /// \param tau the rotation angle of the shape within the plane, in radians.
        std::vector<double> getRadius (const Eigen::VectorXd& params,
                Eigen::Vector2d& centroid, double& tau);

        /// \brief Fit ellipse to a set of points. 2D implementation.
        /// Assumes all points are in a plane with its normal along the Z-axis.
        /// Direct ellipse fit, proposed in article "Direct Least Squares Fitting of Ellipses"
        /// by A. W. Fitzgibbon, M. Pilu, R. B. Fisher in IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
        /// This code is based on the Matlab function DirectEllipseFit by Nikolai Chernov.
        /// It returns ellipses only, even if points are better approximated by a hyperbola.
        /// It is somewhat biased toward smaller ellipses.
        /// \param points set of points in a plane to be approximated.
        Eigen::VectorXd directEllipse(const std::vector<Eigen::Vector3d>& points);

        /// \brief Simple direct method for circle approximation based on a set of points
        /// in a plane. Assumes the plane normal points along the Z-axis.
        /// \param points set of points in a plane to be approximated.
        Eigen::VectorXd directCircle (const std::vector<Eigen::Vector3d>& points);

        /// \brief return normal of plane fitted to set of points.
        /// Modifies the vector of points by replacing the original points with those
        /// projected onto the fitted plane. Also returns the centroid of the plane.
        /// \param points set of points used for plane fitting.
        /// \param planeCentroid centroid of the plane after fitting.
        Eigen::Vector3d projectToPlane (std::vector<Eigen::Vector3d> points, Eigen::Vector3d& planeCentroid);

        /// Create a set of inequalities based on a fcl::CollisionObject. The
        /// returned matrices may be used to find out whether a point is within
        /// the collision object. Returns the inequality matrices as one object (intersect::Inequality).
        /// \param rom fcl:CollisionObject that will be used to create inequalities.
        Inequality fcl2inequalities (const fcl::CollisionObjectPtr_t& rom);

        /// Return true if a point is inside a set of planes described by intersect::Inequality.
        /// \param ineq object comprising the planes that form inequalities.
        /// \param point point to be tested against the inequalities.
        bool is_inside (const Inequality& ineq, const Eigen::Vector3d point);

        /// Get contact points resulting from collision between two fcl::CollisionObjects.
        /// Uses the fcl::collision function to verify collision but for the contact point
        /// computation, a self-implemented triangle-intersection algorithm based on that of
        /// Tomas Möller is used. This function returns a vector of 3D points that form a convex
        /// hull around all computed contact points. Triangle vertices of affordance object within
        /// the rom object are also considered contact points if they form part of the convex hull.
        /// \param rom fcl::CollisionObject that presents the reachability of a robot limb.
        /// \param affordance fcl::CollisionObject presenting the contact surface in collision with a limb.
        std::vector<Eigen::Vector3d> getIntersectionPoints (const fcl::CollisionObjectPtr_t& rom,
               const fcl::CollisionObjectPtr_t& affordance);

    /// \}
    
    } // namespace intersect
} // namespace hpp
    
#endif // HPP_INTERSECT_INTERSECT_HH

