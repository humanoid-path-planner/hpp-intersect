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
#include <hpp/intersect/intersect.hh>
#include <hpp/fcl/collision.h>

namespace hpp {
    namespace intersect {

        /// \brief get major and minor radius of ellipse or raduíus of sphere from given
        /// vector of parameters of the conic function. function modifies parameters centroid
        /// and tau to give centroid of shape and (for ellipses) the rotation angle
        /// within the plane of the shape in radians.
        std::vector<double> getRadius (const Eigen::VectorXd& params,
                Eigen::Vector2d& centroid, double& tau)
        {
          // if number of parameters == 5 -->assume it's a circle ?
          // if number == 6 --> ellipse
          bool ellipse = true;
          if (params.size () < 6) {
              if (params.size () < 5) {
                  // TODO: raise error
                  std::ostringstream oss
              ("getRadius: Wrong number of parameters in conic function!!.");
            throw std::runtime_error (oss.str ());
              }
              ellipse = false;
          }
          
          std::vector<double> radii;
          if (ellipse) {
              double A(params(0)), B(params(1)), C(params(2)), D(params(3)),
                     E(params(4)), F(params(5));

              Eigen::Matrix3d M0;
              M0 << F, D/2.0, E/2.0, D/2.0, A, B/2.0, E/2.0, B/2.0, C;
              Eigen::Matrix2d M;
              M << A, B/2.0, B/2.0, C;
        
               Eigen::EigenSolver<Eigen::Matrix2d> es(M);
               Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType eval = es.eigenvalues ();
               Eigen::Vector2d lambda;

               // make sure eigenvalues are in order for the rest of the computations
               if (fabs(eval(0).real () - A) > fabs(eval(0).real () - C)) {
                  lambda << eval(1).real (), eval(0).real ();   
               } else {
                  lambda << eval(0).real (), eval(1).real ();
               }
               radii.push_back (sqrt (-M0.determinant ()/(M.determinant () * lambda(0))));
               radii.push_back (sqrt (-M0.determinant ()/(M.determinant () * lambda(1))));
               centroid << (B*E - 2*C*D)/(4*A*C - B*B), (B*D - 2*A*E)/(4*A*C - B*B);
               tau = (M_PI/2.0 - atan((A-C)/B))/2.0;
          } else // circle!
          {
            // TODO: Add needed functionality for circles
          }
          return radii;
        }

        /// \brief fit ellipse to a set of points. 2D implementation.
        /// Direct ellipse fit, proposed in article "Direct Least Squares Fitting of Ellipses"
        /// by A. W. Fitzgibbon, M. Pilu, R. B. Fisher in IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
        /// This code is based on the Matlab function DirectEllipseFit by Nikolai Chernov.
        /// It returns ellipses only, even if points are better approximated by a hyperbola.
        /// It is somewhat biased toward smaller ellipses.
        Eigen::VectorXd directEllipse(const std::vector<fcl::Vec3f>& points)
        {
          const size_t nPoints = points.size ();
          // only consider x and y coordinates: supposing points are in a plane
          Eigen::MatrixXd XY(nPoints,2);
          // TODO: optimise
          for (unsigned int i = 0; i < nPoints; ++i) {
              XY (i,0) = points[i][0];
              XY (i,1) = points[i][1];
          }
         
          Eigen::VectorXd centroid (2);
          centroid << XY.block(0,0, nPoints,1).mean (),
                  XY.block(0,1, nPoints, 1).mean ();

          Eigen::MatrixXd D1 (nPoints,3);
          D1 << (XY.block (0,0, nPoints,1).array () - centroid(0)).square (),
             (XY.block (0,0, nPoints,1).array () - centroid(0))*(XY.block (0,1, nPoints,1).array () - centroid(1)),
             (XY.block (0,1, nPoints,1).array () - centroid(1)).square ();
          Eigen::MatrixXd D2 (nPoints,3);
          D2 << XY.block (0,0, nPoints,1).array () - centroid(0),
             XY.block (0,1, nPoints,1).array () - centroid(1),
             Eigen::MatrixXd::Ones (nPoints,1);

          Eigen::Matrix3d S1 = D1.transpose () * D1;
          Eigen::Matrix3d S2 = D1.transpose () * D2;
          Eigen::Matrix3d S3 = D2.transpose () * D2;

          Eigen::Matrix3d T = -S3.inverse () * S2.transpose ();
          Eigen::Matrix3d M_orig = S1 + S2 * T;
          Eigen::Matrix3d M; M.setZero ();
          M << M_orig.block<1,3>(2,0)/2, -M_orig.block<1,3>(1,0), M_orig.block<1,3>(0,0)/2;

          Eigen::EigenSolver<Eigen::Matrix3d> es(M);
          Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType evecCplx = es.eigenvectors ();
          //Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType evalCplx = es.eigenvalues ();

          Eigen::Matrix3d evec = evecCplx.real ();
          //Eigen::Matrix3d eval(Eigen::Matrix3d::Zero());
          //eval.diagonal () = evalCplx.real ();

          Eigen::VectorXd cond (3);
          // The condition has the form 4xz - y^2 > 0 (infinite elliptic cone) for all
          // three eigen vectors. If none of the eigen vectors fulfils the inequality,
          // the direct ellipse method fails.
          cond = (4*evec.block<1,3>(0,0).array () * evec.block<1,3>(2,0).array () -
              evec.block<1,3>(1,0).array ().square ()).transpose ();
          
          Eigen::MatrixXd A0 (0,0);
          // TODO: A0 should always be of size 3x1 --> fix
          for (unsigned int i = 0; i < cond.size (); ++i) {
              if (cond(i) > 0.0) {
                  A0.resize (3,i+1);
                  A0.block(0,i,3,1) = evec.block<3,1>(0,i);
              }
          }
          if (A0.size () < 3) {
              std::ostringstream oss
              ("directEllipse: Could not create ellipse approximation. Maybe try circle instead?");
            throw std::runtime_error (oss.str ());
          }
          // A1.rows () + T.rows () should always be equal to 6!!
          Eigen::MatrixXd A(A0.rows () + T.rows (), A0.cols ());
          A.block(0,0,A0.rows (), A0.cols ()) = A0;
          A.block(A0.rows (), 0, T.rows (), A0.cols ()) = T*A0;

          double A3 = A(3,0) - 2*A(0,0) * centroid(0) - A(1,0) * centroid(1);
          double A4 = A(4,0) - 2*A(2,0) * centroid(1) - A(1,0) * centroid(0);
          double A5 = A(5,0) + A(0,0) * centroid(0)*centroid(0) + A(2,0) * centroid(1)*centroid(1) +
              A(1,0) * centroid(0) * centroid(1) - A(3,0) * centroid(0) - A(4,0) * centroid(1);

          A(3,0) = A3;  A(4,0) = A4;  A(5,0) = A5;
          //A(3,0) = A3/2.0;  A(4,0) = A4/2.0;  A(5,0) = A5;
          //A(1,0) = A(1,0)/2.0;
          A = A/A.norm ();

          return A.block<6,1>(0,0);

        }

        /// \brief return normal of plane fitted to set of points.
        /// Also modifies the vector of points by replacing the points with those
        /// projected onto the fitted plane.
        Eigen::VectorXd projectToPlane (std::vector<fcl::Vec3f> points)
        {
          if (points.size () < 3) {
   					std::ostringstream oss
              ("projectToPlane: Too few input points to create plane.");
            throw std::runtime_error (oss.str ());
          }

          const size_t nPoints = points.size ();
          Eigen::MatrixXd A (nPoints,3);
          Eigen::VectorXd b (nPoints);
          // TODO: optimise
          for (unsigned int i = 0; i < nPoints; ++i) {
              A (i,0) = points[i][0];
              A (i,1) = points[i][1];
              A (i,2) = 1.0;
              b (i) = points[i][2];
          }
          // TODO: which decomposition to choose?
          Eigen::VectorXd params (A.fullPivLu().solve(b));

          if (params.size () < 3) {
            std::ostringstream oss
              ("projectToPlane: Wrong number of parameters for plane function!");
            throw std::runtime_error (oss.str ());
          }
          Eigen::Vector3d normal (params(0), params(1), params(2));
          normal.normalize ();
          // get origin of plane from z = C(0)*x + C(1)*y + C(2);
          Eigen::Vector3d orig(0.0, 0.0, params(2));
          Eigen::Matrix3d origin;
          origin.setZero ();
          origin.diagonal () = orig;
          Eigen::MatrixXd XYZ (nPoints,3);
          XYZ << A.block(0,0,nPoints,2), b;

          Eigen::MatrixXd distance (nPoints,3);

          distance = XYZ - (Eigen::MatrixXd::Ones (nPoints,3)*origin);//orig.array ().transpose ();

          // DEBUG:
          // std::cout << "XYZ: "<< XYZ << std::endl << ", distance :" << distance << std::endl;

          // scalar distance from point to plane along the normal for all points in vector
          Eigen::VectorXd scalarDist (nPoints);
          scalarDist = distance.block(0,0,nPoints,1)*normal(0) + distance.block(0,1,nPoints,1)*normal(1) +
              distance.block(0,2,nPoints,1)*normal(2);
  
          // TODO: optimise
          for (unsigned int i = 0; i > points.size (); ++i) {
            Eigen::Vector3d projectecPoint = XYZ.block(i,0, 1,3).transpose () - scalarDist(i)*normal;
            points[i][0] = projectecPoint (0);
            points[i][1] = projectecPoint (1);
            points[i][2] = projectecPoint (2);
          }

          return params;
        }

        std::vector<fcl::Vec3f> getIntersectionPoints (const fcl::CollisionObjectPtr_t& rom,
               const fcl::CollisionObjectPtr_t& affordance)
        {
            CollisionPair_t col = CollisionPair_t (affordance, rom);
            fcl::CollisionRequest req;
            req.num_max_contacts = 100;
            req.enable_contact = true;
            fcl::CollisionResult res;
            std::vector<fcl::Vec3f> intersectPoints; intersectPoints.clear ();
            // std::size_t collision TODO: should the return value of collision (...) be stored?
            fcl::collide (col.first.get (), col.second.get (), req, res);
          
            //std::cout << "rom trafo: " << rom->getRotation() << rom->getTranslation() << std::endl;
            //std::cout << "aff trafo: " << affordance->getRotation() << affordance->getTranslation() << std::endl;

            if (!res.isCollision ()) {
              // should not happen
              // TODO: make sure only affordance objects in contact are used to avoid this
              // throw std::runtime_error ("Affordance object is out of reach of ROM!!");
                std::cout << "WARNING: Affordance object is out of reach of ROM!! No intersection found."
                    << std::endl;
                return intersectPoints; 
            }

            for (unsigned int i = 0; i < res.numContacts (); ++i) {
                intersectPoints.push_back(res.getContact (i).pos);
            }
            
            //debug
            for (unsigned int i = 0; i < intersectPoints.size (); ++i) {
              std::cout << intersectPoints[i] << std::endl;
            }
            return intersectPoints;
        }





    } // namespace intersect
} // namespace hpp
