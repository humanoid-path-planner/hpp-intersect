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
          cond = (4*evec.block<1,3>(0,0).array () * evec.block<1,3>(2,0).array () -
              evec.block<1,3>(1,0).array ().square ()).transpose ();
          
          Eigen::MatrixXd A0 (3,1);
          A0.setZero ();
          // TODO: A0 shoul always be of size 3x1 --> fix
          for (unsigned int i = 0; i < cond.size (); ++i) {
              if (cond(i) > 0.0) {
                  A0.resize (3,i+1);
                  A0.block(0,i,3,1) = evec.block<3,1>(0,i);
              }
          }
          if (A0.isZero(10e-4)) {
            std::cout << "ERROR: no fitting ellipse found." << std::endl;
          }
          // A1.rows () + T.rows () should always be equal to 6!!
          Eigen::MatrixXd A(A0.rows () + T.rows (), A0.cols ());
          A.block(0,0,A0.rows (), A0.cols ()) = A0;
          A.block(A0.rows (), 0, T.rows (), A0.cols ()) = T*A0;

          double A3 = A(3,0) - 2*A(0,0) * centroid(0) - A(1,0) * centroid(1);
          double A4 = A(4,0) - 2*A(2,0) * centroid(1) - A(1,0) * centroid(0);
          double A5 = A(5,0) + A(0,0) * centroid(0)*centroid(0) + A(2,0) * centroid(1)*centroid(1) +
              A(1,0) * centroid(0) * centroid(1) - A(3,0) * centroid(0) - A(4,0) * centroid(1);

          A(3,0) = A3/2.0;  A(4,0) = A4/2.0;  A(5,0) = A5;
          A(1,0) = A(1,0)/2.0;
          A = A/A.norm ();

          return A.block<6,1>(0,0);

        }

        std::vector<fcl::Vec3f> getIntersection (const fcl::CollisionObjectPtr_t& rom,
               const fcl::CollisionObjectPtr_t& affordance)
        {
            CollisionPair_t col = CollisionPair_t (affordance, rom);
            fcl::CollisionRequest req;
            req.num_max_contacts = 100;
            req.enable_contact = true;
            fcl::CollisionResult res;
            std::vector<fcl::Vec3f> intersectPoints; intersectPoints.clear ();
            std::size_t collision = fcl::collide (col.first.get (),
                    col.second.get (), req, res);
          
            std::cout << "rom trafo: " << rom->getRotation() << rom->getTranslation() << std::endl;
            std::cout << "aff trafo: " << affordance->getRotation() << affordance->getTranslation() << std::endl;


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

            Eigen::VectorXd ellipse = directEllipse(intersectPoints);
            if (ellipse.rows () != 6 || ellipse.cols () != 1) {
              std::cout << "ERROR: wrong number of ellipse coefficients!" << std::endl;
            }
            std::cout << "Ellipse function: " << "(" << ellipse(0) << ")*x^2 + (" << ellipse(1) 
                << ")*x*y + (" << ellipse(2) << ")*y^2 + (" << ellipse(3) << ")*x+ ("
                << ellipse(4) << ")*y + (" << ellipse(5) << ")" << std::endl;


            return intersectPoints;
        }




    } // namespace intersect
} // namespace hpp
