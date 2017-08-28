// created by Aron Fiechter on 2017-07-13.
// This file implements a model of the concept EnvelopeTraits_3.
// It is mostly is a copy of L2_voronoi_traits_2, except that it is for segments
// instead of points.

// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Envelope_3/include/CGAL/Env_plane_traits_3.h $
// $Id: Env_plane_traits_3.h 51989 2009-09-21 10:55:53Z efif $
//
// Author(s)     : Ophir Setter

#ifndef CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H
#define CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/number_utils.h>
#include <CGAL/Envelope_3/Envelope_base.h>

/* to convert from Alg to Rat and viceversa */
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Cartesian.h>

/* to compute bisector of two segments */
#include <CGAL/ch_akl_toussaint.h> // convex hull

/* to format strings */
#include <boost/format.hpp>

namespace CGAL {

template <class Conic_traits_2, class Kernel_>
class L2_segment_voronoi_traits_2 : public Conic_traits_2 {

public:
  typedef Kernel_                                     Kernel;
  typedef Conic_traits_2                              C_traits_2;
  typedef L2_segment_voronoi_traits_2<C_traits_2, Kernel> Self;

  typedef typename C_traits_2::Point_2                Point_2;
  typedef typename C_traits_2::Curve_2                Curve_2;
  typedef typename C_traits_2::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename C_traits_2::Multiplicity           Multiplicity;

  typedef typename C_traits_2::Rat_kernel             Rat_kernel;
  typedef typename C_traits_2::Alg_kernel             Alg_kernel;
  typedef typename C_traits_2::Nt_traits              Nt_traits;

  typedef typename Rat_kernel::FT                     Rational;
  typedef typename Rat_kernel::Point_2                Rat_point_2;
  typedef typename Rat_kernel::Segment_2              Rat_segment_2;
  typedef typename Rat_kernel::Line_2                 Rat_line_2;
  typedef typename Rat_kernel::Ray_2                  Rat_ray_2;
  typedef typename Rat_kernel::Vector_2               Rat_vector_2;
  typedef typename Rat_kernel::Direction_2            Rat_direction_2;
  typedef typename CGAL::Polygon_2<Rat_kernel>        Rat_polygon_2;
  typedef typename Rat_polygon_2::Edge_const_iterator Edge_iterator;

  typedef typename Alg_kernel::FT                     Algebraic;
  typedef typename Alg_kernel::Point_2                Alg_point_2;
  typedef typename Alg_kernel::Segment_2              Alg_segment_2;
  typedef typename Alg_kernel::Line_2                 Alg_line_2;
  typedef typename Alg_kernel::Ray_2                  Alg_ray_2;
  typedef typename Alg_kernel::Direction_2            Alg_direction_2;

  typedef Rational                                    RT;
  typedef Algebraic                                   AT;

  /* Converters */
  typedef CGAL::Cartesian_converter<Alg_kernel, Rat_kernel> AK_to_RK;
  typedef CGAL::Cartesian_converter<Rat_kernel, Alg_kernel> RK_to_AK;
  typedef CGAL::Cartesian<double>                     D_kernel;
  typedef typename D_kernel::Point_2                  D_point_2;
  typedef CGAL::Cartesian_converter<Alg_kernel, D_kernel> AK_to_DK;
  typedef CGAL::Cartesian_converter<D_kernel, Rat_kernel> DK_to_RK;

  /* Define segments as surfaces for simplicity. Each segment should be thought
   * of as a surface representing the function distance, with higher values on
   * the z axis mean farthest points */
  typedef Rat_segment_2                               Surface_3;
  typedef Surface_3                                   Xy_monotone_surface_3;

protected:
  typedef std::pair<X_monotone_curve_2, Multiplicity> Intersection_curve;

private:
  typedef typename std::pair<
    std::pair<Rat_line_2, Rat_line_2>,
    std::pair<Rat_line_2, Rat_line_2>
  >                                                   Rat_delimiter_lines;
  typedef typename std::pair<
    std::pair<Alg_line_2, Alg_line_2>,
    std::pair<Alg_line_2, Alg_line_2>
  >                                                   Alg_delimiter_lines;

  enum Bisector_type {
    PARABOLIC_ARC,
    SUPP_LINE_BISECTOR,
    ENDPOINT_BISECTOR
  };

  /* Parabola class used to provide methods to find intersections with lines,
   * to check point location with respect to a parabola, to create parabolic
   * arcs, to get tangent directions at points. */
  class Parabola {

  private:

    /* Fields */

    /* coefficients of equation: rx^2 + sy^2 + txy + ux + vy + w = 0 */
    RT _r; RT _s; RT _t; RT _u; RT _v; RT _w;

    /* generators of parabola, direction of parabola is the same of directrix;
     * also store orientation for when we build arcs */
    Rat_line_2 _directrix;
    Rat_point_2 _focus;
    Orientation _orientation;

  public:

    /* Empty constructor */
    Parabola() {}

    /* Construct using directrix and focus. Details on the computation of the
     * forumla can be found in doc/parabola.pdf */
    Parabola(Rat_line_2 directrix, Rat_point_2 focus)
      : _directrix(directrix), _focus(focus) {

      /* orientation depends on directrix an focus */
      if (directrix.has_on_positive_side(focus)) {
        this->_orientation = CGAL::COUNTERCLOCKWISE;
      }
      else {
        this->_orientation = CGAL::CLOCKWISE;
      }

      /* extract parameters to compute coefficients */
      RT a = directrix.a();
      RT b = directrix.b();
      RT c = directrix.c();
      RT f_x = focus.x();
      RT f_y = focus.y();
      RT TWO = RT(2);

      /* compute coefficients, see details in doc/parabola.pdf */
      RT r = -CGAL::square(b);
      RT s = -CGAL::square(a);
      RT t = TWO * a * b;
      RT u =
        TWO * a * c +
        TWO * CGAL::square(a) * f_x +
        TWO * CGAL::square(b) * f_x
      ;
      RT v =
        TWO * b * c +
        TWO * CGAL::square(a) * f_y +
        TWO * CGAL::square(b) * f_y
      ;
      RT w =
        CGAL::square(c) -
        CGAL::square(a) * CGAL::square(f_x) -
        CGAL::square(a) * CGAL::square(f_y) -
        CGAL::square(b) * CGAL::square(f_x) -
        CGAL::square(b) * CGAL::square(f_y)
      ;

      this->_r = r;
      this->_s = s;
      this->_t = t;
      this->_u = u;
      this->_v = v;
      this->_w = w;

      //TODO remove
      // std::cout << "Constructed parabola: "
      //           << _r << "x^2 + "
      //           << _s << "y^2 + "
      //           << _t << "xy + "
      //           << _u << "x + "
      //           << _v << "y + "
      //           << _w
      //           << std::endl
      // ;

      RK_to_AK to_alg;
      CGAL_assertion(CGAL::square(_t) - 4 * _r * _s == 0);  // curve is parabola
      CGAL_assertion(this->has_on(to_alg(
        CGAL::midpoint(focus, directrix.projection(focus))  // origin is on p
      )));
    }

    /* Getters */
    RT r() { return _r; }
    RT s() { return _s; }
    RT t() { return _t; }
    RT u() { return _u; }
    RT v() { return _v; }
    RT w() { return _w; }
    Rat_line_2 directrix() { return _directrix; }
    Rat_point_2 focus() { return _focus; }
    Orientation orientation() { return _orientation; }

    /* Methods */

    /* Evaluate the equation of the parabola rx^2 + sy^2 + txy + ux + vy + w = 0
     * using the x and y of the point. If the result is 0, the point is on the
     * parabola, if it is positive the point lies on the positive side, if it
     * is negative the point lies on the negative side. */
    Algebraic evaluate(Point_2 point) {
      Algebraic x = point.x();
      Algebraic y = point.y();
      Algebraic result(
        this->r() * CGAL::square(x) +
        this->s() * CGAL::square(y) +
        this->t() * x * y +
        this->u() * x +
        this->v() * y +
        this->w()
      );
      //TODO remove
      // std::cout << "Evaluated point " << point
      //           << " and the result is " << result
      //           << std::endl
      // ;
      return result;
    }

    /* Check if a given point lies on the parabola by checking if the values of
     * x and y (the point's coordinates) satisfy the equation of the parabola */
    bool has_on(Point_2 point) {
      return CGAL::is_zero(this->evaluate(point));
    }
    /* Check if a given point lies on the positive side of the parabola. The
     * positive side is the one on the left when traveling on the curve in the
     * same direction as the directrix (by construction they have the "same"
     * oriented side) */
    bool has_on_positive_side(Point_2 point) {
      return CGAL::is_positive(this->evaluate(point));
    }
    /* Check if a given point lies on the negative side of the parabola. The
     * negative side is the one on the right when traveling on the curve in the
     * same direction as the directrix (by construction they have the "same"
     * oriented side) */
    bool has_on_negative_side(Point_2 point) {
      return CGAL::is_negative(this->evaluate(point));
    }

    /* Given a point, return the tangent line at that point. The tangent can be
     * found by looking at the line that is perpendicular to the directrix and
     * goes through this point and the line that goes through the focus and this
     * point. The bisector of these two lines is the tangent at this point.
     * Precondition (checked): the point is on the parabola. */
    Alg_line_2 tangent_at_point(Alg_point_2 point) {
      /* check precondition */
      CGAL_precondition_msg(
        this->has_on(point),
        "Given point for tangent has to be on the parabola"
      );

      /* get converter and convert */
      RK_to_AK to_alg;
      Alg_line_2 alg_directrix = to_alg(this->directrix());
      Alg_point_2 focus = to_alg(this->focus());
      Alg_point_2 proj_point = alg_directrix.projection(point);

      /* tangent */
      Alg_line_2 tangent;

      /* based on an orientation test, distinguish three cases */
      switch (CGAL::orientation(proj_point, point, focus)) {

        /* special case: the point is the vertex of the parabola. In this case,
         * the tangent at point has the same direction as the directrix */
        case CGAL::COLLINEAR: {
          tangent = Alg_line_2(point, alg_directrix.direction());
          break;
        }

        /* in this case the tangent is directed towards the positive part of the
         * directrix, so we need to have one line directed from the focus to the
         * point and the other one from the directrix to the point */
        case CGAL::LEFT_TURN: {
          tangent = CGAL::bisector(
            Alg_line_2(focus, point),
            Alg_line_2(proj_point, point)
          );
          break;
        }

        /* in this case the tangent is directed towards the negative part of the
         * directrix, so we need to have one line directed from the point to the
         * focus and the other one from the point to the directrix */
        case CGAL::RIGHT_TURN: {
          tangent = CGAL::bisector(
            Alg_line_2(point, focus),
            Alg_line_2(point, proj_point)
          );
          break;
        }

        /* impossible, throw error */
        default: {
          CGAL_error_msg(
            "Point on parabola, its projection on the directrix and the focus"
            " failed to have one of three orientations."
          );
          break;
        }
      }

      /* invert tangent if necessary */
      if (!generally_same_direction(tangent, alg_directrix.direction())) {
        tangent = tangent.opposite();
      }
      return tangent;
    }

    /* Save into the OutputIterator o the intersection(s) of the parabola with
     * a given line l. The type of o must be Alg_point_2.
     * Return a past the end iterator o. */
    template <class OutputIterator>
    OutputIterator get_intersections(Rat_line_2 line, OutputIterator o) {
      /* equation of line:      ax + by + c = 0
       * equation of parabola:  rx^2 + sy^2 + txy + ux + vy + w = 0
       * we can find intersections by substituting line in parabola; described
       * in detail in docs/parabola.pdf, verified using wxMaxima */
      RT a = line.a();
      RT b = line.b();
      RT c = line.c();

      //TODO remove
      // std::cout << "Finding intersections with line: "
      //           << a << "x + "
      //           << b << "y + "
      //           << c << std::endl
      // ;

      /* convert line to algebraic, get nt_traits to solve quadratic equation */
      RK_to_AK to_alg;
      Alg_line_2 alg_line = to_alg(line);
      Nt_traits nt_traits;

      /* in this case the intersection is simpler, since we can substitute the x
       * in the parabola equation with just:    x = -c/a
       * We get a quadratic equation in y, which we can solve using CGAL */
      if (b == 0) {
        /* the quadratic equation in y is:
         *                                    s y^2
         *                     + (v - (ct / a)) y
         *      + ((rc^2 / a^2) - (cu / a) + w)
         *                                          = 0
         */
        RT EQ_A = this->s();
        RT EQ_B = this->v() - ((c * this->t()) / a);
        RT EQ_C = ((this->r() * CGAL::square(c)) / CGAL::square(a)) -
          ((c * this->u()) / a) +
          this->w()
        ;

        //TODO remove
        // std::cout << "Solving quadratic equation: "
        //           << EQ_A << "y^2 + "
        //           << EQ_B << "y + "
        //           << EQ_C
        //           << std::endl
        // ;

        /* to store the 0, 1, or 2 results, indicating the intersections.
         * For all resulting y, find the corresponding x, and add a point to
         * the OutputIterator o */
        Algebraic  ys[2];
        Algebraic * ys_end;
        int n_ys;
        ys_end = nt_traits.solve_quadratic_equation(EQ_A, EQ_B, EQ_C, ys);
        n_ys = ys_end - ys;

        //TODO remove
        // std::cout << "Found " << n_ys << " intersections." << std::endl;

        /* if no intersections return */
        if (n_ys == 0) {
          return o;
        }
        /* else find xs for all ys, add points to iterator */
        else while (--n_ys >= 0) {
          Algebraic current_y = ys[n_ys];
          Algebraic respective_x = alg_line.x_at_y(current_y);
          Alg_point_2 intersection(respective_x, current_y);
          CGAL_assertion(this->has_on(intersection));
          *o++ = intersection;
        }

        return o; // already one past the end, post-incremented when adding
      }
      /* in the general case we substitute the y in the parabola equation with
       * the value:                             y = -c/b + -ax/b
       * We get a quadratic equation in x, which we can solve using CGAL */
      else {
        /* the quadratic equation in x is:
         *                (r + (sa^2 / b^2) - (at / b)) x^2
         *   + ((2acs / b^2) - (ct / b) - (av / b) + u) x
         *              + ((sc^2 / b^2) - (cv / b) + w)
         *                                                  = 0
         */
        RT EQ_A = this->r() +
          (this->s() * CGAL::square(a) / CGAL::square(b)) -
          (a * this->t() / b)
        ;
        RT EQ_B = (2 * a * c * this->s() / CGAL::square(b)) -
          (c * this->t() / b) -
          (a * this->v() / b) +
          this->u()
        ;
        RT EQ_C = (this->s() * CGAL::square(c) / CGAL::square(b)) -
          (c * this->v() / b) +
          this->w()
        ;

        //TODO remove
        // std::cout << "Solving quadratic equation: "
        //           << EQ_A << "x^2 + "
        //           << EQ_B << "x + "
        //           << EQ_C
        //           << std::endl
        // ;

        /* to store the 0, 1, or 2 results, indicating the intersections.
         * For all resulting x, find the corresponding y, and add a point to
         * the OutputIterator o */
        Algebraic  xs[2];
        Algebraic * xs_end;
        int n_xs;
        xs_end = nt_traits.solve_quadratic_equation(EQ_A, EQ_B, EQ_C, xs);
        n_xs = xs_end - xs;

        //TODO remove
        // std::cout << "Found " << n_xs << " intersections." << std::endl;

        /* if no intersections return */
        if (n_xs == 0) {
          return o;
        }
        /* else find xs for all xs, add points to iterator */
        else while (--n_xs >= 0) {
          Algebraic current_x = xs[n_xs];
          Algebraic respective_y = alg_line.y_at_x(current_x);
          Alg_point_2 intersection(current_x, respective_y);
          CGAL_assertion(this->has_on(intersection));
          *o++ = intersection;
        }

        return o; // already one past the end, post-incremented when adding
      }
    }

    /* Given a point on the parabola and a vector of lines, finds all
     * intersections of the parabola with those lines, then looks for the next
     * point on the parabola (in the same direction as the directrix) after the
     * given start point and returns that point. */
    Alg_point_2 next_intersection(
      Alg_point_2 start,
      std::vector<Rat_line_2> delimiters
    ) {
      /* get intersections */
      std::list<Alg_point_2> intersections;
      for (auto& delimiter : delimiters) {
        this->get_intersections(delimiter, std::back_inserter(intersections));
      }
      CGAL_assertion(intersections.size() > 0);

      /* convert directrix and project start point on it */
      RK_to_AK to_alg;
      Alg_line_2 alg_directrix = to_alg(this->directrix());
      Alg_point_2 alg_start_pt = alg_directrix.projection(start);

      /* project intersections on alg_directrix, filter points that are "before"
       * alg_start_pt on the directrix */
      std::list<Alg_point_2> relevant_proj_intersections;
      for (auto& p : intersections) {
        Alg_point_2 proj_p = alg_directrix.projection(p);
        if (
          alg_directrix.direction()
          ==
          Alg_segment_2(alg_start_pt, proj_p).direction()
        ) {
          relevant_proj_intersections.push_back(proj_p);
        }
      };
      CGAL_assertion(relevant_proj_intersections.size() > 0);

      /* find next intersection on directrix, return corresponding point on the
       * parabola. Just look through the intersections and find it. */
      Alg_point_2 next_on_directrix = closest_point<Alg_kernel>(
        alg_start_pt,
        relevant_proj_intersections
      );
      Alg_point_2 result; bool assigned;
      for (auto& intersection : intersections) {
        if (alg_directrix.projection(intersection) == next_on_directrix) {
          result = intersection;
          assigned = true;
        }
      }

      CGAL_assertion_msg(assigned, "Could not find closest point");
      //TODO remove
      // std::cout << "This is the closest intersection: " << result << std::endl;
      return result;
    }

    /* Construct a parabolic arc on the parabola from point p1 to point p2.
     * Precondition (checked): p1 and p2 are on the parabola */
    Curve_2 construct_parabolic_arc(Point_2 p1, Point_2 p2) {
      /* check precondition: both points lie on the parabola */
      CGAL_assertion(this->has_on(p1));
      CGAL_assertion(this->has_on(p2));

      /* construct the curve using the parameters and the endpoints */
      Curve_2 arc(_r,_s,_t,_u,_v,_w, this->orientation(), p1, p2);
      CGAL_assertion(arc.is_valid()); // valid arc
      return arc;
    }
  }; // end of class Parabola

  /* Returns the squared distance between two points in L2 metric. */
  static Algebraic sqdistance(const Point_2& p1, const Point_2& p2) {
    return CGAL::squared_distance(p1, p2);
  }
  static Rational sqdistance(const Rat_point_2& p1, const Rat_point_2& p2) {
    return CGAL::squared_distance(p1, p2);
  }

  /* Returns the squared distance between a point and a segment in L2 metric. */
  static Algebraic sqdistance(const Point_2& p, const Rat_segment_2& s) {
    /* if segment is degenerate (a point), call other function */
    if (s.is_degenerate()) {
      RK_to_AK to_alg;
      return sqdistance(p, to_alg(s.source()));
    }

    /* find projection of p on supporting line of s */
    Alg_segment_2 alg_s(
      Alg_point_2(s.source().x(), s.source().y()),
      Alg_point_2(s.target().x(), s.target().y())
    );
    Alg_line_2 l = alg_s.supporting_line();
    Alg_point_2 proj = l.projection(p);

    /* if the projection is on s, the distance is d(p,proj) */
    if (alg_s.has_on(proj)) {
      return sqdistance(p, proj);
    }
    /* otherwise, the distance is min(d(p,s1),d(p,s2)), where s1 and s2 are the
    endpoints of the segment s */
    else {
      return CGAL::min(
        sqdistance(p, Alg_point_2(s.source().x(), s.source().y())),
        sqdistance(p, Alg_point_2(s.target().x(), s.target().y()))
      );
    }
  }


  /* Given a point p and a list of points, return the closest point.
   * Precondition (checked): the list of points is not empty */
  template <class K>
  static typename K::Point_2 closest_point(
    typename K::Point_2 p,
    std::list<typename K::Point_2> points
  ) {
    CGAL_precondition_msg(points.size() > 0, "List of points cannot be empty.");
    typename K::Point_2 result;
    typename K::FT smaller_sqdistance = -1;
    for (auto& q : points) {
      typename K::FT sqdist_pq = sqdistance(p, q);
      if (smaller_sqdistance < 0 || smaller_sqdistance > sqdist_pq) {
        result = q;
        smaller_sqdistance = sqdist_pq;
      }
    }

    CGAL_assertion(smaller_sqdistance >= 0);
    return result;
  }

  /* Construct a point in the middle of the curve cv. This function is copied
   * from Env_sphere_traits_3.h */
  static Point_2 construct_middle_point(const X_monotone_curve_2& cv) {
    Point_2 result;

    /* get the x-value of the middle point */
    Alg_point_2 mid_x = CGAL::midpoint(cv.source(),cv.target());

    /* if cv is vertical, it is just a segment */
    if (cv.is_vertical()) result = Point_2(mid_x);
    /* otherwise take the point with the same x coordinate but on cv */
    else result = Point_2(cv.point_at_x(mid_x));

    return result;
  }

  /* Convert the Curve_2 cv into multiple X_monotone_curve_2 using the provided
   * make_x_monotone function. Store the results into the list x_mono_curves.
   * Precondition (checked): cv is a valid curve. */
  static void make_curve_2_into_many_x_monotone_curve_2(Curve_2& cv,
    std::vector<X_monotone_curve_2>& x_mono_curves) {
    /* check precondition */
    CGAL_precondition_msg(
      cv.is_valid(),
      "The given curve cannot be converted to X_monotone parts because it is "
      "not valid"
    );

    /* instantiate traits, we need the provided function */
    C_traits_2 c_traits;
    typename C_traits_2::Make_x_monotone_2 make_x_monotone =
      c_traits.make_x_monotone_2_object();

    /* call the provided function */
    std::vector<CGAL::Object> pre_x_mono_curves;
    make_x_monotone(cv, std::back_inserter(pre_x_mono_curves));

    /* cast all CGAL::Objects into X_monotone_curve_2 and add to list */
    for(size_t i = 0; i < pre_x_mono_curves.size(); i++ ) {
      X_monotone_curve_2 curr;
      bool check = CGAL::assign(curr, pre_x_mono_curves[i]);
      assert(check); CGAL_USE(check);
      x_mono_curves.push_back(curr);
    }
  }

  /* Check if a given segment called edge connects two the other segments s1 and
   * s2 by any of their endpoints.
   * Return true if this is the case, return false if edge is actyally just
   * connecting s1's endpoints (or s2's). We cannot just check for equality
   * because edge could be just the same as one segment but in the other
   * direction. */
  static bool edge_connects_segments(Rat_segment_2 edge, Rat_segment_2 s1,
    Rat_segment_2 s2) {
    /* create a copy of edge but in the other direction, then check equality for
     * both versions of edge */
    Rat_segment_2 rev_edge(edge.target(), edge.source());
    return !(edge == s1 || edge == s2 || rev_edge == s1 || rev_edge == s2);
  }

  /* Given a bisector finds the point that is the furthest intersection
   * (following the direction of the bisector) of the bisector with the four
   * lines saved in delimiters.
   * Return the unbounded ray. */
  static Rat_ray_2 find_unbounded_ray(
    Rat_line_2 bisector,
    Rat_delimiter_lines delimiters
  ) {
    /* get the four intersection points, add them to two lists, sort one by x
     * and the other one by y */
    Rat_point_2 p1; // intersection between bisector and delimiters[0][0]
    Rat_point_2 p2; // intersection between bisector and delimiters[0][1]
    Rat_point_2 p3; // intersection between bisector and delimiters[1][0]
    Rat_point_2 p4; // intersection between bisector and delimiters[1][1]
    CGAL::assign(p1, CGAL::intersection(bisector, delimiters.first.first));
    CGAL::assign(p2, CGAL::intersection(bisector, delimiters.first.second));
    CGAL::assign(p3, CGAL::intersection(bisector, delimiters.second.first));
    CGAL::assign(p4, CGAL::intersection(bisector, delimiters.second.second));
    std::vector<Rat_ray_2> intersections_x = {
      Rat_ray_2(p1, bisector.direction()),
      Rat_ray_2(p2, bisector.direction()),
      Rat_ray_2(p3, bisector.direction()),
      Rat_ray_2(p4, bisector.direction())
    };
    std::vector<Rat_ray_2> intersections_y;
    std::copy(intersections_x.begin(), intersections_x.end(),
              std::back_inserter(intersections_y));
    std::sort(intersections_x.begin(), intersections_x.end(),
      [](Rat_ray_2 a, Rat_ray_2 b) {
      return a.source().x() < b.source().x();
    });
    std::sort(intersections_y.begin(), intersections_y.end(),
      [](Rat_ray_2 a, Rat_ray_2 b) {
      return a.source().y() < b.source().y();
    });

    /* find the farthest point according to the direction of bisector */
    Rat_direction_2 dir = bisector.direction();
    if (dir.dx() == 0) {
      return (dir.dy() > 0) ? intersections_y.back() : intersections_y.front();
    }
    else {
      return (dir.dx() > 0) ? intersections_x.back() : intersections_x.front();
    }
  }


  /* Given a direction and a start point finds the point that is the first
   * intersection (after this start point with the four lines saved in the
   * vector of delimiter lines.
   * Return the found intersection point. */
  static Alg_point_2 find_next_intersection(
    Alg_direction_2 direction,
    Alg_point_2 start_pt,
    std::vector<Rat_line_2> delimiters
  ) {
    /* converter */
    RK_to_AK to_alg;
    /* list to store intersections */
    std::list<Alg_point_2> intersections;

    /* for each delimiter add the intersection with ray, if it's not start_pt */
    Alg_ray_2 ray(start_pt, direction);
    for (auto& delimiter : delimiters) {
      if (CGAL::do_intersect(to_alg(delimiter), ray)) {
        Alg_point_2 intersection;
        CGAL::assign(intersection, CGAL::intersection(to_alg(delimiter), ray));
        if (intersection != start_pt) {
          intersections.push_back(intersection);
        }
      }
    }

    CGAL_assertion_msg(intersections.size() > 0, "There must be intersections");

    /* all intersections are in the correct direction because we used a ray
     * starting from start_pt, so return the closest one */
    return closest_point<Alg_kernel>(start_pt, intersections);
  }

  /* Determine the position of the point p relative to the segments s1 and s2.
   * We do not have the segments though: we have four lines, each orthogonal to
   * one endpoint of one of the two segments. Every line is oriented so to have
   * the inner part of the segment on their negative side (right side).
   * There are three main cases (described by enum Bisector_type):
   * - PARABOLIC_ARC: when p is closer to one segment's inner part and to one of
   *   the other segment's endpoints. In this case, save the supporting_line of
   *   the first segment and the endpoint of the second segment.
   * - SUPP_LINE_BISECTOR: when p is closer to both inner parts of both the two
   *   segments. In this case save the two supporting_lines of the two segments.
   * - ENDPOINT_BISECTOR: when p is closer to two endpoints of the two segments.
   *   In this case save those two endpoints.
   * In all three cases we save in o1 the correct endpoint or supporting_line of
   * s1, and in o2 the same for s2.
   * Precondition (checked): the point p is not on any of the four delimiters
   */
  static Bisector_type find_position(
    Alg_point_2 p,
    Alg_delimiter_lines delimiter_lines,
    Rat_segment_2 s1,
    Rat_segment_2 s2,
    Object& o1,  // to store [directrix1 or focus1]/line1/point1
    Object& o2   // to store [directrix2 or focus2]/line2/point2
  ) {
    /* check precondition */
    CGAL_precondition(!delimiter_lines.first.first.has_on(p));
    CGAL_precondition(!delimiter_lines.first.second.has_on(p));
    CGAL_precondition(!delimiter_lines.second.first.has_on(p));
    CGAL_precondition(!delimiter_lines.second.second.has_on(p));

    /* assume point is not on any delimiter, consider all other cases. To do so,
     * first determine what must be stored in o1, then in o2. Save in two flags
     * information about the case. In the end, determine the case. */
    bool o1_is_line = false;
    bool o2_is_line = false;

    /* determine o1 */
    if (delimiter_lines.first.first.has_on_positive_side(p)) {
      o1 = CGAL::make_object(s1.source());
    }
    else if (delimiter_lines.first.second.has_on_positive_side(p)) {
      o1 = CGAL::make_object(s1.target());
    }
    else {
      CGAL_assertion(
        delimiter_lines.first.first.has_on_negative_side(p)
        &&
        delimiter_lines.first.second.has_on_negative_side(p)
      );
      o1 = CGAL::make_object(s1.supporting_line());
      o1_is_line = true;
    }

    /* determine o2 */
    if (delimiter_lines.second.first.has_on_positive_side(p)) {
      o2 = CGAL::make_object(s2.source());
    }
    else if (delimiter_lines.second.second.has_on_positive_side(p)) {
      o2 = CGAL::make_object(s2.target());
    }
    else {
      CGAL_assertion(
        delimiter_lines.second.first.has_on_negative_side(p)
        &&
        delimiter_lines.second.second.has_on_negative_side(p)
      );
      o2 = CGAL::make_object(s2.supporting_line());
      o2_is_line = true;
    }

    /* determine case using flags */
    if (o1_is_line) {
      return (o2_is_line) ? SUPP_LINE_BISECTOR : PARABOLIC_ARC;
    }
    else {
      return (o2_is_line) ? PARABOLIC_ARC : ENDPOINT_BISECTOR;
    }
  }

  /* Given a line and a direction determine whether the line is oriented in that
   * genreal direction, that is in a range of [-90˚, 90˚] around the Given
   * direction. */
  static bool generally_same_direction(Alg_line_2 line, Alg_direction_2 dir) {
    return line.direction().counterclockwise_in_between(
      dir.vector().perpendicular(CGAL::CLOCKWISE).direction(),
      dir.vector().perpendicular(CGAL::COUNTERCLOCKWISE).direction()
    ); //TODO what if they are perpendicular? which case is it?
  }



public:

  class Make_xy_monotone_3 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Surface_3& s,
                                bool /* is_lower */,
                                OutputIterator o) const {
      /* the surfaces we are considering are distance functions from line
       * segments and are already xy_monotone because there is only one possible
       * distance value for any point on the plane */
      *o++ = s; // just insert the surface in o, return o one past the end
      return o;
    }
  };

  Make_xy_monotone_3 make_xy_monotone_3_object() const {
    return Make_xy_monotone_3();
  }



  class Construct_projected_boundary_2 {
  public:
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const {
      // /* the surfaces we are considering are distance functions of line
      //  * segments and are infinite, so they have no projected boundary */
      // return o; // the iterator remains empty

      /* save boundary for intersection with it */
      RT far_l = 10000;
      std::vector<Rat_segment_2> border = {
        Rat_segment_2(Rat_point_2(-far_l, -far_l), Rat_point_2(far_l, -far_l)),
        Rat_segment_2(Rat_point_2(far_l, -far_l), Rat_point_2(far_l, far_l)),
        Rat_segment_2(Rat_point_2(far_l, far_l), Rat_point_2(-far_l, far_l)),
        Rat_segment_2(Rat_point_2(-far_l, far_l), Rat_point_2(-far_l, -far_l))
      };

      /* The surfaces representing distance functions are infinite, but to make
       * this work with bounded bisectors we need a boundary. We use the same
       * for all segments */
      for (auto& seg : border) {
        X_monotone_curve_2 x_seg = X_monotone_curve_2(seg);
        *o++ = CGAL::make_object(std::make_pair(x_seg, CGAL::ON_NEGATIVE_SIDE));
      }
      return o;
    }
  };

  Construct_projected_boundary_2
  construct_projected_boundary_2_object() const
  {
    return Construct_projected_boundary_2();
  }


/* ########################################################################## */
/* ###            MAIN PART: COMPUTE BISECTOR OF TWO SEGMENTS             ### */
/* ########################################################################## */

  class Construct_projected_intersections_2 {
  public:
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const {

      //TODO remove
      std::cout << "Finding bisector of s1 = "
                << s1 << " and s2 = " << s2 << " --- "
      ;

      /* create converter functors to convert from:
       * - Rational to Algebraic
       * - Algebraic to Cartesian<double>
       * - Cartesian<double> to Rational
       * The last two are used together to convert by approximation from
       * Algebraic to Rational */
      RK_to_AK to_alg;
      AK_to_DK to_dbl;
      DK_to_RK to_rat;

      /* this is a very bad solution, and will need to be changed, for example
       * by using Arr_algebraic_segment_traits_2 that supports both unbounded
       * and bounded curves (also left-/right-unbounded) */
      /* save boundary for intersection with it */
      RT far_l = 10000;
      std::vector<Rat_segment_2> border = {
        Rat_segment_2(Rat_point_2(-far_l, -far_l), Rat_point_2(far_l, -far_l)),
        Rat_segment_2(Rat_point_2(far_l, -far_l), Rat_point_2(far_l, far_l)),
        Rat_segment_2(Rat_point_2(far_l, far_l), Rat_point_2(-far_l, far_l)),
        Rat_segment_2(Rat_point_2(-far_l, far_l), Rat_point_2(-far_l, -far_l))
      };

      /* if the two segments are the same (also if one is just the other but
       * reversed), their distance function is the same, so there is no
       * intersection */
      if (s1 == s2 || s1 == Rat_segment_2(s2.target(), s2.source())) {
        return o;
      }
      /* if one of the segments is degenerate, the bisector is a parabola and
       * two rays, if instead they are both degenerate (that is, they are two
       * two points) the bisector is a line */
      else if (s1.is_degenerate() || s2.is_degenerate()) {
        /* line */
        if (s1.is_degenerate() && s2.is_degenerate()) {
          Rat_line_2 bisector = CGAL::bisector(s1.source(), s2.source());
          Rat_point_2 start, end;
          bool assigned_start = false, assigned_end = false;
          for (auto& segment : border) {
            if (CGAL::do_intersect(bisector, segment)) {
              if (assigned_start && !assigned_end) {
                assigned_end = true;
                CGAL_assertion_msg(
                  CGAL::assign(end, CGAL::intersection(bisector, segment)),
                  "Could not assing end."
                );
              }
              else if (!assigned_start && !assigned_end) {
                assigned_start = true;
                CGAL_assertion_msg(
                  CGAL::assign(start, CGAL::intersection(bisector, segment)),
                  "Could not assing start."
                );
              }
            }
          }
          Rat_segment_2 seg_bisector(start, end);
          X_monotone_curve_2 x_seg_bisector = X_monotone_curve_2(seg_bisector);
          *o++ = CGAL::make_object(
            Intersection_curve(x_seg_bisector, 0)
          );
        }
        /* parabolic arc and two rays (or in other special cases just lines) */
        else {
          /* determine which one is the non-degenerate segment */
          bool s1_is_degenerate = s1.is_degenerate();
          CGAL_assertion_msg(
            (s1_is_degenerate != s2.is_degenerate()),
            "One segment should be degenerate and the other one not."
          );

          /* get directrix and focus */
          Rat_segment_2 directrix_generator = s1_is_degenerate ? s2 : s1;
          Rat_line_2 directrix = directrix_generator.supporting_line();
          Rat_point_2 focus = s1_is_degenerate ? s1.source() : s2.source();

          /* if segment and point collinear, the bisector is just the bisector
           * of the closest endpoint and the point */
          if (CGAL::collinear(
            focus,
            directrix_generator.source(),
            directrix_generator.target())
          ) {
            Rat_line_2 bisector;

            /* if the point is not on the segment just get the bisector */
            if (!directrix_generator.has_on(focus)) {

              /* get bisector */
              if (
                sqdistance(focus, directrix_generator.source())
                <
                sqdistance(focus, directrix_generator.target())
              ) {
                bisector = CGAL::bisector(focus, directrix_generator.source());
              }
              else {
                bisector = CGAL::bisector(focus, directrix_generator.target());
              }
            }
            /* otherwise the bisector is the line on orthogonal to the segment
             * that goes through the point on it */
            else {
              bisector = directrix.perpendicular(focus);
            }

            /* find endpoints */ //TODO abstract this part
            Rat_point_2 start, end;
            bool assigned_start = false, assigned_end = false;
            for (auto& segment : border) {
              if (CGAL::do_intersect(bisector, segment)) {
                if (assigned_start && !assigned_end) {
                  assigned_end = true;
                  CGAL_assertion_msg(
                    CGAL::assign(end, CGAL::intersection(bisector, segment)),
                    "Could not assing end."
                  );
                }
                else if (!assigned_start && !assigned_end) {
                  assigned_start = true;
                  CGAL_assertion_msg(
                    CGAL::assign(
                      start,
                      CGAL::intersection(bisector, segment)
                    ),
                    "Could not assing start."
                  );
                }
              }
            }
            Rat_segment_2 sg_bisector(start, end);
            X_monotone_curve_2 x_sg_bisector = X_monotone_curve_2(sg_bisector);
            *o++ = CGAL::make_object(
              Intersection_curve(x_sg_bisector, 0)
            );
          }
          /* otherwise, we have the actual parabola and the two rays */
          else {
            /* ... */
          }

          /* make parabola */
          // *o++ = CGAL::make_object(
          //   Intersection_curve(curve_seg, 0)
          // );
        }
      }
      /* if the two segments are not the same, compute all parts of their plane
       * bisector */
      else {
        /* first of all, for each segment create the two lines that divide the
         * plane in three areas: one of all points closest to the inner part of
         * the segment, the other two of all points closest to the two endpoints
         * of the segment.
         * The lines are saved with an orientation such that they both have the
         * inner part of the segment on their right side (negative side).
         * Note: a vector constructed using a segment is oriented from source to
         * target of that segment, so to build a line such that the segment lies
         * on the right side of it, we need to use:
         * - the source of the segment and as direction the vector oriented 90
         *   degrees counterclockwise from the segment vector
         * - the target of the segment but as direction the vector oriented 90
         *   degrees clockwise. */
        Rat_delimiter_lines delimiter_lines = {
          {
            Rat_line_2(
              s1.source(),
              Rat_vector_2(s1).perpendicular(CGAL::COUNTERCLOCKWISE)
            ),
            Rat_line_2(
              s1.target(),
              Rat_vector_2(s1).perpendicular(CGAL::CLOCKWISE)
            )
          },
          {
            Rat_line_2(
              s2.source(),
              Rat_vector_2(s2).perpendicular(CGAL::COUNTERCLOCKWISE)
            ),
            Rat_line_2(
              s2.target(),
              Rat_vector_2(s2).perpendicular(CGAL::CLOCKWISE)
            )
          }
        };
        /* also save the lines in alg form for convenience */
        Alg_delimiter_lines alg_delimiter_lines = {
          {
            to_alg(delimiter_lines.first.first),
            to_alg(delimiter_lines.first.second)
          },
          {
            to_alg(delimiter_lines.second.first),
            to_alg(delimiter_lines.second.second)
          }
        };
        /* also save the lines in a vector for convenience */
        std::vector<Rat_line_2> delimiter_lines_vector = {
          delimiter_lines.first.first, delimiter_lines.first.second,
          delimiter_lines.second.first, delimiter_lines.second.second
        };
        /* also save segment endpoints "generating" these lines */
        std::vector<Rat_point_2> segment_endpoints = {
          s1.source(), s1.target(), s2.source(), s2.target()
        };

        /* then compute the 2 or 4 unbounded edges of the bisector.
         * To do this, first compute the convex hull of the endpoints of the
         * segments. The pairs of vertices of the hull that are not of the same
         * segment are the pairs of which the bisector lines contain the
         * unbounded rays that are the unbounded rays of the plane bisector of
         * the two segments */

        /* compute hull of endpoints; the hull points will be stored inside
         * ch_points in counterclockwise order */
        std::list<Rat_point_2> ch_points;
        std::list<Rat_point_2> points = {
          s1.source(), s1.target(), s2.source(), s2.target()
        };
        CGAL::ch_akl_toussaint(
          points.begin(), points.end(),
          std::back_insert_iterator<std::list<Rat_point_2>>(ch_points)
        );

        /* make a polygon out of the hull points, iterate over vertices to find
         * pairs to make rays, directed towards outside of polygon */
        Rat_polygon_2 ch_polygon(ch_points.begin(), ch_points.end());
        CGAL_assertion(ch_polygon.is_convex()); // it is a hull
        CGAL_assertion(ch_polygon.area() >= 0); // it is counterclockwise

        /* list to save the unbounded rays of the bisector */
        std::list<Rat_ray_2> unbounded_ray_list;

        for ( // for all edges
          Edge_iterator eit = ch_polygon.edges_begin();
          eit != ch_polygon.edges_end();
          ++eit
        ) {
          if (edge_connects_segments(*eit, s1, s2)) { // if it's not s1 or s2
            /* create line that bisects the segment, orient it outside */
            Rat_line_2 bisector_line = CGAL::bisector(
              eit->target(), eit->source()
            );
            /* find "farthest" intersection with delimiters, ray starts there */
            Rat_ray_2 unbounded_ray = find_unbounded_ray(
              bisector_line, delimiter_lines
            );
            unbounded_ray_list.push_back(unbounded_ray);

            /* make very long segment to represent an unbounded ray, so that it
             * can be saved as an X_monotone_curve_2, because the Conic_traits
             * require that curves are bounded */
            Rat_point_2 start_point = unbounded_ray.source();
            Rat_point_2 end_point;
            bool assigned = false;
            for (auto& seg : border) {
              if (assigned) break;
              else if (CGAL::do_intersect(unbounded_ray, seg) && !assigned) {
                assigned = true;
                CGAL_assertion_msg(
                  CGAL::assign(
                    end_point,
                    CGAL::intersection(unbounded_ray, seg)
                  ),
                  "Could not assign end."
                );
              }
            }

            CGAL_assertion_msg(assigned, "Could not find ray end_point.");
            Rat_segment_2 segment(start_point, end_point);
            X_monotone_curve_2 segment_curve(segment); // it's just straight
            *o++ = CGAL::make_object(
              Intersection_curve(segment_curve, 0)
            );
          }
        }

        /* if the two segments do NOT intersect, construct the bisector starting
         * from one unbounded edge, finding the correct intersection points
         * using the delimiter_lines.
         * In this case, the ray start points should be only two. */
        if (!CGAL::do_intersect(s1, s2)) { // segments do not intersect
          CGAL_assertion(unbounded_ray_list.size() == 2);

          /* starting from the source of one unbounded ray and finishing at the
           * source of the other, compute the rest of the bisector, consisting
           * of:
           * - parabolic arcs: when we are in the "area of influence" of the
           *   interior of a segment and of one endpoint of the other
           * - segments: when we are in the "area of influence" of the interiors
           *   of the two segments or of two endpoints of the two segments */
          Rat_ray_2 start_ray = unbounded_ray_list.front();
          Rat_ray_2 end_ray = unbounded_ray_list.back();

          Alg_point_2 start_pt = to_alg(start_ray.source());
          Alg_point_2 end_pt = to_alg(end_ray.source());

          Alg_direction_2 curr_direction = to_alg(
            - start_ray.direction()
          );

          /* call big private function that iteratively constructs the parts of
           * the plane bisector of the segments s1 and s2 starting from a start
           * point going in a given direction and finishing at an end point */
          o = this->construct_bisector_from_point_to_point(
            s1, s2,                 // the two segments
            o,                      // OutputIterator
            start_pt, end_pt,       // construct bisector from start to end
            curr_direction,         // initial direction, updated
            alg_delimiter_lines,    // delimiter lines of s1 and s2
            delimiter_lines_vector  // same but as vector and in rational
          );
        } // end of segments do not intersect

        /* if instead they do intersect, assert it, then proceed to computing
         * the bisector in this case.
         * In this case, the ray start points should be four, but only if the
         * intersection is not by one or two endpoints (weak intersection) */
        else {
          CGAL_assertion(CGAL::do_intersect(s1, s2)); // they HAVE to intersect

          //TODO add check for touching segments, maybe even before (also add
          // check for three collinear points among the segments)

          /* starting from the source of one unbounded ray and finishing at the
           * source of the other, compute the rest of the bisector, consisting
           * of:
           * - parabolic arcs: when we are in the "area of influence" of the
           *   interior of a segment and of one endpoint of the other
           * - segments: when we are in the "area of influence" of the interiors
           *   of the two segments or of two endpoints of the two segments
           * Since in this case we have four rays, do this process twice. */
          CGAL_assertion(unbounded_ray_list.size() == 4);
          Rat_ray_2 start_ray_one = unbounded_ray_list.front();
          unbounded_ray_list.pop_front();
          Rat_ray_2 start_ray_two = unbounded_ray_list.front();
          unbounded_ray_list.pop_front();
          Rat_ray_2 end_ray_one = unbounded_ray_list.front();
          unbounded_ray_list.pop_front();
          Rat_ray_2 end_ray_two = unbounded_ray_list.front();
          unbounded_ray_list.pop_front();

          CGAL_assertion_msg(
            unbounded_ray_list.empty(),
            "There should only be 4 rays, and all were popped."
          );

          /* get start point, end point and direction for both cases */
          Alg_point_2 start_pt_one = to_alg(start_ray_one.source());
          Alg_point_2 end_pt_one = to_alg(end_ray_one.source());
          Alg_direction_2 curr_direction_one = to_alg(
            - start_ray_one.direction()
          );
          Alg_point_2 start_pt_two = to_alg(start_ray_two.source());
          Alg_point_2 end_pt_two = to_alg(end_ray_two.source());
          Alg_direction_2 curr_direction_two = to_alg(
            - start_ray_two.direction()
          );

          /* call big private function that iteratively constructs the parts of
           * the bisector of s1 and s2 starting from a start point going in a
           * given direction and finishing at an end point. Do this two times,
           * for both bisectors */
          o = this->construct_bisector_from_point_to_point(
            s1, s2,                   // the two segments
            o,                        // OutputIterator
            start_pt_one, end_pt_one, // construct bisector from start to end
            curr_direction_one,        // initial direction, updated
            alg_delimiter_lines,      // delimiter lines of s1 and s2
            delimiter_lines_vector    // same but as vector and in rational
          );
          o = this->construct_bisector_from_point_to_point(
            s1, s2,                   // the two segments
            o,                        // OutputIterator
            start_pt_two, end_pt_two, // construct bisector from start to end
            curr_direction_two,        // initial direction, updated
            alg_delimiter_lines,      // delimiter lines of s1 and s2
            delimiter_lines_vector    // same but as vector and in rational
          );
        } // end of segments intersect


        //TODO remove
        std::cout << "Created bisector of s1 = "
                  << s1 << " and s2 = " << s2 << std::endl
        ;

      } // end of segments are not the same

      /* return one past the end iterator */
      return o;
    }

  private:

    /* Given a two Curve_2 parabolic arcs and an algebraic segment that
     * connects them, find an approximated segment supported by a line with
     * rational coefficients. The supporting MUST intersect both arcs.
     * The arcs are tangent to the segment to approximate at its endpoints.
     * Precondition (not checked): the arcs and the segment are oriented in the
     * correct way, so that they are directed from the source of prev_arc to the
     * target of next_arc.
     * Return this supporting line.
     */
    Rat_line_2 get_approximated_inner_segment_supporting_line(
      Alg_segment_2& segment,
      Curve_2& prev_arc,
      Curve_2& next_arc
    ) const {

      std::cout << "\n\n____________________________________________________\n"
                << "Approximating algebraic segment: ("
                << segment.source().x().toString() << ", "
                << segment.source().y().toString()
                << ") --> ("
                << segment.target().x().toString() << ", "
                << segment.target().y().toString()
                << ")" << ":"
      ;

      //TODO this is very fake but might work in some cases
      AK_to_DK to_dbl;
      DK_to_RK to_rat;

      /* determine if we have to rotate or translate: if the orientation of the
       * previous and of the next arc is the same, we have to translate the
       * algebraic segment up or down, if the orientations are different instead
       * we have to slightly rotate the segment.
       * This is to ensure that the supporting_line line of the approximated
       * Rational segment still intersects both curves.
       *
       * There are two main cases:
       * - the curves have the same orientation. In this case, the approximated
       *   segment needs to be slightly moved towards the "interior" of both
       *   curves, so that it definitely intersects still both of them
       * - the curves have a different orietnat. In this case, the approximated
       *   segment needs to be slightly rotated so that the source is in the
       *   "interior" of the prev_arc, and the target of the next_arc
       *
       * To approximate, we work on the single coordinates of the points, so we
       * need to create four Rational numbers.
       * According to which direction each point needs to be moved, we multiply
       * the Algebraic coordinate by a large multiple of 2 (say, 2^16), then we
       * take the ceiling or floor of this number so to get an Integer, then
       * create a Rational that is this integer divided by the large multiple of
       * 2 (again, for example 2^16).
       */
      Rat_segment_2 approximated_segment;

      /* rotate */
      if (prev_arc.orientation() != next_arc.orientation()) {
        std::cout << "[by rotation]:\n";
      }
      /* move up or down */
      else {
        std::cout << "[by translation]:\n";
      }

      approximated_segment = to_rat(to_dbl(segment));

      std::cout << "Approximated to: " << approximated_segment;
      std::cout << "\n____________________________________________________\n\n";

      return approximated_segment.supporting_line();
    }

    /* Given two Curve_2 objects, return true if they are the same */
    bool same_curves(Curve_2 cv1, Curve_2 cv2) const {
      return
        (cv1.r() == cv2.r()) &&
        (cv1.s() == cv2.s()) &&
        (cv1.t() == cv2.t()) &&
        (cv1.u() == cv2.u()) &&
        (cv1.v() == cv2.v()) &&
        (cv1.w() == cv2.w()) &&
        (cv1.orientation() == cv2.orientation()) &&
        (cv1.source() == cv2.source()) &&
        (cv1.target() == cv2.target())
      ;
    }

    /* Given three Curve_2 in sequence, update the endpoints to make themselves
     * connected, if they are not already.
     * Keep curr like it is, update the others.
     * Precondition (not checked): the three curved arcs are oriented in the
     * correct way, so that they are directed from the source of prev to the
     * target of next.
     * Precondition (not checked): the source and target of curr lie on prev and
     * on next respectively.
     */
    void update_endpoints(Curve_2& prev, Curve_2& curr, Curve_2 next) const {
      if (prev.target() != curr.source()) prev.set_target(curr.source());
      if (next.source() != curr.target()) next.set_source(curr.target());
      return;
    }

    /* Helper function for the construction of the plane bisecotr of two line
     * segments.
     * Given a start point and and end point, iteratively constructs all pieces
     * of the bisector. Other objects passed to the function are the segments
     * themselves, the initial direction where the bisector must continue from
     * start_pt, the iterator where to store the pieces of the bisector as
     * Intersection_curve objects (containing X_monotone_curve_2), and two
     * different instances of the delimiter lines, that are the lines orthogonal
     * to the segments' endpoints, such that each segment is inside the negative
     * part of his couple of delimiter lines.
     * Returns a one past the end iterator of the list of Intersection_curve.
     */
    template <class OutputIterator>
    OutputIterator construct_bisector_from_point_to_point(
      Rat_segment_2 s1, Rat_segment_2 s2,
      OutputIterator o,
      Alg_point_2 start_pt, Alg_point_2 end_pt,
      Alg_direction_2 curr_direction,
      Alg_delimiter_lines alg_delimiter_lines,
      std::vector<Rat_line_2> delimiter_lines_vector
    ) const {
      /* create converter functors to convert from:
       * - Rational to Algebraic
       * - Algebraic to Cartesian<double>
       * - Cartesian<double> to Rational
       * The last two are used together to convert by approximation from
       * Algebraic to Rational */
      RK_to_AK to_alg;
      AK_to_DK to_dbl;
      DK_to_RK to_rat;

      /* list to store all Curve_2 bisector part. They all will be converted to
       * X_monotone_curve_2 and inserted into OutputIterator o at the end of
       * this method */
      std::list<Curve_2> bisector_parts;

      /* when computing the straight parts of the bisector between two conic
       * arcs, it is likely that the starting and ending points of these
       * straight parts (segments) are algebraic points. These cannot be added
       * as Curve_2 objects because their supporting conic (a line) could have
       * algebraic coefficients, and this is not supported by the arrangement
       * traits class (Arr_conic_traits_2).
       * As a solution, keep the coefficients of the parabola of the previous
       * arc and wait until the coefficients of the parabola of the next arc are
       * available; then, just use the Curve_2 constructor that does not require
       * the exact endpoints, but the three conic coefficients lists and two
       * approximate endpoints: it computes the actual endpoints as the
       * intersections of the three curves closest to the approximate endpoints.
       * Then, get the computed endpoints of the new curve and update the
       * end point of the previous Curve_2 and the start point of the successive
       * Curve_2
       *
      //  * Keep two parts to approximate for the case in which the segment is the
      //  * part of the bisector that intersects the two segments at their
      //  * intersection (so this happens only when the two segments are
      //  * intersecting, of course).
       */
      Curve_2 prev_arc;
      bool part_to_approximate_exists = false;
      Alg_segment_2 part_to_approximate;
      // Alg_segment_2 part_to_approximate_1;
      // Alg_segment_2 part_to_approximate_2;

      /* rename start point */
      Alg_point_2 curr_pt = start_pt;

      //TODO remove
      int iteration = 0;
      std::cout << "\n####################################\n"
                << "#        Starting while loop       #"
                << "\n####################################\n"
      ;

      /* "walk" through the bisector to find all parts until every piece has
       * been created and added to the OutputIterator o */
      while (curr_pt != end_pt) {
        //TODO remove
        std::cout << iteration++ << ": ";

        /* find next intersection with delimiter_lines when going in the
         * direction saved in "curr_direction", then find a middle point
         * between curr_pt and that intersection */
        Alg_point_2 approximate_next_intersection = find_next_intersection(
          curr_direction, curr_pt, delimiter_lines_vector
        );
        Alg_point_2 midpoint = CGAL::midpoint(
          curr_pt,
          approximate_next_intersection
        );

        /* to store the true next intersection and the next direction */
        Alg_point_2 actual_next_intersection;
        Alg_direction_2 next_direction;

        /* determine where this middle point is relative to the two segments
         * s1 and s2, and create the correct piece of the bisector. The
         * objects o1 and o2 that are passed will store in the cases:
         * - PARABOLIC_ARC:       o1 = focus/directrix  o2 = focus/directrix
         * - SUPP_LINE_BISECTOR:  o1 = supp_line1,      o2 = supp_line2
         * - ENDPOINT_BISECTOR:   o1 = endpoint_1,      o2 = endpoint_2   */
        Object o1, o2;
        switch (find_position(midpoint, alg_delimiter_lines, s1, s2, o1, o2)) {

          case PARABOLIC_ARC: {
            std::cout << "- PARABOLIC_ARC - ";
            /* extract directrix and focus */
            Rat_line_2 directrix; Rat_point_2 focus;
            if (CGAL::assign(directrix, o1)) {
              CGAL_assertion(CGAL::assign(focus, o2));
            }
            else {
              CGAL_assertion(CGAL::assign(focus, o1));
              CGAL_assertion(CGAL::assign(directrix, o2));
            }

            /* keep or invert directrix based on curr_direction */
            if (!generally_same_direction(
              to_alg(directrix), curr_direction
            )) {
              directrix = directrix.opposite();
            }

            /* create parabola */
            Parabola supporting_conic(directrix, focus);
            CGAL_assertion(supporting_conic.has_on(curr_pt));

            /* find actual next intersection of parabola */
            actual_next_intersection = supporting_conic.next_intersection(
              curr_pt, delimiter_lines_vector
            );

            /* get tangent with correctly oriented direction */
            Alg_line_2 tangent = supporting_conic.tangent_at_point(
              actual_next_intersection
            );
            next_direction = tangent.direction();

            /* get parabolic arc */
            Curve_2 this_arc = supporting_conic.construct_parabolic_arc(
              curr_pt,
              actual_next_intersection
            );

            //TODO remove
            std::cout << this_arc << '\n';

            /* deal with approximation of previous segment if necessary */
            if (part_to_approximate_exists) {
              /* get the approximated segment supporting_conic (a line) */
              Rat_line_2 approx_last_segment_line =
                get_approximated_inner_segment_supporting_line(
                part_to_approximate,
                prev_arc,
                this_arc
              );

              /* the prev_arc must be the last one in the list */
              Curve_2 prev_arc_in_list = bisector_parts.back();
              CGAL_assertion_msg(
                same_curves(prev_arc_in_list, prev_arc),
                "prev_arc should be last element in list of parts of bisector."
              );
              bisector_parts.pop_back(); // remove last curve

              /* create new segment curve */
              Curve_2 approx_last_segment_curve(
                0,
               	0,  // supporting conic is a line, so it's linear
               	0,
               	approx_last_segment_line.a(),
               	approx_last_segment_line.b(),
               	approx_last_segment_line.c(),
                CGAL::COLLINEAR,
                part_to_approximate.source(),
               	prev_arc.r(),
               	prev_arc.s(),
               	prev_arc.t(),
               	prev_arc.u(),
               	prev_arc.v(),
               	prev_arc.w(),
                part_to_approximate.target(),
               	this_arc.r(),
               	this_arc.s(),
               	this_arc.t(),
               	this_arc.u(),
               	this_arc.v(),
               	this_arc.w()
              );
              CGAL_assertion_msg(
                approx_last_segment_curve.is_valid(),
                "Created approximated segment curve is not valid"
              );

              //TODO remove
              std::cout << " ---> (approximated segment from before to: "
                        << approx_last_segment_curve << ")\n"
              ;

              /* if needed (because the endpoints might just be the same ones,
               * for example if the approximation was exact), update prev_arc
               * and this_arc end and start point to coincide with start and end
               * of this new approximated segment curve */
              update_endpoints(prev_arc, approx_last_segment_curve, this_arc);

              /* push in list of bisector parts the updated prev_arc and the
               * now approximated segment curve */
              bisector_parts.push_back(prev_arc);
              bisector_parts.push_back(approx_last_segment_curve);
            }

            /* save as Curve_2 in list of bisector parts, save this curve */
            bisector_parts.push_back(this_arc);
            prev_arc = this_arc;

            break;
          }

          case SUPP_LINE_BISECTOR: {
            std::cout << "- SUPP_LINE_BISECTOR - ";
            /* extract two supporting lines */
            Rat_line_2 supp_line1; Rat_line_2 supp_line2;
            CGAL_assertion(CGAL::assign(supp_line1, o1));
            CGAL_assertion(CGAL::assign(supp_line2, o2));

            /* orient supporting lines according to curr_direction, get
             * bisector, assert that curr_pt is on it and get direction */
            if (!generally_same_direction(
              to_alg(supp_line1), curr_direction
            )) {
              supp_line1 = supp_line1.opposite();
            }
            if (!generally_same_direction(
              to_alg(supp_line2), curr_direction
            )) {
              supp_line2 = supp_line2.opposite();
            }
            Alg_line_2 supp_line_bisector = CGAL::bisector(
              to_alg(supp_line1),
              to_alg(supp_line2)
            );
            CGAL_assertion_msg(
              supp_line_bisector.has_on(curr_pt),
              "The point curr_pt should be on the bisector, but it is not"
            );

            /* save next_direction, find actual next intersection */
            next_direction = supp_line_bisector.direction();
            actual_next_intersection = find_next_intersection(
              next_direction, curr_pt, delimiter_lines_vector
            );

            /* get bisector part, check if it intersects the two segments (this
             * happens when the segments intersect) and if yes split it in two
             * parts */
            Alg_segment_2 bisector_part(curr_pt, actual_next_intersection);
            if (CGAL::do_intersect(bisector_part, to_alg(s1))) {
              /* must also intersect s2 */
              CGAL_assertion(CGAL::do_intersect(bisector_part, to_alg(s2)));

              /* should also intersect with s2 at exactly the same point */
              Alg_point_2 segments_intersection_1, segments_intersection_2;
              CGAL_assertion(CGAL::assign(
                segments_intersection_1,
                CGAL::intersection(bisector_part, to_alg(s1))
              ));
              CGAL_assertion(CGAL::assign(
                segments_intersection_2,
                CGAL::intersection(bisector_part, to_alg(s2))
              ));
              CGAL_assertion_msg(
                (segments_intersection_1 == segments_intersection_2),
                "The bisector should intersect the two segments exactly on the "
                "intersection of the two segments."
              );

              //TODO update to smart approximation of segment version
              /* save first part as Curve_2 in list of bisector parts */
              bisector_parts.push_back(Curve_2(to_rat(to_dbl(
                Alg_segment_2(curr_pt, segments_intersection_1)
              ))));
              /* save second part as Curve_2 in list of bisector parts */
              bisector_parts.push_back(Curve_2(to_rat(to_dbl(
                Alg_segment_2(segments_intersection_1, actual_next_intersection)
              ))));
            }
            else {
              /* save this algebraic segment in the `part_to_approximate`
               * variable; when the next bisector part will be constructed (and
               * it should be an arc) or at the end of the loop (if this is the
               * last part) this segment will be approximated and added to the
               * list of pieces of the bisector
               */
              part_to_approximate = bisector_part;
              part_to_approximate_exists = true;
            }

            //TODO remove
            std::cout << bisector_part << " (not added, to approximate)" << '\n';

            break;
          }

          case ENDPOINT_BISECTOR: {
            std::cout << "- ENDPOINT_BISECTOR - ";
            /* extract two endpoints */
            Rat_point_2 endpoint1; Rat_point_2 endpoint2;
            CGAL_assertion(CGAL::assign(endpoint1, o1));
            CGAL_assertion(CGAL::assign(endpoint2, o2));

            /* create bisector, orient it according to curr_direction */
            Rat_line_2 endpoint_bisector = CGAL::bisector(endpoint1, endpoint2);
            if (!generally_same_direction(
              to_alg(endpoint_bisector), curr_direction
            )) {
              endpoint_bisector = endpoint_bisector.opposite();
            }

            /* assert curr_pt is on bisector */
            CGAL_assertion(to_alg(endpoint_bisector).has_on(curr_pt));

            /* save next_direction, find actual next intersection */
            next_direction = to_alg(endpoint_bisector).direction();
            actual_next_intersection = find_next_intersection(
              next_direction, curr_pt, delimiter_lines_vector
            );

            /* no need to approximate the segment in this case: the supporting
             * conic (a line) has rational coefficients; just save curve with
             * rational coefficients and the two alg points as endpoints,
             * orientation COLLINEAR.
             * save as Curve_2 in list of bisector parts */
            bisector_parts.push_back(Curve_2(
              0,
              0,  // supporting conic is a line, so it's linear
              0,
              endpoint_bisector.a(),
              endpoint_bisector.b(),
              endpoint_bisector.c(),
              CGAL::COLLINEAR,
              curr_pt,
              actual_next_intersection
            ));

            //TODO remove
            std::cout << bisector_parts.back() << '\n';

            break;
          }

          default: break; // should never happen
        }

        /* update current starting point and current direction of the next piece
         * of the bisector */
        curr_pt = actual_next_intersection;
        curr_direction = next_direction;
      }

      /* convert all Curve_2 to X_monotone_curve_2, add them all to the
       * OutputIterator o */
      for (auto& current_cv : bisector_parts) {
        std::vector<X_monotone_curve_2> arc_x_mono_parts;
        make_curve_2_into_many_x_monotone_curve_2(
          current_cv,
          arc_x_mono_parts
        );
        for (auto& x_mono_curve : arc_x_mono_parts) {
          *o++ = CGAL::make_object(
            Intersection_curve(x_mono_curve, 0)
            //TODO multiplicity? would need to save it in list bisector_parts
          );
        }
      }

      /* return one past the end iterator */
      return o;
    }
  }; // end of class Construct_projected_intersections_2

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const {
    return Construct_projected_intersections_2();
  }

/* ########################################################################## */
/* ###               END: COMPUTE BISECTOR OF TWO SEGMENTS                ### */
/* ########################################################################## */



  class Compare_z_at_xy_3 {
  public:
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      printf("\n ---> Compare at point\n");
      return CGAL::compare(sqdistance(p, s1), sqdistance(p, s2));
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      printf("\n ---> Compare at cv\n");
      /* compare using the middle point */
      Point_2 p = construct_middle_point(cv);
      return this->operator()(p, s1, s2);
    }

    Comparison_result operator()(const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      printf("\n ---> Compare not intersecting\n");
      /* if the two unbounded surfaces do not intersect, then they must
       * represent the same segment's distance function */
      CGAL_assertion_msg(
        (s1 == s2 || s1 == Rat_segment_2(s2.target(), s2.source())),
        "Distance function surfaces do not intersect but they are not the same"
      );
      return CGAL::EQUAL; // they are literally the same surface
    }
  };

  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  {
    printf("\n#################################\nCreated Compare_z_at_xy_3 obj#################################\n");
    return Compare_z_at_xy_3();
  }


  /* Call helper function with flag set to true */
  class Compare_z_at_xy_above_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      return compare_z_at_xy_3_helper(cv, s1, s2, true);
    }

  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3();
  }


  /* Call helper function with flag set to false */
  class Compare_z_at_xy_below_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      return compare_z_at_xy_3_helper(cv, s1, s2, false);
    }
  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const {
    return Compare_z_at_xy_below_3();
  }

  /* Helper function for Compare_z_at_xy_above_3 and Compare_z_at_xy_below_3 */
  static Comparison_result compare_z_at_xy_3_helper(
    const X_monotone_curve_2& cv,
    const Xy_monotone_surface_3& s1,
    const Xy_monotone_surface_3& s2,
    bool compare_above
  ) {
    Algebraic move_by = 10;

    /* construct a point on the curve cv, assert equidistant from s1 and s2 */
    Alg_point_2 midpoint = construct_middle_point(cv);
    Algebraic difference = sqdistance(midpoint, s1) - sqdistance(midpoint, s2);

    std::cout << std::endl << "##############################" << std::endl;
    std::cout << "# Compare s1[" << s1 << "] and s2[" << s2 << "] "
              << (compare_above ? "above" : "below") << " cv=[ "
              << cv.r() << "x^2 + "
              << cv.s() << "y^2 + "
              << cv.t() << "xy + "
              << cv.u() << "x + "
              << cv.v() << "y + "
              << cv.w() << std::endl
    ;

    std::cout << "# midpoint(" << midpoint << ") , ";

    /* print warning if necessary */
    /* colour */
    #define RESET   "\033[0m"
    #define RED     "\033[31m"      /* Red */
    #define GREEN   "\033[32m"      /* Green */
    #define YELLOW  "\033[33m"      /* Yellow */
    #define BLUE    "\033[34m"      /* Blue */
    char message[100];
    sprintf(
      message,
      RED "Warning: " YELLOW "s1 and s2 are not equidistant, difference: %lf" RESET,
      CGAL::to_double(difference)
    );
    CGAL_warning_msg((difference == 0), message);
    std::cout << "the difference is: " << difference << ".";

    /* get converter and convert */
    RK_to_AK to_alg;
    Alg_segment_2 alg_s1 = to_alg(s1);
    Alg_segment_2 alg_s2 = to_alg(s2);

    Alg_point_2 moved_point;
    Algebraic displacement = compare_above ? move_by : -move_by;
    if (!cv.is_vertical()) {
      moved_point = Alg_point_2(midpoint.x(), midpoint.y() + displacement);
    }
    else {
      std::cout << " -- CV IS VERTICAL --";
      moved_point = Alg_point_2(midpoint.x() - displacement, midpoint.y());
    }

    std::cout << " Moved midpoint to pt(" << moved_point << ")" << std::endl;

    if (sqdistance(moved_point, s1) < sqdistance(moved_point, s2)) {
      std::cout << "# Returning CGAL::SMALLER." << std::endl
                << "##############################" << std::endl
      ;
      return CGAL::SMALLER;
    }
    else {
      std::cout << "# Returning CGAL::LARGER." << std::endl
                << "##############################" << std::endl
      ;
      return CGAL::LARGER;
    }

    CGAL_error_msg("This function is not working properly");
    return CGAL::EQUAL;
  }

}; // class L2_segment_voronoi_traits_2
} // namespace CGAL

#endif // CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H
