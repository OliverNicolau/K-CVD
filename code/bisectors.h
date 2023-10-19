#ifndef BISECTORS_H
#define BISECTORS_H

#define CGAL_SDG_VERBOSE
#undef CGAL_SDG_VERBOSE

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

/* We need assertions */
#define CGAL_DEBUG
#if defined(CGAL_NO_ASSERTIONS)
#undef CGAL_NO_ASSERTIONS
#endif
#include <CGAL/assertions.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CGAL_Ipelet_base.h>
#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Object.h>
#include <CGAL/intersections.h>
#include <CGAL/Apollonius_site_2.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>
#include <CGAL/Segment_Delaunay_graph_site_2.h>

#include <CGAL/Polychain_2.h>

#include <CGAL/Linf2D_voronoi_traits_2.h>
#include <CGAL/L2_voronoi_traits_2.h>
#include <CGAL/L2_coarse_HVD_traits_2.h>
#include <CGAL/L2_coarse_FCVD_traits_2.h>
#include <CGAL/envelope_3.h>

// added for L_2 SVD/FSVD, fn == 10, fn == 11
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/L2_segment_voronoi_traits_2.h>
// end added for L_2 SVD/FSVD

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/clustergrabber.h>

// added for L_2 k-CVD, fn == 13
//#include "k_color_VD.cpp"
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace CGAL_bisectors {

typedef CGAL::Cartesian<double>                           Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>            Delaunay;
//typedef CGAL::Regular_triangulation_euclidean_traits_2
//  <Kernel,Kernel::FT>                                     RGt;
typedef CGAL::Regular_triangulation_2<Kernel>             Regular;
//Apollonius
typedef CGAL::Apollonius_graph_traits_2<Kernel>           AT;
typedef CGAL::Apollonius_graph_2<AT>                      Apollonius;
typedef Apollonius::Site_2                                ASite;
// L_infinity bisector
struct Gt_inf
  : public CGAL::Segment_Delaunay_graph_Linf_traits_2<
    Kernel,CGAL::Field_with_sqrt_tag> {};
typedef Gt_inf::Site_2                                    Site_2;
typedef
  CGAL::SegmentDelaunayGraphLinf_2::Bisector_Linf<Gt_inf> Inf_bis;
// --------------------------------------------------------------------

Inf_bis bisector_linf;

/* Voronoi diagrams */
typedef CGAL::Exact_predicates_exact_constructions_kernel VD_Kernel;
typedef VD_Kernel::FT                                   Number_type;
typedef VD_Kernel::Iso_rectangle_2                      Iso_rectangle_2;
typedef VD_Kernel::Point_2                              VD_Point_2;
typedef VD_Kernel::Segment_2                            VD_Segment_2;
typedef std::vector<VD_Point_2>                         Points;

typedef CGAL::Linf2D_voronoi_traits_2<VD_Kernel>        VD_Traits_3;
typedef VD_Traits_3::Surface_3                          VD_Surface_3;
typedef CGAL::Envelope_diagram_2<VD_Traits_3>           VD_Envelope_diagram_2;

typedef CGAL::L2_voronoi_traits_2<VD_Kernel>            L2_VD_Traits_3;
typedef L2_VD_Traits_3::Surface_3                       L2_VD_Surface_3;
typedef CGAL::Envelope_diagram_2<L2_VD_Traits_3>        L2_VD_Envelope_diagram_2;

/* Added for L_2 SVD/FSVD, fn == 10, fn == 11 */
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef typename Rat_kernel::Segment_2                  Rat_segment_2;
typedef typename Rat_kernel::Point_2                    Rat_point_2;
typedef typename Alg_kernel::Point_2                    Alg_point_2;

typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits> Conic_traits_2;
typedef CGAL::L2_segment_voronoi_traits_2<Conic_traits_2, VD_Kernel> L2_FSVD_Traits_3;
typedef CGAL::Envelope_diagram_2<L2_FSVD_Traits_3>      L2_FSVD_Envelope_diagram_2;
typedef Conic_traits_2::X_monotone_curve_2              X_monotone_curve_2;
/* end added for L_2 SVD/FSVD */

typedef CGAL::L2_HVD_traits_2<VD_Kernel>                HVD_Traits_3;
typedef HVD_Traits_3::Surface_3                         HVD_Surface_3;
typedef CGAL::Envelope_diagram_2<HVD_Traits_3>          HVD_Envelope_diagram_2;

typedef CGAL::L2_FCVD_traits_2<VD_Kernel>               FCVD_Traits_3;
typedef FCVD_Traits_3::Surface_3                        FCVD_Surface_3;
typedef CGAL::Envelope_diagram_2<FCVD_Traits_3>         FCVD_Envelope_diagram_2;

typedef HVD_Envelope_diagram_2                          Envelope;

// --------------------------------------------------------------------
//k CVD vertex colored data structure
typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
typedef CGAL::Regular_triangulation_vertex_base_2<Kernel>              Vbase;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::vector<bool>,Kernel,Vbase> Vb;
typedef CGAL::Regular_triangulation_face_base_2<Kernel>                Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    Tds;
typedef CGAL::Regular_triangulation_2<Kernel,Tds>                      Regular_color;

typedef CGAL::Polygon_2<Kernel>                         Cluster_2;
typedef CGAL::Polygon_2<VD_Kernel>                      VD_Cluster_2;

/* Ipelet labels */
const unsigned int num_entries = 15;
const std::string sublabel[] = {
  "two points euclidean bisector",            // 0
  "two points L_inf bisector",                // 1
  "point/segment L_inf-parabola",             // 2
  "two sites L_inf bisector",                 // 3
  "Linf 2D Voronoi Diagram",                  // 4
  "L2 farthest Voronoi Diagram",              // 5
  "Hausdorff Voronoi Diagram",                // 6
  "L2 FVD with polygonal input",              // 7
  "L2 farthest color Voronoi diagram (FCVD)", // 8
  "L2 NVD with polygonal input",              // 9
  "L2 segment Voronoi Diagram",               // 10
  "L2 farthest segment Voronoi Diagram",      // 11
  "L2 FCVDstar",                              // 12
  "L2 k-CVD (with polygonal input?)",         // 13
  "Help"
};
const std::string helpmsg[] = {
  "Draw the L2 bisector of two points",
  "Draw the L_inf bisector of two points",
  "Draw the L_inf parabola for a point and a segment",
  "Draw the L_inf bisector for two sites (point or segment)",
  "Draw the L_inf Voronoi diagram for points with 2D bisectors",
  "Draw the L2 farthest Voronoi diagram for points",
  "Draw the Hausdorff Voronoi diagram for points",
  "Draw the L2 FVD for points of given clusters",
  "Draw the farthest color Voronoi diagram for points",
  "Draw the L2 NVD for points of given clusters",
  "Draw the L2 Voronoi diagram for segments",
  "Draw the L2 farthest Voronoi diagram for segments",
  "Draw the L2 FCVDstar for points of given clusters",
  "Draw the k-order color Voronoi diagram for points of given clusters"
};

class bisectorIpelet
  : public CGAL::Ipelet_base<Kernel,num_entries> {
public:
  bisectorIpelet()
    :CGAL::Ipelet_base<Kernel,num_entries>("Bisectors",sublabel,helpmsg){}
  void protected_run(int);

  // cluster grabber
  template <class output_iterator>
  struct Cluster_grabber
    : public CGAL::internal::Cluster_grabber<Kernel,output_iterator>{
    Cluster_grabber(output_iterator it):
      CGAL::internal::Cluster_grabber<Kernel,output_iterator>(it){}
  };

  template<class output_iterator>
  boost::function_output_iterator<Cluster_grabber<output_iterator> >
  cluster_grabber(output_iterator it){
    return boost::make_function_output_iterator(
      Cluster_grabber<output_iterator>(it)
    );
  }

private:
  Algebraic PRECISION = Algebraic(5); // values under 0.5 take too long
  Rational BOUNDARY = Rational(10000);

  /* Find a number of points on the curve cv such that they are all close enough
   * to approximate the curve. Then create segments between these points and add
   * them to the provided list. */
  void arc_to_segments(const X_monotone_curve_2& cv, std::list<Segment_2>& segments) {

    /* create points on the arc, create segments between them and push segments */
    Algebraic current_x = cv.left().x();
    Algebraic end_x = cv.right().x();

    for (; current_x + PRECISION < end_x; current_x += PRECISION) {
      Alg_point_2 current = Alg_point_2(current_x, 0);
      Alg_point_2 next = Alg_point_2(current_x + PRECISION, 0);
      /* project points on arc */
      Alg_point_2 arc_segment_start_pt = cv.point_at_x(current);
      Alg_point_2 arc_segment_end_pt = cv.point_at_x(next);
      /* create segment between two points */
      Segment_2 * arc_segment = new Segment_2(
        Point_2(
          CGAL::to_double(arc_segment_start_pt.x()),
          CGAL::to_double(arc_segment_start_pt.y())
        ),
        Point_2(
          CGAL::to_double(arc_segment_end_pt.x()),
          CGAL::to_double(arc_segment_end_pt.y())
        )
      );
      segments.push_back(*arc_segment);
    }

    /* push last segment */
    Alg_point_2 current = Alg_point_2(current_x, 0);
    Alg_point_2 last_arc_segment_start_pt = cv.point_at_x(current);
    Alg_point_2 last_arc_segment_end_pt = cv.right();
    Segment_2 * last_arc_segment = new Segment_2(
      Point_2(
        CGAL::to_double(last_arc_segment_start_pt.x()),
        CGAL::to_double(last_arc_segment_start_pt.y())
      ),
      Point_2(
        CGAL::to_double(last_arc_segment_end_pt.x()),
        CGAL::to_double(last_arc_segment_end_pt.y())
      )
    );
    segments.push_back(*last_arc_segment);

    return;
  }

  /* Check if two points are both on the boundary of the diagram. This boundary
   * is added artificially in FSVD (fn: 10). */
  //TODO correct: what if this the diagram is a single long line that touches
  // two bounds?
  bool on_boundary(Point_2 p1, Point_2 p2) {
    return
      false && // comment this line to activate the method
      ((CGAL::abs(p1.x()) == BOUNDARY && CGAL::abs(p2.x()) == BOUNDARY) ||
      (CGAL::abs(p1.y()) == BOUNDARY && CGAL::abs(p2.y()) == BOUNDARY))
    ;
  }

  template <class T,class output_iterator>
  bool cast_obj_into_seg(T& obj, Iso_rectangle_2& bbox,output_iterator out_it);

  template<class iterator,class output_iterator>
  void cast_list_into_seg(iterator first, iterator end,
			              Iso_rectangle_2& bbox, output_iterator out_it);
  void kCVD_draw_dual_in_ipe(Regular_color& rt, Iso_rectangle_2& bbox);

}; // end of class bisectorIpelet
// --------------------------------------------------------------------

} // end of namespace CGAL_bisectors

CGAL_IPELET(CGAL_bisectors::bisectorIpelet)

#endif
