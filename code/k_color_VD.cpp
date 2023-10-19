// Copyright (c) 2005-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6/CGAL_ipelets/demo/CGAL_ipelets/include/CGAL_ipelets/k_delaunay.h $
// $Id: k_delaunay.h 7a62583 2022-11-14T19:14:33+01:00 albert-github
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Nicolas Carrez

// Heavily eddited by : Nico

#ifndef K_COLOR_VD_H
#define K_COLOR_VD_H

//inludes for the typedefs, can we avoid this?
#include "bisectors.h"

using namespace CGAL_bisectors;

double give_weight(const Kernel::Weighted_point_2& P){return CGAL::to_double(P.weight());}
double give_x(const Kernel::Weighted_point_2& P){return CGAL::to_double(P.point().x());}
double give_y(const Kernel::Weighted_point_2& P){return CGAL::to_double(P.point().y());}

double give_weight(const Kernel::Point_2& ){return 0;}
double give_x(const Kernel::Point_2& P){return CGAL::to_double(P.x());}
double give_y(const Kernel::Point_2& P){return CGAL::to_double(P.y());}

// Create the weighted point associated to the current selection of k points
Kernel::Weighted_point_2 compute_weighted_barycenter(std::vector<int>& Current_tuple,
													 std::vector<std::list<Cluster_2>::iterator>& Current_col_it,
													 int order){
  typedef Kernel::Point_2 Point_2;
  typedef Kernel::Weighted_point_2  Weighted_point_2;
  double weight = 0;
  double pt_x = 0;
  double pt_y = 0;
  for(int color_pt_i = 0; color_pt_i < Current_tuple.size(); color_pt_i++){
	int vertex_index_i = Current_tuple[color_pt_i];
	Point_2 pt_i = (Current_col_it[color_pt_i])->vertex(vertex_index_i);
    pt_x = pt_x + give_x(pt_i);
    pt_y = pt_y + give_y(pt_i);
    weight = weight + order * give_weight(pt_i);
    //subtract form the weight the sum of the squared distances between each pair of wpoints selected
    for(int color_pt_j = color_pt_i + 1; color_pt_j < Current_tuple.size(); color_pt_j++){
	  int vertex_index_j = Current_tuple[color_pt_j];
      Point_2 pt_j = Current_col_it[color_pt_j]->vertex(vertex_index_j);
      weight = weight - CGAL::to_double(CGAL::squared_distance(
                                          Kernel::Construct_point_2()(pt_i),
                                          Kernel::Construct_point_2()(pt_j)));
    }
  }
  weight = weight / (double) (order*order);
  pt_x = pt_x / (double) order;
  pt_y = pt_y / (double) order;
  return Weighted_point_2(Point_2(pt_x,pt_y),weight);
}

// Compute the next tuple of points from the current one
// If input Current_sel is the last one, return False (else True)
bool next_tuple_Selection(std::vector<int>& Current_tuple,
						  std::vector<int>& Current_col_sizes,
						  int order){
  
  bool carry = true;
  int i = order - 1; //int order = Current_sel.size();
  while(carry and i >= 0){
	 //iterate tuple to next point
    Current_tuple[i] = Current_tuple[i] + 1;
    // if all points of this color have been visited
    if(Current_tuple[i] >= Current_col_sizes[i]){
	  // set back to the first point from this color
      Current_tuple[i] = 0;
      // and carry the increment to the next color
      carry = true;
    }else{carry = false;}
    i--;
  }
  // If carry after incrementing the first point in vector
  // Then this was the last tuple of points to visit
  return not carry;
}

// Input should be a list of clusters/lists of Weighted/unweighted 2 dimensional points
// TODO :  extend code to higher dimensions
void k_color_VD(Regular_color& rtc, std::list<Cluster_2>& input_clusters, int order){
  typedef Kernel::Point_2 Point_2;
  typedef Kernel::Weighted_point_2  Weighted_point_2;
  
  int N_colors = input_clusters.size(); // number of cluster/colors
  
  // Vector of number of vertices for all clusters/colors
  std::vector<int> Colors_sizes;
  // Vector of pointers to all clusters/colors
  std::vector<std::list<Cluster_2>::iterator> Colors_it;
  
  //                        i.e. it_clusters
  std::list<Cluster_2>::iterator it_colors = input_clusters.begin(); 
  for(; it_colors != input_clusters.end(); ++it_colors){
	  Colors_sizes.push_back(it_colors->size());
	  std::list<Cluster_2>::iterator it = it_colors;
	  Colors_it.push_back(it);
  }
  
  // The bitmask of colors currently selected
  std::vector<bool> Colors_bitmask(N_colors, false);
  std::fill(Colors_bitmask.end() - order, Colors_bitmask.end(), true);
  
  do{
	//for (auto x : Colors_bitmask) std::cout << x;
    //std:: cout  << " colors" << std::endl;
    // The tuple of points from colors currently selected
    std::vector<int> Current_tuple (order, 0);
    // Pointers to the first and last element of the selected colors
    std::vector<int> Current_col_sizes;
    std::vector<std::list<Cluster_2>::iterator> Current_col_it;
    for(int i = 0; i < N_colors; i++){
      // If the i-th color is selected by the bitmask
      if(Colors_bitmask[i]){
		// Initialize the point selection at the first point of the i-th color
        Current_col_sizes.push_back(Colors_sizes[i]);
        Current_col_it.push_back(Colors_it[i]);
      }
    }
    
    do{
      //for (auto x : Current_tuple) std::cout << x;
      //std:: cout << std::endl;
      //insert weighted barycenter to the VD
      Regular_color::Vertex_handle vh;
      vh = rtc.insert(compute_weighted_barycenter(Current_tuple,Current_col_it,order));
      //add color
      vh->info() = Colors_bitmask;
      //select the next tuple of points
    }while(next_tuple_Selection(Current_tuple,Current_col_sizes,order));
    
    //select the next combination of colors (advance the bitmask permutation)
  }while(std::next_permutation(Colors_bitmask.begin(),Colors_bitmask.end()));
  
} // end k_color_VD

typedef typename Kernel::Segment_2                                        Segment_2;
typedef typename Kernel::Ray_2                                            Ray_2;
typedef typename Kernel::Line_2                                           Line_2;

template <class T,class output_iterator>
bool bisectorIpelet::cast_obj_into_seg(T& obj, Iso_rectangle_2& bbox,output_iterator out_it){
  CGAL::Object obj_cgal = CGAL::intersection(obj,bbox);
  Segment_2 s;
  bool ret=CGAL::assign(s, obj_cgal);
  if (ret) *out_it++=s;
  return ret;
}

//Convert infinite objects into drawable segments
template<class iterator,class output_iterator>
void bisectorIpelet::cast_list_into_seg(iterator first, iterator end,
                                     Iso_rectangle_2& bbox, output_iterator out_it){
  for (iterator it=first;it!=end;++it)
	cast_obj_into_seg(*it,bbox,out_it);
}

void bisectorIpelet::kCVD_draw_dual_in_ipe(Regular_color& rtc, Iso_rectangle_2& bbox){
	bool internal_edges=false;
	
	std::list<Ray_2> ray_list;
	std::list<Line_2> line_list;
	std::list<Segment_2> seg_list;
	Regular_color::Edge_iterator eit = rtc.finite_edges_begin();
	int count1 = 0; int count2 = 0;
	for(; eit != rtc.finite_edges_end(); ++eit) {
		Regular_color::Edge e = *eit;
		Regular_color::Face_handle f = e.first;
		int i = e.second;
		Regular_color::Vertex_handle v1 = f->vertex(Regular_color::cw(i));
		Regular_color::Vertex_handle v2 = f->vertex(Regular_color::ccw(i));
		//std::cout << "V1col: ";
		std::vector<bool> color_v1 = v1->info();
		//for (bool x : color_v1) std::cout << x;
		//std::cout << " V2col: ";
		std::vector<bool> color_v2 = v2->info();
		//for (bool x : color_v2) std::cout << x;
		//std::cout << std::endl;
		//evaluate edge, if monocromatic then skip it
		if (color_v1 == color_v2) continue;
		CGAL::Object o = rtc.dual(e);
		Line_2  l;
		Ray_2   r;
		Segment_2 s;
		if(CGAL::assign(s,o)) seg_list.push_back(s);
		if(CGAL::assign(r,o)) ray_list.push_back(r);
		if(CGAL::assign(l,o)) line_list.push_back(l);
	}
	std::vector<Segment_2> seg_cont;
	//filter degenerate segments
	for(typename std::list<Segment_2>::iterator iteS = seg_list.begin();iteS!=seg_list.end();){
		typename std::list<Segment_2>::iterator itc=iteS++;
		if (itc->is_degenerate()) seg_list.erase(itc);
	}

	cast_list_into_seg(ray_list.begin(), ray_list.end(), bbox, std::back_inserter(seg_cont));//cast rays into segments in bbox
	cast_list_into_seg(line_list.begin(),line_list.end(),bbox,std::back_inserter(seg_cont));//cast lines into segments in bbox
	cast_list_into_seg(seg_list.begin(),seg_list.end(),bbox,std::back_inserter(seg_cont));//cast lines into segments in bbox
	
	bool makegrp=true;
	bool deselect_all=false;
	
	draw_in_ipe(seg_cont.begin(), seg_cont.end(),makegrp);
	
	if (deselect_all) get_IpePage()->deselectAll();
}

#endif
