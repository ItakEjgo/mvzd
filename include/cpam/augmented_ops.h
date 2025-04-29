#pragma once
#include "utils.h"
#include "map_ops.h"

// *******************************************
//   AUGMENTED MAP OPERATIONS
// *******************************************

extern size_t visited_leaf;
extern size_t visited_inte;
extern double leaf_time;
extern double inte_time;

namespace cpam {

template<class Map>
struct augmented_ops : Map {
  using Entry = typename Map::Entry;
  using node = typename Map::node;
  using Seq = typename Map::_Seq;
  using ET = typename Map::ET;
  using GC = typename Map::GC;
  using K = typename Map::K;
  using aug_t = typename Entry::aug_t;
  using ptr = typename GC::ptr;
  using Map::B;
  using Map::kBaseCaseSize;
  using Map::kNodeLimit;

  static inline aug_t aug_val(node* b) {
    return Seq::aug_val(b);
  }

  struct aug_sum_t {
    aug_t result;
    aug_sum_t() : result(Entry::get_empty()) {}
    void add_entry(ET e) {
      result = Entry::combine(result, Entry::from_entry(e));
    }
    void add_aug_val(aug_t av) {
      result = Entry::combine(result, av);
    }
  };

  // the sum left of or at key
  template<class aug>
  static void aug_sum_left(node* b, const K& key, aug& a) {
    while (b) {
      if (Map::is_compressed(b)) {
        auto fn = [&] (const auto& et) {
          if (Map::comp(Entry::get_key(et), key)) {
            a.add_entry(et);
            return true;
          }
          return false;
        };
        Map::iterate_cond(b, fn);
        return;
      }
      auto rb = Map::cast_to_regular(b);
      if (!Map::comp(key, Map::get_key(rb))) {
        a.add_entry(Map::get_entry(rb));
        if (rb->lc) a.add_aug_val(Map::aug_val(rb->lc));
        b = rb->rc;
      } else {
        b = rb->lc;
      }
    }
  }

  // the sum right of or at key
  template<class aug>
  static void aug_sum_right(node* b, const K& key, aug& a) {
    while (b) {
      if (Map::is_compressed(b)) {
        auto fn = [&] (const auto& et) {
          if (!Map::comp(Entry::get_key(et), key)) {
            a.add_entry(et);
          }
        };
        Map::iterate_seq(b, fn);
        return;
      }
      auto rb = Map::cast_to_regular(b);
      if (!Map::comp(Map::get_key(rb), key)) {
	a.add_entry(Map::get_entry(rb));
	if (rb->rc) a.add_aug_val(Map::aug_val(rb->rc));
	b = rb->lc;
      } else b = rb->rc;
    }
  }

  template<class aug>
  static void aug_sum_range(node* b, const K& key_left, const K& key_right, aug& a) {
    node* r = Map::range_root_2(b, key_left, key_right);
    if (r) {
      if (Map::is_compressed(r)) {
        auto fn = [&] (const auto& et) {
          if (Map::comp(key_left, Entry::get_key(et)) &&
              Map::comp(Entry::get_key(et), key_right)) {
            a.add_entry(et);
          }
        };
        Map::iterate_seq(r, fn);
      } else {
        auto rr = Map::cast_to_regular(r);
        // add in left side (right of or at key_left)
        aug_sum_right(rr->lc, key_left, a);
        // add in middle
        a.add_entry(Map::get_entry(rr));
        // add in right side (left of or at key_right)
        aug_sum_left(rr->rc, key_right, a);
      }
    }
  }

  template<typename Func>
  static std::optional<ET> aug_select(node* b, const Func& f) {
    if (!b) return {};
    if (Map::is_compressed(b)) {
      std::optional<ET> ret;
      auto fn = [&] (const auto& et) {
        if (!f(Entry::from_entry(et))) {
          ret = et;
          return false;  // stop
        }
        return true;  // keep iterating
      };
      Map::iterate_cond(b, fn);
      return ret;
    }
    auto rb = Map::cast_to_regular(b);
    if (f(Map::aug_val(rb->lc))) {
      if (f(Entry::from_entry(Map::get_entry(rb)))) {
        return aug_select(rb->rc, f);
      }
      return Map::get_entry(rb);
    }
    return aug_select(rb->lc, f);
  }

  template <class Func>
  static node* aug_filter_bc(ptr b1, const Func& f) {
    assert(b1.size() > 0);
    ET stack[kBaseCaseSize + 1];

    auto b1_node = b1.node_ptr();
    size_t offset = 0;
    aug_t cur = Entry::get_empty();
    auto copy_f = [&] (ET a) {  // has to be a copy since we move
      cur = Entry::combine(cur, Entry::from_entry(a));
      if (f(cur)) {
        parlay::move_uninitialized(stack[offset++], a);
      }
    };
    Map::iterate_seq(b1_node, copy_f);
    assert(offset <= kBaseCaseSize);

    Map::decrement_recursive(b1_node);

    if (offset < B) {
      return Map::to_tree_impl((ET*)stack, offset);
    } else {
      return Map::make_compressed(stack, offset);
    }
  }

  template <class Func>
  static node* aug_filter_bc_mid(ptr b1, const Func& f) {
    assert(b1.size() > 0);
    ET stack[kBaseCaseSize + 1];

    auto b1_node = b1.node_ptr();
    size_t offset = 0;
    aug_t cur = Entry::get_empty();
    auto copy_f = [&] (ET a) {  // has to be a copy since we move
      cur = Entry::combine(cur, Entry::from_entry(a));
      if (f(cur)) {
        parlay::move_uninitialized(stack[offset++], a);
      }
    };
    Map::iterate_seq_mid(b1_node, copy_f);
    assert(offset <= kBaseCaseSize);

    Map::decrement_recursive(b1_node);

    if (offset < B) {
      return Map::to_tree_impl((ET*)stack, offset);
    } else {
      return Map::make_compressed(stack, offset);
    }
  }

  template<class F, typename F2>
  static size_t range_count_filter(ptr b, const F &f, const F2 &f2, size_t granularity=kNodeLimit) {
    if (b.empty()) return 0;
    auto cur_aug = aug_val(b.unsafe_ptr());
    auto [lc, e, rc, root] = Map::expose(std::move(b));

    auto cur_pt = std::get<1>(e);
    // auto cur_pt2 = std::get<1>(e);

    auto flag = f(cur_aug.first);
    // auto flag2 = f(cur_aug.first);
  
    if (flag < 0) {
      GC::decrement(root);
      return 0;
    }
    if (flag == 1) {
      GC::decrement(root);
      return cur_aug.second;
    }

    auto cur_pt_inside = f2(cur_pt) > 0 ? 1 : 0;
    // auto cur_pt_inside2 = f2(cur_pt) > 0 ? 1 : 0;

    auto l = range_count_filter(std::move(lc), f, f2, granularity); 
    auto r = range_count_filter(std::move(rc), f, f2, granularity); 

    GC::decrement(root);

    return l + r + cur_pt_inside;
  }

  // // F check bounding box.
  // template<typename F, typename P>
  // auto spatial_diff_filter(node *lhs, node*rhs, F &f, size_t granularity=kNodeLimit){
  //   P add = {}, remove = {}; 
  //   if (lhs == rhs){
  //     return make_tuple(add, remove);
  //   }
  //   if (!lhs){
  //     remove = range_report_filter2();
  //   }
  // }

  // F for box check, F2 for point check
  template<typename F, typename F2> 
  static size_t range_count_filter2(node* b, const F &f, const F2 &f2, size_t granularity=kNodeLimit){
    if (!b) return 0;
    auto cur_aug = aug_val(b);
    auto flag = f(cur_aug.first);
    if (flag < 0) return 0; //exclude
    if (flag == 1){
      return cur_aug.second;  // fully contained
    }

    if (Map::is_compressed(b)){ // leaf node
      auto ret = 0;
      auto f_filter = [&](const auto& et){
        auto cur_pt = std::get<1>(et);
        if (f2(cur_pt) == 1){ ret++; }
      };
      Map::iterate_seq(b, f_filter);
      return ret;
    }

    auto rb = Map::cast_to_regular(b);
    auto cur_pt = Map::get_val(rb);
    auto flag2 = f2(cur_pt) == 1 ? 1 : 0;

    auto l = range_count_filter2(rb->lc, f, f2, granularity); 
    auto r = range_count_filter2(rb->rc, f, f2, granularity); 
    
    return l + r + flag2;
  }

  template<class F, typename Out>
  static void range_report_filter(ptr b, const F &f, int64_t &cnt, Out &out, size_t granularity=kNodeLimit) {
    if (b.empty()) return;
    // auto b_ptr = b.node_ptr();

    // if (Map::is_compressed(b_ptr)){ //  touch leaf node
    //   std::cout << "touch a leaf" << std::endl;
    //   auto cur_aug = aug_val(b.unsafe_ptr());
    //   auto flag = f(cur_aug.first);
    //   if (flag < 0) return; //  exclude
    //   if (flag == 1){ //  fully covered
    //     auto f_collect = [&](const auto& et){ out[cnt++] = std::get<1>(et); };
    //     Map::iterate_seq(b_ptr, f_collect);
    //   }
    //   else{ //  partially overlapped
    //     auto f_filter = [&](const auto& et){
    //       auto cur_pt = std::get<1>(et);
    //       auto cur_box = std::make_pair(cur_pt, cur_pt);
    //       if (f(cur_box) == 1){ out[cnt++] = cur_pt; }
    //     };
    //     Map::iterate_seq(b_ptr, f_filter);
    //   }
    //   return;
    // }

    // std::cout << "touch an interior" << std::endl;
    auto cur_aug = aug_val(b.unsafe_ptr());
    
    auto [lc, e, rc, root] = Map::expose(std::move(b));

    auto flag = f(cur_aug.first);

    if (flag < 0) {
      GC::decrement(root);
      return; //  exclude
    }

    auto cur_pt = std::get<1>(e);
    auto pt_box = std::make_pair(cur_pt, cur_pt);
    auto cur_pt_inside = f(pt_box) > 0 ? 1 : 0;

    range_report_filter(std::move(lc), f, cnt, out);
    if (cur_pt_inside) {
      out[cnt++] = cur_pt;
    }
    range_report_filter(std::move(rc), f, cnt, out);

    GC::decrement(root);
    return;
  }


  template<class F, typename Out>
  static void range_report_filter2(node* b, const F &f, int64_t &cnt, Out &out, size_t granularity=kNodeLimit) {
    if (!b) return;
    auto cur_aug = aug_val(b);
    auto flag = f(cur_aug.first);
    
    if (flag < 0) return; //exclude

    if (Map::is_compressed(b)){ // leaf node
      // visited_leaf++;
      
      // auto mbr = cur_aug.first;
      // std::cout << std::fixed << std::setprecision(6) << "[(" << mbr.first.x << ", " << mbr.first.y << "), (" << mbr.second.x << ", " << mbr.second.y << ")" << "]" << std::endl;

      if (flag == 1){
        auto f_filter = [&](const auto &et){
          // out[cnt++] = std::get<1>(et);
          parlay::assign_uninitialized(out[cnt++], std::get<1>(et));
        };
        Map::iterate_seq(b, f_filter);
        // Map::iterate_seq(b, f_filter);
      }
      else{
        auto f_filter = [&](const auto& et){
          auto cur_pt = std::get<1>(et);
          auto pt_box = std::make_pair(cur_pt, cur_pt);
          if (f(pt_box) == 1){ 
            parlay::assign_uninitialized(out[cnt++], cur_pt);
            // out[cnt++] = cur_pt; 
          }
        };
        Map::iterate_seq(b, f_filter);
      }
      return;
    }

    // visited_inte++;

    auto rb = Map::cast_to_regular(b);
    auto cur_pt = Map::get_val(rb);
    
    auto pt_box = std::make_pair(cur_pt, cur_pt);
    auto flag2 = f(pt_box) == 1 ? 1 : 0;
    range_report_filter2(rb->lc, f, cnt, out, granularity); 
    if (flag2) {
      parlay::assign_uninitialized(out[cnt++], cur_pt);
      // out[cnt++] = cur_pt;
    }
    range_report_filter2(rb->rc, f, cnt, out, granularity); 
  }

  //  F is point-point dis, F2 is point-mbr dis
  template<typename F, typename F2, typename Out>
  static void knn_filter(node* b, const F &f, const F2 &f2, size_t &k, Out &out) {
    if (!b) return;

    auto pt_check = [&](const auto &cur_pt){
      auto cur_dis = f(cur_pt);
      if (out.size() < k) out.push(std::make_pair(cur_pt, cur_dis));
      else if (cur_dis < out.top().second){
        out.pop();
        out.push(std::make_pair(cur_pt, cur_dis));
      }
    };

    if (Map::is_compressed(b)){ // leaf node
      auto f_filter = [&](const auto &et){
        auto cur_pt = std::get<1>(et);
        pt_check(cur_pt);
      };
      Map::iterate_seq(b, f_filter);
      return; 
    }

    auto rb = Map::cast_to_regular(b);
    auto cur_pt = Map::get_val(rb);
    pt_check(cur_pt);

    auto l_dis = std::numeric_limits<double>::max();
    auto r_dis = std::numeric_limits<double>::max();
    if (rb->lc){
      auto cur_aug = aug_val(rb->lc);
      l_dis = f2(cur_aug.first);
    }
    if (rb->rc){
      auto cur_aug = aug_val(rb->rc);
      r_dis = f2(cur_aug.first);
    }
    auto go_left = [&](){
      if (out.size() < k || l_dis < out.top().second){
        knn_filter(rb->lc, f, f2, k, out);
      }
    };
    auto go_right = [&](){
      if (out.size() < k || r_dis < out.top().second){
        knn_filter(rb->rc, f, f2, k, out);
      }
    };

    if (l_dis <= r_dis){  //  go left first
      go_left();
      go_right();
    }
    else{
      go_right();
      go_left();
    }
  }


  template<class Func>
  static node* aug_filter_mid(ptr b, const Func& f, size_t granularity=kNodeLimit) {
    if (b.empty()) return NULL;
    // if (b.size() <= kBaseCaseSize) {
    //   return aug_filter_bc_mid(std::move(b), f);
    // }
    // TODO: better functionality for getting aug_val from b
    // std::cout << "My aug_val = " << aug_val(b.unsafe_ptr()) << std::endl;
    if (!f(aug_val(b.unsafe_ptr()))) return NULL;

    // size_t n = b.size();
    auto [lc, e, rc, root] = Map::expose(std::move(b));

    // auto [l, r] = utils::fork<node*>(n >= granularity,
    //   [&]() {return aug_filter_mid(std::move(lc), f, granularity);},
    //   [&]() {return aug_filter_mid(std::move(rc), f, granularity);});
    
    auto l = aug_filter_mid(std::move(lc), f, granularity);
    auto r = aug_filter_mid(std::move(rc), f, granularity);

    if (f(Entry::from_entry(e))) {
      return Map::join(l, e, r, root);
    } else {
      GC::decrement(root);
      return Map::join2(l, r);
    }
  }

  template <class Func>
  static node* aug_filter_bc2(ptr b1, const Func& f) {
    assert(b1.size() > 0);
    ET stack[kBaseCaseSize + 1];

    auto b1_node = b1.node_ptr();
    size_t offset = 0;
    // aug_t cur = Entry::get_empty();
    auto copy_f = [&] (ET a) {  // has to be a copy since we move
      // cur = Entry::combine(cur, Entry::from_entry(a));
      if (f(a)) {
        parlay::move_uninitialized(stack[offset++], a);
      }
    };
    Map::iterate_seq(b1_node, copy_f);
    assert(offset <= kBaseCaseSize);

    Map::decrement_recursive(b1_node);

    if (offset < B) {
      return Map::to_tree_impl((ET*)stack, offset);
    } else {
      return Map::make_compressed(stack, offset);
    }
  }

  template<class Fpt, class Fbb>
  static node* aug_filter2(ptr b, const Fpt& fpt, const Fbb &fbb, size_t granularity=kNodeLimit) {
    if (b.empty()) return NULL;
    if (b.size() <= kBaseCaseSize) {
      return aug_filter_bc2(std::move(b), fpt);
    }
    // TODO: better functionality for getting aug_val from b
    // std::cout << "My aug_val = " << aug_val(b.unsafe_ptr()) << std::endl;
    if (!fbb(aug_val(b.unsafe_ptr()))) return NULL;

    size_t n = b.size();
    auto [lc, e, rc, root] = Map::expose(std::move(b));

    auto [l, r] = utils::fork<node*>(n >= granularity,
      [&]() {return aug_filter2(std::move(lc), fpt, fbb, granularity);},
      [&]() {return aug_filter2(std::move(rc), fpt, fbb, granularity);});

    if (fpt(e)) {
      return Map::join(l, e, r, root);
    } else {
      GC::decrement(root);
      return Map::join2(l, r);
    }
  }

  template<class Func>
  static node* aug_filter(ptr b, const Func& f, size_t granularity=kNodeLimit) {
    if (b.empty()) return NULL;
    if (b.size() <= kBaseCaseSize) {
      return aug_filter_bc(std::move(b), f);
    }
    // TODO: better functionality for getting aug_val from b
    // std::cout << "My aug_val = " << aug_val(b.unsafe_ptr()) << std::endl;
    if (!f(aug_val(b.unsafe_ptr()))) return NULL;

    size_t n = b.size();
    auto [lc, e, rc, root] = Map::expose(std::move(b));

    auto [l, r] = utils::fork<node*>(n >= granularity,
      [&]() {return aug_filter(std::move(lc), f, granularity);},
      [&]() {return aug_filter(std::move(rc), f, granularity);});

    if (f(Entry::from_entry(e))) {
      return Map::join(l, e, r, root);
    } else {
      GC::decrement(root);
      return Map::join2(l, r);
    }
  }

  template <class Func>
  static node* insert_lazy(node* b, const ET& e, const Func& f) {
    aug_t av = Entry::from_entry(e);
    auto g = [&] (const aug_t& a) { return Entry::combine(av,a);};

    auto lazy_join = [&] (node* l, node* r, node* _m) -> node* {
      auto m = Map::cast_to_regular(_m);
      m->rc = r; m->lc = l;
      if (Map::is_balanced(m)) {
	Map::lazy_update(m,g);
	return m;
      } else return Map::node_join(l,r,m);
    };

    return Map::template insert_tmpl<Func, decltype(lazy_join), false>(b, e, f, lazy_join);
  }

};

}  // namespace cpam
