#pragma once
#include "geobase.h"

using namespace std;

namespace ZDTree{

    using parlay::sequence;
    using geobase::Point;

    struct ZDNode{
        using node = void;
        using node_size_t = size_t;

        struct BaseNode{
		    BaseNode(){}

		    virtual ~BaseNode() = default;
		    virtual bool is_leaf(){ return false; }
		    virtual node_size_t get_num_points(){ return 0; }
		    virtual node_size_t get_ref_count(){ return 1; }
        };

        struct InteNode: BaseNode{  //  interior node, store two pointers
            node_size_t ref_cnt;
            node_size_t num_pts;
            BaseNode* l_son;
            BaseNode* r_son;

            InteNode(): ref_cnt(1), num_pts(0), l_son(nullptr), r_son(nullptr){} // default constructor, ref_cnt should be 1

		    virtual bool is_leaf(){ return false; }
		    virtual node_size_t get_num_points(){ return num_pts; }
		    virtual node_size_t get_ref_count(){ return ref_cnt; }
        };

        struct LeafNode: BaseNode{
            node_size_t ref_cnt;
            sequence<Point*> records;

            template<typename Records>
            LeafNode(Records r): ref_cnt(1){       //  constructor, store records pointers
			    records.resize(r.size());
			    for (size_t i = 0; i < r.size(); i++){
				    records[i] = r[i];
			    }
            }

		    virtual bool is_leaf(){ return true; }
		    virtual node_size_t get_num_points(){ return records.size(); }
		    virtual node_size_t get_ref_count(){ return ref_cnt; }

            void print_records(){
			    cout << records.size() << endl;
			    for (size_t i = 0; i < records.size(); i++){
				    cout << "(" << records[i]->x << ", " << records[i]->y << ")" << endl; 
			    }
            }
        };

        using inte_node_allocator = parlay::type_allocator<InteNode>;   //  parlay memory management
        using leaf_node_allocator = parlay::type_allocator<LeafNode>;   //  parlay memory management

        static InteNode* cast_to_inte(BaseNode* x){
            assert(!x->is_leaf());
            return (InteNode*)x;
        }

        static LeafNode* cast_to_leaf(BaseNode* x){
            assert(x->is_leaf());
            return (LeafNode*)x;
        }

        template<typename P_set>
        static LeafNode* alloc_leaf(P_set &P){  //  return a new leaf, call constructor for a point set
            LeafNode* new_leaf = leaf_node_allocator::create(P);
            return new_leaf;
        }

        static InteNode* alloc_inte(){  //  return a new interior node, call default constructor
            InteNode* new_inte = inte_node_allocator::create();
            return new_inte;
        }

        static void free_leaf(LeafNode* x){
            leaf_node_allocator::retire(x);
            x = nullptr;
        }

        static void free_inte(InteNode* x){
            inte_node_allocator::retire(x);
            x = nullptr;
        }

        // free a tree node x
        static void free_node(BaseNode* x){
            if (x->is_leaf()){
                auto a = cast_to_leaf(x);
                // leaf_node_allocator::free(a);
                leaf_node_allocator::retire(a);
            }
            else{
                auto a = cast_to_inte(x);
                // free_inte(a);
                // inte_node_allocator::free(a);
                inte_node_allocator::retire(a);
            }
        }
    };
}