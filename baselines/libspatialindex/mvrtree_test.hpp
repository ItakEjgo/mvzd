#pragma once

#include <iostream>
#include <cstring>
#include <spatialindex/SpatialIndex.h>
#include <chrono>
#include <malloc.h>

#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "helper/time_loop.h"
#include "helper/geobase.h"

extern size_t m_mvrtree_size;

namespace MVRTest{
    using namespace SpatialIndex;
    using namespace std;

    class CollectPointVisitor : public IVisitor{
    public:
        parlay::sequence<geobase::Point> pts;
    
        auto retrive_point(const IData &d){
            auto id = d.getIdentifier();
            IShape* pS;
            d.getShape(&pS);
            Point p;
            pS->getCenter(p);
            auto x = p.m_pCoords[0];
            auto y = p.m_pCoords[1];
            delete pS;
            return geobase::Point(id, x, y);
        }

        void visitNode(const INode&) override{}
        void visitData(const IData &d) override{
            // cout << "ok" << endl;
            pts.push_back(retrive_point(d));
        }
        void visitData(std::vector<const IData*> &v) override{
            /* Bulk Query Processing for Rtree Queries */
            // cout << "ok" << endl;
            for (auto ptr: v){
                pts.push_back(retrive_point(*ptr));
            }
        }
    };

    double memory_in_MB(size_t num_of_bytes){
         return 1.0 * num_of_bytes / 1024.0 / 1024.0;
    }

    template<typename Tree, typename RQ>
    auto spatial_diff(Tree &tree, RQ &q, size_t d = 2){
        double plow[2], phigh[2];    
        plow[0] = q.first.x, plow[1] = q.first.y;
        phigh[0] = q.second.x, phigh[1] = q.second.y;
        CollectPointVisitor vis1, vis2;
        auto r = TimeRegion(plow, phigh, 0, 0.9, d);
        auto r2 = TimeRegion(plow, phigh, 1, 2.1, d);

        tree->intersectsWithQuery(r, vis1);
        // geobase::print_Pset_info(vis1.pts, "test");
        tree->intersectsWithQuery(r2, vis2);
        // geobase::print_Pset_info(vis2.pts, "test");
        // auto [added, removed] = geobase::diff_by_id(vis1.pts, vis2.pts);
        // cout << "two ok" << endl;
        // auto [added, removed] = geobase::diff_by_id(vis1.pts, vis2.pts);
        // cout << "ready to return" << endl;
        return geobase::diff_by_id(vis1.pts, vis2.pts);
    }
    
    template<typename Tree, typename PT>
    auto apply_insert(Tree &tree, PT &P, double st = 0.0, double ed = 0.0, size_t d = 2){
        double plow[2], phigh[2];
        for (auto &pt: P){
            plow[0] = pt.x, plow[1] = pt.y;
            phigh[0] = pt.x, phigh[1] = pt.y;
            auto r = TimeRegion(plow, phigh, st, ed, d);
            tree->insertData(0, nullptr, r, pt.id);
        }           
    }

    template<typename Tree, typename PT>
    auto apply_delete(Tree &tree, PT &P, double st = 0.0, double ed = 1.0, size_t d = 2){
        double plow[2], phigh[2];
        for (auto &pt: P){
            plow[0] = pt.x, plow[1] = pt.y;
            phigh[0] = pt.x, phigh[1] = pt.y;
            auto r = TimeRegion(plow, phigh, st, ed, d);
            tree->deleteData(r, pt.id);
        }           
    }

    /* 
        batch size: total number of modified points, include both insert and delete;
        insert ratio: from 0 to 10, determine the ratio of insert and delete points.
    */
    template<typename PT, typename RQ>
    void spatial_diff_test(PT &P, RQ &range_queries, size_t &batch_size, size_t &insert_ratio, size_t d = 2){
        uint32_t capacity = 32;
        id_type indexIdentifier;
        IStorageManager* storageManager = StorageManager::createNewMemoryStorageManager();
        ISpatialIndex* tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);

        /*  build tree */
        apply_insert(tree, P, 0, 0);

        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, batch_size);
        auto [P_insert, P_delete] = geobase::split_insert_delete(P_test, insert_ratio, P.size());
        // print_Pset_info(P_delete, "delete", 10);
		// print_Pset_info(P_insert, "insert", 10);
        
        /* apply delete, insert */
        apply_delete(tree, P_delete, 0, 1);
        apply_insert(tree, P_insert, 1, 2);

        parlay::sequence<size_t> addCnt(range_queries.size());
        parlay::sequence<size_t> removeCnt(range_queries.size());

        // auto [added, removed] = spatial_diff(tree, range_queries[0]);
        // geobase::print_mbr(range_queries[0]);
        // geobase::print_Pset_info(added, "added");
        // geobase::print_Pset_info(removed, "removed");


        auto avg_time = time_loop(
            3, 1.0,
            [&](){},
            [&](){
                for (size_t i = 0; i < range_queries.size(); i++){
                    auto [added, removed] = spatial_diff(tree, range_queries[i]);
                    // cout << i <<  " ok." << endl;
                    addCnt[i] = added.size();
                    removeCnt[i] = removed.size();
                    // cout << i <<  " finished: " << addCnt[i] << ", " << removeCnt[i] << endl;
                }
            },
            [&]{}
        );
        cout << fixed << setprecision(6) << "[MVRtree] spatial-diff time (avg): " << avg_time << endl;
        
        // string file_name = "../../src/output/mvr_spatial_diff-" + to_string(batch_size); 
        // ofstream spatialDiffOut(file_name);
        // for (size_t i = 0; i < range_queries.size(); i++){
        //     spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
        // }

        delete tree;
        delete storageManager;
    }

    void build_test(parlay::sequence<geobase::Point> &P, ISpatialIndex* tree, size_t d = 2){
        id_type id = 0;
        double plow[2], phigh[2];

        // string output_name = "mvrtree.idx";

        IStorageManager* storageManager = nullptr;
        uint32_t capacity = 32;
        id_type indexIdentifier;
        // SpatialIndex::TimeRegion r;
        // SpatialIndex::MVRTree::MVRTree* test;

        auto tree_build_avg = time_loop(
            1, 1.0, [&]() {
                delete tree;
                delete storageManager;
                storageManager = StorageManager::createNewMemoryStorageManager();
                // storageManager = StorageManager::createNewDiskStorageManager(output_name, 4096);
                tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);
                id = 0;
                m_mvrtree_size = 0;
            },
            [&]() {
                for (auto &pt: P){
                    // plow[0] = pt.first, plow[1] = pt.second;
                    // phigh[0] = pt.first, phigh[1] = pt.second;
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, plow[1] = pt.y;

                    auto r = TimeRegion(plow, phigh, 0, 0, d);
                    tree->insertData(0, nullptr, r, pt.id);
                }
                // cout << *tree << endl;
                // cout << "MVRtree used: " <<  memory_in_MB(m_mvrtree_size) << " MB" << endl;
            },
        [&](){
            // delete tree;
            // delete storageManager;
        } );

        cout << *tree << endl;
        cout << "MVRtree used: " <<  memory_in_MB(m_mvrtree_size) << " MB" << endl;

        cout << "[MVRtree-Build]: " << tree_build_avg << endl;
        delete tree;
        delete storageManager;
    }
    
    void batch_insert_test(parlay::sequence<geobase::Point> &P, int batch_percent = 10, size_t d = 2){
        /* build tree */
        IStorageManager* storageManager = nullptr;
        uint32_t capacity = 32;
        id_type indexIdentifier;
        ISpatialIndex* tree = nullptr;
        double plow[2], phigh[2];

        
        /* calculate inserted points */
        auto n = P.size();
        auto num_processed = batch_percent * n / 10000; 
        auto P_insert = P.substr(0, num_processed);
        for (size_t i = 0; i < P_insert.size(); i++){
            P_insert[i].id += n;
        }
        IStatistics* stats = nullptr;

        auto avg_time = time_loop(
            1, 1.0, 
            [&]() {
                storageManager= StorageManager::createNewMemoryStorageManager();
                tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);
                for (auto &pt: P){
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, phigh[1] = pt.y;    
                    auto r = TimeRegion(plow, phigh, 0, 0, d);
                    tree->insertData(0, nullptr, r, pt.id); 
                }
            },
            [&]() {    
                for (auto &pt: P_insert){
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, phigh[1] = pt.y;
                    auto r = TimeRegion(plow, phigh, 1.0, 1.0, d);
                    tree->insertData(0, nullptr, r, pt.id);
                }
            },
            [&](){
                tree->getStatistics(&stats);
                delete tree;
                delete storageManager;
            }
        );

        cout << *stats << endl;
        delete stats;
        cout << "Batch size: " << setprecision(2) << 1.0 * batch_percent / 100 << "% "; 
        cout << "[MVRtree-Insert]: " << avg_time << endl;
    }

    void batch_delete_test(parlay::sequence<geobase::Point> &P, int batch_percent, size_t d = 2){
        /* build tree */
        IStorageManager* storageManager = nullptr;
        uint32_t capacity = 32;
        id_type indexIdentifier;
        ISpatialIndex* tree = nullptr;
        double plow[2], phigh[2];
        
        /* calculate deleted points */
        auto n = P.size();
        auto num_processed = batch_percent * n / 10000; 
        auto P_delete = P.substr(0, num_processed);
        IStatistics* stats = nullptr;

        /* Note: for deletions, id start from 1 */
        auto avg_time = time_loop(
            1, 1.0, 
            [&]() {
                storageManager= StorageManager::createNewMemoryStorageManager();
                tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);
                for (auto &pt: P){
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, phigh[1] = pt.y;    
                    auto r = TimeRegion(plow, phigh, 0, 0, d);
                    tree->insertData(0, nullptr, r, pt.id); 
                }
            },
            [&]() {    
                for (auto &pt: P_delete){
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, phigh[1] = pt.y;
                    auto r = TimeRegion(plow, phigh, 0, 0, d);
                    tree->deleteData(r, pt.id);
                    

                }
            },
            [&](){
                tree->getStatistics(&stats);
                delete tree;
                delete storageManager;
            }
        );

        cout << *stats << endl;
        delete stats;
        cout << "Batch size: " << setprecision(2) << 1.0 * batch_percent / 100 << "% "; 
        cout << "[MVRtree-Delete]: " << avg_time << endl;
    }

    template <typename Vistor>
    void range_count_test(parlay::sequence<geobase::Point> &P, parlay::sequence<geobase::Bounding_Box> &query, parlay::sequence<size_t> &cnt, size_t d = 2){
        /* build tree */
        uint32_t capacity = 32;
        id_type indexIdentifier;
        IStorageManager* storageManager = StorageManager::createNewMemoryStorageManager();
        ISpatialIndex* tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);

        double plow[2], phigh[2];

        // IStatistics* stats = nullptr;

        for (auto &pt: P){
            plow[0] = pt.x, plow[1] = pt.y;
            phigh[0] = pt.x, phigh[1] = pt.y;
            auto r = TimeRegion(plow, phigh, 0, 0, d);
            tree->insertData(0, nullptr, r, pt.id);
        } 
        cout << *tree << endl;
        
        parlay::sequence<size_t> rangeCnt(query.size());

        // Vistor vis;
        auto avg_time = time_loop(
            3, 1.0, 
            [&]() {},
            [&]() {    
                for (size_t i = 0; i < query.size(); i++){
                    plow[0] = query[i].first.x, plow[1] = query[i].first.y;
                    phigh[0] = query[i].second.x, phigh[1] = query[i].second.y;
                    Vistor vis;
                    auto r = TimeRegion(plow, phigh, 0, 1, d);
                    cout << r << endl;
                    tree->intersectsWithQuery(r, vis);
                    rangeCnt[i] = vis.range_count_res;
                }
            },
            [&](){
                // tree->getStatistics(&stats);
            }
        );
        for (size_t i = 0; i < query.size(); i++){
            if (rangeCnt[i] != cnt[i]){
                cout << "[ERROR] incorrect result " << rangeCnt[i] << "-" << cnt[i] << endl;
            }
        }
        cout << fixed << setprecision(6) << "[MVRtree-Range Count]: " << avg_time << endl;
        delete tree;
        delete storageManager;
        
        /*          range query example        */
        // plow[0] = 19; plow[1] = 52;
        // phigh[0] = 22; phigh[1] = 54;
        // MyVisitor vis;
        // auto r = TimeRegion(plow, phigh, 0, 1, 2);
        // tree->intersectsWithQuery(r, vis);
    }

    template <typename Vistor>
    void knn_test(parlay::sequence<geobase::Point> &P, size_t k = 10, size_t q_num = 50000, size_t d = 2){
        /* build tree */
        uint32_t capacity = 32;
        id_type indexIdentifier;
        IStorageManager* storageManager = StorageManager::createNewMemoryStorageManager();
        ISpatialIndex* tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);

        double plow[2], phigh[2];

        // IStatistics* stats = nullptr;

        for (auto &pt: P){
            plow[0] = pt.x, plow[1] = pt.y;
            phigh[0] = pt.x, phigh[1] = pt.y;
            auto r = TimeRegion(plow, phigh, 0, 0, d);
            tree->insertData(0, nullptr, r, pt.id);
        } 
        // cout << *tree << endl;
        parlay::sequence<size_t> knn_sqrdis(q_num);

        auto avg_time = time_loop(
            1, 1.0, 
            [&]() {},
            [&]() {    
                for (size_t i = 0; i < q_num; i++){
                    plow[0] = P[i].x, plow[1] = P[i].y;
                    phigh[0] = P[i].x, phigh[1] = P[i].y;
                    Vistor vis;
                    auto r = TimeRegion(plow, phigh, -0.0, 1.0, d);
                    tree->nearestNeighborQuery(k, r, vis);
                    knn_sqrdis[i] = std::numeric_limits<double>::min();
                    // cout << "[INFO] res size = " << vis.ids.size() << endl;
                    for (auto &id: vis.ids){
                        auto cur_dis = geobase::point_point_sqrdis(P[id], P[i]);
                        if (cur_dis > knn_sqrdis[i]) knn_sqrdis[i] = cur_dis;
                    }
                    // int xx; cin >> xx;
                }
            },
            [&](){
            }
        );

        cout << fixed << setprecision(6) << "[MVRtree-KNN]: " << avg_time << endl;

        auto output_file = "mvr-knn.res";
		ofstream RES(output_file);
		for (size_t i = 0; i < q_num; i++){
			RES << knn_sqrdis[i] << endl;
		}

    }


    void test(){
        cout << "ww" << endl;
    }
}