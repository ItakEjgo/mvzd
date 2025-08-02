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

    // Get resident set size (RSS) in MB on Linux
    double getProcessRSS_MB() {
        ifstream stat_stream("/proc/self/statm");
        size_t size = 0, resident = 0;
        stat_stream >> size >> resident;
        stat_stream.close();

        long page_size = sysconf(_SC_PAGESIZE); // bytes per page
        size_t rss_bytes = resident * page_size;

        return static_cast<double>(rss_bytes) / (1024.0 * 1024.0); // convert to MB
    }

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

    template<typename PSeq>
    auto build_mvrtree(IStorageManager* &storageManager, ISpatialIndex* &tree, PSeq &P, size_t d = 2){
        /* build tree */
        uint32_t capacity = 32;
        id_type indexIdentifier;

        storageManager = StorageManager::createNewMemoryStorageManager();
        tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);

        double plow[2], phigh[2];
        for (auto &pt: P){
            plow[0] = pt.x, plow[1] = pt.y;
            phigh[0] = pt.x, phigh[1] = pt.y;
            auto r = TimeRegion(plow, phigh, 0, 0, d);
            tree->insertData(0, nullptr, r, pt.id);
        } 
        // cout << *tree << endl;
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

    void build_test(parlay::sequence<geobase::Point> &P, size_t d = 2){

        IStorageManager* storageManager = nullptr;
        ISpatialIndex* mvrtree = nullptr;
        double before_mem = 0.0, after_mem = 0.0;        
        /* Build initial version */
        auto build_avg = time_loop(
            1, -1.0, [&]() {
                delete mvrtree;
                delete storageManager;
            },
            [&]() {
                before_mem = getProcessRSS_MB();
                build_mvrtree(storageManager, mvrtree, P);
                after_mem = getProcessRSS_MB();
            },
            [&](){});

        // cout << *mvrtree << endl;

        cout << "[MVRtree build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;
        cout << "[MVRtree build memory]: " << fixed << setprecision(6) << after_mem - before_mem << " MB" << endl;
        cout << fixed << setprecision(2) << "[Before/After Mem]: " << before_mem << ", " << after_mem << endl;

        delete mvrtree;
        delete storageManager;
    }
    
    void batch_insert_test(parlay::sequence<geobase::Point> &P, parlay::sequence<size_t> &batch_sizes, size_t d = 2){
        /* build tree */
        IStorageManager* storageManager = nullptr;
        ISpatialIndex* tree = nullptr;

        auto rand_p = geobase::shuffle_point(P);
        auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        if (max_batch_size > P.size()) max_batch_size = P.size();
        
        /* calculate inserted points */

        auto P_insert = rand_p.substr(0, max_batch_size);
        for (size_t i = 0; i < P_insert.size(); i++){
            P_insert[i].id += P.size();
        }

        build_mvrtree(storageManager, tree, P);
        cout << *tree << endl;

        size_t l = 0, r = 0;
        double tot_t = 0;
        for (auto &batch_size: batch_sizes){
            if (batch_size > P.size()) break;
            cout << ">---------------------------------------------<" << endl;
            r = batch_size;
            auto cur_insert = P_insert.substr(l, r - l);

            auto avg_time = time_loop(
                1, -1.0, 
                [&]() {},
                [&]() {    
                    apply_insert(tree, cur_insert);
                },
                [&](){}
            );
            tot_t += avg_time;
            cout << "[INFO] batch_size: " << batch_size << endl;
            cout << *tree << endl;
            cout << "[MVRTree] batch insert time = " << fixed << setprecision(6) << tot_t << endl;

            l = r;
        }

        delete tree;
        delete storageManager;
    }

    void batch_delete_test(parlay::sequence<geobase::Point> &P, parlay::sequence<size_t> &batch_sizes, size_t d = 2){
          /* build tree */
        IStorageManager* storageManager = nullptr;
        ISpatialIndex* tree = nullptr;

        auto rand_p = geobase::shuffle_point(P);
        auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        if (max_batch_size > P.size()) max_batch_size = P.size();
        
        /* calculate inserted points */

        auto P_delete = rand_p.substr(0, max_batch_size);

        build_mvrtree(storageManager, tree, P);
        cout << *tree << endl;

        size_t l = 0, r = 0;
        double tot_t = 0;
        for (auto &batch_size: batch_sizes){
            if (batch_size > P.size()) break;
            cout << ">---------------------------------------------<" << endl;
            r = batch_size;
            auto cur_insert = P_delete.substr(l, r - l);

            auto avg_time = time_loop(
                1, -1.0, 
                [&]() {},
                [&]() {    
                    apply_delete(tree, cur_insert);
                },
                [&](){}
            );
            tot_t += avg_time;
            cout << "[INFO] batch_size: " << batch_size << endl;
            cout << *tree << endl;
            cout << "[MVRTree] batch delete time = " << fixed << setprecision(6) << tot_t << endl;

            l = r;
        }

        delete tree;
        delete storageManager;

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
        // cout << *tree << endl;
        
        parlay::sequence<size_t> rangeCnt(query.size());

        // Vistor vis;
        for (size_t i = 0; i < query.size(); i++){
            auto avg_time = time_loop(
                3, 1.0, 
                [&]() {},
                [&]() {    
                        plow[0] = query[i].first.x, plow[1] = query[i].first.y;
                        phigh[0] = query[i].second.x, phigh[1] = query[i].second.y;
                        Vistor vis;
                        auto r = TimeRegion(plow, phigh, 0, 1, d);
                        // cout << r << endl;
                        tree->intersectsWithQuery(r, vis);
                        rangeCnt[i] = vis.range_count_res;
                },
                [&](){
                    // tree->getStatistics(&stats);
                }
            );
            if (rangeCnt[i] != cnt[i]){
                cout << "[ERROR] incorrect result " << rangeCnt[i] << "-" << cnt[i] << endl;
            }
            else{
                cout << fixed << setprecision(6) << cnt[i] << " " << avg_time << endl;
            }
        }

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
    void knn_test(ISpatialIndex* tree, parlay::sequence<geobase::Point> &P, size_t k = 10, size_t q_num = 50000, size_t d = 2){
        /* build tree */
        double plow[2], phigh[2];
        parlay::sequence<size_t> knn_sqrdis(q_num);
        // q_num = min(q_num, P.size());
        auto avg_time = time_loop(
            3, 1.0, 
            [&]() {},
            [&]() {    
                for (size_t i = 0; i < q_num; i++){
                    plow[0] = P[i].x, plow[1] = P[i].y;
                    phigh[0] = P[i].x, phigh[1] = P[i].y;
                    Vistor vis;
                    auto r = TimeRegion(plow, phigh, -0.0, 1.0, d);
                    tree->nearestNeighborQuery(k, r, vis);
                    /* Correctness Check */
                    // knn_sqrdis[i] = std::numeric_limits<double>::min();
                    // for (auto &id: vis.ids){
                    //     auto cur_dis = geobase::point_point_sqrdis(P[id], P[i]);
                    //     if (cur_dis > knn_sqrdis[i]) knn_sqrdis[i] = cur_dis;
                    // }
                }
            },
            [&](){
            }
        );

        cout << fixed << setprecision(6) << "[MVRtree-KNN]: " << avg_time << endl;
        // auto output_file = "mvr-knn.res";
		// ofstream RES(output_file);
		// for (size_t i = 0; i < q_num; i++){
		// 	RES << knn_sqrdis[i] << endl;
		// }
    }

    template<typename PT>
    void multi_version_test(PT P, string dir, int start_year = 14, int version_num = 5){
        auto cur_year = start_year;
    
        parlay::sequence<geobase::Point> P_delete[version_num], P_insert[version_num], P_update[version_num], P_updove[version_num];

        for (auto i = 0; i != version_num; i++){
            auto delete_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-delete.txt";
            auto insert_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-insert.txt";
            auto update_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update.txt";
            auto updove_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update_remove.txt";

            ifstream fin_delete(delete_file_name);
            ifstream fin_insert(insert_file_name);
            ifstream fin_update(update_file_name);
            ifstream fin_updove(updove_file_name);

            auto delete_mbr = read_pts(P_delete[i], fin_delete, 1);
            auto insert_mbr = read_pts(P_insert[i], fin_insert, 1);
            auto update_mbr = read_pts(P_update[i], fin_update, 1);
            auto updove_mbr = read_pts(P_updove[i], fin_updove, 1);

            P_delete[i].append(P_updove[i]);
            P_insert[i].append(P_update[i]);

            cur_year += 1;
            delete_mbr = insert_mbr; // useless, just remove warning
            update_mbr = updove_mbr;
        }


        IStorageManager* storageManager = nullptr;
        ISpatialIndex* mvrtree = nullptr;
        double before_mem = 0.0, after_mem = 0.0;        
        /* Build initial version */
        auto build_avg = time_loop(
            1, -1.0, [&]() {
                delete mvrtree;
                delete storageManager;
            },
            [&]() {
                before_mem = getProcessRSS_MB();
                build_mvrtree(storageManager, mvrtree, P);
                after_mem = getProcessRSS_MB();
            },
            [&](){});

        // cout << *mvrtree << endl;

        cout << "[MVRtree init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;
        cout << "[MVRtree init build memory]: " << fixed << setprecision(6) << after_mem - before_mem << " MB" << endl;

        for (auto i = 0; i < version_num; i++){
            cout << "dealing with version " << i + 1 << ":" << endl;

            auto commit_avg = time_loop(
                1, -1.0, 
                [&]() {},
                [&]() {
                    before_mem = getProcessRSS_MB();
                    apply_delete(mvrtree, P_delete[i], 1.0 * i);
                    apply_insert(mvrtree, P_insert[i], 1.0 * i + 1.0);
                    after_mem = getProcessRSS_MB();
                },
                [&](){
                });
            // cout << *mvrtree << endl;
            cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;
            cout << "[Rtree new version memory]: " << fixed << setprecision(6) << after_mem - before_mem << " MB" << endl;
        }
        
        delete mvrtree;
        delete storageManager;
    }

}