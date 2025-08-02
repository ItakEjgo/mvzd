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

namespace MV3RTest{
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

    double memory_in_MB(size_t num_of_bytes){
         return static_cast<double>(num_of_bytes) / 1048576.0;
    }

    template<typename MVRTree, typename RTree, typename PT>
    auto apply_insert(MVRTree &mvrtree, RTree &rtree, PT &P, double st = 0.0, double ed = 0.0){
        double plow[2], phigh[2];
        double plow_3d[3], phigh_3d[3];
        double inf_time = 1000;

        for (auto &pt: P){
            plow[0] = pt.x, plow[1] = pt.y;
            phigh[0] = pt.x, phigh[1] = pt.y;

            plow_3d[0] = pt.x, plow_3d[1] = pt.y, plow_3d[2] = st;
            phigh_3d[0] = pt.x, phigh_3d[1] = pt.y, phigh_3d[2] = inf_time;

            auto r = TimeRegion(plow, phigh, st, ed, 2);
            auto r3 = Region(plow_3d, phigh_3d, 3);

            mvrtree->insertData(0, nullptr, r, pt.id);  //  mvrtree is 2 dimensional
            rtree->insertData(0, nullptr, r3, pt.id);   //  rtree is 3 dimensional
        }           
    }

    template<typename MVRTree, typename RTree, typename PT>
    auto apply_delete(MVRTree &mvrtree, RTree &rtree, PT &P, double st = 0.0, double ed = 1.0){
        double plow[2], phigh[2];
        double plow_3d[3], phigh_3d[3];
        double inf_time = 1000;

        for (auto &pt: P){
            plow[0] = pt.x, plow[1] = pt.y;
            phigh[0] = pt.x, phigh[1] = pt.y;
            plow_3d[0] = pt.x, plow_3d[1] = pt.y, plow_3d[2] = st;
            phigh_3d[0] = pt.x, phigh_3d[1] = pt.y, phigh_3d[2] = inf_time;

            auto r = TimeRegion(plow, phigh, st, ed, 2);
            auto r3 = Region(plow_3d, phigh_3d, 3);
            
            mvrtree->deleteData(r, pt.id);
            rtree->deleteData(r3, pt.id);
        }           
    }

    template<typename PSeq>
    auto build_mv3rtree(IStorageManager* &storageManager_mvrtree, ISpatialIndex* &mvrtree,
                       IStorageManager* &storageManager_rtree, ISpatialIndex* &rtree, PSeq &P){
        /* build tree */
        uint32_t capacity = 32;

        id_type indexIdentifier_mvrtree, indexIdentifier_rtree;

        storageManager_mvrtree = StorageManager::createNewMemoryStorageManager();
        storageManager_rtree = StorageManager::createNewMemoryStorageManager();
        mvrtree = MVRTree::createNewMVRTree(*storageManager_mvrtree, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier_mvrtree);
        rtree = RTree::createNewRTree(*storageManager_rtree, 0.7, capacity, capacity, 3, SpatialIndex::RTree::RV_RSTAR, indexIdentifier_rtree);

        apply_insert(mvrtree, rtree, P);
    }

    //  mvrtree + 3d rtree
    void build_test(parlay::sequence<geobase::Point> &P, size_t d = 2){

        IStorageManager* storageManager_mvrtree = nullptr;
        ISpatialIndex* mvrtree = nullptr;
        IStorageManager* storageManager_rtree = nullptr;
        ISpatialIndex* rtree = nullptr;

        double before_mem = 0.0, after_mem = 0.0;        
        /* Build initial version */
        auto build_avg = time_loop(
            1, -1.0, [&]() {},
            [&]() {
                before_mem = getProcessRSS_MB();
                build_mv3rtree(storageManager_mvrtree, mvrtree, storageManager_rtree, rtree, P);
                after_mem = getProcessRSS_MB();
            },
            [&](){});

        // cout << *mvrtree << endl;
        // cout << *rtree << endl;

        cout << "[MV3Rtree build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;
        cout << "[MV3Rtree build memory]: " << fixed << setprecision(6) << after_mem - before_mem << " MB" << endl;
        cout << fixed << setprecision(2) << "[Before/After Mem]: " << before_mem << ", " << after_mem << endl;

        delete mvrtree;
        delete rtree;
        delete storageManager_mvrtree;
        delete storageManager_rtree;
    }
    
    void batch_insert_test(parlay::sequence<geobase::Point> &P, parlay::sequence<size_t> batch_sizes, size_t d = 2){
        /* build tree */
        IStorageManager* storageManager_mvrtree = nullptr;
        IStorageManager* storageManager_rtree = nullptr;
        ISpatialIndex* mvrtree = nullptr;
        ISpatialIndex* rtree = nullptr;
        
        auto rand_p = geobase::shuffle_point(P);
        auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        if (max_batch_size > P.size()) max_batch_size = P.size();
        
        /* calculate inserted points */

        auto P_insert = rand_p.substr(0, max_batch_size);
        for (size_t i = 0; i < P_insert.size(); i++){
            P_insert[i].id += P.size();
        }

        build_mv3rtree(storageManager_mvrtree, mvrtree, storageManager_rtree, rtree, P);
        cout << *mvrtree << endl;
        cout << *rtree << endl;

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
                    apply_insert(mvrtree, rtree, cur_insert);
                },
                [&](){}
            );
            tot_t += avg_time;
            cout << "[INFO] batch_size: " << batch_size << endl;
            cout << *mvrtree << endl;
            cout << *rtree << endl;
            cout << "[MV3RTree] batch insert time = " << fixed << setprecision(6) << tot_t << endl;

            l = r;
        }
        delete mvrtree;
        delete rtree;
        delete storageManager_mvrtree;
        delete storageManager_rtree;
    }

    void batch_delete_test(parlay::sequence<geobase::Point> &P, parlay::sequence<size_t> batch_sizes, size_t d = 2){
              /* build tree */
        IStorageManager* storageManager_mvrtree = nullptr;
        IStorageManager* storageManager_rtree = nullptr;
        ISpatialIndex* mvrtree = nullptr;
        ISpatialIndex* rtree = nullptr;
        
        auto rand_p = geobase::shuffle_point(P);
        auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        if (max_batch_size > P.size()) max_batch_size = P.size();
        
        /* calculate inserted points */

        auto P_delete = rand_p.substr(0, max_batch_size);

        build_mv3rtree(storageManager_mvrtree, mvrtree, storageManager_rtree, rtree, P);
        cout << *mvrtree << endl;
        cout << *rtree << endl;

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
                    apply_delete(mvrtree, rtree, cur_insert);
                },
                [&](){}
            );
            tot_t += avg_time;
            cout << "[INFO] batch_size: " << batch_size << endl;
            cout << *mvrtree << endl;
            cout << *rtree << endl;
            cout << "[MV3RTree] batch delete time = " << fixed << setprecision(6) << tot_t << endl;

            l = r;
        }
        delete mvrtree;
        delete rtree;
        delete storageManager_mvrtree;
        delete storageManager_rtree;
    }

    template <typename Vistor>
    void range_count_test(parlay::sequence<geobase::Point> &P, parlay::sequence<geobase::Bounding_Box> &query, parlay::sequence<size_t> &cnt, size_t d = 2){
        /* build tree */
        uint32_t capacity = 32;
        id_type indexIdentifier;
        IStorageManager* storageManager = StorageManager::createNewMemoryStorageManager();
        ISpatialIndex* tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);
        const size_t n = P.size();

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
                // parlay::parallel_for(0, query.size(), [&](int i){
                //     double plow_cur[2], phigh_cur[2];
                //     plow_cur[0] = query[i].first.x, plow_cur[1] = query[i].first.y;
                //     phigh_cur[0] = query[i].second.x, phigh_cur[1] = query[i].second.y;

                //     // plow[0] = query[i].first.x, plow[1] = query[i].first.y;
                //     // phigh[0] = query[i].second.x, phigh[1] = query[i].second.y;

                //     Vistor vis;
                //     auto r = TimeRegion(plow_cur, phigh_cur, 0, 1, d);
                //     // auto r = TimeRegion(plow, phigh, 0, 1, d);
                //     tree->intersectsWithQuery(r, vis);
                //     rangeCnt[i] = vis.range_count_res;
                // });

                for (auto i = 0; i < query.size(); i++){
                    plow[0] = query[i].first.x, plow[1] = query[i].first.y;
                    phigh[0] = query[i].second.x, phigh[1] = query[i].second.y;
                    Vistor vis;
                    auto r = TimeRegion(plow, phigh, 0, 1, d);
                    tree->intersectsWithQuery(r, vis);
                    rangeCnt[i] = vis.range_count_res;
                }
            },
            [&](){
                // tree->getStatistics(&stats);
            }
        );
        for (auto i = 0; i < query.size(); i++){
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


        IStorageManager* storageManager_mvrtree = nullptr;
        ISpatialIndex* mvrtree = nullptr;
        IStorageManager* storageManager_rtree = nullptr;
        ISpatialIndex* rtree = nullptr;

        double before_mem = 0.0, after_mem = 0.0;        
        /* Build initial version */
        auto build_avg = time_loop(
            1, -1.0, [&]() {},
            [&]() {
                before_mem = getProcessRSS_MB();
                build_mv3rtree(storageManager_mvrtree, mvrtree, storageManager_rtree, rtree, P);
                after_mem = getProcessRSS_MB();
            },
            [&](){});

        cout << *mvrtree << endl;
        cout << ">------------------------------------------<" << endl;
        cout << *rtree << endl;

        cout << "[MV3Rtree init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;
        cout << "[MV3Rtree init build memory]: " << fixed << setprecision(6) << after_mem - before_mem << " MB" << endl;

        for (auto i = 0; i < version_num; i++){
            cout << "dealing with version " << i + 1 << ":" << endl;

            auto commit_avg = time_loop(
                1, -1.0, 
                [&]() {},
                [&]() {
                    before_mem = getProcessRSS_MB();
                    apply_delete(mvrtree, rtree, P_delete[i], 1.0 * i);
                    apply_insert(mvrtree, rtree, P_insert[i], 1.0 * i + 1.0);
                    after_mem = getProcessRSS_MB();
                },
                [&](){
                });
            // cout << *mvrtree << endl;

            cout << *mvrtree << endl;
            cout << ">------------------------------------------<" << endl;
            cout << *rtree << endl;
            
            cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;
            cout << "[Rtree new version memory]: " << fixed << setprecision(6) << after_mem - before_mem << " MB" << endl;
        }

        delete mvrtree;
        delete rtree;
        delete storageManager_mvrtree;
        delete storageManager_rtree;
    }


    void test(){
        cout << "ww" << endl;
    }
}