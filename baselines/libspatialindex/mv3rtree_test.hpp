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

    double memory_in_MB(size_t num_of_bytes){
         return static_cast<double>(num_of_bytes) / 1048576.0;
    }

    //  mvrtree + 3d rtree
    void build_test(parlay::sequence<geobase::Point> &P, ISpatialIndex* mvrtree, ISpatialIndex* rtree, size_t d = 2){
        double plow[2], phigh[2], plow_3d[3], phigh_3d[3];

        // string output_name = "mvrtree.idx";

        uint32_t capacity = 32;
        id_type indexIdentifier_mvrtree, indexIdentifier_rtree;
        IStorageManager* storageManager_mvrtree = nullptr;
        IStorageManager* storageManager_rtree = nullptr;
        auto inf_time = 1000; 

        auto tree_build_avg = time_loop(
            1, 1.0, [&]() {
                delete mvrtree;
                delete rtree;
                delete storageManager_mvrtree;
                delete storageManager_rtree;
                storageManager_mvrtree = StorageManager::createNewMemoryStorageManager();
                storageManager_rtree = StorageManager::createNewMemoryStorageManager();
                // storageManager = StorageManager::createNewDiskStorageManager(output_name, 4096);
                mvrtree = MVRTree::createNewMVRTree(*storageManager_mvrtree, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier_mvrtree);
                rtree = RTree::createNewRTree(*storageManager_rtree, 0.7, capacity, capacity, 3, SpatialIndex::RTree::RV_RSTAR, indexIdentifier_rtree);
                m_mvrtree_size = 0;
            },
            [&]() {
                for (auto &pt: P){
                    // plow[0] = pt.first, plow[1] = pt.second;
                    // phigh[0] = pt.first, phigh[1] = pt.second;
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, plow[1] = pt.y;
                    plow_3d[0] = pt.x, plow_3d[1] = pt.y, plow_3d[2] = 0;
                    phigh_3d[0] = pt.x, phigh_3d[1] = pt.y, phigh_3d[2] = inf_time;
                    auto r = TimeRegion(plow, phigh, 0, 0, 2);
                    auto r3 = Region(plow_3d, phigh_3d, 3);
                    mvrtree->insertData(0, nullptr, r, pt.id);  //  mvrtree is 2 dimensional
                    rtree->insertData(0, nullptr, r3, pt.id);   //  rtree is 3 dimensional
                }
            },
        [&](){
        } );

        cout << "[MVRtree info]: " << endl;
        cout << *mvrtree << endl;
        cout << "[3DRtree info]: " << endl;
        cout << *rtree << endl;
        cout << "[MV3Rtree-Build]: " << tree_build_avg << endl;

        delete mvrtree;
        delete rtree;
        delete storageManager_mvrtree;
        delete storageManager_rtree;
    }
    
    void batch_insert_test(parlay::sequence<geobase::Point> &P, int batch_percent = 10, size_t d = 2){
        /* build tree */
        double plow[2], phigh[2], plow_3d[3], phigh_3d[3];
        // string output_name = "mvrtree.idx";
        uint32_t capacity = 32;
        id_type indexIdentifier_mvrtree, indexIdentifier_rtree;
        IStorageManager* storageManager_mvrtree = nullptr;
        IStorageManager* storageManager_rtree = nullptr;
        ISpatialIndex* mvrtree = nullptr;
        ISpatialIndex* rtree = nullptr;
        auto inf_time = 1000; 
        
        /* calculate inserted points */
        auto n = P.size();
        auto num_processed = batch_percent * n / 10000; 
        auto P_insert = P.substr(0, num_processed);
        for (size_t i = 0; i < P_insert.size(); i++){
            P_insert[i].id += n;
        }

        auto avg_time = time_loop(
            1, 1.0, 
            [&]() {
                delete mvrtree;
                delete rtree;
                delete storageManager_mvrtree;
                delete storageManager_rtree;
                storageManager_mvrtree = StorageManager::createNewMemoryStorageManager();
                storageManager_rtree = StorageManager::createNewMemoryStorageManager();
                mvrtree = MVRTree::createNewMVRTree(*storageManager_mvrtree, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier_mvrtree);
                rtree = RTree::createNewRTree(*storageManager_mvrtree, 0.7, capacity, capacity, 3, SpatialIndex::RTree::RV_RSTAR, indexIdentifier_rtree);

                for (auto &pt: P){
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, phigh[1] = pt.y;    
                    plow_3d[0] = pt.x, plow_3d[1] = pt.y, plow_3d[2] = 0;
                    phigh_3d[0] = pt.x, phigh_3d[1] = pt.y, phigh_3d[2] = inf_time;
                    auto r = TimeRegion(plow, phigh, 0, 0, d);
                    auto r3 = Region(plow_3d, phigh_3d, 3);
                    mvrtree->insertData(0, nullptr, r, pt.id); 
                    rtree->insertData(0, nullptr, r3, pt.id);
                }

            },
            [&]() {    
                for (auto &pt: P_insert){
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, phigh[1] = pt.y;
                    plow_3d[0] = pt.x, plow_3d[1] = pt.y, plow_3d[2] = 1;
                    phigh_3d[0] = pt.x, phigh_3d[1] = pt.y, phigh_3d[2] = inf_time;
                    auto r = TimeRegion(plow, phigh, 1.0, 1.0, 2);
                    auto r3 = Region(plow_3d, phigh_3d, 3);
                    mvrtree->insertData(0, nullptr, r, pt.id);
                    rtree->insertData(0, nullptr, r3, pt.id);
                }
            },
            [&](){
                // tree->getStatistics(&stats);
            }
        );

        cout << "[MVRtree info]: " << endl;
        cout << *mvrtree << endl;
        cout << "[3DRtree info]: " << endl;
        cout << *rtree << endl;

        cout << "Batch size: " << setprecision(2) << 1.0 * batch_percent / 100 << "% "; 
        cout << "[MV3Rtree-Insert]: " << avg_time << endl;

        delete mvrtree;
        delete rtree;
        delete storageManager_mvrtree;
        delete storageManager_rtree;
    }

    void batch_delete_test(parlay::sequence<geobase::Point> &P, int batch_percent, size_t d = 2){
        /* build tree */
        double plow[2], phigh[2], plow_3d[3], phigh_3d[3];
        // string output_name = "mvrtree.idx";
        uint32_t capacity = 32;
        id_type indexIdentifier_mvrtree, indexIdentifier_rtree;
        IStorageManager* storageManager_mvrtree = nullptr;
        IStorageManager* storageManager_rtree = nullptr;
        ISpatialIndex* mvrtree = nullptr;
        ISpatialIndex* rtree = nullptr;
        auto inf_time = 1000; 
        
        /* calculate inserted points */
        auto n = P.size();
        auto num_processed = batch_percent * n / 10000; 
        auto P_delete = P.substr(0, num_processed);

        /* Note: for deletions, id start from 1 */
        auto avg_time = time_loop(
            1, 1.0, 
            [&]() {
                delete mvrtree;
                delete rtree;
                delete storageManager_mvrtree;
                delete storageManager_rtree;
                storageManager_mvrtree = StorageManager::createNewMemoryStorageManager();
                storageManager_rtree = StorageManager::createNewMemoryStorageManager();
                mvrtree = MVRTree::createNewMVRTree(*storageManager_mvrtree, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier_mvrtree);
                rtree = RTree::createNewRTree(*storageManager_mvrtree, 0.7, capacity, capacity, 3, SpatialIndex::RTree::RV_RSTAR, indexIdentifier_rtree);

                for (auto &pt: P){
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, phigh[1] = pt.y;    
                    plow_3d[0] = pt.x, plow_3d[1] = pt.y, plow_3d[2] = 0;
                    phigh_3d[0] = pt.x, phigh_3d[1] = pt.y, phigh_3d[2] = inf_time;
                    auto r = TimeRegion(plow, phigh, 0, 0, 2);
                    auto r3 = Region(plow_3d, phigh_3d, 3);
                    mvrtree->insertData(0, nullptr, r, pt.id); 
                    rtree->insertData(0, nullptr, r3, pt.id);
                }
            },
            [&]() {    
                for (auto &pt: P_delete){
                    plow[0] = pt.x, plow[1] = pt.y;
                    phigh[0] = pt.x, phigh[1] = pt.y;
                    plow_3d[0] = pt.x, plow_3d[1] = pt.y, plow_3d[2] = 0;
                    phigh_3d[0] = pt.x, phigh_3d[1] = pt.y, phigh_3d[2] = inf_time;
                    auto r = TimeRegion(plow, phigh, 0, 0, 2);
                    auto r3 = Region(plow_3d, phigh_3d, 3);
                    mvrtree->deleteData(r, pt.id);
                    rtree->deleteData(r3, pt.id);
                }
            },
            [&](){
            }
        );

        cout << "[MVRtree info]: " << endl;
        cout << *mvrtree << endl;
        cout << "[3DRtree info]: " << endl;
        cout << *rtree << endl;

        cout << "Batch size: " << setprecision(2) << 1.0 * batch_percent / 100 << "% "; 
        cout << "[MV3Rtree-Delete]: " << avg_time << endl;

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


    void test(){
        cout << "ww" << endl;
    }
}