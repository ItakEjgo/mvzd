#include <iostream>
#include <cstring>
#include <spatialindex/SpatialIndex.h>
#include <chrono>
#include "mvrtree_test.hpp"
#include "mv3rtree_test.hpp"
#include "helper/parse_command_line.h"
#include "helper/geobase.h"
#include <typeinfo>


using namespace SpatialIndex;
using namespace std;

size_t maxSize = 0;
vector<pair<TimeRegion, id_type> > removed, added;

// example of a Visitor pattern.
// see RTreeQuery for a more elaborate example.
// class MyVisitor : public IVisitor
// {
// public:
// 	void visitNode(const INode& /* n */) override {}

// 	void visitData(const IData& d) override
// 	{
// 		std::cout << d.getIdentifier() << std::endl;
// 			// the ID of this data entry is an answer to the query. I will just print it to stdout.
// 	}

// 	void visitData(std::vector<const IData*>& /* v */) override {}
// };

class MyVisitor : public IVisitor {
public:
    size_t range_count_res = 0;

    void visitNode(const INode& /* n */) override {}
    void visitData(const IData& d) override {
        range_count_res += 1;
        // cout << "Data ID: " << d.getIdentifier() << endl;
		// IShape* pS;
		// d.getShape(&pS);
        // Point p;
        // pS->getCenter(p);
        // cout << p << endl;
        // delete pS;
    }
    void visitData(std::vector<const IData*>& v) override {

        range_count_res += v.size();
        // for (const auto& d : v) {
        //     cout << "Data ID: " << d->getIdentifier() << endl;
        // }
    }
};

class NNVisitor : public IVisitor{
public:
    parlay::sequence<size_t> ids;

    void visitNode(const INode&) override{}
    void visitData(const IData &d) override{
        ids.push_back(d.getIdentifier());
    }
    void visitData(std::vector<const IData*> &v) override{
        for (auto d: v){
            ids.push_back(d->getIdentifier());
        }
    }
};

void build_test(char** argv){
    IStorageManager* storageManager = StorageManager::createNewMemoryStorageManager();
    uint32_t capacity = 32;
    id_type indexIdentifier;
    ISpatialIndex* tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);
    
    id_type id = 0;
    double plow[2], phigh[2];
    ifstream fin(argv[1]);
    int n, d;
    
    fin >> n >> d;
    vector<pair<double, double> > pts(n);

    for (auto i = 0; i < n; i++){
        fin >> pts[i].first >> pts[i].second;
    }

    auto t_s = chrono::high_resolution_clock::now();
    for (auto i = 0; i < n; i++){
        plow[0] = pts[i].first, plow[1] = pts[i].second;
        phigh[0] = pts[i].first, phigh[1] = pts[i].second;
        if (i < 2){
            auto r = TimeRegion(plow, phigh, 0, 0, d);
            tree->insertData(0, nullptr, r, id++);
        }
        else{
            auto r = TimeRegion(plow, phigh, 1, 1, d);
            tree->insertData(0, nullptr, r, id++);
        }
        // tree->insertData_MV3Rtree(0, nullptr, r, id++, removed, added);
    }

    MyVisitor vis;
    // auto qp = TimeRegion(plow, phigh, 0, 0.1, 2);
    // Point qp = Point(plow, 2);
    TimeRegion r = TimeRegion(plow, phigh, (double)1.0, (double)1.1, 2);
    tree->nearestNeighborQuery(3, r, vis);

    // plow[0] = 18, plow[1] = 51;
    // phigh[0] = 24, phigh[1] = 54;
    // TimeRegion r = TimeRegion(plow, phigh, (double)0, (double)0.1, 2);
    // tree->intersectsWithQuery(r, vis);

    auto t_e = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t_e - t_s);

    // cout << "build time: " << 1.0 * duration.count() / 1000 << endl;
    // cout << *tree << endl;

    delete tree;
    delete storageManager;
}

void insert_test(char** argv){
    IStorageManager* storageManager = StorageManager::createNewMemoryStorageManager();
    uint32_t capacity = 20;
    id_type indexIdentifier;
    ISpatialIndex* tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);
    
    double x, y;
    double plow[2], phigh[2];
    ifstream fin(argv[1]);
    int n, d;
    
    fin >> n >> d;
    vector<pair<double, double> > pts1, pts2;
    auto num_processed = 10 * n / 100;
    for (auto i = 0; i < n; i++){
        fin >> x >> y;
        if (i < num_processed){
            pts1.emplace_back(x, y);
        }
        else{
            pts2.emplace_back(x, y);
        }
        // fin >> pts[i].first >> pts[i].second;
    }

    for (size_t i = 0; i < pts2.size(); i++){
        plow[0] = pts2[i].first, plow[1] = pts2[i].second;
        phigh[0] = pts2[i].first, phigh[1] = pts2[i].second;
        auto r = TimeRegion(plow, phigh, 0, 0, d);
        tree->insertData(0, nullptr, r, num_processed + i);
    }   //  build for (1-x)%

    auto t_s = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < pts1.size(); i++){
        plow[0] = pts1[i].first, plow[1] = pts1[i].second;
        phigh[0] = pts1[i].first, phigh[1] = pts1[i].second;
        auto r = TimeRegion(plow, phigh, 0, 0, d);
        tree->insertData(0, nullptr, r, num_processed);
    }
    auto t_e = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t_e - t_s);
    cout << "num inserted: " << pts1.size() << endl;
    cout << "insert time: " << 1.0 * duration.count() / 1000 << endl;

    /*          range query example        */
    // plow[0] = 19; plow[1] = 52;
    // phigh[0] = 22; phigh[1] = 54;
    // MyVisitor vis;
    // auto r = TimeRegion(plow, phigh, 0, 1, 2);
    // tree->intersectsWithQuery(r, vis);

    delete tree;
    delete storageManager;
}

void delete_test(char** argv){
    IStorageManager* storageManager = StorageManager::createNewMemoryStorageManager();
    uint32_t capacity = 20;
    id_type indexIdentifier;
    ISpatialIndex* tree = MVRTree::createNewMVRTree(*storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);
    
    double x, y;
    double plow[2], phigh[2];
    ifstream fin(argv[1]);
    int n, d;
    
    fin >> n >> d;
    vector<pair<double, double> > pts1, pts2;
    auto num_processed = 10 * n / 100;
    for (auto i = 0; i < n; i++){
        fin >> x >> y;
        if (i < num_processed){
            // pts1.emplace_back(x, y);
            pts2.emplace_back(x, y);
        }
        else{
            pts2.emplace_back(x, y);
        }
        // fin >> pts[i].first >> pts[i].second;
    }
    for (size_t i = 0; i < pts1.size(); i++){
        pts2.emplace_back(pts1[i]);
    }

    for (size_t i = 0; i < pts2.size(); i++){
        plow[0] = pts2[i].first, plow[1] = pts2[i].second;
        phigh[0] = pts2[i].first, phigh[1] = pts2[i].second;
        auto r = TimeRegion(plow, phigh, 0, 0, d);
        tree->insertData(0, nullptr, r, i);
    }   //  build for (1+x)%

    auto t_s = chrono::high_resolution_clock::now();
    for (auto i = 0; i < num_processed; i++){
        plow[0] = pts2[i].first, plow[1] = pts2[i].second;
        phigh[0] = pts2[i].first, phigh[1] = pts2[i].second;
        auto r = TimeRegion(plow, phigh, 0, 0, d);
        tree->deleteData(r, i);
    }
    // for (size_t i = 0; i < pts1.size(); i++){
    //     plow[0] = pts1[i].first, plow[1] = pts1[i].second;
    //     phigh[0] = pts1[i].first, phigh[1] = pts1[i].second;
    //     auto r = TimeRegion(plow, phigh, 0, 0, d);
    //     tree->deleteData(r, n + i);
    // }
    auto t_e = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t_e - t_s);
    cout << *tree << endl;
    cout << "num deleted: " << pts1.size() << endl;
    cout << "delete time: " << 1.0 * duration.count() / 1000 << endl;


    delete tree;
    delete storageManager;
}

/* 3D R-tree Test*/
void R_tree_test(char** argv){
    IStorageManager* storageManager = StorageManager::createNewMemoryStorageManager();
    uint32_t capacity = 20;
    id_type indexIdentifier;
    auto dim = 2;
    ISpatialIndex* tree = RTree::createNewRTree(*storageManager, 0.7, capacity, capacity, dim, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
    
    double x, y;
    double plow[2], phigh[2];
    ifstream fin(argv[1]);
    int n, d;
    
    fin >> n >> d;
    vector<pair<double, double> > pts1, pts2;

    for (auto i = 0; i < n; i++){
        fin >> x >> y;
        pts2.emplace_back(x, y);
    }
    cout << pts2.size() << endl;
    for (size_t i = 0; i < pts2.size(); i++){
        plow[0] = pts2[i].first, plow[1] = pts2[i].second;
        phigh[0] = pts2[i].first, phigh[1] = pts2[i].second;
        auto r = Region(plow, phigh, dim);
        // cout << r << endl;
        tree->insertData(0, nullptr, r, i);
    } 
    
    MyVisitor vis;

    Point qp = Point(plow, 2);
    tree->nearestNeighborQuery(3, qp, vis);

    // plow[0] = 19, plow[1] = 51;
    // phigh[0] = 23, phigh[1] = 54;
    // auto r = Region(plow, phigh, 2);
    // tree->intersectsWithQuery(r, vis);

    // auto num_processed = 10 * n / 100;
    // for (auto i = 0; i < n; i++){
    //     fin >> x >> y;
    //     if (i < num_processed){
    //         pts1.emplace_back(x, y);
    //     }
    //     else{
    //         pts2.emplace_back(x, y);
    //     }
    // }
    // for (size_t i = 0; i < pts2.size(); i++){
    //     plow[0] = pts2[i].first, plow[1] = pts2[i].second, plow[2] = 1;
    //     phigh[0] = pts2[i].first, phigh[1] = pts2[i].second, phigh[2] = 2;
    //     auto r = Region(plow, phigh, dim);
    //     cout << r << endl;
    //     tree->insertData(0, nullptr, r, num_processed + i);
    // }   //  build for (1-x)%
    // cout << *tree << endl;
}

void MV3Rtree_build(char** argv){
    IStorageManager* MVRtree_storageManager = StorageManager::createNewMemoryStorageManager();
    uint32_t capacity = 32;
    id_type indexIdentifier;
    ISpatialIndex* mvrtree = MVRTree::createNewMVRTree(*MVRtree_storageManager, 0.7, capacity, capacity, 2, SpatialIndex::MVRTree::RV_RSTAR, indexIdentifier);

    IStorageManager* Rtree_storageManager = StorageManager::createNewMemoryStorageManager();
    ISpatialIndex* rtree_3d = RTree::createNewRTree(*Rtree_storageManager, 0.7, capacity, capacity, 3, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);


    
    id_type id = 0;
    double plow[2], phigh[2], plow_3d[3], phigh_3d[3];
    ifstream fin(argv[1]);
    int n, d;
    
    fin >> n >> d;
    vector<pair<double, double> > pts(n);

    auto num_processed = 10 * n / 100;
    for (auto i = 0; i < n; i++){
        fin >> pts[i].first >> pts[i].second;
    }

    auto t_ss = chrono::high_resolution_clock::now();
    for (auto i = num_processed; i < n; i++){
        plow[0] = pts[i].first, plow[1] = pts[i].second;
        phigh[0] = pts[i].first, phigh[1] = pts[i].second;
        auto r = TimeRegion(plow, phigh, 0, 0, d);
        vector<pair<TimeRegion, id_type> > removed = {}, added = {};
        mvrtree->insertData_MV3Rtree(0, nullptr, r, id++, removed, added);

        for (auto &pr: removed){
            plow_3d[0] = pr.first.m_pLow[0], plow_3d[1] = pr.first.m_pLow[1], plow_3d[2] = pr.first.m_startTime;
            phigh_3d[0] = pr.first.m_pHigh[0], phigh_3d[1] = pr.first.m_pHigh[1], phigh_3d[2] = 1;
            auto r = Region(plow_3d, phigh_3d, 3);
            // cout << r << " " << pr.second << endl;
            rtree_3d->deleteData(r, pr.second);
            // cout << "delete finished:" << flag << endl;
        }

        for (auto &pr: added){
            plow_3d[0] = pr.first.m_pLow[0], plow_3d[1] = pr.first.m_pLow[1], plow_3d[2] = pr.first.m_startTime;
            phigh_3d[0] = pr.first.m_pHigh[0], phigh_3d[1] = pr.first.m_pHigh[1], phigh_3d[2] = 1;
            auto r = Region(plow_3d, phigh_3d, 3);
            // cout << r << " " << pr.second << endl;
            rtree_3d->insertData(0, nullptr, r, pr.second);
            // cout << "insert finished" << endl;
        }
        // tree->insertData(0, nullptr, r, id++);
    }

    double remove_time = 0.0, add_time = 0.0;
    auto t_s = chrono::high_resolution_clock::now();
    for (auto i = 0; i < num_processed; i++){
        plow[0] = pts[i].first, plow[1] = pts[i].second;
        phigh[0] = pts[i].first, phigh[1] = pts[i].second;
        auto r = TimeRegion(plow, phigh, 0, 0, d);
        vector<pair<TimeRegion, id_type> > removed = {}, added = {};
        mvrtree->insertData_MV3Rtree(0, nullptr, r, id++, removed, added);

        auto removed_s = chrono::high_resolution_clock::now();
        for (auto &pr: removed){
            plow_3d[0] = pr.first.m_pLow[0], plow_3d[1] = pr.first.m_pLow[1], plow_3d[2] = pr.first.m_startTime;
            phigh_3d[0] = pr.first.m_pHigh[0], phigh_3d[1] = pr.first.m_pHigh[1], phigh_3d[2] = 1;
            auto r = Region(plow_3d, phigh_3d, 3);
            // cout << r << " " << pr.second << endl;
            rtree_3d->deleteData(r, pr.second);
            // cout << "delete finished:" << flag << endl;
        }
        auto removed_t = chrono::high_resolution_clock::now();
        auto removed_duration = chrono::duration_cast<chrono::microseconds>(removed_t - removed_s);
        remove_time += 1.0 * removed_duration.count();

        auto added_s = chrono::high_resolution_clock::now();
        for (auto &pr: added){
            plow_3d[0] = pr.first.m_pLow[0], plow_3d[1] = pr.first.m_pLow[1], plow_3d[2] = pr.first.m_startTime;
            phigh_3d[0] = pr.first.m_pHigh[0], phigh_3d[1] = pr.first.m_pHigh[1], phigh_3d[2] = 1;
            auto r = Region(plow_3d, phigh_3d, 3);
            // cout << r << " " << pr.second << endl;
            rtree_3d->insertData(0, nullptr, r, pr.second);
            // cout << "insert finished" << endl;
        }
        auto added_t = chrono::high_resolution_clock::now();
        auto added_duration = chrono::duration_cast<chrono::microseconds>(added_t - added_s);
        add_time += 1.0 * added_duration.count();
        // tree->insertData(0, nullptr, r, id++);
    }
    auto t_e = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t_e - t_ss);
    cout << "build time: " << 1.0 * duration.count() / 1000 << endl;
    duration = chrono::duration_cast<chrono::milliseconds>(t_e - t_s);
    cout << "insert time: " << 1.0 * duration.count() / 1000 << endl;
    cout << "remove time: " << 1.0 * remove_time / 1000000 << endl;
    cout << "add time: " << 1.0 * add_time / 1000000 << endl;

    cout << mvrtree->isIndexValid() << endl;
    cout << rtree_3d->isIndexValid() << endl;

    delete mvrtree;
    delete MVRtree_storageManager;
    delete rtree_3d;
    delete Rtree_storageManager;
}

void run(int argc, char** argv){
	cpam::commandLine cmd(argc, argv, "[-i <Path-to-Input>] [-t <Task-Name>] [-b <Path-to-Batch-File>] [-bf <batch-fraction>] "
                                      "[-r <Path-to-Range-Query>] [-real <Is-Real-Dataset?>] [-a <Algorithm-Name>]"
									  );

    if (!cmd.getOption("-i")){
        cout << "[ERROR] <Path-to-Input> is not specified" << endl;
        return;
    }
    string input_file = cmd.getOptionValue("-i");
    if (!cmd.getOption("-t")){
        cout << "[ERROR] <Task-Name> is not specified" << endl;
        return;
    }
    string algo = cmd.getOptionValue("-a");
    string task = cmd.getOptionValue("-t");
    int is_real = cmd.getOptionIntValue("-real", 0);

    ifstream fin(input_file);
    parlay::sequence<geobase::Point> P;
    // auto largest_mbr = geobase::read_pts(P, fin, is_real);
    geobase::read_pts(P, fin, is_real);

    if (task == "build"){
        auto build_mvrtree = [&](){
            MVRTest::build_test(P);
        };
        auto build_mv3rtree = [&](){
            MV3RTest::build_test(P);
        };
        if (algo == "mvrtree"){
            build_mvrtree();
        }
        else if (algo == "mv3rtree"){
            build_mv3rtree();
        }
        else if (algo == "combined"){
            build_mvrtree();
            malloc_trim(0);  
            build_mv3rtree();
        }
        else{
            cout << "not supported algorithm." << endl;
        }
    }

    if (task == "batch-insert"){
        parlay::sequence<size_t> batch_sizes = {
			10000,
			20000,
			50000,
			100000,
			200000,
			500000,
			1000000,
			2000000,
			5000000,
			10000000,
			20000000,
			50000000,
			100000000
		};
        if (algo == "mvrtree"){
            MVRTest::batch_insert_test(P, batch_sizes);
        }
        if (algo == "mv3rtree"){
            MV3RTest::batch_insert_test(P, batch_sizes);
        }
    }

    if (task == "batch-delete"){
         parlay::sequence<size_t> batch_sizes = {
			10000,
			20000,
			50000,
			100000,
			200000,
			500000,
			1000000,
			2000000,
			5000000,
			10000000,
			20000000,
			50000000,
			100000000
		};

        if (algo == "mvrtree"){
            MVRTest::batch_delete_test(P, batch_sizes);
        }
        if (algo == "mv3rtree"){
            MV3RTest::batch_delete_test(P, batch_sizes);
        }
    }

    IStorageManager* storageManager = nullptr;
    ISpatialIndex* mvrtree = nullptr;
    MVRTest::build_mvrtree(storageManager, mvrtree, P);
    // cout << mvrtree << endl;

    if (task == "range-report"){
        string query_file = cmd.getOptionValue("-r");
        auto [cnt, query] = geobase::read_range_query(query_file, 8, maxSize);
        MVRTest::range_count_test<MyVisitor>(P, query, cnt);
    }
    if (task == "knn"){
        parlay::sequence<size_t> k_vals = {
            1,
            2,
            5,
            10,
            20,
            50,
            100
        };
        for (auto &k: k_vals){
            cout << "[INFO] k = " << k << endl;
            MVRTest::knn_test<NNVisitor>(mvrtree, P, k);
        }    
    }
    if (task == "spatial-diff"){
        string query_file = cmd.getOptionValue("-r");
        auto [cnt, query] = geobase::read_range_query(query_file, 8, maxSize);
        size_t batch_size = 10000;
        parlay::sequence<size_t> batch_sizes = {
            10000,
            20000,
            50000,
            100000,
            200000,
            500000,
            1000000,
            2000000,
            5000000,
            10000000,
            20000000, 
            50000000,
            100000000
        };
        parlay::sequence<size_t> ratios = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        
        size_t ratio = 5;
        cout << "[INFO]: Exp for Batch Size" << endl;
        for (auto &batch_size: batch_sizes){
            if (batch_size > P.size()) break;
            cout << "[INFO] Dealing with Batch Size: " << batch_size << endl;
            MVRTest::spatial_diff_test(P, query, batch_size, ratio);
        }

        batch_size = P.size() / 10 < 10000000 ? P.size() / 10 : 10000000;
        cout << "[INFO]: Exp for Insert/Delete Ratio" << endl;
        for (auto &ratio: ratios){
            cout << "[INFO] Dealing with Ratio: " << ratio << endl;
            MVRTest::spatial_diff_test(P, query, batch_size, ratio);
        }
    }

    if (task == "multi-version-test"){
		if (!cmd.getOption("-mv")){
			cout << "[ERROR]: <Dir-to-Multi-Version> is not specified." << endl;
		}
		else{
            string mv_dir = cmd.getOptionValue("-mv");
            if (algo == "mvrtree"){
                // MVRTest::multi_version_test(P, mv_dir, 18, 7);
                MVRTest::multi_version_test(P, mv_dir, 14, 11);   // Japan
            }
            else if (algo == "mv3rtree"){
                // MV3RTest::multi_version_test(P, mv_dir, 18, 7);
                MV3RTest::multi_version_test(P, mv_dir, 14, 11);   // Japan
            }
            else{
                MVRTest::multi_version_test(P, mv_dir, 18, 7);
                malloc_trim(0);
                MV3RTest::multi_version_test(P, mv_dir, 18, 7);
            }
		}
		return;
	}

    delete mvrtree;
    delete storageManager;
}

int main(int argc, char** argv) {
    run(argc, argv);
    // build_test(argv);
    // insert_test(argv);
    // delete_test(argv);
    // R_tree_test(argv);
    // MV3Rtree_build(argv);

    return 0;
}
