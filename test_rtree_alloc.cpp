#include <iostream>
#include <vector>
#include <chrono>
#include <atomic>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef pair<BoostPoint, size_t> Value;

inline std::atomic<size_t> boost_live_mem(0);

template <typename T>
class TrackingAllocator {
public:
    typedef T value_type;
    TrackingAllocator() = default;
    template <typename U> TrackingAllocator(const TrackingAllocator<U>&) {}
    
    T* allocate(std::size_t n) {
        boost_live_mem.fetch_add(n * sizeof(T), std::memory_order_relaxed);
        return static_cast<T*>(::operator new(n * sizeof(T)));
    }
    void deallocate(T* p, std::size_t n) {
        boost_live_mem.fetch_sub(n * sizeof(T), std::memory_order_relaxed);
        ::operator delete(p);
    }
};

template <typename T, typename U>
bool operator==(const TrackingAllocator<T>&, const TrackingAllocator<U>&) { return true; }
template <typename T, typename U>
bool operator!=(const TrackingAllocator<T>&, const TrackingAllocator<U>&) { return false; }

typedef bgi::rtree<Value, bgi::quadratic<32>, bgi::indexable<Value>, bgi::equal_to<Value>, TrackingAllocator<Value>> RTree;

int main() {
    vector<Value> pts;
    for(size_t i=0; i<6500000; i++) {
        pts.push_back(make_pair(BoostPoint(i*1.0, i*1.0), i));
    }
    auto t1 = chrono::high_resolution_clock::now();
    cout << "Bulk loading..." << endl;
    RTree tree(pts.begin(), pts.end());
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Bulk load took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " ms" << endl;
    
    cout << "Live Mem: " << boost_live_mem.load() / (1024.0*1024.0) << " MB" << endl;
    return 0;
}
