#include <iostream>
#include <vector>
#include <queue>

using namespace std;

struct Point {
    size_t id;
    double x, y;
    Point(size_t i, double x, double y) : id(i), x(x), y(y) {}
};

typedef pair<Point, double> nn_pair;

struct nn_pair_cmp {
    bool operator()(const nn_pair &lhs, const nn_pair &rhs) const {
        return lhs.second < rhs.second ||
               (lhs.second == rhs.second && lhs.first.id > rhs.first.id);
    }
};

int main() {
    nn_pair p1 = {Point(1, 0, 0), 10.0};
    nn_pair p2 = {Point(2, 0, 0), 10.0};
    nn_pair p3 = {Point(3, 0, 0), 5.0};
    
    priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> pq;
    pq.push(p1);
    pq.push(p2);
    pq.push(p3);
    
    while(!pq.empty()) {
        cout << "ID: " << pq.top().first.id << " Dist: " << pq.top().second << endl;
        pq.pop();
    }
    
    return 0;
}
