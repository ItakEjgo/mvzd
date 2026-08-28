# Multiversion Spatial Index Benchmark Analysis

## 1. Overview of the Current Benchmark
The current benchmark (`run_experiments.py` executing `verify_bench.cpp`) performs an **event-driven simulation** using real-world OpenStreetMap (OSM) data from Bhutan (2018–2026). 
- **Initialization:** It loads a base snapshot from 2017 to initialize the indices (`MVZD`, `CPAMBB`, `RlogTree`).
- **Event Loop:** It reads all historical OSM commits/changesets and queues them as `START` and `END` events. It applies updates (insertions and deletions) sequentially to simulate real-world database progression.
- **Checkpointing & Querying:** Every `q_step` commits (default 1000), the benchmark pauses to run a suite of spatial queries against the current state of the trees.
- **Query Types:** It generates dynamically sized Range Queries (`Range_Small` [0.05%], `Range_Med` [0.5%], `Range_Large` [5%]) and K-Nearest Neighbor Queries (`KNN_1`, `KNN_10`, `KNN_100`).
- **Metrics Collected:** Query execution time, structural hash (for correctness verification), query result counts, and memory consumption.

## 2. Source Code Implementation Issues & Bugs

### A. Critical Bug: Unfair Memory Tracking for RlogTree
The memory profiling methodology compares "apples to oranges," severely penalizing the `RlogTree` baseline.
- **MVZD & CPAMBB:** Track memory precisely by overriding allocators or traversing tree structures to count exactly the bytes used by index nodes (e.g., `mvq::global_live_mem` and `size_in_bytes_all`).
- **RlogTree:** Uses `mem_boost()`, which calls `get_rss_bytes()` to read the **total physical memory (RSS) of the entire process**. 
- **The Issue:** The `initial_rss` variable is declared but never initialized to the baseline memory. The process loads the entire Bhutan OSM commit history (years 2018-2026) into `std::vector<Event>` and massive hash maps (`active_nodes`, `changeset_ops`). This test-runner overhead (hundreds of MBs/GBs) is falsely attributed entirely to `RlogTree`'s memory footprint, making `RlogTree` look terribly inefficient compared to `MVZD`/`CPAMBB`.
- **Missed Intent:** The developer wrote a `TrackingAllocator` to track Boost's memory, but **forgot to pass it** as the 5th template argument to `bgi::rtree`. 

### B. Critical Logic Flaw: RlogTree Compaction Thresholds
The `RlogTree` variants (`Rlog_1yr`, `Rlog_2yr`, etc.) are intended to simulate rebuilding the RTree every 1, 2, or N years.
- **The Issue:** The thresholds are hardcoded to synthetic values (e.g., `20` batches for `1yr`). This was likely copied from an older synthetic workload where 100 total batches represented 5 years.
- **Impact:** In the real-world Bhutan dataset, 2018 alone has ~1,868 changesets. By using a threshold of `20`, `Rlog_1yr` actually rebuilds the entire RTree every **few days** (20 commits) instead of every year. 
- **Distorted Results:** This will cause `RlogTree` to have astronomically high update latencies (rebuilding 90+ times a year instead of 1 time) and artificially fast query latencies (because the delta log is always nearly empty, avoiding the need to scan a large log).

### C. KNN Search Implementation in RlogTree (Observation)
The implementation of `RlogTree::knn_report` is mathematically correct but complex. It correctly fetches the top-$k$ elements from the base RTree and merges them with the `insert_log` using a max-heap. 
- It uses a custom comparator `MaxHeapCmp`.
- **Optimization 5 (Early Termination):** It correctly breaks the RTree traversal once the distance of the elements from the RTree exceeds the largest distance in the full max-heap. This is sound because `bgi::nearest` guarantees non-decreasing distance order.

## 3. Information Gathered & Unreasonable Aspects

1. **Unreasonable Semantic Mismatch:** The naming convention `Rlog_1yr` directly contradicts the hardcoded code logic (`threshold = 20` commits). The benchmark does not measure 1 year of data accumulation for the log-structured baseline.
2. **CPAMBB Global State Dependency:** The memory tracking for `CPAMBB` relies on a global variable `cpambb_global_mem` and an `unordered_map` (`cpambb_mmp`) to deduplicate shared nodes. While mathematically correct for tracking the *delta* of newly allocated nodes across versions in a single sequential run, it relies heavily on the state not being corrupted and prohibits concurrent evaluation.
3. **Data Pre-loading Overhead:** The `verify_bench.cpp` program pre-loads all historical data into memory before running. For larger datasets (like the entire world OSM), this architecture will immediately OOM. The data should ideally be streamed from disk as the simulation progresses.
4. **Boost RTree Deletions (Optimization 1/2):** The author correctly avoided using Boost's native `remove()` (which triggers expensive tree rebalancing) and instead uses a `removed_ids` hash set (lazy deletion). However, because of the frequent compaction bug mentioned above, this optimization's true benefit is never fully realized.
