# SILVA: Spatial Index Library with Version Access

## Overview
**SILVA** is a high-performance library for managing multi-version spatial data (points). It enables storing multiple versions of spatial indexes and supports real-time spatial queries on any stored version, including range counting/reporting and k-nearest neighbors (KNN). 
In addition to *standard spatial queries*, it also supports
*version control* operations (commit, merge, spatial diff, purge) are optimized for parallel execution. 
Our experimental results show that SILVA can achieve  **100x faster operations** and **70% lower memory usage** than state-of-the-art multi-version spatial index MV(3)R-tree from *libspatialindex* library.

## Key Features
- **Multi-Version Control**  
   versioned operations:  
   - `commit`
   - `merge`
   - `spatial diff`, and 
   - `purge`
- **Standrad Single Version Spatial Queries**  
  Execute on any historical version:
  - `range_count` / `range_report`
  - `knn_search`
- **Parallel Updates**
  - Update to the multi-version index can be highly in parallel.  
- **High Efficiency**  
  - Single-version queries match traditional spatial indexes (Rtrees in boost) with no extra overhead.

---

## Implemented Indexes
### 1. PaC-Z-Code & PaC-Z-BB
- **High-level Idea**  
  Maps spatial points to 1D via Z-curve, then extends the **PaC Tree** (1D functional binary search tree) for versioned operations and spatial queries.

### 2. MVZD Tree
- **High-level Idea**  
  Leverages **spatial invariance** to maximize sharing of unchanged regions across versions.
  Achieves minimal memory footprint as well as fast versioned operations.
---

## Performance Comparison
- Please refer to Experiment section in our paper.
---

## Code Structure
.

├── baselines/ # Comparison implementations

│ ├── libspatialindex/ # MVRtree & MV3Rtree

│ └── boostRtree/ # Rtree

├── src/ # SILVA core

├── include/ # PaC-Tree libraries

├── parlaylib/ # parallel primitives


## Compilation
bash

cd src/ # or baselines/

make # Generates binaries

---

## Command-Line Parameters
| Parameter              | Description                          | Example                   |
|------------------------|--------------------------------------|---------------------------|
| `-i <Path-to-Input>`   | Input data file path                 | `-i data/points.in`      |
| `-o <Path-to-Output>`  | Output results path                  | `-o results.log`      |
| `-t <Task-Name>`       | Operation type (`build`, `insert`,`range-count` etc.)| `-t build`          |
| `-a <Algorithm-Name>` | Index type (`mvzd`, `paczz`, `paczbb`, etc.) | `-a MVZD`                |
| `-b <Batch-file>`      | Batch operations file                | `-b batch/ops.in`        |
| `-bf <batch-fraction>` | Percent of batch to process | `-bf 10`               |
| `-real <0/1>`   | Real-world dataset flag              | `-real 1`              |

### Example Command
bash

./main -i data/uniform_100M.in -o build.log -t build -a mvzd -real 0
