# Evaluation Report: Multi-Version Spatial Indexing

This report consolidates the experimental results evaluating the performance of pure persistent trees (MVZD, CPAMBB) against log-structured approaches (RlogTree variants) over a simulated 5-year data evolution period (2021-2025). 
**Note:** "Update Time" has been formally renamed to **Branching & Merging Time** throughout this report, as the core operation evaluates the cost of creating a branch and merging it into the master snapshot.

---

## 0. Dataset Overview

To properly contextualize the performance trends (especially the severe spikes in 2021 and the flattening in 2025), it is crucial to understand the absolute volume of operations injected into the database each year. The table below outlines the exact number of incremental update points (used symmetrically for both Forward Insertions and Backward Deletions) derived from the Bhutan Evolution Dataset.

![Dataset Operations](Paper_Dataset_Stats.png)

---

## 1. Comprehensive Performance Matrix

Table 1 presents the macro-level performance across all tested granularities, from BP=100% (1 batch/year) to BP=3.125% (32 batches/year). 

**Formatting Notes:**
* All execution times are reported in milliseconds (ms).
* The fastest time in each column per year is bolded and underlined, with relative performance degradation color-coded via a 5-level logarithmic scale (< 3x, < 10x, < 35x, < 100x, > 100x).
* Algorithms triggering full tree compactions during a specific phase are denoted with an asterisk (`*`) and their corresponding update latencies are highlighted with a bounding box.

![Table 1: All Years Combined](Paper_Table_All_Years_Combined.png)

---

## 2. Memory Footprint and Spatial Evolution

Figure 1 illustrates the Resident Set Size (RSS) evolution during data ingestion. The timeline is dynamically scaled to align intra-year batches while compressing inter-year gaps.

**Key Observations:**
* **RSS Bloat in Log-Structured Trees:** `Rlog_1yr` exhibits the highest peak memory (~1187 MB). This is driven by heap fragmentation and unreleased `glibc` free-lists resulting from the node allocations/deallocations during its periodic `compact()` operations.
* **Contiguous Allocation Efficiency:** Despite the internal capacity fragmentation inherent to fat-leaves, `MVZD` maintains a stable memory footprint (~628 MB). Its reliance on large, contiguous array allocations effectively eliminates external fragmentation.

### Memory Lifecycle (BP = 100%)
![Memory BP 100](plots_and_tables/bp_100/Global_bp100_Fig1_Memory.png)

### Memory Lifecycle (BP = 50%)
![Memory BP 50](plots_and_tables/bp_50/Global_bp50_Fig1_Memory.png)

### Memory Lifecycle (BP = 25%)
![Memory BP 25](plots_and_tables/bp_25/Global_bp25_Fig1_Memory.png)

### Memory Lifecycle (BP = 12.5%)
![Memory BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig1_Memory.png)

### Memory Lifecycle (BP = 6.25%)
![Memory BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig1_Memory.png)

### Memory Lifecycle (BP = 3.125%)
![Memory BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig1_Memory.png)

---

## 3. Branching & Merging Time Evolution

Figure 2 visualizes the write-side overhead, contrasting the periodic latency spikes of log-structured merging with the continuous structural maintenance of persistent trees.
This section is broken down into Incremental Insertion, Incremental Deletion, and Cumulative Cost.

**Amortized Maintenance Cost Explanation:** "Amortized" refers to the total accumulated time divided by the number of operations. While `RlogTree` experiences massive latency spikes (e.g., hundreds of milliseconds during compaction), its everyday log appending is nearly 0ms. When accumulated over 5 years (the Cumulative plot), its *total* cost is lower than persistent trees without bulk-merge.
**Bulk-Merge Optimization:** By implementing the multi-branch bulk-merge operator, `MVZD` aggregates intermediate deltas and eliminates redundant path-copying traversals, keeping its update throughput highly competitive.

### 3.1 Incremental Insertion Cost
#### BP = 100%
![Branching Forward BP 100](plots_and_tables/bp_100/Global_bp100_Fig2_Branching_Forward.png)

#### BP = 50%
![Branching Forward BP 50](plots_and_tables/bp_50/Global_bp50_Fig2_Branching_Forward.png)

#### BP = 25%
![Branching Forward BP 25](plots_and_tables/bp_25/Global_bp25_Fig2_Branching_Forward.png)

#### BP = 12.5%
![Branching Forward BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig2_Branching_Forward.png)

#### BP = 6.25%
![Branching Forward BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig2_Branching_Forward.png)

#### BP = 3.125%
![Branching Forward BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig2_Branching_Forward.png)

### 3.2 Incremental Deletion Cost
#### BP = 100%
![Branching Backward BP 100](plots_and_tables/bp_100/Global_bp100_Fig2_Branching_Backward.png)

#### BP = 50%
![Branching Backward BP 50](plots_and_tables/bp_50/Global_bp50_Fig2_Branching_Backward.png)

#### BP = 25%
![Branching Backward BP 25](plots_and_tables/bp_25/Global_bp25_Fig2_Branching_Backward.png)

#### BP = 12.5%
![Branching Backward BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig2_Branching_Backward.png)

#### BP = 6.25%
![Branching Backward BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig2_Branching_Backward.png)

#### BP = 3.125%
![Branching Backward BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig2_Branching_Backward.png)

### 3.3 Cumulative Branching & Merging Time
#### BP = 100%
![Branching Cumulative BP 100](plots_and_tables/bp_100/Global_bp100_Fig2_Branching_Cumulative.png)

#### BP = 50%
![Branching Cumulative BP 50](plots_and_tables/bp_50/Global_bp50_Fig2_Branching_Cumulative.png)

#### BP = 25%
![Branching Cumulative BP 25](plots_and_tables/bp_25/Global_bp25_Fig2_Branching_Cumulative.png)

#### BP = 12.5%
![Branching Cumulative BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig2_Branching_Cumulative.png)

#### BP = 6.25%
![Branching Cumulative BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig2_Branching_Cumulative.png)

#### BP = 3.125%
![Branching Cumulative BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig2_Branching_Cumulative.png)

---

## 4. Query Scalability and Read Penalty

Query performances are broken down into individual charts. The X-axis aligns with the batch-completion snapshots from the memory timeline.

**Key Observations:**
* **The Read Penalty (Sawtooth Effect):** `RlogTree` variants experience exponential performance degradation as queries scale up. The accumulation of redundant logs forces expensive filtering traversals, causing query latency to climb steadily until a compaction resets the tree.
* **Strict O(1) Guarantees:** Because `MVZD` and `CPAMBB` pay the structural update cost upfront, they maintain pristine tree structures at all versions. This guarantees optimal, flat-line query latencies.

### Query Latency: Range_Small
#### BP = 100%
![Query Range_Small BP 100](plots_and_tables/bp_100/Global_bp100_Fig3_Query_Range_Small.png)

#### BP = 50%
![Query Range_Small BP 50](plots_and_tables/bp_50/Global_bp50_Fig3_Query_Range_Small.png)

#### BP = 25%
![Query Range_Small BP 25](plots_and_tables/bp_25/Global_bp25_Fig3_Query_Range_Small.png)

#### BP = 12.5%
![Query Range_Small BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig3_Query_Range_Small.png)

#### BP = 6.25%
![Query Range_Small BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig3_Query_Range_Small.png)

#### BP = 3.125%
![Query Range_Small BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig3_Query_Range_Small.png)

### Query Latency: Range_Med
#### BP = 100%
![Query Range_Med BP 100](plots_and_tables/bp_100/Global_bp100_Fig3_Query_Range_Med.png)

#### BP = 50%
![Query Range_Med BP 50](plots_and_tables/bp_50/Global_bp50_Fig3_Query_Range_Med.png)

#### BP = 25%
![Query Range_Med BP 25](plots_and_tables/bp_25/Global_bp25_Fig3_Query_Range_Med.png)

#### BP = 12.5%
![Query Range_Med BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig3_Query_Range_Med.png)

#### BP = 6.25%
![Query Range_Med BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig3_Query_Range_Med.png)

#### BP = 3.125%
![Query Range_Med BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig3_Query_Range_Med.png)

### Query Latency: Range_Large
#### BP = 100%
![Query Range_Large BP 100](plots_and_tables/bp_100/Global_bp100_Fig3_Query_Range_Large.png)

#### BP = 50%
![Query Range_Large BP 50](plots_and_tables/bp_50/Global_bp50_Fig3_Query_Range_Large.png)

#### BP = 25%
![Query Range_Large BP 25](plots_and_tables/bp_25/Global_bp25_Fig3_Query_Range_Large.png)

#### BP = 12.5%
![Query Range_Large BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig3_Query_Range_Large.png)

#### BP = 6.25%
![Query Range_Large BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig3_Query_Range_Large.png)

#### BP = 3.125%
![Query Range_Large BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig3_Query_Range_Large.png)

### Query Latency: KNN_1
#### BP = 100%
![Query KNN_1 BP 100](plots_and_tables/bp_100/Global_bp100_Fig3_Query_KNN_1.png)

#### BP = 50%
![Query KNN_1 BP 50](plots_and_tables/bp_50/Global_bp50_Fig3_Query_KNN_1.png)

#### BP = 25%
![Query KNN_1 BP 25](plots_and_tables/bp_25/Global_bp25_Fig3_Query_KNN_1.png)

#### BP = 12.5%
![Query KNN_1 BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig3_Query_KNN_1.png)

#### BP = 6.25%
![Query KNN_1 BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig3_Query_KNN_1.png)

#### BP = 3.125%
![Query KNN_1 BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig3_Query_KNN_1.png)

### Query Latency: KNN_10
#### BP = 100%
![Query KNN_10 BP 100](plots_and_tables/bp_100/Global_bp100_Fig3_Query_KNN_10.png)

#### BP = 50%
![Query KNN_10 BP 50](plots_and_tables/bp_50/Global_bp50_Fig3_Query_KNN_10.png)

#### BP = 25%
![Query KNN_10 BP 25](plots_and_tables/bp_25/Global_bp25_Fig3_Query_KNN_10.png)

#### BP = 12.5%
![Query KNN_10 BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig3_Query_KNN_10.png)

#### BP = 6.25%
![Query KNN_10 BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig3_Query_KNN_10.png)

#### BP = 3.125%
![Query KNN_10 BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig3_Query_KNN_10.png)

### Query Latency: KNN_100
#### BP = 100%
![Query KNN_100 BP 100](plots_and_tables/bp_100/Global_bp100_Fig3_Query_KNN_100.png)

#### BP = 50%
![Query KNN_100 BP 50](plots_and_tables/bp_50/Global_bp50_Fig3_Query_KNN_100.png)

#### BP = 25%
![Query KNN_100 BP 25](plots_and_tables/bp_25/Global_bp25_Fig3_Query_KNN_100.png)

#### BP = 12.5%
![Query KNN_100 BP 12.5](plots_and_tables/bp_12.5/Global_bp12.5_Fig3_Query_KNN_100.png)

#### BP = 6.25%
![Query KNN_100 BP 6.25](plots_and_tables/bp_6.25/Global_bp6.25_Fig3_Query_KNN_100.png)

#### BP = 3.125%
![Query KNN_100 BP 3.125](plots_and_tables/bp_3.125/Global_bp3.125_Fig3_Query_KNN_100.png)

