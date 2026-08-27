import os

bps = ["100", "50", "25", "12.5", "6.25", "3.125"]
queries = ["Range_Small", "Range_Med", "Range_Large", "KNN_1", "KNN_10", "KNN_100"]

content = """# Evaluation Report: Multi-Version Spatial Indexing

This report consolidates the experimental results evaluating the performance of pure persistent trees (MVZD, CPAMBB) against log-structured approaches (RlogTree variants) over a simulated 5-year data evolution period (2021-2025). 
**Note:** "Update Time" has been formally renamed to **Branching & Merging Time** throughout this report, as the core operation evaluates the cost of creating a branch and merging it into the master snapshot.

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

"""

for bp in bps:
    content += f"### Memory Lifecycle (BP = {bp}%)\n"
    content += f"![Memory BP {bp}](plots_and_tables/bp_{bp}/Global_bp{bp}_Fig1_Memory.png)\n\n"

content += """---

## 3. Branching & Merging Time Evolution

Figure 2 visualizes the write-side overhead, contrasting the periodic latency spikes of log-structured merging with the continuous structural maintenance of persistent trees.
This section is broken down into Incremental Insertion, Incremental Deletion, and Cumulative Cost.

**Amortized Maintenance Cost Explanation:** "Amortized" refers to the total accumulated time divided by the number of operations. While `RlogTree` experiences massive latency spikes (e.g., hundreds of milliseconds during compaction), its everyday log appending is nearly 0ms. When accumulated over 5 years (the Cumulative plot), its *total* cost is lower than persistent trees without bulk-merge.
**Bulk-Merge Optimization:** By implementing the multi-branch bulk-merge operator, `MVZD` aggregates intermediate deltas and eliminates redundant path-copying traversals, keeping its update throughput highly competitive.

### 3.1 Incremental Insertion Cost
"""

for bp in bps:
    content += f"#### BP = {bp}%\n"
    content += f"![Branching Forward BP {bp}](plots_and_tables/bp_{bp}/Global_bp{bp}_Fig2_Branching_Forward.png)\n\n"

content += """### 3.2 Incremental Deletion Cost
"""

for bp in bps:
    content += f"#### BP = {bp}%\n"
    content += f"![Branching Backward BP {bp}](plots_and_tables/bp_{bp}/Global_bp{bp}_Fig2_Branching_Backward.png)\n\n"

content += """### 3.3 Cumulative Branching & Merging Time
"""

for bp in bps:
    content += f"#### BP = {bp}%\n"
    content += f"![Branching Cumulative BP {bp}](plots_and_tables/bp_{bp}/Global_bp{bp}_Fig2_Branching_Cumulative.png)\n\n"

content += """---

## 4. Query Scalability and Read Penalty

Query performances are broken down into individual charts. The X-axis aligns with the batch-completion snapshots from the memory timeline.

**Key Observations:**
* **The Read Penalty (Sawtooth Effect):** `RlogTree` variants experience exponential performance degradation as queries scale up. The accumulation of redundant logs forces expensive filtering traversals, causing query latency to climb steadily until a compaction resets the tree.
* **Strict O(1) Guarantees:** Because `MVZD` and `CPAMBB` pay the structural update cost upfront, they maintain pristine tree structures at all versions. This guarantees optimal, flat-line query latencies.

"""

for qt in queries:
    content += f"### Query Latency: {qt}\n"
    for bp in bps:
        content += f"#### BP = {bp}%\n"
        content += f"![Query {qt} BP {bp}](plots_and_tables/bp_{bp}/Global_bp{bp}_Fig3_Query_{qt}.png)\n\n"

with open("Evaluation_Report.md", "w", encoding="utf-8") as f:
    f.write(content)

