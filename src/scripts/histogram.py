import matplotlib.pyplot as plt
import numpy as np

def read_numbers(filename):
    """Read integers (one per line) from a text file."""
    with open(filename, "r") as f:
        numbers = [int(line.strip()) for line in f if line.strip().isdigit()]
    return numbers

def plot_histogram(data_list, labels, bins=50, overlay=False):
    """Plot one or multiple histograms."""
    plt.figure(figsize=(10, 6))

    if overlay:
        # Overlay multiple datasets
        for data, label in zip(data_list, labels):
            plt.hist(data, bins=bins, alpha=0.5, label=label, edgecolor='black')
        plt.legend()
    else:
        # Plot each dataset separately
        for i, (data, label) in enumerate(zip(data_list, labels)):
            plt.figure(figsize=(10, 6))
            plt.hist(data, bins=bins, edgecolor='black')
            plt.title(f"Histogram of {label}")
            plt.xlabel("Value")
            plt.ylabel("Frequency")
            plt.grid(True, alpha=0.4)
            plt.tight_layout()
            plt.show()
        return

    plt.title("Overlayed Histogram")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()

def describe_data(data, name="Dataset"):
    """Print basic statistics."""
    data_np = np.array(data)
    print(f"\n{name} Statistics:")
    print(f"  Count: {len(data_np)}")
    print(f"  Min: {data_np.min():,}")
    print(f"  Max: {data_np.max():,}")
    print(f"  Mean: {data_np.mean():,.2f}")
    print(f"  Median: {np.median(data_np):,.2f}")
    print(f"  Std Dev: {data_np.std():,.2f}")

if __name__ == "__main__":
    # Example: replace with your file paths
    file1 = "dataset1.txt"
    file2 = "dataset2.txt"

    data1 = read_numbers(file1)
    data2 = read_numbers(file2)

    # Print summary stats
    describe_data(data1, "Dataset 1")
    describe_data(data2, "Dataset 2")

    # Plot histograms separately
    plot_histogram([data1, data2], ["Dataset 1", "Dataset 2"], bins=50)

    # Plot overlay comparison
    plot_histogram([data1, data2], ["Dataset 1", "Dataset 2"], bins=50, overlay=True)
