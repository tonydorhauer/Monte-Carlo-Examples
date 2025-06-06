import matplotlib.pyplot as plt
import csv

def read_samples(filename):
    samples = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row:  # Make sure row is not empty
                samples.append(float(row[0]))
    return samples

def plot_histogram(filename, bins=100):
    samples = read_samples(filename)
    plt.hist(samples, bins=bins, edgecolor='black')
    plt.xlabel('Sample Value')
    plt.ylabel('Frequency')
    plt.grid(False)
    plt.show()
    
def read_2d_samples(filename):
    """Reads 2D samples (x, y pairs) from a CSV."""
    samples = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) == 2:
                x, y = float(row[0]), float(row[1])
                samples.append((x, y))
    return samples

def plot_2d_histogram(filename, bins=100):
    """Plots a 2D histogram from a CSV file of x, y samples."""
    samples = read_2d_samples(filename)
    x_vals = [s[0] for s in samples]
    y_vals = [s[1] for s in samples]

    plt.figure(figsize=(6, 5))
    plt.hist2d(x_vals, y_vals, bins=bins, density=True, cmap='RdYlBu_r')
    plt.colorbar(label='Density')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('2D Metropolis Samples')
    plt.tight_layout()
    plt.show()
    

if __name__ == "__main__":
    #plot_histogram("samples.csv")
    plot_2d_histogram("samples2D.csv")