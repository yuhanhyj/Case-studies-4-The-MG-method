import matplotlib.pyplot as plt
import glob
import os

def load_residuals(filename):
    iterations = []
    residuals = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                iterations.append(int(parts[0]))
                residuals.append(float(parts[1]))
    return iterations, residuals

def main():
    # Find all residual_N*.txt files
    files = sorted(glob.glob("residual_N*.txt"))

    if not files:
        print("No residual_N*.txt files found.")
        return

    plt.figure(figsize=(10, 6))
    for file in files:
        iterations, residuals = load_residuals(file)
        label = os.path.splitext(os.path.basename(file))[0]  # e.g., residual_N64_lmax2
        plt.semilogy(iterations, residuals, label=label)

    plt.xlabel("Iteration")
    plt.ylabel("Residual (log scale)")
    plt.title("Multigrid Convergence: Residual vs Iteration")
    plt.legend()
    plt.grid(True, which="both", ls="--", lw=0.5)
    plt.tight_layout()
    plt.savefig("residual_convergence.png")
    plt.show()

if __name__ == "__main__":
    main()
