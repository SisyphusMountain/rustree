# Lineages Through Time (LTT) Plots

## Overview

The `rustree` Python API provides built-in support for creating Lineages Through Time (LTT) plots, a common visualization in phylogenetics that shows how the number of lineages changes over time.

## Quick Start

```python
import rustree
import matplotlib.pyplot as plt

# Simulate a species tree
tree = rustree.simulate_species_tree(n=100, lambda_=1.0, mu=0.5, seed=42)

# Option 1: Get data and plot manually
ltt = tree.get_ltt_data()
plt.plot(ltt['times'], ltt['lineages'])
plt.xlabel('Time (past to present)')
plt.ylabel('Number of lineages')
plt.show()

# Option 2: Use built-in plotting
tree.plot_ltt()  # Display interactively
tree.plot_ltt(filepath='ltt.png')  # Save to file
```

## API Reference

### `get_ltt_data()`

Returns a dictionary with LTT data suitable for plotting.

**Returns:**
- `dict` with keys:
  - `'times'`: List of time points (from present at 0 going backward)
  - `'lineages'`: Number of lineages alive at each time point

**Example:**
```python
tree = rustree.simulate_species_tree(n=50, lambda_=1.0, mu=0.3, seed=42)
ltt = tree.get_ltt_data()

print(f"Time points: {len(ltt['times'])}")
print(f"Time range: {min(ltt['times']):.2f} to {max(ltt['times']):.2f}")
print(f"Max lineages: {max(ltt['lineages'])}")
```

### `plot_ltt(filepath=None, title=None, xlabel=None, ylabel=None)`

Creates and displays/saves an LTT plot using matplotlib.

**Parameters:**
- `filepath` (str, optional): Path to save the plot. If `None`, displays interactively.
- `title` (str, optional): Plot title. Default: `"Lineages Through Time"`
- `xlabel` (str, optional): X-axis label. Default: `"Time (past to present)"`
- `ylabel` (str, optional): Y-axis label. Default: `"Number of lineages"`

**Example:**
```python
tree = rustree.simulate_species_tree(n=100, lambda_=1.0, mu=0.7, seed=123)

# Display interactively
tree.plot_ltt()

# Save with custom labels
tree.plot_ltt(
    filepath='my_ltt.png',
    title='Birth-Death Tree (λ=1.0, μ=0.7)',
    xlabel='Time before present',
    ylabel='Lineage count'
)
```

## Advanced Usage

### Log-Scale LTT Plots

For trees with many lineages, log-scale plots are often more informative:

```python
tree = rustree.simulate_species_tree(n=500, lambda_=1.0, mu=0.5, seed=42)
ltt = tree.get_ltt_data()

import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Linear scale
ax1.plot(ltt['times'], ltt['lineages'])
ax1.set_ylabel('Lineages')
ax1.set_title('Linear Scale')

# Log scale
ax2.semilogy(ltt['times'], ltt['lineages'])
ax2.set_ylabel('Lineages (log)')
ax2.set_title('Log Scale')

plt.savefig('ltt_comparison.png')
```

### Comparing Multiple Trees

```python
import matplotlib.pyplot as plt

scenarios = [
    {"n": 50, "mu": 0.0, "label": "No extinction"},
    {"n": 50, "mu": 0.5, "label": "Moderate extinction"},
    {"n": 50, "mu": 0.9, "label": "High extinction"},
]

plt.figure(figsize=(10, 6))

for scenario in scenarios:
    tree = rustree.simulate_species_tree(
        n=scenario['n'],
        lambda_=1.0,
        mu=scenario['mu'],
        seed=42
    )
    ltt = tree.get_ltt_data()
    plt.plot(ltt['times'], ltt['lineages'], label=scenario['label'], linewidth=2)

plt.xlabel('Time (past to present)')
plt.ylabel('Number of lineages')
plt.title('Effect of Extinction Rate on LTT')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('ltt_extinction_comparison.png')
```

### Analyzing LTT Statistics

```python
tree = rustree.simulate_species_tree(n=100, lambda_=1.0, mu=0.5, seed=42)
ltt = tree.get_ltt_data()

# Summary statistics
max_lineages = max(ltt['lineages'])
time_to_mrca = max(ltt['times'])
avg_lineages = sum(ltt['lineages']) / len(ltt['lineages'])

print(f"Maximum lineages: {max_lineages}")
print(f"Time to MRCA: {time_to_mrca:.3f}")
print(f"Average lineages: {avg_lineages:.1f}")

# Find time when lineages peaked
max_idx = ltt['lineages'].index(max_lineages)
peak_time = ltt['times'][max_idx]
print(f"Lineages peaked at time {peak_time:.3f}")
```

### Smoothing LTT Plots

For noisy trees, you can apply smoothing:

```python
import numpy as np
from scipy.ndimage import gaussian_filter1d

tree = rustree.simulate_species_tree(n=200, lambda_=1.0, mu=0.6, seed=42)
ltt = tree.get_ltt_data()

# Apply Gaussian smoothing
smoothed_lineages = gaussian_filter1d(ltt['lineages'], sigma=2)

plt.figure(figsize=(10, 6))
plt.plot(ltt['times'], ltt['lineages'], alpha=0.3, label='Raw')
plt.plot(ltt['times'], smoothed_lineages, linewidth=2, label='Smoothed')
plt.xlabel('Time (past to present)')
plt.ylabel('Number of lineages')
plt.legend()
plt.savefig('ltt_smoothed.png')
```

## Understanding LTT Plots

### Interpretation

- **Slope**: Steep increases indicate rapid speciation
- **Plateaus**: Periods of balanced speciation/extinction
- **Decreases**: Occur when going backward through extinction events
- **Pure birth trees**: Always decrease monotonically going backward
- **Trees with extinction**: Can show complex patterns

### Key Features

1. **Present (time=0)**: Starts with *n* extant species
2. **Past**: Number of lineages changes at each birth/death event
3. **MRCA**: Ends at 1 lineage (most recent common ancestor)

### Example Patterns

```python
import rustree
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. Pure birth - smooth exponential
tree = rustree.simulate_species_tree(n=50, lambda_=1.0, mu=0.0, seed=1)
ltt = tree.get_ltt_data()
axes[0, 0].plot(ltt['times'], ltt['lineages'])
axes[0, 0].set_title('Pure Birth: Smooth exponential decay')

# 2. Low extinction - gentle fluctuations
tree = rustree.simulate_species_tree(n=50, lambda_=1.0, mu=0.2, seed=2)
ltt = tree.get_ltt_data()
axes[0, 1].plot(ltt['times'], ltt['lineages'])
axes[0, 1].set_title('Low Extinction: Gentle fluctuations')

# 3. High extinction - large fluctuations
tree = rustree.simulate_species_tree(n=50, lambda_=1.0, mu=0.8, seed=3)
ltt = tree.get_ltt_data()
axes[1, 0].plot(ltt['times'], ltt['lineages'])
axes[1, 0].set_title('High Extinction: Large fluctuations')

# 4. Log scale view
tree = rustree.simulate_species_tree(n=100, lambda_=1.0, mu=0.5, seed=4)
ltt = tree.get_ltt_data()
axes[1, 1].semilogy(ltt['times'], ltt['lineages'])
axes[1, 1].set_title('Log Scale: Reveals early diversification')

for ax in axes.flat:
    ax.set_xlabel('Time (past to present)')
    ax.set_ylabel('Lineages')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('ltt_patterns.png')
```

## Working with Parsed Trees

LTT plots work with both simulated and parsed Newick trees:

```python
# Parse a Newick tree
tree = rustree.parse_species_tree("((A:1,B:1):1,(C:0.5,D:0.5):1.5):0;")

# Generate LTT plot
# For parsed trees, events are inferred from topology
ltt = tree.get_ltt_data()
tree.plot_ltt(title='Parsed Newick Tree')
```

## Notes

- Time flows backward from present (0) to the root (MRCA)
- For trees with extinction, more total nodes exist than extant species
- The algorithm counts lineages by processing events chronologically
- Works with any tree size (tested up to millions of species)
