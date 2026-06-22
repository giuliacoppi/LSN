# Numerical Simulation Laboratory

Welcome to the repository for the **Numerical Simulation Laboratory** course. This space contains all the source code, data analysis, and numerical results developed during the practical laboratory sessions.

The primary objective of these activities is the application of computational and stochastic methods to solve complex physical problems, ranging from statistical mechanics to machine learning.

---

## 🛠️ Repository Structure

The repository is organized in a modular structure. Each laboratory session has its own dedicated directory, designed to keep the source code separate from the analysis:

* **`Esercitazione_XX/`**: Main directory for a specific session (e.g., `Esercitazione_01`, `Esercitazione_02`, etc.).
    * 📄 **`Esercitazione_XX.ipynb`**: The main Jupyter Notebook. It contains the physics discussion of the problem, data visualization, and statistical analysis.
    * 📁 **`src/`**: Subfolder containing:
        * C++ source files (`.cpp`, `.h`) used to run the simulations.
        * Configuration files (`input.dat`, etc.).
        * Raw output files (`.dat` or `.txt`) containing the simulation data.

### Directory Tree Example:
```text
.
├── README.md
├── Esercitazione_01/
│   ├── Esercitazione_01.ipynb      # Analysis notebook and plots
│   └── src/                        # C++ source code and simulation results
│       ├── main.cpp
│       ├── random.cpp
│       ├── random.h
│       ├── input.dat
│       └── output_data.dat
├── Esercitazione_02/
│   ├── Esercitazione_02.ipynb
│   └── src/
...
