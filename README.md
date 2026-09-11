# QuRAMSearch: Quantum-Assisted Genomic Pattern Search and Scar Identification

**QuRAMSearch** is a hybrid bioinformatics and quantum computing tool engineered to scan genomic databases and locate target sequence motifs (scars) with high fidelity[cite: 2, 3]. By combining classical sliding-window alignment with quantum state encoding and Grover amplitude amplification, the software identifies optimal genomic anchoring sites while mitigating false-positive alignments and memory bottlenecks.

---

## 🧬 Key Features

* **Multimodal Data Ingestion:** Flexible ingestion pipeline supporting:
  * FASTA/TXT sequence loading via standard file paths
  * Manual nucleotide sequence entry through terminal prompts
  * Integrated JCVI-Syn3B benchmark genome for rapid local testing
* **Orthogonal Quantum Encoding:** Maps individual nitrogenous bases into discrete, distinguishable Hilbert space vectors:
  * **A:** $|1\rangle$ via Pauli-X gate
  * **C:** $|+\rangle$ via Hadamard gate
  * **G:** $|-\rangle$ via Pauli-X followed by Hadamard gate
  * **T:** $|0\rangle$ ground state
* **Biological Oracle with Uncomputation:** Transforms target scar sequences into the $|1\rangle$ activation basis and applies multi-controlled phase inversion to mark matching states, followed by exact uncomputation to preserve quantum register coherence.
* **Hybrid Search Strategy:** Couples classical Hamming distance computation and sliding-window indexing with quantum state mapping to ensure bounds safety and identify candidate loci.
* **Grover Amplification Modeling:** Computes optimal iteration steps $\kappa^{opt} \approx \frac{\pi}{4}\sqrt{N}$ using angular displacement $\theta = \arcsin(1/\sqrt{N})$ to simulate theoretical quantum acceleration curves.
* **Diagnostic Visualization Suite:** Automatically exports high-resolution graphical metrics and structured reports:
  * `graph1_grover.png`: Grover theoretical probability curve
  * `graph2_histogram.png`: Oracle state distribution histogram
  * `graph3_scars.png`: Target scar alignment accuracy chart
  * `analysis_result.txt`: Complete alignment summary with mutation markers and genomic context

---

## 📦 Requirements and Dependencies

* **Python 3.8+**
* [Qiskit](https://qiskit.org/)[cite: 4, 10]
* [Qiskit Aer](https://github.com/Qiskit/qiskit-aer)[cite: 2, 3]
* [NumPy](https://numpy.org/)[cite: 2, 3]
* [Matplotlib](https://matplotlib.org/)[cite: 2, 3]

---

## 💻 Installation

1. Clone this repository:
```bash
git clone [https://github.com/YOUR_USERNAME/QuRAMSearch.git](https://github.com/YOUR_USERNAME/QuRAMSearch.git)
cd QuRAMSearch
