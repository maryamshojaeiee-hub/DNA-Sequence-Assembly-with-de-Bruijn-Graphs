# DNA-Sequence-Assembly-with-de-Bruijn-Graphs

Reconstructing a DNA sequence from short, overlapping fragments, a core problem in genome sequencing. The program cleans noisy fragment data, builds a de Bruijn graph from k-mers, and reconstructs the original sequence by finding an Eulerian path through the graph.

Project for the course "Advanced Programming in Python", Master of Statistics and Data Science, Hasselt University, 2024–2025.

Read the full report → · View the code → · View the tests →

## About the Project

Sequencing machines cannot read a long DNA molecule in one go. Instead, the DNA is cut into short fragments that are read separately and then pieced back together using their overlaps. This project implements that assembly step:

Data cleaning: fragments with duplicated, missing or invalid positions, and duplicated fragments, are removed.
Graph construction: each fragment is split into k-mers, and every k-mer becomes an edge between its left and right (k−1)-mers in a directed multigraph.
Validity check: the graph is tested for connectivity and balanced in- and out-degrees, the conditions for an Eulerian path to exist.
Assembly: an Eulerian path is found with Hierholzer's algorithm, implemented from scratch, and converted back into the DNA sequence.

The program also saves a plot of the de Bruijn graph, and the main functions are covered by unit tests.

*Example:* for the input `DNA_2_5`, the cleaned fragments produced a graph with 26 nodes and 40 edges, which was assembled into the sequence `TTAATTACTCACTACGCACTGGGTCACTGGCTAATTACTCACTG`.

## How to Run
Place the input file in a folder called `csv-files` next to the script. File names follow the pattern `DNA_x_k.csv`, where `x` is the dataset number and k the k-mer size. Set `file_name` in `main_codes.py` (e.g. `file_name = "DNA_2_5"`) and run:

```bash
pip install pandas numpy networkx matplotlib
python main_codes.py
```
The reconstructed sequence is saved as DNA_x.txt and the graph as DNA_x.png. To run the tests:

```bash
pip install pytest
pytest test.py
```

## Tools
Python, pandas, NumPy, NetworkX, Matplotlib, pytest

## Author
Maryam Shojaei Shahrokhabadi · [LinkedIn](https://www.linkedin.com/in/maryam-shojaei-210740250)· [maryam.shojaei.ee@gmail.com]
