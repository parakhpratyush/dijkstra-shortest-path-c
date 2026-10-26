# Dijkstra's Shortest Path Algorithm in C (Binary Min-Heap)

Data Structure (BCS-301) - Project Based Learning (PBL), B.Tech CSE, ITS Engineering College, Greater Noida (2026-27).

Finds the minimum distance (and the actual route) from one source vertex to every other vertex of a weighted graph.

## What is inside
- Graph stored as an **adjacency list** (linked lists of edges)
- **Indexed binary min-heap** (insert, extract-min, decrease-key with a `pos[]` array) - O((V+E) log V)
- Simple **O(V^2) array version** for comparison
- **Bellman-Ford** as an independent checker
- Self-test (306 tests), benchmark, file input, input validation (negative weights and bad vertices are rejected)

## Build
```bash
gcc -O2 -Wall -Wextra -o dijkstra dijkstra.c
```
Windows (MinGW): `gcc -O2 -Wall -Wextra -o dijkstra.exe dijkstra.c`

## Run
```bash
./dijkstra demo                       # built-in 8-place road network, source = ITS College
./dijkstra path 0 7                   # shortest route from vertex 0 to vertex 7
./dijkstra file sample_graph.txt 0    # graph from a text file, source = vertex 0
./dijkstra test                       # heap vs array vs Bellman-Ford on 300 random graphs + edge cases
./dijkstra bench                      # speed comparison up to 16,000 vertices
```

## Input file format
```
V E
name_0 name_1 ... name_{V-1}   (one name per line or space separated, no spaces inside a name)
u v w                          (E lines: undirected edge between vertex u and v with weight w >= 0)
```

## Team
Nikhil Yadav (2502220100516)

Parakh Pratyush (2502220100546)

Niraj Prajapati (2502220100520)

## Guide
Prof. Ashish Kumar Gupta

SDG mapping: SDG 9 (Industry, Innovation and Infrastructure), SDG 11 (Sustainable Cities and Communities).
