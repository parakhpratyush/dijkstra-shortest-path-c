/*
 * Dijkstra's Shortest-Path Algorithm using a Binary Min-Heap (Priority Queue)
 * Course  : Data Structure (BCS-301), B.Tech CSE, ITS Engineering College
 * Features: adjacency list, indexed binary min-heap with decrease-key,
 *           O(V^2) array version for comparison, Bellman-Ford oracle,
 *           random self-test, benchmark, file input.
 * Build   : gcc -O2 -Wall -Wextra -o dijkstra dijkstra.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#define INF INT_MAX
#define NAME_LEN 32
#define MAX_NODES 64

/* ---------------------------------------------------------------- Graph */
typedef struct Edge {
    int to;
    int w;
    struct Edge *next;
} Edge;

typedef struct {
    int V;
    Edge **adj;
    char names[MAX_NODES][NAME_LEN]; /* used only by the demo / file graphs */
} Graph;

static Graph *graph_create(int V) {
    Graph *g = (Graph *)calloc(1, sizeof(Graph));
    if (!g) { perror("calloc"); exit(1); }
    g->V = V;
    g->adj = (Edge **)calloc((size_t)V, sizeof(Edge *));
    if (!g->adj) { perror("calloc"); exit(1); }
    return g;
}

static void graph_free(Graph *g) {
    for (int i = 0; i < g->V; i++) {
        Edge *e = g->adj[i];
        while (e) { Edge *n = e->next; free(e); e = n; }
    }
    free(g->adj);
    free(g);
}

/* returns 0 on success, -1 on invalid input (bad vertex / negative weight) */
static int graph_add_edge(Graph *g, int u, int v, int w, int undirected) {
    if (u < 0 || v < 0 || u >= g->V || v >= g->V) return -1;
    if (w < 0) return -1; /* Dijkstra is only valid for non-negative weights */
    Edge *e = (Edge *)malloc(sizeof(Edge));
    if (!e) { perror("malloc"); exit(1); }
    e->to = v; e->w = w; e->next = g->adj[u]; g->adj[u] = e;
    if (undirected) {
        Edge *r = (Edge *)malloc(sizeof(Edge));
        if (!r) { perror("malloc"); exit(1); }
        r->to = u; r->w = w; r->next = g->adj[v]; g->adj[v] = r;
    }
    return 0;
}

/* ------------------------------------------------------------- Min-Heap */
typedef struct { int v; int dist; } HeapNode;

typedef struct {
    HeapNode *a;   /* 0-based array representation of the complete tree   */
    int *pos;      /* pos[v] = index of vertex v inside a[], -1 if absent */
    int size;
} MinHeap;

static MinHeap *heap_create(int cap) {
    MinHeap *h = (MinHeap *)malloc(sizeof(MinHeap));
    if (!h) { perror("malloc"); exit(1); }
    h->a = (HeapNode *)malloc((size_t)cap * sizeof(HeapNode));
    h->pos = (int *)malloc((size_t)cap * sizeof(int));
    if (!h->a || !h->pos) { perror("malloc"); exit(1); }
    h->size = 0;
    for (int i = 0; i < cap; i++) h->pos[i] = -1;
    return h;
}

static void heap_free(MinHeap *h) { free(h->a); free(h->pos); free(h); }

static void heap_swap(MinHeap *h, int i, int j) {
    HeapNode t = h->a[i]; h->a[i] = h->a[j]; h->a[j] = t;
    h->pos[h->a[i].v] = i;
    h->pos[h->a[j].v] = j;
}

static void sift_up(MinHeap *h, int i) {
    while (i > 0) {
        int p = (i - 1) / 2;
        if (h->a[p].dist <= h->a[i].dist) break;
        heap_swap(h, p, i);
        i = p;
    }
}

static void sift_down(MinHeap *h, int i) {
    for (;;) {
        int l = 2 * i + 1, r = 2 * i + 2, s = i;
        if (l < h->size && h->a[l].dist < h->a[s].dist) s = l;
        if (r < h->size && h->a[r].dist < h->a[s].dist) s = r;
        if (s == i) break;
        heap_swap(h, i, s);
        i = s;
    }
}

static void heap_insert(MinHeap *h, int v, int dist) {
    int i = h->size++;
    h->a[i].v = v; h->a[i].dist = dist; h->pos[v] = i;
    sift_up(h, i);
}

static HeapNode heap_extract_min(MinHeap *h) {
    HeapNode top = h->a[0];
    h->size--;
    h->pos[top.v] = -1;
    if (h->size > 0) {
        h->a[0] = h->a[h->size];
        h->pos[h->a[0].v] = 0;
        sift_down(h, 0);
    }
    return top;
}

static void heap_decrease_key(MinHeap *h, int v, int dist) {
    int i = h->pos[v];
    h->a[i].dist = dist;
    sift_up(h, i);
}

/* ------------------------------------------------- Dijkstra (heap) O(E log V) */
static void dijkstra_heap(const Graph *g, int src, int *dist, int *parent) {
    MinHeap *h = heap_create(g->V);
    for (int v = 0; v < g->V; v++) { dist[v] = INF; parent[v] = -1; }
    dist[src] = 0;
    for (int v = 0; v < g->V; v++) heap_insert(h, v, dist[v]);
    while (h->size > 0) {
        HeapNode u = heap_extract_min(h);
        if (u.dist == INF) break;               /* rest are unreachable */
        for (Edge *e = g->adj[u.v]; e; e = e->next) {
            if (h->pos[e->to] == -1) continue;  /* already finalised */
            long long nd = (long long)u.dist + e->w;
            if (nd < dist[e->to]) {             /* relaxation */
                dist[e->to] = (int)nd;
                parent[e->to] = u.v;
                heap_decrease_key(h, e->to, (int)nd);
            }
        }
    }
    heap_free(h);
}

/* ------------------------------------------------- Dijkstra (array) O(V^2) */
static void dijkstra_array(const Graph *g, int src, int *dist, int *parent) {
    char *done = (char *)calloc((size_t)g->V, 1);
    for (int v = 0; v < g->V; v++) { dist[v] = INF; parent[v] = -1; }
    dist[src] = 0;
    for (int it = 0; it < g->V; it++) {
        int u = -1;
        for (int v = 0; v < g->V; v++)
            if (!done[v] && (u == -1 || dist[v] < dist[u])) u = v;
        if (u == -1 || dist[u] == INF) break;
        done[u] = 1;
        for (Edge *e = g->adj[u]; e; e = e->next) {
            long long nd = (long long)dist[u] + e->w;
            if (!done[e->to] && nd < dist[e->to]) { dist[e->to] = (int)nd; parent[e->to] = u; }
        }
    }
    free(done);
}

/* --------------------------------------- Bellman-Ford (independent oracle) */
static void bellman_ford(const Graph *g, int src, int *dist) {
    for (int v = 0; v < g->V; v++) dist[v] = INF;
    dist[src] = 0;
    for (int pass = 0; pass < g->V - 1; pass++) {
        int changed = 0;
        for (int u = 0; u < g->V; u++) {
            if (dist[u] == INF) continue;
            for (Edge *e = g->adj[u]; e; e = e->next) {
                long long nd = (long long)dist[u] + e->w;
                if (nd < dist[e->to]) { dist[e->to] = (int)nd; changed = 1; }
            }
        }
        if (!changed) break;
    }
}

/* ------------------------------------------------------------- Printing */
static void print_path_rec(const Graph *g, const int *parent, int v) {
    if (parent[v] != -1) { print_path_rec(g, parent, parent[v]); printf("->"); }
    printf("%s", g->names[v]);
}

static void print_table(const Graph *g, int src, const int *dist, const int *parent) {
    printf("\nShortest distances from source: %s\n", g->names[src]);
    printf("%-3s %-9s %s\n", "ID", "Dist(km)", "Shortest path (ends at the destination)");
    printf("-----------------------------------------------------------------------------\n");
    for (int v = 0; v < g->V; v++) {
        printf("%-3d ", v);
        if (dist[v] == INF) { printf("%-9s %s\n", "INF", "unreachable"); continue; }
        printf("%-9d ", dist[v]);
        print_path_rec(g, parent, v);
        printf("\n");
    }
}

static void print_adjacency(const Graph *g) {
    printf("\nAdjacency list (%d vertices):\n", g->V);
    for (int u = 0; u < g->V; u++) {
        printf("  [%d] %-14s:", u, g->names[u]);
        for (Edge *e = g->adj[u]; e; e = e->next) printf(" (%d,w=%d)", e->to, e->w);
        printf("\n");
    }
}

/* ---------------------------------------------------------- Demo graph */
static Graph *build_demo_graph(void) {
    static const char *nm[] = {
        "ITS College", "Pari Chowk", "Knowledge Pk", "Alpha-1",
        "Botanical Gdn", "Sector 62", "City Centre", "Sector 18"
    };
    Graph *g = graph_create(8);
    for (int i = 0; i < 8; i++) strncpy(g->names[i], nm[i], NAME_LEN - 1);
    /* illustrative road distances in km (sample data, not surveyed values) */
    graph_add_edge(g, 0, 1, 7, 1);
    graph_add_edge(g, 0, 2, 3, 1);
    graph_add_edge(g, 1, 3, 4, 1);
    graph_add_edge(g, 2, 3, 5, 1);
    graph_add_edge(g, 2, 4, 14, 1);
    graph_add_edge(g, 3, 4, 8, 1);
    graph_add_edge(g, 3, 5, 12, 1);
    graph_add_edge(g, 4, 6, 6, 1);
    graph_add_edge(g, 5, 6, 3, 1);
    graph_add_edge(g, 5, 7, 9, 1);
    graph_add_edge(g, 6, 7, 5, 1);
    graph_add_edge(g, 1, 4, 10, 1);
    return g;
}

/* ------------------------------------------------------- File-based graph
 * Format:  V E
 *          name_0 ... (one name per line, V lines)
 *          u v w  (E lines, undirected, 0-based ids)                        */
static Graph *load_graph_file(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return NULL; }
    int V, E;
    if (fscanf(f, "%d %d", &V, &E) != 2 || V <= 0 || V > MAX_NODES) {
        fprintf(stderr, "Bad header (1 <= V <= %d)\n", MAX_NODES); fclose(f); return NULL;
    }
    Graph *g = graph_create(V);
    for (int i = 0; i < V; i++) {
        if (fscanf(f, "%31s", g->names[i]) != 1) { fprintf(stderr, "Bad name list\n"); fclose(f); graph_free(g); return NULL; }
    }
    for (int i = 0; i < E; i++) {
        int u, v, w;
        if (fscanf(f, "%d %d %d", &u, &v, &w) != 3 || graph_add_edge(g, u, v, w, 1) != 0) {
            fprintf(stderr, "Invalid edge #%d (bad vertex or negative weight)\n", i + 1);
            fclose(f); graph_free(g); return NULL;
        }
    }
    fclose(f);
    return g;
}

/* ---------------------------------------------------- Random graph tools */
static Graph *random_graph(int V, int E, int maxw) {
    Graph *g = graph_create(V);
    for (int v = 1; v < V; v++)                 /* spanning chain keeps it connected */
        graph_add_edge(g, v - 1, v, 1 + rand() % maxw, 1);
    for (int i = V - 1; i < E; i++)
        graph_add_edge(g, rand() % V, rand() % V, 1 + rand() % maxw, 1);
    return g;
}

static double now_ms(void) { return 1000.0 * (double)clock() / CLOCKS_PER_SEC; }

/* ------------------------------------------------------------ Self test */
static int run_self_test(void) {
    srand(2026);
    int total = 0, pass = 0;
    /* Test 1: 300 random graphs, heap == array == Bellman-Ford */
    for (int t = 0; t < 300; t++) {
        int V = 2 + rand() % 60, E = V + rand() % (4 * V);
        Graph *g = random_graph(V, E, 50);
        int *d1 = malloc(V * sizeof(int)), *p1 = malloc(V * sizeof(int));
        int *d2 = malloc(V * sizeof(int)), *p2 = malloc(V * sizeof(int));
        int *d3 = malloc(V * sizeof(int));
        int src = rand() % V;
        dijkstra_heap(g, src, d1, p1);
        dijkstra_array(g, src, d2, p2);
        bellman_ford(g, src, d3);
        int ok = 1;
        for (int v = 0; v < V; v++) if (d1[v] != d2[v] || d1[v] != d3[v]) ok = 0;
        total++; pass += ok;
        free(d1); free(p1); free(d2); free(p2); free(d3); graph_free(g);
    }
    printf("[T1] 300 random graphs: heap == array == Bellman-Ford : %d/%d passed\n", pass, total);

    int local_total = 0, local_pass = 0;
    /* Test 2: single vertex */
    { Graph *g = graph_create(1); int d[1], p[1]; dijkstra_heap(g, 0, d, p);
      local_total++; local_pass += (d[0] == 0); graph_free(g); }
    /* Test 3: disconnected vertex stays INF */
    { Graph *g = graph_create(3); graph_add_edge(g, 0, 1, 4, 1);
      int d[3], p[3]; dijkstra_heap(g, 0, d, p);
      local_total++; local_pass += (d[1] == 4 && d[2] == INF); graph_free(g); }
    /* Test 4: negative weight rejected */
    { Graph *g = graph_create(2); local_total++;
      local_pass += (graph_add_edge(g, 0, 1, -5, 1) == -1); graph_free(g); }
    /* Test 5: invalid vertex rejected */
    { Graph *g = graph_create(2); local_total++;
      local_pass += (graph_add_edge(g, 0, 9, 1, 1) == -1); graph_free(g); }
    /* Test 6: zero-weight edges + parallel edges (shorter one wins) */
    { Graph *g = graph_create(3);
      graph_add_edge(g, 0, 1, 9, 1); graph_add_edge(g, 0, 1, 2, 1); graph_add_edge(g, 1, 2, 0, 1);
      int d[3], p[3]; dijkstra_heap(g, 0, d, p);
      local_total++; local_pass += (d[1] == 2 && d[2] == 2); graph_free(g); }
    /* Test 7: demo graph known answer (ITS College -> Sector 18 = 27 km, hand-verified) */
    { Graph *g = build_demo_graph(); int d[8], p[8]; dijkstra_heap(g, 0, d, p);
      local_total++; local_pass += (d[7] == 27 && d[3] == 8 && d[6] == 22); graph_free(g); }
    printf("[T2-T7] 6 edge-case tests (see Table 3.3)  : %d/%d passed\n",
           local_pass, local_total);
    total += local_total; pass += local_pass;
    printf("\nOVERALL: %d/%d tests passed -> %s\n", pass, total, pass == total ? "ALL OK" : "FAILURE");
    return pass == total ? 0 : 1;
}

/* ------------------------------------------------------------ Benchmark */
static void run_benchmark(void) {
    srand(7);
    int sizes[] = {500, 1000, 2000, 4000, 8000, 16000};
    printf("%-8s %-9s %-16s %-16s %-8s\n", "V", "E", "Heap O(ElogV) ms", "Array O(V^2) ms", "Speedup");
    printf("--------------------------------------------------------------\n");
    for (int k = 0; k < 6; k++) {
        int V = sizes[k], E = 5 * V;
        Graph *g = random_graph(V, E, 100);
        int *d = malloc(V * sizeof(int)), *p = malloc(V * sizeof(int));
        int reps = 5;
        double t0 = now_ms();
        for (int r = 0; r < reps; r++) dijkstra_heap(g, 0, d, p);
        double th = (now_ms() - t0) / reps;
        t0 = now_ms();
        for (int r = 0; r < reps; r++) dijkstra_array(g, 0, d, p);
        double ta = (now_ms() - t0) / reps;
        char sp[16]; snprintf(sp, sizeof sp, "%.1fx", th > 0 ? ta / th : 0.0);
        printf("%-8d %-9d %-16.3f %-16.3f %-8s\n", V, E, th, ta, sp);
        free(d); free(p); graph_free(g);
    }
}

/* ----------------------------------------------------------------- main */
static void usage(const char *p) {
    printf("Usage:\n  %s demo                 run on built-in city graph from ITS College\n"
           "  %s path <src> <dst>     shortest path between two vertex ids on demo graph\n"
           "  %s file <graph.txt> <src>   run on graph from file\n"
           "  %s test                 self-test (heap vs array vs Bellman-Ford)\n"
           "  %s bench                benchmark heap vs array\n", p, p, p, p, p);
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    if (!strcmp(argv[1], "test")) return run_self_test();
    if (!strcmp(argv[1], "bench")) { run_benchmark(); return 0; }

    Graph *g = NULL; int src = 0, dst = -1;
    if (!strcmp(argv[1], "demo")) { g = build_demo_graph(); }
    else if (!strcmp(argv[1], "path") && argc == 4) { g = build_demo_graph(); src = atoi(argv[2]); dst = atoi(argv[3]); }
    else if (!strcmp(argv[1], "file") && argc == 4) { g = load_graph_file(argv[2]); src = atoi(argv[3]); }
    else { usage(argv[0]); return 1; }
    if (!g) return 1;
    if (src < 0 || src >= g->V || dst >= g->V) { fprintf(stderr, "Vertex id out of range (0..%d)\n", g->V - 1); graph_free(g); return 1; }

    int *dist = malloc(g->V * sizeof(int)), *parent = malloc(g->V * sizeof(int));
    dijkstra_heap(g, src, dist, parent);
    if (!strcmp(argv[1], "demo")) print_adjacency(g);
    if (dst >= 0) {
        printf("\nQuery: %s  ->  %s\n", g->names[src], g->names[dst]);
        if (dist[dst] == INF) printf("No path exists.\n");
        else { printf("Minimum distance = %d km\nRoute            = ", dist[dst]); print_path_rec(g, parent, dst); printf("\n"); }
    } else print_table(g, src, dist, parent);
    free(dist); free(parent); graph_free(g);
    return 0;
}
