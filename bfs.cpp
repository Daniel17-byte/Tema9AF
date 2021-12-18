#include <stdlib.h>
#include <string.h>
#include "bfs.h"
#include <iostream>

using namespace std;

int get_neighbors(const Grid* grid, Point p, Point neighb[])
{
    // TODO: fill the array neighb with the neighbors of the point p and return the number of neighbors
    // the point p will have at most 4 neighbors (up, down, left, right) => 4 if-uri
    // avoid the neighbors that are outside the grid limits or fall into a wall
    // note: the size of the array neighb is guaranteed to be at least 4 
    int  a = p.row, b = p.col;

    if (a > 0 && b > 0 && a < grid->rows && b < grid->cols && grid->mat[a][b] == 0) {

        int n = 0;
        //sus - N
        a = p.row - 1;
        b = p.col;
        if (a > 0 && b > 0 && a < grid->rows && b < grid->cols && grid->mat[a][b] == 0)
        {
            neighb[n].row = a;
            neighb[n].col = b;
            n++;
        }

        //jos - S
        a = p.row + 1;
        if (a > 0 && b > 0 && a < grid->rows && b < grid->cols && grid->mat[a][b] == 0)
        {
            neighb[n].row = a;
            neighb[n].col = b;
            n++;
        }

        //stanga - V
        a = p.row;
        b = p.col - 1;
        if (a > 0 && b > 0 && a < grid->rows && b < grid->cols && grid->mat[a][b] == 0)
        {
            neighb[n].row = a;
            neighb[n].col = b;
            n++;
        }

        //dreapta - E
        b = p.col + 1;
        if (a > 0 && b > 0 && a < grid->rows && b < grid->cols && grid->mat[a][b] == 0)
        {
            neighb[n].row = a;
            neighb[n].col = b;
            n++;
        }
        return n;
    }
    else
        return 0;

}

void grid_to_graph(const Grid* grid, Graph* graph)
{
    //we need to keep the nodes in a matrix, so we can easily refer to a position in the grid
    Node* nodes[MAX_ROWS][MAX_COLS];
    int i, j, k;
    Point neighb[4];

    //compute how many nodes we have and allocate each node
    graph->nrNodes = 0;
    for (i = 0; i < grid->rows; ++i) {
        for (j = 0; j < grid->cols; ++j) {
            if (grid->mat[i][j] == 0) {
                nodes[i][j] = (Node*)malloc(sizeof(Node));
                memset(nodes[i][j], 0, sizeof(Node)); //initialize all fields with 0/NULL
                nodes[i][j]->position.row = i;
                nodes[i][j]->position.col = j;
                ++graph->nrNodes;
            }
            else {
                nodes[i][j] = NULL;
            }
        }
    }
    graph->v = (Node**)malloc(graph->nrNodes * sizeof(Node*));
    k = 0;
    for (i = 0; i < grid->rows; ++i) {
        for (j = 0; j < grid->cols; ++j) {
            if (nodes[i][j] != NULL) {
                graph->v[k++] = nodes[i][j];
            }
        }
    }

    //compute the adjacency list for each node
    for (i = 0; i < graph->nrNodes; ++i) {
        graph->v[i]->adjSize = get_neighbors(grid, graph->v[i]->position, neighb);
        if (graph->v[i]->adjSize != 0) {
            graph->v[i]->adj = (Node**)malloc(graph->v[i]->adjSize * sizeof(Node*));
            k = 0;
            for (j = 0; j < graph->v[i]->adjSize; ++j) {
                if (neighb[j].row >= 0 && neighb[j].row < grid->rows &&
                    neighb[j].col >= 0 && neighb[j].col < grid->cols &&
                    grid->mat[neighb[j].row][neighb[j].col] == 0) {
                    graph->v[i]->adj[k++] = nodes[neighb[j].row][neighb[j].col];
                }
            }
            if (k < graph->v[i]->adjSize) {
                //get_neighbors returned some invalid neighbors
                graph->v[i]->adjSize = k;
                graph->v[i]->adj = (Node**)realloc(graph->v[i]->adj, k * sizeof(Node*));
            }
        }
    }
}

void free_graph(Graph* graph)
{
    if (graph->v != NULL) {
        for (int i = 0; i < graph->nrNodes; ++i) {
            if (graph->v[i] != NULL) {
                if (graph->v[i]->adj != NULL) {
                    free(graph->v[i]->adj);
                    graph->v[i]->adj = NULL;
                }
                graph->v[i]->adjSize = 0;
                free(graph->v[i]);
                graph->v[i] = NULL;
            }
        }
        free(graph->v);
        graph->v = NULL;
    }
    graph->nrNodes = 0;
}

typedef struct nQ {
    Node* n;
    struct nQ* next;
}NodeQ;

typedef struct q {
    NodeQ* head;
    NodeQ* tail;
}queue;

void enque(queue* Q, Node* ss) {
    NodeQ* s = (NodeQ*)malloc(sizeof(NodeQ));
    s->next = NULL;
    s->n = ss;


    if (Q->head == NULL) {
        Q->head = s;
        Q->tail = s;
    }
    else {
        Q->tail->next = s;
        Q->tail = s;
    }
}

Node* deque(queue* Q) {
    if (Q->head != NULL) {
        Node* s = Q->head->n;
        Q->head = Q->head->next;
        return s;
    }
    else {
        Q->head = NULL;
        Q->tail = Q->head;
        return NULL;
    }
}

void bfs(Graph* graph, Node* s, Operation* op)
{
    // TOOD: implement the BFS algorithm on the graph, starting from the node s
    // for counting the number of operations, the optional op parameter is received
    // since op can be NULL (when we are calling the bfs for display purposes), you should check it before counting:
    // if(op != NULL) op->count();

    for (int i = 0; i < graph->nrNodes; i++)
    {

        graph->v[i]->color = COLOR_WHITE;
        graph->v[i]->parent = NULL;
        graph->v[i]->dist = 0;

    }

    s->color = COLOR_GRAY;
    s->dist = 0;
    s->parent = NULL;
    queue* Q = (queue*)malloc(sizeof(queue));
    Q->head = NULL;
    Q->tail = Q->head;
    enque(Q, s);

    while (Q->head != NULL) {
        Node* aux = deque(Q);

        for (int ii = 0; ii < aux->adjSize; ii++) {
            if (op != NULL)
                op->count();
            if (aux->adj[ii]->color == COLOR_WHITE) {
                aux->adj[ii]->dist = aux->dist + 1;  // for all the visited nodes, the minimum distance from s (dist) and the parent in the BFS tree should be set
                aux->adj[ii]->color = COLOR_GRAY;
                aux->adj[ii]->parent = aux;
                if(op != NULL)
                    op->count(3);
                enque(Q, aux->adj[ii]);

            }
        }
        if (op != NULL)
            op->count();
        aux->color = COLOR_BLACK;
        // at the end of the algorithm, every node reachable from s should have the color BLACK
    }

}




void PP(int p[], Point repr[], int nivel, int root, int nrc)
{


    cout << endl;
    for (int j = 0; j < nivel; j++)
    {
        cout << "     ";
    }

    cout << "{" << repr[root].row << " x " << repr[root].col << "}" << endl;
    for (int j = 0; j < nrc; j++)
    {
        if (p[j] == root)
            PP(p, repr, nivel + 1, j, nrc);
    }



}

void print_bfs_tree(Graph* graph)
{
    //first, we will represent the BFS tree as a parent array
    int n = 0; //the number of nodes
    int* p = NULL; //the parent array
    Point* repr = NULL; //the representation for each element in p

    //some of the nodes in graph->v may not have been reached by BFS
    //p and repr will contain only the reachable nodes
    int* transf = (int*)malloc(graph->nrNodes * sizeof(int));
    for (int i = 0; i < graph->nrNodes; ++i) {
        if (graph->v[i]->color == COLOR_BLACK) {
            transf[i] = n;
            ++n;
        }
        else {
            transf[i] = -1;
        }
    }
    if (n == 0) {
        //no BFS tree
        free(transf);
        return;
    }

    int err = 0;
    p = (int*)malloc(n * sizeof(int));
    repr = (Point*)malloc(n * sizeof(Node));
    for (int i = 0; i < graph->nrNodes && !err; ++i) {
        if (graph->v[i]->color == COLOR_BLACK) {
            if (transf[i] < 0 || transf[i] >= n) {
                err = 1;
            }
            else {
                repr[transf[i]] = graph->v[i]->position;
                if (graph->v[i]->parent == NULL) {
                    p[transf[i]] = -1;
                }
                else {
                    err = 1;
                    for (int j = 0; j < graph->nrNodes; ++j) {
                        if (graph->v[i]->parent == graph->v[j]) {
                            if (transf[j] >= 0 && transf[j] < n) {
                                p[transf[i]] = transf[j];
                                err = 0;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    free(transf);
    transf = NULL;

    if (!err) {
        // TODO: pretty print the BFS tree
        // the parrent array is p (p[k] is the parent for node k or -1 if k is the root)
        // when printing the node k, print repr[k] (it contains the row and column for that point)
        // you can adapt the code for transforming and printing multi-way trees from the previous labs
        int root;
        for (int f = 0;f < n;f++)
        {
            if (p[f] == -1) {
                root = f; break;
            }
        }
        PP(p, repr, 0, root, n);

    }


    if (p != NULL) {
        free(p);
        p = NULL;
    }
    if (repr != NULL) {
        free(repr);
        repr = NULL;
    }
}

int shortest_path(Graph* graph, Node* start, Node* end, Node* path[])
{
    // TODO: compute the shortest path between the nodes start and end in the given graph
    // note: the size of the array path is guaranteed to be at least 1000

    bfs(graph, start, NULL); // pt a afla distanta de la start la end
    if (end->parent == NULL)
    {
        cout << "unreachable";
        return -1;
    }    // if end is not reachable from start, return -1

    else {
        int distanta = end->dist; //initializata prin parcurgere bfs a grafului incepand cu start
        for (int nod = end->dist - 1; nod >= 0; nod--)
        {
            path[nod] = end;
            end = end->parent;
        }
        // the nodes from the path, should be filled, in order, in the array path
        return distanta;    // the number of nodes filled in the path array should be returned

    }


}

void performance()
{
    int n, i;
    Profiler p("bfs");

    // vary the number of edges
    for (n = 1000; n <= 4500; n += 100) {
        Operation op = p.createOperation("bfs-edges", n);
        Graph graph;
        graph.nrNodes = 100;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node));
        for (i = 0; i < graph.nrNodes; ++i) {
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }
        // TODO: generate n random edges
        // make sure the generated graph is connected
        for (int i = 0; i < graph.nrNodes; i++)
        {
            graph.v[i]->adj = (Node**)malloc(graph.nrNodes * graph.nrNodes * sizeof(Node));
            graph.v[i]->adjSize = 0;
        }
        for (int i = 0; i < graph.nrNodes - 1; i++)
        {
            graph.v[i]->adj[graph.v[i]->adjSize] = graph.v[i + 1];
            graph.v[i]->adjSize++;
            graph.v[i + 1]->adj[graph.v[i + 1]->adjSize] = graph.v[i];
            graph.v[i + 1]->adjSize++;
        }
        for (int i = graph.nrNodes - 1; i < n; i++)
        {
            int nr1 = rand() % graph.nrNodes;
            int nr2 = rand() % graph.nrNodes;
            graph.v[nr1]->adjSize++;
            graph.v[nr1]->adj[graph.v[nr1]->adjSize - 1] = graph.v[nr2];
            graph.v[nr2]->adjSize++;
            graph.v[nr2]->adj[graph.v[nr2]->adjSize - 1] = graph.v[nr1];
        }

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    // vary the number of vertices
    for (n = 100; n <= 200; n += 10) {
        Operation op = p.createOperation("bfs-vertices", n);
        Graph graph;
        graph.nrNodes = n;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node));
        for (i = 0; i < graph.nrNodes; ++i) {
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }
        // TODO: generate 4500 random edges
        // make sure the generated graph is connected
        for (int i = 0; i < graph.nrNodes; i++)
        {
            graph.v[i]->adj = (Node**)malloc(graph.nrNodes * graph.nrNodes * sizeof(Node));
            graph.v[i]->adjSize = 0;
        }
        for (int i = 0; i < graph.nrNodes - 1; i++)
        {
            graph.v[i]->adj[graph.v[i]->adjSize] = graph.v[i + 1];
            graph.v[i]->adjSize++;
            graph.v[i + 1]->adj[graph.v[i + 1]->adjSize] = graph.v[i];
            graph.v[i + 1]->adjSize++;
        }
        for (int i = graph.nrNodes - 1; i < 4500; i++)
        {
            int nr1 = rand() % graph.nrNodes;
            int nr2 = rand() % graph.nrNodes;
            graph.v[nr1]->adjSize++;
            graph.v[nr1]->adj[graph.v[nr1]->adjSize - 1] = graph.v[nr2];
            graph.v[nr2]->adjSize++;
            graph.v[nr2]->adj[graph.v[nr2]->adjSize - 1] = graph.v[nr1];
        }

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    p.showReport();
}
