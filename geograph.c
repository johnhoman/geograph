#include <Python.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

/*
pair-wise distances will have to be counted for each
set of nodes added....is this feasible?

to do
----
I need to figure out if I want this to be pure c or
run in python, because I can't use the python code
base without running the interpreter
*/

#ifdef ALLOCATE_PyHEAP
#    define GeoGraph_ALLOC(x) PyMem_RawMalloc((x))
#    define GeoGraph_FREE(x) PyMem_RawFree((x))
#    define GeoGraph_REALLOC(ob, x) PyMem_RawRealloc((ob), (x))
#else
#    define GeoGraph_ALLOC(x) malloc((x))
#    define GeoGraph_FREE(x) free((x))
#    define GeoGraph_REALLOC(ob, x) realloc((ob), (x))
#endif

#ifndef GEONODE_DEFAULT_ALLOC
#  define GEONODE_DEFAULT_ALLOC 10
#endif

#define deallocator void

int tests_run = 0;

#define TEST(x) test_##x
#define FAIL() printf("%s -> failed (line %d)\n", __func__, __LINE__)
#define PASS(test) tests_run++; printf("%s -> passed - %d tests passed \n", #test, tests_run)
#define _assert(test) \
    do { \
        if (!(test)) { \
           FAIL(); \
           return 1; \
        } \
    } while(0)

#define _verify(test) \
    do { \
        int r = test(); \
        if(!r) { \
            PASS(test); \
        } \
    } while(0)

typedef struct _geonode {
    PyObject_HEAD
    PyObject *index;
    struct _geonode *next;
    struct _geonode *prev;
    double lat;
    double lon;
    struct _geonode **leafs;
    size_t _alloc;
    size_t _used;
    int8_t _visited;
} GeoNode;

/* double linked list */
typedef struct _geograph {
    PyObject_HEAD
    struct _geonode *head;
    struct _geonode *tail;
    double (* isedge)(GeoNode *, GeoNode *);
    double thresh;
} GeoGraph;

typedef struct _geograph_iter {
    PyObject_HEAD
    struct _geonode *_current;
    size_t _rank;
} GeoGraphIter;

typedef struct _geonode_iter {
    PyObject_HEAD
    GeoNode *_node;
    size_t _cur;
    size_t _size;
} GeoNodeIter;

typedef struct {
    PyObject_HEAD
    struct _geonode **nodes;
    size_t size;
} ConnectedEdges;

typedef struct {
    /* edge node */
    GeoNode *lv;
    GeoNode *rv;
} Edge;


static int
geonode_allocate_connections(GeoNode *nd)
{
    GeoNode ***con = &nd->leafs;
    size_t bytes = GEONODE_DEFAULT_ALLOC * sizeof(GeoNode *);
    *con = (GeoNode **)GeoGraph_ALLOC(bytes);
    if (con == NULL)
        return 0;
    else {
        nd->_used = 0;
        nd->_alloc = GEONODE_DEFAULT_ALLOC;
        return 1;
    }
    printf("internal error, should not reach code (LINE %d)", __LINE__);
    exit(0);
}

static GeoNode *
geonode_new(PyObject *index, double lat, double lon)
{
    GeoNode *node;

    node = (GeoNode *)GeoGraph_ALLOC(sizeof(GeoNode));
    if (node == NULL)
        return NULL;
    node->index = index;
    node->next = NULL;
    node->prev = NULL;
    node->lat = lat;
    node->lon = lon;
    if (!geonode_allocate_connections(node)) {
        GeoGraph_FREE(node);
        return NULL;
    } else {
        node->_visited = 0;
    }
    return node;
}

/* leafs are all connected via linked list,
   so they can be taken care of during the graph dealloc */
deallocator
geonode_dealloc(GeoNode *self)
{
    if (self != NULL)
        GeoGraph_FREE(self);
}


static GeoNodeIter *
geonodeiter_new(GeoNode *node)
{
    GeoNodeIter *it;

    if (node->_used == 0) {
        return NULL;
    } else {
        it = GeoGraph_ALLOC(sizeof(GeoNodeIter));
        it->_size = node->_used;
        it->_cur = 0;
        it->_node = node;
        return it;
    }
}

static int
geonode_richcompare(GeoNode *self, GeoNode *other, int opid)
{
    int result;
    if (self && other && (self == other || self->index == other->index))
        if (self->index && other->index)
            result = 1;
        else
            result = 0;
    else if (self->index == NULL || other->index == NULL)
        result = 0;
    else {
        result = PyObject_RichCompareBool(self->index, other->index, opid);
    }

    return result;
}

static GeoNode *
geonodeiter_next(GeoNodeIter *it)
{
    if (it == NULL || it->_cur == it->_size) {
        return NULL;
    } else {
        return it->_node->leafs[it->_cur++];
    }
}

/* check if node already exists in the leafs */
static int
geonode_exists(GeoNode *self, GeoNode *other)
{
    GeoNodeIter *it;
    GeoNode *nd;

    it = geonodeiter_new(self);
    if (it == NULL)
        return 0;

    while ((nd = geonodeiter_next(it)) != NULL) {
        /* points to same object or the index is the same
           and both cannot be null*/
        if (nd && other && (nd == other ||
                PyObject_RichCompareBool(other->index, nd->index, Py_EQ))) {
            GeoGraph_FREE(it);
            return 1;
        }
    }

    // does not exist
    return 0;
}

static int
geonode_appendnode(GeoNode *self, GeoNode *node)
{
    size_t sz;

    if (self != node && !geonode_exists(self, node)) {
        if (self->_alloc <= self->_used) {
            sz = self->_alloc * 2;
            void *ptr = GeoGraph_REALLOC(self->leafs, sz*sizeof(GeoNode *));
            if (ptr == NULL) {
                fprintf(stderr, "could not reallocate space for leafs\n");
                return 0;
            }
            self->_alloc = sz;
            self->leafs = (GeoNode **)ptr;
        }
        self->leafs[self->_used++] = node;
        self->leafs[self->_used] = NULL;
    }

    return 1;
}

static int
geonode_connectnodes(GeoNode *a, GeoNode *b)
{
    if (geonode_appendnode(a, b) && geonode_appendnode(b, a))
        return 1;
    return 0;
}

static int
geonode_isconnected(GeoNode *self, GeoNode *other,
        int (*func)(GeoNode *, GeoNode *))
{
    return func(self, other);
}

static double
geonode_distance(GeoNode *self, GeoNode *other)
{
    // for now just skip the ecef transformation and used lat/lon as 2d distances

    double delta_lat, delta_lon;

    delta_lat = self->lat - other->lat;
    delta_lon = self->lon - other->lon;

    return sqrt(delta_lat * delta_lat + delta_lon * delta_lon);
}

static GeoGraph *
geograph_new(void)
{
    GeoGraph *graph;

    graph = (GeoGraph *)GeoGraph_ALLOC(sizeof(GeoGraph));
    if (graph == NULL)
       return NULL;
    else {
       graph->head = NULL;
       graph->tail = NULL;
       graph->isedge = &geonode_distance;
       graph->thresh = 11;
    }
    return graph;
}

/* forward declare */
static GeoGraphIter *geographiter_new(GeoGraph *);
static GeoNode *geographiter_next(GeoGraphIter *);

deallocator
geograph_dealloc(GeoGraph *self)
{
    GeoGraphIter *it;
    GeoNode *nd;

    it = geographiter_new(self);
    if (it == NULL) // TODO: this will cause a memory leak
        return;

    while ((nd = geographiter_next(it)) != NULL)
        geonode_dealloc(nd);
}

static int
geograph_isempty(GeoGraph *graph)
{
    if (graph->head == NULL)
        return 1;
    else
        return 0;
}

static int
geograph_isedge(GeoGraph *graph, GeoNode *a, GeoNode *b, double d)
{
    /* geograph_isedge(graph, a, b, d) - determine based on relative distance
       whether or not a and b form and edge. */
    if (graph->isedge == NULL) {
        return 0;
    } else {
        return (*graph->isedge)(a, b) < d ? 1: 0;
    }
}

static GeoNode *geographiter_next(GeoGraphIter *graph);
static GeoGraphIter *geographiter_new(GeoGraph *graph);

static size_t
geograph_count(GeoGraph *graph)
{
    GeoGraphIter *it;

    it = geographiter_new(graph);
    while (geographiter_next(it) != NULL);

    return it->_rank + 1;
}

/* forward declare  */
static int _geograph_findedges(GeoGraph *, GeoNode *,
      double (*isedge)(GeoNode *, GeoNode *));

static GeoGraphIter *geographiter_new(GeoGraph *graph);

static int
geograph_contains(GeoGraph *graph, GeoNode *node)
{
    GeoGraphIter *it;
    GeoNode *nd;

    it = geographiter_new(graph);
    if (it == NULL)
        return 0;
    else {
        // TODO: not ideal. linear lookup
        while ((nd = geographiter_next(it)) != NULL) {
            if (geonode_richcompare(nd, node, Py_EQ))
                return 1;
        }
        return 0;
    }

}

static int
geograph_append(GeoGraph *graph, GeoNode *node)
{
    // graph the tail and append to that chain
    if (geograph_isempty(graph)) {
        graph->head = node;
        graph->tail = node;
    // don't append nodes with an index already in the graph
    } else if (!geograph_contains(graph, node)) {

        struct _geonode *last = graph->tail;
        /* don't want to drop any.
           the last link in the chain should have a NULL next pointer
        */
        if (last == NULL || last->next) {
            printf("error: last should not be null unless the graph is emtpy.");
            return 0;
        }
        last->next = node;
        node->prev = last;
        graph->tail = node;
        if (!_geograph_findedges(graph, node, graph->isedge)) {
            return 0;
        }
    }
    return 1;
}

static GeoGraphIter *
geographiter_new(GeoGraph *graph)
{
    GeoGraphIter *it;

    it = (GeoGraphIter *)GeoGraph_ALLOC(sizeof(GeoGraphIter));
    if (it == NULL)
        return NULL;
    else {
        it->_rank = -1; // pointing to head
        it->_current = graph->head;
    }
    return it;
}

static GeoNode *
geographiter_next(GeoGraphIter *it)
{
    GeoNode *next;
    if (it->_current == NULL)
        return NULL;
    else {
        next = it->_current;
        it->_current = it->_current->next;
        it->_rank++;
        return next;
    }
}

/* search the graph for all near neighbors within a distance threshold d */
static int
_geograph_findedges(GeoGraph *graph, GeoNode *node,
                   double (*isedge)(GeoNode *, GeoNode *))
{
    /* TODO: fix this. Super dangerous, should not mutate if "finding"
       something. Should return the edges that are connected and then
       assign them as leafs or something like that.
    */
    GeoGraphIter *it;
    GeoNode *nd;
    double distance, d = graph->thresh;

    if (isedge == NULL) {
        isedge = graph->isedge;
    }

    it = geographiter_new(graph);
    while ((nd = geographiter_next(it)) != NULL) {
         // not an edge if it's itself
         // nd == node should be the tail of the list
         if (nd == node) {
             continue;
         }
         // calculate distance
         distance = (*isedge)(node, nd);
         if (distance <= d) {
             if (geonode_connectnodes(node, nd))
                 return 1;
             else {
                 /* probably couldn't allocate memory. error message
                    should be handled in geonode_connectnodes. Just return 0
                 */
                 return 0;
             }
         }
    }
    return 1;
}

int TEST(geograph_append_01)(void)
{
    GeoGraph *gg = geograph_new();
    GeoNode *gn = geonode_new(NULL, 77.5, -87.3);

    geograph_append(gg, gn);
    _assert(gg->head == gn);
    _assert(gg->tail == gn);
    return 0;
}

int TEST(geograph_append_02)(void)
{
    GeoGraph *gg = geograph_new();
    GeoNode *gn[] = {
        geonode_new(NULL, 77.7, -87.3),
        geonode_new(NULL, 77.5, -87.4),
        geonode_new(NULL, 77.6, -87.5)
    };

    geograph_append(gg, gn[0]);
    geograph_append(gg, gn[1]);
    geograph_append(gg, gn[2]);

    _assert(gg->head == gn[0]);
    _assert(gg->tail == gn[2]);
    _assert(gn[0]->next == gn[1]);
    _assert(gn[1]->prev == gn[0]);
    _assert(gn[1]->next == gn[2]);
    _assert(gn[2]->prev == gn[1]);
    return 0;
}

int TEST(geonode_new_01)()
{
    GeoNode nd = *geonode_new(NULL, 88.8, -99.0);

    _assert(nd._alloc == GEONODE_DEFAULT_ALLOC);
    _assert(nd._used == 0);
    return 0;
}

int TEST(geonode_appendnode_01)()
{
    GeoNode *self = geonode_new(NULL, 88.3, -33.4 );
    GeoNode *other = geonode_new(NULL, 89.3, -23.4);

    geonode_appendnode(self, other);

    _assert(self->_used == 1);
    _assert(self->leafs[0] == other);
    _assert(self->leafs[1] == NULL); /* keep null terminated*/
    return 0;
}

int TEST(geonode_appendnode_realloc_01)()
{
    GeoNode *self = geonode_new(PyLong_FromLong(0), 88.3, -33.4);

    size_t i;
    for (i = 0; i < GEONODE_DEFAULT_ALLOC * 2 + 1; ++i) {
        _assert(geonode_appendnode(self,
            geonode_new(PyLong_FromLong(i + 1), 88.3 + (int)i, -33.4 - (int)i)
        ));
    }

    _assert(self->_used == GEONODE_DEFAULT_ALLOC * 2 + 1);
    _assert(self->_alloc == 40);
    return 0;
}

int TEST(geonode_connectnodes_01)()
{
    GeoNode *a = geonode_new(PyLong_FromLong(1), 88.3, -33.4);
    GeoNode *b = geonode_new(PyLong_FromLong(10), 89.3, -23.4);
    GeoNode *c = geonode_new(PyLong_FromLong(101), 87.3, -25.4);

    _assert(geonode_connectnodes(a, b));
    _assert(geonode_connectnodes(a, c));

    _assert(b->_used == 1);
    _assert(a->_used == 2);
    _assert(c->_used == 1);

    _assert(a->leafs[0] == b);
    _assert(a->leafs[1] == c);
    _assert(a->leafs[2] == NULL);

    _assert(b->leafs[0] == a); /* keep null terminated*/
    _assert(b->leafs[1] == NULL); /* keep null terminated*/
    return 0;
}

int TEST(geograph_findedges_01)()
{
    GeoGraph *gg = geograph_new();
    gg->thresh = 2.0;
    GeoNode *node;
    GeoNode *nodes[] = {
        geonode_new(PyLong_FromLong(0), 3., -4.),
        geonode_new(PyLong_FromLong(1), 3., -4.),
        geonode_new(PyLong_FromLong(2), 3., -4.),
        geonode_new(PyLong_FromLong(3), -3.0, 4.0),
        geonode_new(PyLong_FromLong(4), -3.0, 4.0),
        geonode_new(PyLong_FromLong(5), -3.0, 4.0)
    };

    size_t i;
    for (i = 0; i < 6; ++i) {
        _assert(geograph_append(gg, nodes[i]));
    }

    node = nodes[0];
    _assert(node->_used == 2);
    return 0;
}

int TEST(geograph_count_01)()
{
    GeoGraph *gg = geograph_new();
    GeoNode *nodes[] = {
        geonode_new(PyLong_FromLong(1), 88.3, -33.4),
        geonode_new(PyLong_FromLong(1), 89.3, -23.4),
        geonode_new(PyLong_FromLong(2), 87.3, -25.4),
        geonode_new(NULL, 98.3, -13.4),
        geonode_new(NULL, 99.3, -15.4),
        geonode_new(NULL, 97.3, -15.4)
    };

    size_t i;
    for (i = 0; i < 6; ++i) {
        _assert(geograph_append(gg, nodes[i]));
    }

    _assert(geograph_count(gg) == 5);

    geograph_append(gg, geonode_new(PyLong_FromLong(2), 66.4, 55.6)); // index already exists
    _assert(geograph_count(gg) == 5);
    GeoNode *nd = geonode_new(NULL, 55.5, 44.5);
    geograph_append(gg, nd);
    _assert(geograph_count(gg) == 6);
    return 0;
}

int TEST(geograph_connectedges)()
{
    GeoNode *a = geonode_new(PyLong_FromLong(1), 0, 0);
    GeoNode *b = geonode_new(PyLong_FromLong(2), 1, 0);
    GeoNode *c = geonode_new(PyLong_FromLong(3), 1, 0);

    _assert(geonode_connectnodes(a, b));

    _assert(a->_used == 1);
    _assert(a->leafs[0] == b);
    _assert(b->_used == 1);
    _assert(b->leafs[0] == a);
    _assert(geonode_connectnodes(a, c));
    _assert(a->_used == 2);
    _assert(a->leafs[1] == c);
    return 0;
}


int geograph_run_all_tests()
{
    _verify(TEST(geograph_append_01));
    _verify(TEST(geograph_append_02));
    _verify(TEST(geonode_new_01));
    _verify(TEST(geonode_appendnode_01));
    _verify(TEST(geonode_appendnode_realloc_01));
    _verify(TEST(geonode_connectnodes_01));
    _verify(TEST(geograph_connectedges));
    _verify(TEST(geograph_findedges_01));
    _verify(TEST(geograph_count_01));
    return 0;
}

int main(int argc, const char **argv)
{
    Py_Initialize();
    return geograph_run_all_tests();
}
