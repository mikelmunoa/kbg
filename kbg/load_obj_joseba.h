#ifndef WAVEFRONT_READER_H
#define WAVEFRONT_READER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#define MAXLINE 200

typedef struct {
    double x, y, z;
} coordinate;

typedef struct {
    int r, g, b;
} color;

typedef struct {
    coordinate coord;
    int num_faces;
} vertex;

typedef struct {
    int num_vertices;
    int* vertex_ind_table;
} face;

typedef struct {
    vertex* vertex_table;
    face* face_table;
    int num_vertices;
    int num_faces;
    coordinate min;
    coordinate max;
    color rgb;
} object3d;

int read_wavefront(char* file_name, object3d* object_ptr);

#endif /* WAVEFRONT_READER_H */
