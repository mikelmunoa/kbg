//	Program developed by
//	
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut
//
// 
//programa honek 3Dko objektuak marrazten ditu, eta testura bat erabiltzen du.


#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
// #include "cargar-triangulo.h"
#include <math.h>
#include "obj.h"

typedef struct erreferentzi_sistema
    {
    double xc[4];
    double yc[4];
    double zc[4];
    } erreferentzi_sistema;

// texture information

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);
unsigned char *bufferra;
int dimx,dimy;

int indexx;
int i,j,k;
face *triangulosptr;
object3d *foptr;
object3d *sel_ptr;
object3d *aux;
object3d *kamera;
object3d *berria;
light eguzkia;

double angelua;



int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;
double* E, *at, *up;
double *u, *v, *w;
double *view_matrix, result[16];
char mode, backulling, bektoreak, perspektiba,gouraud;



char fitxiz[100];




void objektuari_aldaketa_sartu_ezk(double m[16])
{
}



void objektuari_aldaketa_sartu_esk(double m[16])
{
}



/*********************************MATRIZE ETA BEKTORE ERAGIKETAK*********************************/



//Matrizeen arteko biderkete kalkulatu eta emaitza matrize batean gordeko du
void matrize_biderkaketa(double m1[16] , double m2[16], double emaitza[16]){
    for ( i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            double count = 0.0;
            for (k = 0; k < 4; k++) count += m1[i*4+k]*m2[k*4+j];
            emaitza[i*4+j] = count;
        }   
    }
}

//funtzio honek bi bektoreen arteko biderkaketa kalkulatzen du eta emaitza bektore batean gordeko du
void bektoreen_biderkaketa(double v1[4], double v2[4], double emaitza[4])
{
    emaitza[0] = v1[1]*v2[2] - v1[2]*v2[1];
    emaitza[1] = v1[2]*v2[0] - v1[0]*v2[2];
    emaitza[2] = v1[0]*v2[1] - v1[1]*v2[0];
    emaitza[3] = 0.0;
}


//funtzio honek bektore bat normalizatzen du eta emaitza bektore batean gordeko du
void bektore_normalizatu(double v[4], double emaitza[4])
{
    int i;
    double norma = 0.0;
    for ( i = 0; i < 3; i++)
    {
        norma += v[i]*v[i];
    }
    norma = sqrt(norma);
    for ( i = 0; i < 3; i++)
    {
        emaitza[i] = v[i]/norma;
    }
}


//matrize bat pantailaratzen du
void print_matrizea(char *str)
{
int i;

printf("%s\n",str);
for (i = 0;i<4;i++)
   printf("%lf, %lf, %lf, %lf\n",sel_ptr->mptr->m[i*4],sel_ptr->mptr->m[i*4+1],sel_ptr->mptr->m[i*4+2],
                                 sel_ptr->mptr->m[i*4+3]);
}


//matrizearen alderintzezkoa kalkulatzen du eta emaitza matrize batean gordeko du
void matrizearen_alderantzizkoa(double m[16], double inb[16]){
    for ( i = 0; i < 4; i++)for (j = 0; j < 4; j++)inb[i*4+j] = m[j*4+i];
}

void vertx_to_punto(vertex *vptr, punto *pptr)
{
    pptr->x = vptr->coord.x;
    pptr->y = vptr->coord.y;
    pptr->z = vptr->coord.z;
}


//hiruki baten bektore normala kalulatzen du eta emaitza bektore batean gordeko du
void hiruki_baten_bektore_normala(punto p1, punto p2, punto p3, double emaitza[4])
{
    double v1[4], v2[4];
    v1[0] = p2.x - p1.x;
    v1[1] = p2.y - p1.y;
    v1[2] = p2.z - p1.z;
    v1[3] = 0.0;
    v2[0] = p3.x - p1.x;
    v2[1] = p3.y - p1.y;
    v2[2] = p3.z - p1.z;
    v2[3] = 0.0;
    bektoreen_biderkaketa(v1,v2,emaitza);
    bektore_normalizatu(emaitza,emaitza);
}

void mxv(vertex *vptr, double m[16], vertex p)
{
    vptr->coord.x = p.coord.x * m[0] +p.coord.y * m[1] +p.coord.z * m[2] +  m[3];
    vptr->coord.y = p.coord.x * m[4] +p.coord.y * m[5] +p.coord.z * m[6] +  m[7];
    vptr->coord.z = p.coord.x * m[8] +p.coord.y * m[9] +p.coord.z * m[10] +  m[11];

}

double kosinua_kalkulatu(double *v1, double *v2){
    double biderketa = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    double norma1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    double norma2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    return biderketa/(norma1*norma2);
}

/*********************************OBJEKTUAREKIN ERAGIKETAK*********************************/

void hasieratu_kamara(){
    E = (double *)malloc(4*sizeof(double));
    at = (double *)malloc(4*sizeof(double));
    up = (double *)malloc(4*sizeof(double));

    u = (double *)malloc(4*sizeof(double));
    v = (double *)malloc(4*sizeof(double));
    w = (double *)malloc(4*sizeof(double));

    view_matrix = (double*)malloc(16*sizeof(double));

    u[0] = 1.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;

    v[0] = 0.0;
    v[1] = 1.0;
    v[2] = 0.0;
    v[3] = 0.0;

    w[0] = 0.0;
    w[1] = 0.0;
    w[2] = 1.0;
    w[3] = 0.0;

    E[0] = 0.0;
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = 1.0;

    at[0] = 0.0;
    at[1] = 0.0;
    at[2] = 0.0;
    at[3] = 1.0;

    up[0] = 0.0;
    up[1] = 1.0;
    up[2] = 0.0;
    up[3] = 1.0;

    kamera = (object3d *)malloc(sizeof(object3d));
    kamera->mptr = (mlist *)malloc(sizeof(mlist));
    for(i=0;i<16;i++)kamera->mptr->m[i]=0.0;
    read_wavefront("cam.obj", kamera);
    for(i=0;i<4;i++){
        kamera->mptr->m[i*4+i] = 1.0;
    }



}
//puntu bat matrizearekin biderkatzen du
void mxp(punto *pptr, double m[16], punto p)
{
pptr->x = p.x * m[0] +p.y * m[1] +p.z * m[2] +  m[3];
pptr->y = p.x * m[4] +p.y * m[5] +p.z * m[6] +  m[7];
pptr->z = p.x * m[8] +p.y * m[9] +p.z * m[10] +  m[11];
pptr->u = p.u;
pptr->v = p.v;
}

void calculate_camera_position(){
    punto lag;
    lag.x = 0;
    lag.y = 0;
    lag.z = 100;
    mxp(&lag,kamera->mptr->m,lag);
    E[0] = lag.x;
    E[1] = lag.y;
    E[2] = lag.z;
}


void calculate_camera_coordenates(double *E, double *at, double *up, double *u, double *v, double *w)
{
    double *aux;
    aux = (double *)malloc(4*sizeof(double));
    for(i=0;i<4;i++)aux[i]=0.0;
    for(i=0;i<4;i++)u[i]=0.0;
    for(i=0;i<4;i++)v[i]=0.0;
    for(i=0;i<4;i++)w[i]=0.0;
    for(i=0;i<4;i++)aux[i] = at[i] - E[i];
    bektore_normalizatu(aux,w);
    bektoreen_biderkaketa(w,up,u);
    bektoreen_biderkaketa(u,w,v);
    free(aux);
}

void calculate_projection_matrix(double *V, double *P)
{
    double n = 1.0;
    double f = 1000.0;
    double t = 1.0;
    double b = -1.0;
    double r = 1.0;
    double l = -1.0;
    P[0] = 2*n/(r-l);
    P[1] = 0.0;
    P[2] = 0.0;
    P[3] = 0.0;
    P[4] = 0.0;
    P[5] = 2*n/(t-b);
    P[6] = 0.0;
    P[7] = 0.0;
    P[8] = (r+l)/(r-l);
    P[9] = (t+b)/(t-b);
    P[10] = -(f+n)/(f-n);
    P[11] = -1.0;
    P[12] = 0.0;
    P[13] = 0.0;
    P[14] = -2*f*n/(f-n);
    P[15] = 0.0;
}


void kalkulatu_Mesa(object3d cam,double *emaitza){
    double *Mo;
    double *V;
    Mo = (double *)malloc(16*sizeof(double));
    V = (double *)malloc(16*sizeof(double));
    calculate_camera_coordenates(E,at,up,u,v,w);
    Mo[0] = u[0];
    Mo[1] = u[1];
    Mo[2] = u[2];
    Mo[3] = 0.0;
    Mo[4] = v[0];
    Mo[5] = v[1];
    Mo[6] = v[2];
    Mo[7] = 0.0;
    Mo[8] = w[0];
    Mo[9] = w[1];
    Mo[10] = w[2];
    Mo[11] = 0.0;
    Mo[12] = 0.0;
    Mo[13] = 0.0;
    Mo[14] = 0.0;
    Mo[15] = 1.0;

    V[0] = 1.0;
    V[1] = 0.0;
    V[2] = 0.0;
    V[3] = -E[0];
    V[4] = 0.0;
    V[5] = 1.0;
    V[6] = 0.0;
    V[7] = -E[1];
    V[8] = 0.0;
    V[9] = 0.0;
    V[10] = 1.0;
    V[11] = -E[2];
    V[12] = 0.0;
    V[13] = 0.0;
    V[14] = 0.0;
    V[15] = 1.0;

    matrize_biderkaketa(Mo,V,emaitza);
    free(Mo);
    free(V);
}


//x-ren gaineko eragiketak
void x_aldaketa(int dir)
{
    mlist *matrize = (mlist *)malloc(sizeof(mlist));
    matrize->hptr = sel_ptr->mptr;
    double ald_mat[16];
    for (int i = 0; i < 16; i++) {
        ald_mat[i] = 0.0;
    }

    if (aldaketa == 't') {
        ald_mat[0] = 1.0;
        ald_mat[5] = 1.0;
        ald_mat[10] = 1.0;
        ald_mat[15] = 1.0;
        ald_mat[3] = 20 * (dir - 0.5);
        matrize_biderkaketa(ald_mat, sel_ptr->mptr->m, matrize->m);
    }
    else if (aldaketa == 'r') {
        ald_mat[0] = 1.0;
        ald_mat[5] = 0.968912421710;
        ald_mat[6] = -0.247403959254;
        ald_mat[9] = 0.247403959254;
        ald_mat[10] = 0.968912421710;
        ald_mat[15] = 1.0;
        if(ald_lokala == 1){
            matrize_biderkaketa( sel_ptr->mptr->m, ald_mat,matrize->m);
        }
        else{
            matrize_biderkaketa(ald_mat, sel_ptr->mptr->m, matrize->m);
            }
    }
    sel_ptr->mptr = matrize;
}


//y-ren gaineko eragiketak
void y_aldaketa(int dir)
{
    mlist *matrize = (mlist *)malloc(sizeof(mlist));
    matrize->hptr = sel_ptr->mptr;
    double ald_mat[16];
    for(i = 0;i<16;i++)ald_mat[i]=0.0;
    if (aldaketa == 't'){
            ald_mat[0] = 1.0;
            ald_mat[5] = 1.0;
            ald_mat[10] = 1.0;
            ald_mat[15] = 1.0;
            ald_mat[7] = 20*(dir-0.5);
            matrize_biderkaketa(ald_mat,sel_ptr->mptr->m,matrize->m);
    }
    else if(aldaketa == 'r'){
        ald_mat[0] = 0.968912421710;
        ald_mat[2] = 0.247403959254;
        ald_mat[5] = 1.0;
        ald_mat[8] = -0.247403959254;
        ald_mat[10] = 0.968912421710;
        ald_mat[15] = 1.0;
        if(ald_lokala ==1){
            matrize_biderkaketa(sel_ptr->mptr->m,ald_mat,matrize->m);
        }
        else{
            matrize_biderkaketa(ald_mat,sel_ptr->mptr->m,matrize->m);
        }

    }
        sel_ptr->mptr = matrize;
}

//z-ren gaineko eragiketak
void z_aldaketa(int dir)
{
    print_matrizea;
    mlist *matrize = (mlist *)malloc(sizeof(mlist));
    matrize->hptr = sel_ptr->mptr;
    double ald_mat[16];
    for(i = 0;i<16;i++)ald_mat[i]=0.0;

    if (aldaketa == 't'){
        ald_mat[0] = 1.0;
        ald_mat[5] = 1.0;
        ald_mat[10] = 1.0;
        ald_mat[15] = 1.0;
        ald_mat[11] = 20*(dir-0.5);
        matrize_biderkaketa(ald_mat,sel_ptr->mptr->m,matrize->m);
    }
    else if(aldaketa == 'r'){
        ald_mat[0] = 0.968912421710;
        ald_mat[1] = -0.247403959254;
        ald_mat[4] = 0.247403959254;
        ald_mat[5] = 0.968912421710;
        ald_mat[10] = 1;
        if(ald_lokala ==1){
            matrize_biderkaketa(sel_ptr->mptr->m,ald_mat,matrize->m);
        }
        else{
            matrize_biderkaketa(ald_mat,sel_ptr->mptr->m,matrize->m);
        }
    }
    sel_ptr->mptr = matrize;
}

//undo funtzioa
void undo()
{
    if(sel_ptr->mptr->hptr != 0){
        sel_ptr->mptr = sel_ptr->mptr->hptr;
    }
}


/**********************************ARGIAREN FUNTZIOAK**********************************/

void argiak_hasieratu(){
    eguzkia.onoff = 1;
    eguzkia.type = 0;
    eguzkia.I.r=0.5;
    eguzkia.I.g=0.5;
    eguzkia.I.b=0.5;
    eguzkia.dir[0] = 0.0;
    eguzkia.dir[1] = 0.0;
    eguzkia.dir[2] = 1.0;
    eguzkia.dir[3] = 0.0;
}

void intentsitatea_kalkulatu(unsigned char *colorv, object3d *optr, int aurpegi_index, int erpin_index) {
    double NL;
    // Calculate the ambient component
    colorv[0] = (unsigned char)(optr->Ka.r * 255);
    colorv[1] = (unsigned char)(optr->Ka.g * 255);
    colorv[2] = (unsigned char)(optr->Ka.b * 255);

    if(gouraud == 'y'){
        NL = kosinua_kalkulatu(eguzkia.dir, optr->vertex_table[erpin_index].N);
    }else{
        NL = angelua;
    }
    if(NL <0) NL=0;

    // Calculate the diffuse component
    colorv[0] += (unsigned char)(eguzkia.I.r * eguzkia.onoff*255 * optr->kd.r * NL);
    colorv[1] += (unsigned char)(eguzkia.I.g * eguzkia.onoff*255 * optr->kd.g * NL);
    colorv[2] += (unsigned char)(eguzkia.I.b * eguzkia.onoff*255 * optr->kd.b * NL);

    // Calculate the specular component
    colorv[0] += (unsigned char)(eguzkia.I.r * 255 * optr->ks.r);
    colorv[1] += (unsigned char)(eguzkia.I.g * 255 * optr->ks.g);
    colorv[2] += (unsigned char)(eguzkia.I.b * 255 * optr->ks.b);
}

/*********************************TESTURA ETA MARRAZTE ERAGIKETAK*********************************/



//texturako pixelak kalkulatzen eta bueltatzen ditu
unsigned char * color_textura(float u, float v,unsigned char *colorv ,object3d *optr)
{
    //kalkulatu objektuaren kolorea
    int x = (int)(u*dimx);
    int y = (int)(v*dimy);




    return colorv;
}


// lerro bat marrazten du
void  dibujar_linea_z(int linea,float c1x, float c1z, float c1u,float c1v,float c2x,float c2z,float c2u,float c2v,object3d *optr,int aurpegi)
{
float xkoord,zkoord;
float u,v;
int erpin;
unsigned char r,g,b;
unsigned char *colorv;
float malda_z;
float koef;
float bektore_normala;
unsigned char *erpin_koloreak;
if(gouraud == 'n'){
    colorv = (unsigned char *)malloc(3*sizeof(unsigned char));
    erpin = 0;
    intentsitatea_kalkulatu(colorv,optr,aurpegi,erpin);
}
else{
    colorv = (unsigned char *)malloc(3*optr->face_table[aurpegi].num_vertices*sizeof(unsigned char));
    for(i=0;i<optr->face_table[aurpegi].num_vertices;i++){
        erpin = optr->face_table[aurpegi].vertex_ind_table[i];
        intentsitatea_kalkulatu(colorv + 3*i,optr,aurpegi,erpin);
    }
    
}

glBegin( GL_POINTS );
for (xkoord = c1x,zkoord =c1z, u = c1u, v=c1v; xkoord <= c2x; xkoord ++)
    {
    colorv=  color_textura(u, v,colorv,optr); 
    r= colorv[0];
    g=colorv[1];
    b=colorv[2];    
    glColor3ub(r,g,b);
    glVertex3f(xkoord, linea, zkoord );
    zkoord = c1z + (c2z-c1z)*(xkoord-c1x)/(c2x-c1x);
    u = c1u + (c2u-c1u)*(xkoord-c1x)/(c2x-c1x);
    v = c1v + (c2v-c1v)*(xkoord-c1x)/(c2x-c1x);
    }
glEnd();
}


//Ebaki puntuak kalkulatzen ditu
void ebaketaKalkulatu(vertex *gptr,vertex *bptr,int h, punto *ebakipuntu){
    float dif_y = gptr->coord.y - bptr->coord.y;
    float dif_x = gptr->coord.x - bptr->coord.x;
    float dif_z = gptr->coord.z - bptr->coord.z;
    // float dif_u = gptr->coord.u - bptr->coord.u;
    // float dif_v = gptr->coord.v - bptr->coord.v;
    float altuera = gptr->coord.y - h;
    ebakipuntu->y = h;
    if(dif_y != 0){
    ebakipuntu->x = gptr->coord.x - altuera*dif_x/dif_y;
    ebakipuntu->z = gptr->coord.z - altuera*dif_z/dif_y;
    // ebakipuntu->u = gptr->coord.u - altuera*dif_u/dif_y;
    // ebakipuntu->v = gptr->coord.v - altuera*dif_v/dif_y;
    }

}


//Triangelua marrazten duen funtzioa
void dibujar_triangulo(object3d *optr, int i)
{
face *tptr;
double normal_vector[4];
vertex *pgoiptr, *pbeheptr, *perdiptr, *lag;
float x1,h1,z1,u1,v1,x2,h2,z2,u2,v2,x3,h3,z3,u3,v3;
float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v;
int linea;
int h;
float cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
vertex p1,p2,p3;
punto l[3];


if (i >= optr->num_faces) return;
tptr = optr->face_table + i;

for(j=0;j<tptr->num_vertices;j++){
    p1 = optr->vertex_table[tptr->vertex_ind_table[j]];
    mxv(&p1,optr->mptr->m,p1);
    vertx_to_punto(&p1,&l[j]);
}

hiruki_baten_bektore_normala(l[0],l[1],l[2],tptr->N);
for(j=0;j<tptr->num_vertices;j++){
    optr->vertex_table[tptr->vertex_ind_table[j]].N[0] += tptr->N[0];
    optr->vertex_table[tptr->vertex_ind_table[j]].N[1] += tptr->N[1];
    optr->vertex_table[tptr->vertex_ind_table[j]].N[2] += tptr->N[2];
    bektore_normalizatu(optr->vertex_table[tptr->vertex_ind_table[j]].N,optr->vertex_table[tptr->vertex_ind_table[j]].N);
}







if (backulling == 'b'){
    double kosinu;
    kosinu = kosinua_kalkulatu(tptr->N,w);
    if(kosinu>0)return;

    
} 

if(gouraud == 'n'){
    angelua = kosinua_kalkulatu(tptr->N,eguzkia.dir);
}





p1 = optr->vertex_table[tptr->vertex_ind_table[0]];
p2 = optr->vertex_table[tptr->vertex_ind_table[1]];
p3 = optr->vertex_table[tptr->vertex_ind_table[2]];

kalkulatu_Mesa(*kamera,view_matrix);
matrize_biderkaketa(view_matrix,optr->mptr->m,result);
if(perspektiba == 'p'){
    double *P;
    P = (double *)malloc(16*sizeof(double));
    calculate_projection_matrix(view_matrix,P);
    matrize_biderkaketa(P,result,result);
    free(P);
}
if (lineak == 1)
{
    glBegin(GL_POLYGON);
    for(i=0;i<tptr->num_vertices;i++){
        p1 = optr->vertex_table[tptr->vertex_ind_table[i]];
        mxv(&p1,result,p1);
        vertx_to_punto(&p1,&l[i]);
        glVertex3d(p1.coord.x,p1.coord.y,p1.coord.z);
    }
    glEnd();
    if(bektoreak == 'v'){
        glBegin(GL_LINES);
        glColor3f(1.0,0.0,0.0);
        glVertex3d(l[0].x,l[0].y,l[0].z);
        glVertex3d(l[0].x+50*tptr->N[0],l[0].y+50*tptr->N[1],l[0].z+50*tptr->N[2]);
        glEnd();
    }
    return;
}


mxv(&p1,result,p1);
mxv(&p2,result,p2);
mxv(&p3,result,p3);

if(p1.coord.y>p2.coord.y){
    pgoiptr = &(p1);
    pbeheptr = &(p2);
}else{
    pgoiptr = &(p2);
    pbeheptr = &(p1);
}
if(p3.coord.y>pgoiptr->coord.y){
    perdiptr = pgoiptr;
    pgoiptr = &(p3);
}else if(p3.coord.y<pbeheptr->coord.y){
    perdiptr = pbeheptr;
    pbeheptr = &(p3);
}else{
    perdiptr = &(p3);
}

//set triangle vertices in order of x

if(pgoiptr->coord.y == perdiptr->coord.y && pgoiptr->coord.x<perdiptr->coord.x){
    lag = pgoiptr;
    pgoiptr = perdiptr;
    perdiptr = lag;
}
if(pbeheptr->coord.y == perdiptr->coord.y && pbeheptr->coord.x<perdiptr->coord.x){
    lag = pbeheptr;
    pbeheptr = perdiptr;
    perdiptr = lag;
}
if(pbeheptr->coord.y == pgoiptr->coord.y && pbeheptr->coord.x<pgoiptr->coord.x){
    lag = pbeheptr;
    pbeheptr = pgoiptr;
    pgoiptr = lag;
}

punto e1, e2;




//triangelua marrazten du
for(h =pgoiptr->coord.y;h>perdiptr->coord.y;h--){
    ebaketaKalkulatu(pgoiptr,perdiptr,h,&e1);
    ebaketaKalkulatu(pgoiptr,pbeheptr,h,&e2);
    if(e1.x<e2.x){
        dibujar_linea_z(h,e1.x,e1.z,e1.u,e1.v,e2.x,e2.z,e2.u,e2.v,optr,i);
    }
    else{
        dibujar_linea_z(h,e2.x,e2.z,e2.u,e2.v,e1.x,e1.z,e1.u,e1.v,optr,i);

    }
    
}
for(h=perdiptr->coord.y;h>pbeheptr->coord.y;h--){
    ebaketaKalkulatu(perdiptr,pbeheptr,h,&e1);
    ebaketaKalkulatu(pgoiptr,pbeheptr,h,&e2);
    if(e1.x<e2.x){
        dibujar_linea_z(h,e1.x,e1.z,e1.u,e1.v,e2.x,e2.z,e2.u,e2.v,optr,i);
    }
    else{
        dibujar_linea_z(h,e2.x,e2.z,e2.u,e2.v,e1.x,e1.z,e1.u,e1.v,optr,i);
    }
    
}

}


// This function will be called whenever the screen has to be drawn
static void marraztu(void)
{
float u,v;
int i,j;
object3d *auxptr;
/*
unsigned char* colorv;
unsigned char r,g,b;
*/

  // marrazteko objektuak behar dira
  // no se puede dibujar sin objetos
if (foptr ==0) return;

// clear viewport...
if (objektuak == 1) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
    else 
      {
      if (denak == 0) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
      }

glMatrixMode(GL_PROJECTION);
glLoadIdentity();
glOrtho(-5000.0, 5000.0, -5000.0, 5000.0,-5000.0, 5000.0);

calculate_camera_position();
triangulosptr = sel_ptr->face_table;
if (objektuak == 1)
    {
    if (denak == 1)
        {
        for (auxptr =foptr; auxptr != 0; auxptr = auxptr->hptr)
            {
            for (i =0; i < auxptr->num_faces; i++)
                {
                dibujar_triangulo(auxptr,i);
                }
            }
        }
      else
        {
        for (i =0; i < sel_ptr->num_faces; i++)
            {
            dibujar_triangulo(sel_ptr,i);
            }
        }
    }
  else
    {
     dibujar_triangulo(sel_ptr,indexx);
    }
glFlush();





}

/*********************************USER INPUT*********************************/



static void teklatua (unsigned char key, int x, int y)
{
int retval;
int i;
FILE *obj_file;

switch(key)
	{
	case 13: 
	        if (foptr != 0)  // objekturik ez badago ezer ez du egin behar
	                         // si no hay objeto que no haga nada
	            {
	            indexx ++;  // azkena bada lehenengoa bihurtu
		                // pero si es el último? hay que controlarlo!
		    if (indexx == sel_ptr->num_faces) 
		        {
		        indexx = 0;
		        if ((denak == 1) && (objektuak == 0))
		            {
		            glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
		            glFlush();
		            }
		        }
		    }
		break;
	case 'd':
		if (denak == 1) denak = 0;
		    else denak = 1;
		break;
	case 'o':
		if (objektuak == 1) objektuak = 0;
		    else objektuak = 1;
		break;
	case 'l':
		if (lineak == 1) lineak = 0;
		    else lineak = 1;
		break;
	case 't':
	        aldaketa = 't';
		break;
	case 'r':
		aldaketa = 'r';
		break;
    case 'b':
        if (backulling == 'b') backulling = 'f';
            else backulling = 'b';
        break;
	case 'g':
		if (ald_lokala == 1) ald_lokala = 0;
		    else ald_lokala = 1;
		break;
        case 'x':
                x_aldaketa(1);
                break;
        case 'y':
                y_aldaketa(1);
                break;
        case 'z':
                z_aldaketa(1);
                break;
        case 'X':
                x_aldaketa(0);
                break;
        case 'Y':
                y_aldaketa(0);
                break;
        case 'Z':
                z_aldaketa(0);
                break;
        case 'u':
                undo();
                break;
    case 'c':
            mode = 'n';
            aux = sel_ptr;
            sel_ptr = kamera;
            print_matrizea("kamera mugitzen");
        break;
    case 'C':
        //change camera
            mode = 'C';
            aldaketa = 'r';
        break;
    case 'v':
            if(bektoreak == 'v'){
                bektoreak = 'e';
            }
            else{
                bektoreak = 'v';
            }
        break;
    case 'p':
        if(perspektiba == 'p'){
            perspektiba = 'o';
        }
        else{
            perspektiba = 'p';
        }
        aux = sel_ptr;
        sel_ptr = kamera;

        
        break;
	case 'f':
    /*Ask for file*/
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
            berria = (object3d*)malloc(sizeof(object3d));
            berria->mptr = (mlist *)malloc(sizeof(mlist));
	        printf("read wavefront with output: %d\n",read_wavefront(fitxiz,berria));
            for(i=0;i<16;i++)berria->mptr->m[i]=0.0;
            berria->mptr->m[0]=1.0;
            berria->mptr->m[5]=1.0;
            berria->mptr->m[10]=1.0;
            berria->mptr->m[15]=1.0;  

	        indexx = 0;
            foptr = berria;
            berria->hptr = sel_ptr;

            
            for(aux =foptr;aux!=0;aux=aux->hptr){
                printf("aux: %p\t",aux);
                printf("berria: %p\t",berria);
                printf("sel_ptr: %p\n",sel_ptr);
            }
    break;
    case 'j':
    if(gouraud == 'n'){
        gouraud = 'y';

    }else{
        gouraud = 'n';
    }
    
    break;
    case '0':
        if(eguzkia.onoff == 1){
            eguzkia.onoff = 0;
        }else{
            eguzkia.onoff = 1;
        }

    break;
       /* case 'S':  // save to file
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
                if ((obj_file = fopen(fitxiz, "w")) == NULL)
                         {
                         printf("ezin fitxategia ireki\n");
                         }
                     else
                         {
                         for (i =0; i < sel_ptr->num_triangles; i++)
                            {
                            fprintf(obj_file,"t %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                                 sel_ptr->triptr[i].p1.x-250, sel_ptr->triptr[i].p1.y-250, sel_ptr->triptr[i].p1.z, 
                                 sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v,
                                 sel_ptr->triptr[i].p2.x-250, sel_ptr->triptr[i].p2.y-250, sel_ptr->triptr[i].p2.z, 
                                 sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v,
                                 sel_ptr->triptr[i].p3.x-250, sel_ptr->triptr[i].p3.y-250, sel_ptr->triptr[i].p3.z, 
                                 sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v );
                            }
                         fclose(obj_file);
                         }
                break; */
        case 9: /* <TAB> */
            if (foptr != 0) // objekturik gabe ez du ezer egin behar
                            // si no hay objeto no hace nada
                {
                sel_ptr = sel_ptr->hptr;
                /*The selection is circular, thus if we move out of the list we go back to the first element*/
                if (sel_ptr == 0) sel_ptr = foptr;
                indexx =0; // the selected polygon is the first one
                }
            break;
	case 27:  // <ESC>
		exit( 0 );
		break;
	default:
		printf("%d %c\n", key, key );
	}

// The screen must be drawn to show the new triangle
glutPostRedisplay();
}




//main funtzioa
int main(int argc, char** argv)
{
int retval;

	printf("Press <ESC> to finish\n");
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_RGB|GLUT_DEPTH );
	glutInitWindowSize ( 1000, 1000 );
	glutInitWindowPosition ( 0, 0 );
	glutCreateWindow( "KBG/GO praktika" );
	
	glutDisplayFunc( marraztu );
	glutKeyboardFunc( teklatua );
	/* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */ 
        retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
        if (retval !=1) 
            {
            printf("Ez dago texturaren fitxategia (testura.ppm)\n");
            exit(-1);
            }
        
	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        denak = 1;
        lineak =0;
        objektuak = 1;
        foptr = 0;
        sel_ptr = 0;
        aldaketa = 'r';
        ald_lokala = 1;
        backulling = 'b';
        gouraud = 'n';
        foptr = (object3d *)malloc(sizeof(object3d));
        foptr->mptr = (mlist *)malloc(sizeof(mlist));
        sel_ptr = (object3d *)malloc(sizeof(object3d));
        sel_ptr->mptr = (mlist *)malloc(sizeof(mlist));
        for(i=0;i<16;i++)foptr->mptr->m[i]=0.0;
        foptr->mptr->m[0]=1.0;
        foptr->mptr->m[5]=1.0;
        foptr->mptr->m[10]=1.0;
        foptr->mptr->m[15]=1.0;
        if (argc>1) {
            read_wavefront(argv[1],foptr);
            sel_ptr = foptr;
        }
        else{
            read_wavefront("k.obj", foptr);
            sel_ptr = foptr;
        } 
        sel_ptr->Ka.r = 0.24725f;
        sel_ptr->Ka.g = 0.1995f;
        sel_ptr->Ka.b = 0.0745f;
        sel_ptr->kd.r = 0.75164f;
        sel_ptr->kd.g = 0.60648f;
        sel_ptr->kd.b = 0.22648f;
        sel_ptr->ks.b = 0.628281f;
        sel_ptr->ks.g = 0.555802f;
        sel_ptr->ks.r = 0.366065f;


                
            
    hasieratu_kamara();
    argiak_hasieratu();
    glScalef(10,10,10);
	glutMainLoop();

	return 0;   
}