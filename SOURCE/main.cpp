#include <windows.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <cmath>
#include <bits/stdc++.h>
#include "BmpLoader.h"
#include <string>
#include <sstream>

using namespace std;

//int font=(int)GLUT_BITMAP_8_BY_13;
// zz is font position x_look is side position
unsigned int ID;
//BmpLoader bl;'
float rot=0;

bool pause = false;
bool start= true;

static float zz = 50;
static float yy=2;
static float x_look = 0;
double kon=0.1;
double factor=0.2;
double eyex=0.0+ x_look;
double eyey=1.0;
double eyez=7.5 + zz;
double posx=x_look;
double posy=0.0;
double posz=-10.0+zz;

static float score = 0;
static float final_score = 0;

GLfloat alpha = 0.0, theta = 0.0, axis_x=0.0, axis_y=0.0;
double Txval=0,Tyval=0,Tzval=0;
double moving=0;

///from curve code///

const double PI = 3.14159265389;

/* GLUT callback Handlers */

int anglex= 0, angley = 0, anglez = 0;          //rotation angles
int window;
int wired=0;
int shcpt=1;
int animat = 0;
const int L=3;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;				//number of slices along x-direction
const int ntheta = 20;
int counter=0;

GLfloat ctrlpoints[L+1][3] =
{
    // {2,3,0},{0,3,0},{0,-3,0},{2,-3,0}
    //{-7,2,0},{-6,2,0},{-5,2,0},{-4,2,0},{-3,3,0},{-2,4,0},{-1,5,0},{0,6,0},{1,5,0},{2,4,0},{3,3,0},{4,2,0},{5,2,0},{6,2,0},{7,2,0}
    {-2,2,0},{-1,2,0},{0,2,0},{1,2,0}
};


double ex=0, ey=0, ez=15, lx=0,ly=0,lz=0, hx=0,hy=1,hz=0;

float wcsClkDn[3],wcsClkUp[3];
///////////////////////////////
class point1
{
public:
    point1()
    {
        x=0;
        y=0;
    }
    int x;
    int y;
} clkpt[2];
int flag=0;
GLint viewport[4]; //var to hold the viewport info
GLdouble modelview[16]; //var to hold the modelview info
GLdouble projection[16]; //var to hold the projection matrix info

//////////////////////////
void scsToWcs(float sx,float sy, float wcsv[3] );
void processMouse(int button, int state, int x, int y);
void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back=0, float ambFactor=1.0, float specFactor=1.0);

void scsToWcs(float sx,float sy, float wcsv[3] )
{

    GLfloat winX, winY, winZ; //variables to hold screen x,y,z coordinates
    GLdouble worldX, worldY, worldZ; //variables to hold world x,y,z coordinates

    //glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
    glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
    glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

    winX = sx;
    winY = (float)viewport[3] - (float)sy;
    winZ = 0;

    //get the world coordinates from the screen coordinates
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);
    wcsv[0]=worldX;
    wcsv[1]=worldY;
    wcsv[2]=worldZ;
}


//control points
long long nCr(int n, int r)
{
    if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);
}

///////////////////////
void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}

void bottleBezier()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta/2; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glTexCoord2d(theta,t-dt);
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glTexCoord2d(theta,t);
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}

//curve code end//

static void resize(int width, int height)
{
    const float ar = (float) width / (float) height;
    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float d=1;
    glFrustum(-ar*d, ar*d, -1.0*d, 1.0*d, 2.0, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;
}

static GLfloat v_cube[8][3] =
{
    {0,0,0},
    {0,0,1},
    {0,1,0},
    {0,1,1},
    {1,0,0},
    {1,0,1},
    {1,1,0},
    {1,1,1}
};
static GLubyte c_ind[6][4] =
{
    {3,1,5,7},
    {2,0,1,3},
    {7,5,4,6},
    {2,3,7,6},
    {1,0,4,5},
    {6,4,0,2}
};


static void getNormal3p(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat x3, GLfloat y3, GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

void drawcube(float cr, float cg, float cb,int n=1, bool e=false)
{
    GLfloat m_no[] = {0, 0, 0, 1.0};
    GLfloat m_amb[] = {cr,cg,cb,1};
    GLfloat m_diff[] = {cr,cg,cb,1};
    GLfloat m_spec[] = {1,1,1,1};
    GLfloat m_sh[] = {90};

    GLfloat m_em[] = {1,1,1,1};

    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);

    if(e)
    {
        glMaterialfv(GL_FRONT,GL_EMISSION, m_diff);
    }
    else
    {
        glMaterialfv(GL_FRONT,GL_EMISSION, m_no);
    }



    glBegin(GL_QUADS);
    for(GLint i = 0; i<6; i++)
    {
        //glColor3f(cr,cg,cb);
        getNormal3p(v_cube[c_ind[i][0]][0], v_cube[c_ind[i][0]][1], v_cube[c_ind[i][0]][2],
                    v_cube[c_ind[i][1]][0], v_cube[c_ind[i][1]][1], v_cube[c_ind[i][1]][2],
                    v_cube[c_ind[i][2]][0], v_cube[c_ind[i][2]][1], v_cube[c_ind[i][2]][2]);

        glVertex3fv(&v_cube[c_ind[i][0]][0]);
        glTexCoord2f(0,0);
        glVertex3fv(&v_cube[c_ind[i][1]][0]);
        glTexCoord2f(n,0);
        glVertex3fv(&v_cube[c_ind[i][2]][0]);
        glTexCoord2f(n,n);
        glVertex3fv(&v_cube[c_ind[i][3]][0]);
        glTexCoord2f(0,n);
    }
    glEnd();
}

void LoadTexture(const char*filename)
{
    glGenTextures(1, &ID);
    glBindTexture(GL_TEXTURE_2D, ID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, ID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}

void rotation(void)
{
    double e1,e2,e3,e4,e5,degree;
    e1=eyex-posx;
    e3=eyez-posz;
    degree = kon*3.14159/180;//convert degree to radian
    //result = sin(x);
    e4=(e1*cos(degree))+(e3*sin(degree));
    e5=(e3*cos(degree))-(e1*sin(degree));
    eyex=e4+posx;
    eyez=e5+posz;
}


void Cylinder3D(double a1,double b1,double c1)
{

    GLfloat mat_ambient[] = { 1, 1, 1, 1.0 };
    GLfloat mat_diffuse[] = { 1, 1, 1, 1.0 };
    GLfloat mat_specular[] = { 1,1,1, 1.0 };
    GLfloat mat_shininess[] = {90};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    GLUquadricObj *quadratic;
    
    quadratic = gluNewQuadric();
    gluQuadricTexture(quadratic, GL_TRUE);
    glRotatef(-90.0f, 1.0f,0.0f, 0.0f);
    gluCylinder(quadratic,a1,b1,c1,32,32);
    //gluQuadricTexture(quadratic, TRUE);
    gluDeleteQuadric(quadratic);
}

void circle_3D(GLdouble radius)
{
    GLfloat mat_ambient[] = { 1, 1, 1, 1.0 };
    GLfloat mat_diffuse[] = { 1, 1, 1, 1.0 };
    GLfloat mat_specular[] = { 1,1,1, 1.0 };
    GLfloat mat_shininess[] = {90};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    GLUquadric *qobj = gluNewQuadric();
    gluQuadricTexture(qobj, GL_TRUE);

    glRotatef(270, 1, 0, 0);
    gluSphere(qobj, radius, 20, 20);
    gluDeleteQuadric(qobj);
}


void Tree()
{
    int randm;
    //srand(5);
    randm = (rand() % 9) + 8;

    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D,19);
    glEnable(GL_TEXTURE_2D);
    Cylinder3D(0.7,0.6,randm);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,randm+1,0);
    glBindTexture(GL_TEXTURE_2D,2);
    glEnable(GL_TEXTURE_2D);
    circle_3D(3.5);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

void Tire()
{
    glPushMatrix();
    GLfloat mat_ambient[] = { 1, 1, 1, 1.0 };
    GLfloat mat_diffuse[] = { 1, 1, 1, 1.0 };
    GLfloat mat_specular[] = { 0.1,0.1,0.1, 1.0 };
    GLfloat mat_shininess[] = {90};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);
    glutSolidTorus(0.2, 0.8, 5, 50);
    glPopMatrix();
    glPushMatrix();
    glTranslated(0,-0.8,0);
    glScaled(0.1,1.5,0.1);
    drawcube(1,1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.8,0.0,0);
    glRotatef(90,0,0,1);
    glScaled(0.1,1.5,0.1);
    drawcube(1,1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.1,0.0,0);
    glRotatef(45,0,0,1);
    glPushMatrix();
    glTranslated(0.8,0.0,0);
    glRotatef(90,0,0,1);
    glScaled(0.1,1.5,0.1);
    drawcube(1,1,1,1);
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glRotatef(-45,0,0,1);
    glPushMatrix();
    glTranslated(0.8,0.0,0);
    glRotatef(90,0,0,1);
    glScaled(0.1,1.5,0.1);
    drawcube(1,1,1,1);
    glPopMatrix();
    glPopMatrix();
}

void reff(void)
{
    glPushMatrix();
    glScaled(1,10,1);
    drawcube(0,1,0);
    glPopMatrix();
}

void light()
{
    //light
    GLfloat l_amb[] = {0.1, 0.1, 0.1, 1.0};
    GLfloat l_dif[] = {1,1,1,1};
    GLfloat l_spec[] = {0.2,0.2,0.2,1};
    //GLfloat l_pos1[] = {l1pos,1.0};
    GLfloat l_pos1[] = {0,50,50,1.0};
    GLfloat l_pos2[] = {0,50,-850,1.0};
    //GLfloat l_pos3[] = {l3posx,l3posy,l3posz,1.0};

    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0,GL_AMBIENT,l_amb);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,l_dif);
    glLightfv(GL_LIGHT0,GL_SPECULAR,l_spec);
    glLightfv(GL_LIGHT0,GL_POSITION,l_pos1);

    glEnable(GL_LIGHT1);
    glLightfv(GL_LIGHT1,GL_AMBIENT,l_amb);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,l_dif);
    glLightfv(GL_LIGHT1,GL_SPECULAR,l_spec);
    glLightfv(GL_LIGHT1,GL_POSITION,l_pos2);
}


void axis(void)
{
    glBegin(GL_LINES);
    glColor3f (1.0, 0.0, 0.0); ///red is X
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(50.0, 0.0, 0.0);

    glColor3f (0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 50.0, 0.0); /// green is Y

    glColor3f (0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 50.0); ///blue is Z
    glEnd();
}

void Floor(void)
{
    for(int i=1; i<180; i++)
    {

        if(i<=18||i>=34)
        {
            glBindTexture(GL_TEXTURE_2D,2);
            glEnable(GL_TEXTURE_2D);

            glPushMatrix();
            glTranslatef(-100,-3,i*-10);
            glScaled(200,1,110);
            drawcube(.322,.745,.5,30);
            glPopMatrix();
            glDisable(GL_TEXTURE_2D);
        }
    }
    ///pani
    glBindTexture(GL_TEXTURE_2D,26);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(-10,-3.5,21*-10);
    glScaled(20,1,110);
    drawcube(1,1,1,30);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    
   for(int i=-9; i<11; i++)
    {

        if(i<0||i>1)
        {
            glBindTexture(GL_TEXTURE_2D,2);
            glEnable(GL_TEXTURE_2D);

            glPushMatrix();
            glTranslatef(i*-10,-3,-230);
            glScaled(10,1,115);
            drawcube(0.322,0.745,0.5,30);
            glPopMatrix();
            glDisable(GL_TEXTURE_2D);
        }
    }

}

void track(void)
{

    ///road
    for(int i=1; i<180; i++)
    {
        if(i<=19||i>=22)
        {
            glBindTexture(GL_TEXTURE_2D,1);
            glEnable(GL_TEXTURE_2D);
            glPushMatrix();
            glTranslatef(-7.5,-2.9,-10*i);
            //reff();
            glScaled(17,1,10);
            drawcube(1,1,1,1);
            glPopMatrix();
            glDisable(GL_TEXTURE_2D);
        }
    }
}


void piler()
{
    for(int i=0; i<180; i++)
    {
        glPushMatrix();
        glTranslatef(-8.5,-2.9,-10*i);
        glPushMatrix();
        glTranslatef(-3.7,0,-10*i);
        glScaled(1,1,1);
        Tree();
        glPopMatrix();

        GLfloat mat_ambient[] = { 1, 1, 1, 1.0 };
        GLfloat mat_diffuse[] = { 1, 1, 1, 1.0 };
        GLfloat mat_specular[] = { 0.1,0.1,0.1, 1.0 };
        GLfloat mat_shininess[] = {90};

        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

        glBindTexture(GL_TEXTURE_2D,6);
        glEnable(GL_TEXTURE_2D);

        glPushMatrix();

        GLUquadricObj *quadratic;

        quadratic = gluNewQuadric();
        gluQuadricTexture(quadratic, GL_TRUE);
        glRotatef(-90.0f, 1.0f,0.0f, 0.0f);
        gluCylinder(quadratic,0.2f,0.2f,3.0f,32,32);
        gluDeleteQuadric(quadratic);
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();
    }
    for(int i=0; i<180; i++)
    {

        glPushMatrix();
        glTranslatef(10,-2.9,-10*i);
	        glPushMatrix();
	        glTranslatef(2.8,0,-10*i);
	        //glScaled(1,1,1);
	        Tree();
	        glPopMatrix();
	        //glPopMatrix();
	
	        GLfloat mat_ambient[] = { 1, 1, 1, 1.0 };
	        GLfloat mat_diffuse[] = { 1, 1, 1, 1.0 };
	        GLfloat mat_specular[] = { 0.1,0.1,0.1, 1.0 };
	        GLfloat mat_shininess[] = {80};
	
	        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
	        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
	        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
	        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);
	
	        glBindTexture(GL_TEXTURE_2D,6);
	        glEnable(GL_TEXTURE_2D);
	        glPushMatrix();
	        //glColor3d(1,0,0);
	        GLUquadricObj *quadratic;
	        quadratic = gluNewQuadric();
	        gluQuadricTexture(quadratic, GL_TRUE);
	        glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
	        gluCylinder(quadratic,0.2f,0.2f,3.0f,32,32);
	        //gluQuadricTexture(quadratic, TRUE);
	        gluDeleteQuadric(quadratic);
	        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();
    }
}

void footpath()
{
    for(int i=1; i<180; i++)
    {

        glBindTexture(GL_TEXTURE_2D,10);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(-15,-2.7,-10*i);
        //reff();
        glScaled(7.5,1,10);
        drawcube(1,1,1,3);
        glPopMatrix();

        glDisable(GL_TEXTURE_2D);
    }
    for(int i=1; i<180; i++)
    {
        // glBindTexture(GL_TEXTURE_2D,1);
        //glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,10);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(9,-2.7,-10*i);
        //reff();
        glScaled(6,1,10);
        drawcube(1,1,1,3);
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
    }
}

void building(void)
{
    /// building
    for(int i = 1; i<18; i++)
    {
        //srand(i);
        int rand1,rand2,rand3,rand4;
        srand(i);
        rand1 = (rand() % 13) + 7;
        srand(i+1);
        rand2 = (rand() % 13) + 7;
        srand(i+2);
        rand3 = (rand() % 13) + 7;
        srand(i+3);
        rand4 = (rand() % 13) + 7;
        //printf("%d :rand1 %d rand2 %d rand3 \n",rand1,rand2,rand3);

        ///right buildings
        glBindTexture(GL_TEXTURE_2D,11);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(15,-4,-i*100);
        glScaled(10,rand1,10);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();

        glBindTexture(GL_TEXTURE_2D,12);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(15,-4,(-i*100)+25);
        glScaled(10,rand3,10);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();

        glBindTexture(GL_TEXTURE_2D,5);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(15,-4,(-i*100)+50);
        glScaled(10,rand2,10);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();

        glBindTexture(GL_TEXTURE_2D,4);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(15,-4,(-i*100)+75);
        glScaled(10,rand4,10);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();

        ///left buildings
        glBindTexture(GL_TEXTURE_2D,4);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(-25,-4,-i*100);
        glScaled(10,rand4,10);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();

        glBindTexture(GL_TEXTURE_2D,5);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(-25,-4,(-i*100)+25);
        glScaled(10,rand2,10);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();

        glBindTexture(GL_TEXTURE_2D,12);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(-25,-4,(-i*100)+50);
        glScaled(10,rand1,10);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();

        glBindTexture(GL_TEXTURE_2D,11);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(-25,-4,(-i*100)+75);
        glScaled(10,rand3,10);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();

    }
}

void Reklame(){
	for (int i=1; i<18; i++){
		if (i%5==0 || i==1) {
			//Papan atas nama anggota
			glBindTexture(GL_TEXTURE_2D,21);
        	glEnable(GL_TEXTURE_2D);
       		glPushMatrix();
        	glTranslatef(-10,6,(-i*210)+205);
        	glScaled(22,5,2);
        	drawcube(1,1,1,1);
        	glDisable(GL_TEXTURE_2D);
        	glPopMatrix();
		}
		
		if ((i-1)%5==0 || i==2) {
        	//Papan atas meme 1
			glBindTexture(GL_TEXTURE_2D,22);
    		glEnable(GL_TEXTURE_2D);
       		glPushMatrix();
        	glTranslatef(-10,6,(-i*210)+205);
        	glScaled(22,5,2);
        	drawcube(1,1,1,1);
        	glDisable(GL_TEXTURE_2D);
        	glPopMatrix();
    	}
    	
    	if ((i-2)%5==0 || i==2) {
        	//Papan atas meme 2
			glBindTexture(GL_TEXTURE_2D,23);
    		glEnable(GL_TEXTURE_2D);
       		glPushMatrix();
        	glTranslatef(-10,6,(-i*210)+205);
        	glScaled(22,5,2);
        	drawcube(1,1,1,1);
        	glDisable(GL_TEXTURE_2D);
        	glPopMatrix();
    	}
    	
    	if ((i-3)%5==0 || i==3) {
        	//Papan atas meme 3
			glBindTexture(GL_TEXTURE_2D,24);
    		glEnable(GL_TEXTURE_2D);
       		glPushMatrix();
        	glTranslatef(-10,6,(-i*210)+205);
        	glScaled(22,5,2);
        	drawcube(1,1,1,1);
        	glDisable(GL_TEXTURE_2D);
        	glPopMatrix();
    	}
    	
    	if ((i-4)%5==0 || i==4) {
        	//Papan atas meme 4
			glBindTexture(GL_TEXTURE_2D,25);
    		glEnable(GL_TEXTURE_2D);
       		glPushMatrix();
        	glTranslatef(-10,6,(-i*210)+205);
        	glScaled(22,5,2);
        	drawcube(1,1,1,1);
        	glDisable(GL_TEXTURE_2D);
        	glPopMatrix();
    	}
    	
        //tiang kanan
		glBindTexture(GL_TEXTURE_2D,20);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(-10,-2,(-i*210)+205);
        glScaled(2,8,2);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();
		
		//tiang kiri
		glBindTexture(GL_TEXTURE_2D,20);
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glTranslatef(10,-2,(-i*210)+205);
        glScaled(2,8,2);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();	        
	}
}

void Moon(){
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(-25+x_look,45,-150.6+zz);
    circle_3D(3.3);
    glPopMatrix();
}

void truk(void)
{
    ///truk nicher body
    glBindTexture(GL_TEXTURE_2D,7);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(-0.1-1+x_look,-2,-4.6+zz);
    glScaled(2.2,1,2.6);
    drawcube(1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    ///truk TOP
    glBindTexture(GL_TEXTURE_2D,8);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(-1+.05+x_look,-1,-4+zz);
    glScaled(1.9,.7,2);
    drawcube(1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    
    ///Red light1
    glBindTexture(GL_TEXTURE_2D,9);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(-0.8+x_look,-2,4.65+zz);
    glScaled(0.4,0.17,0.8);
    drawcube(1,1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    ///red light 2
    glBindTexture(GL_TEXTURE_2D,9);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(0.4+x_look,-2,4.65+zz);
    glScaled(0.4,0.17,0.81);
    drawcube(1,1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    
    ///yellow light1
    glBindTexture(GL_TEXTURE_2D,17);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(-0.8+x_look,-1.8,-4.65+zz);
    glScaled(0.2,0.2,0.8);
    drawcube(1,1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    ///yellow light 2
    glBindTexture(GL_TEXTURE_2D,17);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(0.6+x_look,-1.8,-4.65+zz);
    glScaled(0.2,0.2,0.81);
    drawcube(1,1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    
    ///white light1
    glPushMatrix();
    glTranslatef(-0.8+x_look,-1.4,-4.65+zz);
    glScaled(0.2,0.2,0.8);
    drawcube(1,1,1,1,1);
    glPopMatrix();

    ///white light 2
    glPushMatrix();
    glTranslatef(0.6+x_look,-1.4,-4.65+zz);
    glScaled(0.2,0.2,0.81);
    drawcube(1,1,1,1,1);
    glPopMatrix();

    ///numberplate
    glBindTexture(GL_TEXTURE_2D,21);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(-0.3+x_look,-1.8,-4.65+zz);
    glScaled(.7,.4,0.8);
    drawcube(1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    //BUMPER1
    glPushMatrix();
    glTranslatef(-1+x_look,-2,-2.0+zz);
    glScaled(2,.2,3.5);
    drawcube(1,1,1);
    glPopMatrix();
	
    //BUMPER2
    glPushMatrix();
    glTranslatef(-1+x_look,-2,2.0+zz);
    glScaled(2,.2,3.4);
    drawcube(1,1,1);
    glPopMatrix();
	
    ///Muatan
    glBindTexture(GL_TEXTURE_2D,18);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslated(-1+0.05+x_look,-2,-1.7+zz);
    glScaled(1.9,1.9,3.1);
    drawcube(1,1,1,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    ///Muatan2
    glBindTexture(GL_TEXTURE_2D,18);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslated(-1+0.05+x_look,-2,2.1+zz);
    glScaled(1.9,1.9,3.1);
    drawcube(1,1,1,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

	//Gandengan
    glBindTexture(GL_TEXTURE_2D,8);
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glTranslatef(-0.15+x_look,-2,1.28+zz);
    glScaled(0.4,0.17,0.81);
    drawcube(1,1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

	//Roda1
    ///Tyre
    glPushMatrix();
    glTranslatef(0.35-1.5+x_look,-1.8,-3.1+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    ///Tyre
    glPushMatrix();
    glTranslatef(-0.35+1.5+x_look,-1.8,-3.1+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    
	//Roda2
    ///Tyre
    glPushMatrix();
    glTranslatef(0.4-1.5+x_look,-1.8,-1.2+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    ///Tyre
    glPushMatrix();
    glTranslatef(-0.40+1.5+x_look,-1.8,-1.2+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    
	//Roda3
    ///Tyre
    glPushMatrix();
    glTranslatef(0.4-1.5+x_look,-1.8,0.8+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    ///Tyre
    glPushMatrix();
    glTranslatef(-0.4+1.5+x_look,-1.8,0.8+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    
    //Roda4
    ///Tyre
    glPushMatrix();
    glTranslatef(0.4-1.5+x_look,-1.8,2.5+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    ///Tyre
    glPushMatrix();
    glTranslatef(-0.40+1.5+x_look,-1.8,2.5+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    
	//Roda5
    ///Tyre
    glPushMatrix();
    glTranslatef(0.4-1.5+x_look,-1.8,4.5+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
    ///Tyre
    glPushMatrix();
    glTranslatef(-0.4+1.5+x_look,-1.8,4.5+zz);
    glRotatef(90,0,1,0);
    glScaled(0.5,0.5,0.5);
    Tire();
    glPopMatrix();
}

void sphere()
{
    //glTranslated(3,0,3);
    // GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { 1, 0, 1, 1.0 };
    GLfloat mat_diffuse[] = { 1, 0, 1, 1.0 };
    GLfloat mat_specular[] = { 1,1,1, 1.0 };
    GLfloat mat_shininess[] = {10};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);
    glutSolidSphere (1.0, 20, 16);
}

void cone()
{
    glRotatef(-90,1,0,0);
    glEnable(GL_TEXTURE_GEN_S); //enable texture coordinate generation
    glEnable(GL_TEXTURE_GEN_T);
    glBindTexture(GL_TEXTURE_2D,2);
    glEnable(GL_TEXTURE_2D);
    glutSolidCone(9,10,16,16);
    glDisable(GL_TEXTURE_GEN_S); //enable texture coordinate generation
    glDisable(GL_TEXTURE_GEN_T);
    //glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}


void bulbinbuilding()
{
    GLfloat mat_ambient[] = { 1, 1, 1, 1.0 };
    GLfloat mat_diffuse[] = { 1, 1, 1, 1.0 };
    GLfloat mat_specular[] = { 1,1,1, 1.0 };
    GLfloat mat_shininess[] = {10};
    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);
    glPushMatrix();
    glutSolidDodecahedron();
    glPopMatrix();
}

void culvert()
{

    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D,15);
    glEnable(GL_TEXTURE_2D);
    //glRotatef(-30,0,0,1);
    //glTranslatef(-1,0,1.5);
    glScaled(3,0.1,21);
    drawcube(1,1,1,1);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();


    for(int i=0; i<15; i++)
    {
        glPushMatrix();
        glBindTexture(GL_TEXTURE_2D,13);
        glEnable(GL_TEXTURE_2D);
        //glRotatef(-30,0,0,1);
        glTranslatef(0.2,0,i*1.5);
        glScaled(0.1,2,0.1);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();
    }
    for(int i=0; i<15; i++)
    {
        glPushMatrix();
        glBindTexture(GL_TEXTURE_2D,13);
        glEnable(GL_TEXTURE_2D);
        //glRotatef(-30,0,0,1);
        glTranslatef(2.7,0,i*1.5);
        glScaled(0.1,2,0.1);
        drawcube(1,1,1,1);
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();
    }
}


void textDisplay(string str,int x,int y,int z)
{
    GLfloat mat_ambient[] = { 1, 1, 0, 1.0 };
    GLfloat mat_diffuse[] = { 1, 1, 0, 1.0 };
    GLfloat mat_specular[] = { 1,1,1, 1.0 };
    GLfloat mat_shininess[] = {10};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1);
    glPushMatrix();
    glTranslatef(x, y,z);
    glScalef(0.003f,0.002f,1);

    for (int i=0; i<str.size(); i++)
    {
        glutStrokeCharacter(GLUT_STROKE_ROMAN, str[i]);
    }
    glPopMatrix();
}

float scoresave=0;

float rotat;
float life=3;

string str="";
string strlife="";
stringstream stringtext;
stringstream pausetext;
stringstream openingtext;
stringstream lifee;
float speed=1.0;

static void display(void)
{
    // axis();
    //string st;
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*90.0;
    glClearColor(0.1, 0.15, 0.3, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    /// GLUE LOOK AT
    glLoadIdentity();

    ///Looping in a fixed range of track and making it look like spotless.
    if (zz>-850)
    {
        if(pause)
        {
            zz = zz - speed;//1.9
        }
        else
        {
            zz=zz-0;
        }
    }
    if(zz<-850)
    {

        zz=-50;
        speed+=0.3;///increasing speed after each loop
        factor+=0.3;///increasing left right moving ratio
    }
    if(pause)
    {
        scoresave=scoresave+1;
    }
    else
    {
        scoresave=scoresave;
    }


    gluLookAt(0.0+ x_look, yy, 7.5 + zz,
              x_look, 0.0, -10.0 + zz,
              0.0, 1.0, 0.0);
    glRotatef(rot,0,1,0);

    ///showing score in sky,game over,instructions....
    
    glPushMatrix();
    stringtext.str("");
    if(life==0.0)
    {
        pause=false;
        //stringstream lifeoutput;
        string str2;
        str2="Game Over!!";
        textDisplay(str2,x_look-1,3,zz);
        //string str3;
        stringstream lifeoutput;
        lifeoutput<<scoresave;
        string str3;
        str3="Score: " + lifeoutput.str();

        textDisplay(str3,x_look-1,2,zz);
        string str4;
        str4="Press Z to Restart!!";
        textDisplay(str4,x_look-1,1,zz);
    }

    else
    {
        if(start)
        {
            string pausestr;
            pausestr="P = Mulai";
            textDisplay(pausestr,-4,4,50);
            pausestr="E = Pause";
            textDisplay(pausestr,-4,3,50);
            pausestr="A = Kiri";
            textDisplay(pausestr,-4,2,50);
            pausestr="W = NOS";
            textDisplay(pausestr,-4,0,50);
            pausestr="S = Rem";
            textDisplay(pausestr,-4,-1,50);
            pausestr="D = Kanan";
            textDisplay(pausestr,-4,1,50);
        }
        else
        {

            if(pause)
            {
                stringtext<<scoresave;
                str="Score : " + stringtext.str();
            }
            else
            {
                stringtext<<scoresave;
                str="Score : " + stringtext.str();
                strlife="LIFE Remaining"+lifee.str();
                //cout<<strlife<<endl;
                string pausestr;
                pausestr="Press 'E' to Resume..";
                textDisplay(pausestr,x_look-2,2,zz);
            }

            stringstream lifeoutput;
            lifeoutput<<life;
            string str1;
            str1="Life: " + lifeoutput.str()+"/3";
            textDisplay(str1,x_look-3.5,3,zz);
            //string str1;
            
            textDisplay(str,x_look,3,zz);
            str="";
            strlife="";
        }
        glPopMatrix();
    
        glPushMatrix();
        glTranslatef(0.8,-2,-350);
        glRotatef(107,0.0,0.0,1.0);
        glPushMatrix();
        glRotatef(90,0.0,1.0,0.0);
        glScaled(50,5,5);
        //z,

        ///tunnel code and texturing..
        glPushMatrix();
        glBindTexture(GL_TEXTURE_2D,16);
        glEnable(GL_TEXTURE_2D);
        bottleBezier();
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();
        glPopMatrix();

        glPopMatrix();
    }

    light();

    glPushMatrix();
    glTranslatef(-5.4+2+1,-1.9,-2-2-2-5-200);
    culvert();
    glPopMatrix();

    glPushMatrix();
    //cone();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,-20);
    glScaled(5,5,5);
    //pyramid();
    glPopMatrix();

    ///calling objects
    piler();
    footpath();
    Floor();
    track();
    building();
    glPushMatrix();
    glTranslatef(0,0.17,0);
    truk();
    Reklame();
    Moon();
    glPopMatrix();

    glPushMatrix();
    glLoadIdentity();
    glPopMatrix();

    ///collusion detection
    ///Texture cube
    
	///obstacle no 1
    float slide=0,look=0;
    
    slide = 75;
    look = 3;
    glPushMatrix();
    glTranslatef(-9.0+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    
    if (zz>-6.7 - slide&& zz<-2.3- slide)
    {
        if(x_look <1 || x_look > 7)
        {
            zz = zz + 60;
            life--;
        }
    }
    
    ///obstacle no 2
    slide = 135;
    look = 3;
	glPushMatrix();
    glTranslatef(-1.5+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    

    if (zz>-6.7 - slide&& zz<5- slide)
    {
        if(x_look <-10 || x_look > 1)
        {
            zz = zz + 60;
            life--;
        }
    }

    ///obstacle no 3
    slide = 200;
    look = -3;
    if (zz>-6.7 - slide&& zz<10- slide)
    {
        if(x_look <-1.4 || x_look > -.6)
        {
            zz = zz + 60;
            life--;
        }
    }

    ///obstacle no 4
    slide = 275;
    look = 3;
    glPushMatrix();
    glTranslatef(-9.0+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    
    if (zz>-6.7 - slide&& zz<-2.3- slide)
    {
        if(x_look <1 || x_look > 7)
        {
            zz = zz + 60;
            life--;
        }
    }
    
	//obstacle no 5
    slide = 325;
    look = 3;
	glPushMatrix();
    glTranslatef(-1.5+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    

    if (zz>-6.7 - slide&& zz<5- slide)
    {
        if(x_look <-10 || x_look > 1)
        {
            zz = zz + 60;
            life--;
        }
    }
    
	//obstacle no 6
    slide = 375;
    look = 3;
    glPushMatrix();
    glTranslatef(-9.0+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    
    if (zz>-6.7 - slide&& zz<-2.3- slide)
    {
        if(x_look <1 || x_look > 7)
        {
            zz = zz + 60;
            life--;
        }
    }
    
	//obstacle no 7
    slide = 425;
    look = -3;
    glPushMatrix();
    glTranslatef((-5.4+2 +look+moving)*-1,0,-5-slide);
    glBindTexture(GL_TEXTURE_2D,6);
    glEnable(GL_TEXTURE_2D);
    glRotatef(rotat,0,0,1);
    circle_3D(2.3);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    glPopMatrix();


    if (zz>-435 && zz<-421)
    {
        if(abs((x_look)-(-1*(-5.4+2+look+moving)))<=1.7)
        {
            zz = zz + 60;
            life--;
        }
    } 
    
    //obstacle no 8
    slide = 475;
    look = 3;
	glPushMatrix();
    glTranslatef(-1.5+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    

    if (zz>-6.7 - slide&& zz<5- slide)
    {
        if(x_look <-10 || x_look > 1)
        {
            zz = zz + 60;
            life--;
        }
    }
    
    //obstacle no 9
    slide = 525;
    look = 3;
    glPushMatrix();
    glTranslatef(-9.0+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    
    if (zz>-6.7 - slide&& zz<-2.3- slide)
    {
        if(x_look <1 || x_look > 7)
        {
            zz = zz + 60;
            life--;
        }
    }
    
    //obstacle no 10
    slide = 575;
    look = 3;
	glPushMatrix();
    glTranslatef(-1.5+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    

    if (zz>-6.7 - slide&& zz<5- slide)
    {
        if(x_look <-10 || x_look > 1)
        {
            zz = zz + 60;
            life--;
        }
    }
    
    //obstacle no 11
    slide = 625;
    look = 3;
    glPushMatrix();
    glTranslatef(-9.0+look,-2.3,-5-slide);
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glScaled(6.5,4,2);
    drawcube(1,1,1,2);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    
    if (zz>-6.7 - slide&& zz<-2.3- slide)
    {
        if(x_look <1 || x_look > 7)
        {
            zz = zz + 60;
            life--;
        }
    }
    
    
    ///obstacle no 12
    slide = 700;
    look = 3;
    glPushMatrix();
    glTranslatef(-7.4+look+moving,0,-5-slide);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D,15);
    glEnable(GL_TEXTURE_2D);
    glRotatef(rotat,0,0,-1);
    circle_3D(2.3);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    glPopMatrix();


    if (zz>-707 && zz<-702)
    {
        if(abs((x_look)-(-7.4+look+moving))<=1.7)
        {
            zz = zz +60;
            life--;
        }
    }

    ///obstacle no 13
    slide = 800;
    look = -3;
    glPushMatrix();
    glTranslatef((-5.4+2 +look+moving)*-1,0,-5-slide);
    glBindTexture(GL_TEXTURE_2D,6);
    glEnable(GL_TEXTURE_2D);
    glRotatef(rotat,0,0,1);
    circle_3D(2.3);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
    glPopMatrix();


    if (zz>-810 && zz<-796)
    {
        if(abs((x_look)-(-1*(-5.4+2+look+moving)))<=1.7)
        {
            zz = zz + 60;
            life--;
        }
    }


    ///show objects........
    ///obstacle no 9-1
    slide = 900;
    look = 3;
    glPushMatrix();

    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glTranslatef(-7.4+look,-2.3,-5-slide);
    glScaled(6,4,2);
    drawcube(1,1,1,2);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    ///obstacle no 10-2
    slide = 975;
    look = -3;
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D,16);
    glEnable(GL_TEXTURE_2D);
    glTranslatef(-5.4+2 + look,-2.3,-5-slide);
    glScaled(6,4,2);
    drawcube(1,1,1,2);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glutSwapBuffers();

}


static void key(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 27 :
    case 'q':
        exit(0);
        break;
        
    case'e':
    case'E':
        pause=1-pause;
        break;

    case 'w':
        zz = zz - 1.8;
        break;
        
    case 's':
        zz = zz + 0.3;
        break;

    case 'd':
        if ( x_look < 6.7)
            x_look = x_look + factor;
        else
            x_look = x_look;
        break;
        
    case 'a':
        if ( x_look > -5)
            x_look = x_look - factor;
        else
            x_look = x_look;
        break;
        
        break;
    case 'm':
    case 'M':
        //eyex
        //rotation();
        kon+=0.01;
        rot+=5;
        break;

    case 'n':
    case 'N':
        //rotation();
        kon-=0.01;
        rot-=5;
        break;
        
    case 'K':
    case 'k':
        yy=yy+0.1;

    case 'j':
    case 'J':
        yy=yy-0.1;
        
    case 'p':
    case 'P':
        start=1-start;
        pause=1-pause;
        break;
        
    case'z':
    case'Z':
        //pause=true;
        life=3;
        scoresave=0;
        zz=-10;
        start=true;
        x_look=0;
        break;

    }
    glutPostRedisplay();
}
//bool flag=true;
static void idle(void)
{

    rotat+=10;
    if(rotat > 360)
        rotat = 0;

    moving+=0.20;
    if(moving>=15)
    {
        moving=0;
        // flag=false;
    }

    //printf("%lf\n",moving);
    //std::cout << flag << std::endl;

    glutPostRedisplay();
}


int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowSize(640,480);
    glutInitWindowPosition(10,10);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    glutInitWindowPosition (100, 100);
    glutCreateWindow("Car Game");
    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutIdleFunc(idle);


    glShadeModel( GL_SMOOTH );
    glEnable( GL_DEPTH_TEST );
    glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING);
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\jalan.bmp");///1
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\rumput.bmp");///2
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\trotoar.bmp");///3
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\gedung.bmp");///4
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\gedung2.bmp");///5
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\hitam putih.bmp");///6
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\Steel.bmp");///7
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\stainless.bmp");///8
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\redlight.bmp");///9
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\batu jembatan.bmp");///10
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\gedung3.bmp");///11
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\gedung4.bmp");///12
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\batangpohon.bmp");///13
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\ramp.bmp");///14
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\bola.bmp");///15
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\tembok.bmp");///16
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\lampu.bmp");///17
    LoadTexture("D:\\semester 2\\gkv\\Tubes\\badan.bmp");///18
	LoadTexture("D:\\semester 2\\gkv\\Tubes\\batangpohon2.bmp");///19
	LoadTexture("D:\\semester 2\\gkv\\Tubes\\beton.bmp");///20
	LoadTexture("D:\\semester 2\\gkv\\Tubes\\beton1.bmp");///21
	LoadTexture("D:\\semester 2\\gkv\\Tubes\\Foto1.bmp");///22
	LoadTexture("D:\\semester 2\\gkv\\Tubes\\Foto1.bmp");///23
	LoadTexture("D:\\semester 2\\gkv\\Tubes\\Foto1.bmp");///24
	LoadTexture("D:\\semester 2\\gkv\\Tubes\\Foto1.bmp");///25
	LoadTexture("D:\\semester 2\\gkv\\Tubes\\air.bmp");///26
	
    glutMainLoop();

    return EXIT_SUCCESS;
}

