// standard includes
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <utility>

#include "VoronoiCurve.h"

using namespace std;
using namespace Voronoicrv;


#define NENDS 2           /* number of end "points" to draw */

GLdouble width, height;   /* window width and height */
int wd;                   /* GLUT window handle */
int ends1[NENDS][2];       /* array of 2D points */
int minx=999,miny=999,maxx=-999,maxy=-999;

vector<pair<double, double> > pointVec;
vector<pair<int,int> > *boundary;

void init(void)
{
    width  = 1280.0;                 /* initial window width and height, */
    height = 800.0;                  /* within which we draw. */
    ends1[0][0] = (int)(0.25*width);  /* (0,0) is the lower left corner */
    ends1[0][1] = (int)(0.75*height);
    ends1[1][0] = (int)(0.75*width);
    ends1[1][1] = (int)(0.25*height);
    
    return;
}
void drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius)
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    int i,triangleAmount = 20;
    GLfloat twicePi = 2.0f * 3.14159265;
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y);
    //  glColor3f(1,0,0);
    for(i = 0; i <= triangleAmount;i++)
        glVertex2f(x + (radius * cos(i *  twicePi / triangleAmount)),y + (radius * sin(i * twicePi / triangleAmount)));
    glEnd();
}


void pointset(void)
{
    glPointSize(6.0);
    glColor3f(0.0, 0.0, 0.0);
    for(int i=0;i<pointVec.size();i++)
        drawFilledCircle(pointVec[i].first,pointVec[i].second,3.0);
    glColor3f(1.0, 0.2, 0.2);
    glLineWidth(3.0);
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );
    glColor3f(0.18, 0.69, 0.25);
    
    
for(auto & val : *boundary)
    {
    std::pair<int, int> indices=val;
        
            glBegin(GL_LINES);
            glVertex2f(pointVec[indices.first].first,pointVec[indices.first].second);
            glVertex2f(pointVec[indices.second].first,pointVec[indices.second].second);
            glEnd();
            
        
    }
    
}
void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    pointset();
    glFlush();
    return;
}

/* Called when window is resized,
 also when window is first created,
 before the first call to display(). */
void reshape(int w, int h)
{
    /* save new screen dimensions */
    width = (GLdouble) w;
    height = (GLdouble) h;
    
    /* tell OpenGL to use the whole window for drawing */
    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    
    /* do an orthographic parallel projection with the coordinate
     system set to first quadrant, limited by screen/window size */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(minx-20.0,maxx+20.0,miny-20.0,maxy+20.0, -1.f, 1.f);
    return;
}

void kbd(unsigned char key, int x, int y)
{
    switch((char)key) {
        case 'q':
        case 27:    /* ESC */
            glutDestroyWindow(wd);
            exit(0);
        default:
            break;
    }
    
    return;
}

void loadPointSet(string name)
	{
		pointVec.clear();
		ifstream file(name.data());

		if (file)
		{
			while (!file.eof())
			{
				float x, y;
				file >> x >> y;
	if(x>maxx)
            maxx=x;
        if(x<minx)
            minx=x;
        if(y>maxy)
            maxy=y;
        if(y<miny)
	miny=y;


				if (!file.eof())
					pointVec.push_back(pair<double, double>(x, y));
			}
		cout<<pointVec.size();
		}
		else
		{
			cerr << "ERROR: input file " << name << " could not be read." << endl;
			exit(2);
		}

		file.close();

		if (pointVec.size() < 3)
		{
			cerr << "ERROR: input file " << name << " contains less than 3 points." << endl;
			exit(3);
		}
	}

int main(int argc, char **argv)
{
	string filename;
		
	if(argc>=2){
		filename=string(argv[1]);				
	}
	else{
		std::cout<< "Usage" << argv[0] << " filename" << std::endl;
		return 0;
	}
	
    init();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(250, 250);
    wd = glutCreateWindow("Voronoi based Curve Reconstruction");
	
	
    loadPointSet(filename);	
    VoronoiCurve voronoiInstance(pointVec);
    boundary=voronoiInstance.getBoundary();
	
	
    glutReshapeFunc(reshape);
    glutKeyboardFunc(kbd);
    glutDisplayFunc(display);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(3.0);
    glutMainLoop();
    exit(0);
    return 0;
}
