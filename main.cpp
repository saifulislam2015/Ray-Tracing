#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<bits/stdc++.h>
#include <windows.h>
#include <glut.h>
using namespace std;

int recursion_level;

#include "bitmap_image.hpp"
#include "object.hpp"
#define pi (2*acos(0.0))

struct point
{
	double x,y,z;
};

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
point pos,vec_u,vec_l,vec_r;
double a,r;

int image_width;
double VIEW_ANGLE = 80;
double window_height = 500;
double window_width = 500;


void loadTestData()
{

    Object *temp;
    vector3 Center(40,0,10),C(-30,60,20),C2(15,15,45);
    int Radius = 10;

    temp=new Sphere(Center, Radius); // Center(0,0,10), Radius 10
    temp->setColor(0,1,0);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(10);
    objects.push_back(temp);

    temp=new Sphere(C, 20); // Center(0,0,10), Radius 10
    temp->setColor(0,0,1);
    temp->setCoEfficients(0.2,0.2,0.4,0.2);
    temp->setShine(15);
    objects.push_back(temp);

    temp=new Sphere(C2, 15); // Center(0,0,10), Radius 10
    temp->setColor(1,1,0);
    temp->setCoEfficients(0.4,0.3,0.1,0.2);
    temp->setShine(5);
    objects.push_back(temp);

    vector3 a(50,30,0),b(70,60,0),c(50,45,50);
    temp = new Triangle(a,b,c);
    temp->setColor(1,0,0);
    temp->setCoEfficients(0.4,0.2,0.1,0.3);
    temp->setShine(5);
    objects.push_back(temp);

    vector3 d(70,60,0),e(30,60,0),f(50,45,50);
    temp = new Triangle(d,e,f);
    temp->setColor(1,0,0);
    temp->setCoEfficients(0.4,0.2,0.1,0.3);
    temp->setShine(5);
    objects.push_back(temp);

    vector3 g(30,60,0),h(50,30,0),i(50,45,50);
    temp = new Triangle(g,h,i);
    temp->setColor(1,0,0);
    temp->setCoEfficients(0.4,0.2,0.1,0.3);
    temp->setShine(5);
    objects.push_back(temp);

    double coeff[] = {1, 1, 1, 0, 0, 0, 0, 0, 0, -100};
    vector3 reff(0, 0, 0);

    temp = new General_Quadratic(coeff, reff, 0, 0, 20);
    temp->setColor(0, 1, 0);
    temp->setCoEfficients(0.4, 0.2, 0.1, 0.3);
    temp->setShine(10);
    objects.push_back(temp);

    double coeff2[] = {0.0625, 0.04, 0.04, 0, 0, 0, 0, 0, 0, -36};
    vector3 reff2(0, 0, 0);

    temp = new General_Quadratic(coeff2, reff2, 0, 0, 15);
    temp->setColor(1, 0, 0);
    temp->setCoEfficients(0.4, 0.2, 0.1, 0.3);
    temp->setShine(15);
    objects.push_back(temp);

    vector3 light1(-70,70,70);
    lights.push_back(light1);
    vector3 light2(70,70,70);
    lights.push_back(light2);



    temp=new Floor(1000, 20);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);



    image_width = 768;
    recursion_level = 3;
}

void loadActualData()
{
    ifstream input;
    string ob;
    Object *temp;
    int n,l;

    input.open("scene.txt");

    input >> recursion_level;
    input >> image_width;
    input >> n;

    for(int i=0;i<n;i++){
        input >> ob;

        if(ob=="sphere"){
            double a,b,c,r,col[3],am,dif,spec,refl,coff,s;

            input >> a >> b >> c;
            input >> r;
            input >> col[0] >> col[1] >> col[2];
            input >> am >> dif >> spec >> refl;
            input >> s;

            vector3 C(a,b,c);
            temp = new Sphere(C,r);
            temp->setColor(col[0],col[1],col[2]);
            temp->setCoEfficients(am,dif,spec,refl);
            temp->setShine(s);
            objects.push_back(temp);

        }

        if(ob=="triangle"){
            double col[3],am,dif,spec,refl,coff,s,x[3],y[3],z[3];
            for(int i=0;i<3;i++){
                input>> x[i] >> y[i] >> z[i];
                cout << x[i] << y[i] << z[i] << endl;
            }
            input >> col[0] >> col[1] >> col[2];
            input >> am >> dif >> spec >> refl;
            input >> s;

            vector3 a(x[0],y[0],z[0]);
            vector3 b(x[1],y[1],z[1]);
            vector3 c(x[2],y[2],z[2]);

            temp = new Triangle(a,b,c);
            temp->setColor(col[0],col[1],col[2]);
            temp->setCoEfficients(am,dif,spec,refl);
            temp->setShine(s);
            objects.push_back(temp);
        }
        if(ob=="general"){
            double col[3],am,dif,spec,refl,s,coff[10];
            double x,y,z,l,w,h;
            for(int i=0;i<10;i++){
                input>> coff[i] ;
            }
            input >> x >> y >> z >> l >> w >> h;
            vector3 p(x,y,z);
            input >> col[0] >> col[1] >> col[2];
            input >> am >> dif >> spec >> refl;
            input >> s;

            temp = new General_Quadratic(coff,p,l,w,h);
            temp->setColor(col[0],col[1],col[2]);
            temp->setCoEfficients(am,dif,spec,refl);
            temp->setShine(s);
            objects.push_back(temp);
        }
    }

    input >> l;

    for(int i=0;i<l;i++){
        double x,y,z;
        input >> x >> y >> z;
        vector3 lght(x,y,z);
        lights.push_back(lght);
    }
    temp=new Floor(1000, 20);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);
    input.close();

}

void Capture()
{

    vector3** frameBuffer;
    frameBuffer = new vector3* [image_width];
    vector3 black(0,0,0);
    double colorAt[3];
    bitmap_image image(image_width,image_width);

    for(int i=0;i<image_width;i++){
        frameBuffer[i] = new vector3 [image_width];
        for(int j=0;j<image_width;j++){
            frameBuffer[i][j] = black;
        }
    }


    double plane_distance= (window_height/2)/tan(VIEW_ANGLE*pi/360);
    vector3 eye(pos.x,pos.y,pos.z);
    point tl;

    tl.x = pos.x - vec_l.x*plane_distance - vec_r.x*(window_width/2) + vec_u.x*(window_height/2);
    tl.y = pos.y - vec_l.y*plane_distance - vec_r.y*(window_width/2) + vec_u.y*(window_height/2);
    tl.z = pos.z - vec_l.z*plane_distance - vec_r.z*(window_width/2) + vec_u.z*(window_height/2);

    /*tl.x = pos.x + vec_l.x*plane_distance - vec_r.x*(window_width/2) + vec_u.x*(window_height/2);
    tl.y = pos.y + vec_l.y*plane_distance - vec_r.y*(window_width/2) + vec_u.y*(window_height/2);
    tl.z = pos.z + vec_l.z*plane_distance - vec_r.z*(window_width/2) + vec_u.z*(window_height/2);*/

    vector3 topleft(tl.x,tl.y,tl.z);

    double du = window_width/image_width;
    double dv = window_height/image_width;


    for(int i=0;i<image_width;i++){
        for(int j=0;j<image_width;j++){
            double cx = topleft.x + vec_r.x*j*du - vec_u.x*i*dv;
            double cy = topleft.y + vec_r.y*j*du - vec_u.y*i*dv;
            double cz = topleft.z + vec_r.z*j*du - vec_u.z*i*dv;

            vector3 corner(cx,cy,cz);
            vector3 vl = eye.subtract(corner);
            //vector3 vl = corner.subtract(eye);
            vl.normalize();
            Ray ray(eye,vl);

            int nearest = -1;
            double t_min = DBL_MAX;
            for(int k=0;k<objects.size();k++){
                double t = objects[k]->intersect(&ray,colorAt,0);

                if(t<0) continue;
                else if(t<t_min){
                    t_min = t;
                    nearest = k;
                }
            }

            if(nearest!=-1){
                double t = objects[nearest]->intersect(&ray,colorAt,1);
                //vector3 cl(colorAt[0],colorAt[1],colorAt[2]);
                //image.set_pixel(j, i, colorAt[0]*255, colorAt[1]*255, colorAt[2]*255);
                frameBuffer[i][j].x = colorAt[0];
                frameBuffer[i][j].y = colorAt[1];
                frameBuffer[i][j].z = colorAt[2];

                //cout << frameBuffer[i][j].x << " " << frameBuffer[i][j].y << " " << frameBuffer[i][j].z<< endl;
            }

        }
    }


    cout<<"A"<<endl;

    for (int i=0; i<image_width; i++) {
        for (int j=0; j<image_width; j++) {
            image.set_pixel(image_width-j, image_width-i, frameBuffer[i][j].x*255, frameBuffer[i][j].y*255, frameBuffer[i][j].z*255);
            //image.set_pixel(j, i, frameBuffer[i][j].x*255, frameBuffer[i][j].y*255, frameBuffer[i][j].z*255);
            //image.set_pixel(j, i, colorAt[0]*255, colorAt[1]*255, colorAt[2]*255);
            //image.set_pixel(j, i, 1, 0, 0);
        }
    }

    image.save_image("b.bmp");
    cout << "saved" << endl;
}

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    //glColor3f(0.0,1.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,a);
		glVertex3f( a,-a,a);
		glVertex3f(-a,-a,a);
		glVertex3f(-a, a,a);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawCylinder(double radius, double height, int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*(pi/2));
        points[i].y=radius*sin(((double)i/(double)segments)*(pi/2));
    }

    //draw segments using generated points
    for(i=0;i<segments;i++)
    {

        glBegin(GL_QUADS);
        {
			glVertex3f(points[i].x,points[i].y,-height/2);
			glVertex3f(points[i+1].x,points[i+1].y,-height/2);
			glVertex3f(points[i+1].x,points[i+1].y,height/2);
			glVertex3f(points[i].x,points[i].y,height/2);
        }
        glEnd();
    }
}




void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*pi/2);
			points[i][j].y=r*sin(((double)j/(double)slices)*pi/2);
			points[i][j].z=h>0?h:-h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        glColor3f(1,0,0);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                /*glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);*/
			}glEnd();
		}
	}
}


void drawSS()
{
    drawAxes();

}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        case '0':
            {
                Capture();
                break;
            }

		case '1':
			{
			    point px = vec_l,py = vec_r;

                vec_l.x = px.x*cos(3*pi/180) + py.x*sin(3*pi/180);
                vec_l.y = px.y*cos(3*pi/180) + py.y*sin(3*pi/180);
                double r = sqrt(pow(vec_l.x,2)+pow(vec_l.y,2)+pow(vec_l.z,2));
                vec_l.x/=r;
                vec_l.y/=r;
                vec_l.z/=r;

                vec_r.x = py.x*cos(3*pi/180) - px.x*sin(3*pi/180);
                vec_r.y = py.y*cos(3*pi/180) - px.y*sin(3*pi/180);
                double r2 = sqrt(pow(vec_r.x,2)+pow(vec_r.y,2)+pow(vec_r.z,2));
                vec_r.x/=r2;
                vec_r.y/=r2;
                vec_r.z/=r2;

                break;
			}


        case '2':
            {
                point px = vec_l,py = vec_r;

                vec_l.x = px.x*cos(-3*pi/180) + py.x*sin(-3*pi/180);
                vec_l.y = px.y*cos(-3*pi/180) + py.y*sin(-3*pi/180);
                double r = sqrt(pow(vec_l.x,2)+pow(vec_l.y,2)+pow(vec_l.z,2));
                vec_l.x/=r;
                vec_l.y/=r;
                vec_l.z/=r;

                vec_r.x = py.x*cos(-3*pi/180) - px.x*sin(-3*pi/180);
                vec_r.y = py.y*cos(-3*pi/180) - px.y*sin(-3*pi/180);
                double r2 = sqrt(pow(vec_r.x,2)+pow(vec_r.y,2)+pow(vec_r.z,2));
                vec_r.x/=r2;
                vec_r.y/=r2;
                vec_r.z/=r2;

                break;
            }


        case '3':
            {
                point px = vec_l,py = vec_u;

                vec_l.x = px.x*cos(3*pi/180) + py.x*sin(3*pi/180);
                vec_l.y = px.y*cos(3*pi/180) + py.y*sin(3*pi/180);
                vec_l.z = px.z*cos(3*pi/180) + py.z*sin(3*pi/180);
                double r = sqrt(pow(vec_l.x,2)+pow(vec_l.y,2)+pow(vec_l.z,2));
                vec_l.x/=r;
                vec_l.y/=r;
                vec_l.z/=r;

                vec_u.x = py.x*cos(3*pi/180) - px.x*sin(3*pi/180);
                vec_u.y = py.y*cos(3*pi/180) - px.y*sin(3*pi/180);
                vec_u.z = py.z*cos(3*pi/180) - px.z*sin(3*pi/180);
                double r2 = sqrt(pow(vec_u.x,2)+pow(vec_u.y,2)+pow(vec_u.z,2));
                vec_u.x/=r2;
                vec_u.y/=r2;
                vec_u.z/=r2;

                break;
            }

        case '4':
            {
                point px = vec_l,py = vec_u;

                vec_l.x = px.x*cos(-3*pi/180) + py.x*sin(-3*pi/180);
                vec_l.y = px.y*cos(-3*pi/180) + py.y*sin(-3*pi/180);
                vec_l.z = px.z*cos(-3*pi/180) + py.z*sin(-3*pi/180);
                double r = sqrt(pow(vec_l.x,2)+pow(vec_l.y,2)+pow(vec_l.z,2));
                vec_l.x/=r;
                vec_l.y/=r;
                vec_l.z/=r;

                vec_u.x = py.x*cos(-3*pi/180) - px.x*sin(-3*pi/180);
                vec_u.y = py.y*cos(-3*pi/180) - px.y*sin(-3*pi/180);
                vec_u.z = py.z*cos(-3*pi/180) - px.z*sin(-3*pi/180);
                double r2 = sqrt(pow(vec_u.x,2)+pow(vec_u.y,2)+pow(vec_u.z,2));
                vec_u.x/=r2;
                vec_u.y/=r2;
                vec_u.z/=r2;
                break;
            }

        case '5':
            {
                point px = vec_u,py = vec_r;

                vec_u.x = px.x*cos(3*pi/180) + py.x*sin(3*pi/180);
                vec_u.y = px.y*cos(3*pi/180) + py.y*sin(3*pi/180);
                vec_u.z = px.z*cos(3*pi/180) + py.z*sin(3*pi/180);
                double r = sqrt(pow(vec_u.x,2)+pow(vec_u.y,2)+pow(vec_u.z,2));
                vec_u.x/=r;
                vec_u.y/=r;
                vec_u.z/=r;

                vec_r.x = py.x*cos(3*pi/180) - px.x*sin(3*pi/180);
                vec_r.y = py.y*cos(3*pi/180) - px.y*sin(3*pi/180);
                vec_r.z = py.z*cos(3*pi/180) - px.z*sin(3*pi/180);
                double r2 = sqrt(pow(vec_r.x,2)+pow(vec_r.y,2)+pow(vec_r.z,2));
                vec_r.x/=r2;
                vec_r.y/=r2;
                vec_r.z/=r2;

                break;
            }

        case '6':
            {
                point px = vec_u,py = vec_r;

                vec_u.x = px.x*cos(-3*pi/180) + py.x*sin(-3*pi/180);
                vec_u.y = px.y*cos(-3*pi/180) + py.y*sin(-3*pi/180);
                vec_u.z = px.z*cos(-3*pi/180) + py.z*sin(-3*pi/180);
                double r = sqrt(pow(vec_u.x,2)+pow(vec_u.y,2)+pow(vec_u.z,2));
                vec_u.x/=r;
                vec_u.y/=r;
                vec_u.z/=r;

                vec_r.x = py.x*cos(-3*pi/180) - px.x*sin(-3*pi/180);
                vec_r.y = py.y*cos(-3*pi/180) - px.y*sin(-3*pi/180);
                vec_r.z = py.z*cos(-3*pi/180) - px.z*sin(-3*pi/180);
                double r2 = sqrt(pow(vec_r.x,2)+pow(vec_r.y,2)+pow(vec_r.z,2));
                vec_r.x/=r2;
                vec_r.y/=r2;
                vec_r.z/=r2;

                break;
            }

	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraHeight -= 3.0;
			pos.x -= vec_l.x*2;
			pos.y -= vec_l.y*2;
			pos.z -= vec_l.z*2;
			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
			pos.x += vec_l.x*2;
			pos.y += vec_l.y*2;
			pos.z += vec_l.z*2;
			break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;
			pos.x += vec_r.x*2;
			pos.y += vec_r.y*2;
			pos.z += vec_r.z*2;
			break;
		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;
			pos.x -= vec_r.x*2;
			pos.y -= vec_r.y*2;
			pos.z -= vec_r.z*2;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x += vec_u.x*2;
			pos.y += vec_u.y*2;
			pos.z += vec_u.z*2;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos.x -= vec_u.x*2;
			pos.y -= vec_u.y*2;
			pos.z -= vec_u.z*2;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
		    r++;
			break;
		case GLUT_KEY_END:
		    r--;
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
	gluLookAt(pos.x,pos.y,pos.z, pos.x + vec_l.x, pos.y + vec_l.y, pos.z + vec_l.z, vec_u.x,vec_u.y,vec_u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    drawSS();

    for(int i=0;i<objects.size();i++){
        objects[i]->draw();
    }

    for(int i=0;i<lights.size();i++){
        glColor3f(1.0, 1.0, 0);
        glBegin(GL_POINTS);
        {
            glVertex3f(lights[i].x, lights[i].y, lights[i].z);
        }
        glEnd();
    }
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

	pos.x = 100;
	pos.y = 100;
	pos.z = 0;

	vec_u.x = 0;
	vec_u.y = 0;
	vec_u.z = 1;

	vec_l.x = -1/sqrt(2);
	vec_l.y = -1/sqrt(2);
	vec_l.z = 0;

	vec_r.x = -1/sqrt(2);
	vec_r.y = 1/sqrt(2);
	vec_r.z = 0;

	a = 50;
	r = 10;

	glClearColor(0,0,0,0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(VIEW_ANGLE,	1,	1,	1000.0);
}




int main(int argc, char **argv){

    //loadTestData();
    loadActualData();
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
