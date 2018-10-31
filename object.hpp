#ifndef OBJECT_HPP_INCLUDED
#define OBJECT_HPP_INCLUDED

#define pi (2*acos(0.0))
#define AMBIENT 0
#define DIFFUSE 1
#define SPECULAR 2
#define REFLECTION 3

class vector3{
    public:
        double x,y,z;

    vector3(){}
    vector3(double x,double y,double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    void normalize()
    {
        double r = sqrt(x*x+y*y+z*z);

        x/=r;
        y/=r;
        z/=r;
    }

    vector3 subtract(vector3 b)
    {

        vector3 result(x - b.x,y - b.y,z - b.z);

        return result;
    }

};


class Ray{
    public:
        vector3 start,dir;

        Ray(vector3 s,vector3 d){
            this->start = s;
            this->dir = d;
            dir.normalize();
            //cout << dir.x << " " << dir.y << " " << dir.z << endl;
        }

};


class Object{
    public:
        vector3 reference_point;
        double height, width, length;
        double source_factor = 1;
        double eta = 1.5;
        int Shine;
        double color[3];
        double co_efficients[4];
        Object(){ }
        virtual void draw(){}
        void setColor(double r,double g,double b)
        {
            color[0] = r;
            color[1] = g;
            color[2] = b;
        }
        void setShine(int s)
        {
            Shine = s;
        }
        void setCoEfficients(double a,double b,double c,double d){
            co_efficients[0] = a;
            co_efficients[1] = b;
            co_efficients[2] = c;
            co_efficients[3] = d;
        }

        virtual double intersect(Ray *r,double *arr,int level){
            return -1;
        }
        virtual double getIntersectingT(Ray *r){}
        virtual vector3 getNormal(vector3 intersection){}

        double dot(vector3 a,vector3 b){
            return a.x*b.x + a.y*b.y + a.z*b.z ;
        }

        vector3 inverse(vector3 a){
            vector3 res(-a.x,-a.y,-a.z);
            return res;
        }

        vector3 getReflection(Ray* ray, vector3 normal) {
            if(dot(ray->dir, normal)>0) normal=inverse(normal);
            double m = 2*dot(ray->dir,normal);

            double x = ray->dir.x - normal.x*m;
            double y = ray->dir.y - normal.y*m;
            double z = ray->dir.z - normal.z*m;

            vector3 r(x,y,z);
            r.normalize();
            return r;
        }

        vector3 getRefraction(Ray* ray, vector3 normal) {

            double cos = dot(normal,ray->dir);
            double sin = 1.0 - eta*eta*(1-pow(cos,2));

            if (sin >= 0) {
                vector3 refraction = subtract(scalarmultiply(eta,ray->dir),scalarmultiply((eta*cos + sqrt(sin)),normal));
                refraction.normalize();
                return refraction;
            }
            else{
                vector3 r(0,0,0);
                return r;
            }
        }

        vector3 scalarmultiply(double m,vector3 a){
            vector3 result(m*a.x,m*a.y,m*a.z);
            return result;
        }

        vector3 sum(vector3 a,vector3 b){
            vector3 sum(a.x+b.x,a.y+b.y,a.z+b.z);
            return sum;
        }

        vector3 subtract(vector3 a,vector3 b){
            vector3 result(a.x - b.x,a.y - b.y,a.z - b.z);
            return result;
        }

        vector3 cross(vector3 a,vector3 b){
            vector3 res(a.y*b.z - a.z*b.y,a.z*b.x - a.x*b.z,a.x*b.y - a.y*b.x);
            return res;
        }

};

vector <Object*> objects;
vector <vector3> lights;


class Sphere:public Object{
    public:
        Sphere(vector3 Center,double Radius){
            reference_point=Center;
            length=Radius;
        }

        vector3 getNormal(vector3 intersect)
        {
            vector3 vec(intersect.x-reference_point.x,intersect.y-reference_point.y,intersect.z-reference_point.z);
            vec.normalize();
            return vec;
        }



        void draw(){
            glColor3f(color[0], color[1], color[2]);

        vector3 points[100][100];

        double h, r;

        int slices = 24, stacks = 20;

        //generate points
        for(int i=0; i<=stacks; i++) {
            h = length * sin(((double)i/(double)stacks)*(pi/2));
            r = length * cos(((double)i/(double)stacks)*(pi/2));
            for(int j=0; j<=slices; j++) {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }

        //draw quads using generated points
        for(int i=0; i<stacks; i++) {
            for(int j=0; j<slices; j++) {
                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    glVertex3f(points[i][j].x+reference_point.x,points[i][j].y+reference_point.y,points[i][j].z+reference_point.z);
                    glVertex3f(points[i][j+1].x+reference_point.x,points[i][j+1].y+reference_point.y,points[i][j+1].z+reference_point.z);
                    glVertex3f(points[i+1][j+1].x+reference_point.x,points[i+1][j+1].y+reference_point.y,points[i+1][j+1].z+reference_point.z);
                    glVertex3f(points[i+1][j].x+reference_point.x,points[i+1][j].y+reference_point.y,points[i+1][j].z+reference_point.z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x+reference_point.x,points[i][j].y+reference_point.y,-points[i][j].z+reference_point.z);
                    glVertex3f(points[i][j+1].x+reference_point.x,points[i][j+1].y+reference_point.y,-points[i][j+1].z+reference_point.z);
                    glVertex3f(points[i+1][j+1].x+reference_point.x,points[i+1][j+1].y+reference_point.y,-points[i+1][j+1].z+reference_point.z);
                    glVertex3f(points[i+1][j].x+reference_point.x,points[i+1][j].y+reference_point.y,-points[i+1][j].z+reference_point.z);
                }
                glEnd();
            }
        }
        }

        double getIntersectingT(Ray *r)
        {
            vector3 R0 = r->start.subtract(reference_point);


            /*vector3 p(0,100,10);
            vector3 d(0,1,0);
            Ray custom(p,d);
            vector3 R0 = custom.start.subtract(reference_point);*/

            double A = dot(r->dir,r->dir);
            double B = 2*dot(r->dir,R0);
            //double A = dot(custom.dir,custom.dir);
            //double B = 2*dot(custom.dir,R0);
            double C = dot(R0,R0) - length*length;

            double D = B*B - 4*A*C;

            if(D<0) return -1;

            double t1 = (-B+sqrt(D))/2*A;
            double t2 = (-B-sqrt(D))/2*A;
            //cout << t1 << "  " << t2 << endl;

            double t = t1>t2?t2:t1;
            return t;
        }

        double intersect(Ray *r, double *current_color, int level){
            double t= getIntersectingT(r);
            if(t<=0) return -1;
            if(level==0) return t;

            vector3 intersectionPoint = sum(r->start,scalarmultiply(t,r->dir));
            //colorAt = getColorAt(intersectionPoint)
            for (int i=0; i<3; i++) {
                current_color[i] = color[i]*co_efficients[AMBIENT];
            }

            vector3 normal = getNormal(intersectionPoint);
            vector3 reflection = getReflection(r, normal);
            vector3 refraction = getRefraction(r,normal);

            for(int i=0;i<lights.size();i++){
                    vector3 direction = subtract(lights[i],intersectionPoint);
                    double r1 = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);
                    direction.normalize();

                    vector3 start = sum(intersectionPoint,scalarmultiply(1,direction));

                    Ray L(start, direction);
                    double arr[3];
                    int nearest = -1;
                    double t_min = DBL_MAX;
                    bool d = false;

                    for(int j=0;j<objects.size();j++){
                        double t = objects[j]->intersect(&L,arr,0);
                        if(t>0 && t<=r1){
                            d = true;
                            break;
                        }
                    }
                    if(!d){
                        double lambert = dot(L.dir,normal);
                        double phong = pow(dot(reflection,r->dir),Shine);

                        if(lambert<0) lambert = 0;
                        if(phong<0) phong = 0;

                        for (int k=0; k<3; k++) {
                            current_color[k] += source_factor*lambert*co_efficients[DIFFUSE]*color[k];
                            current_color[k] += source_factor*phong*co_efficients[SPECULAR]*color[k];
                        }
                    }

                    if(level<recursion_level){
                        start = sum(intersectionPoint,scalarmultiply(1,reflection));
                        Ray reflectionRay(start, reflection);

                        double reflected_color[3];
                        int nearest = -1;
                        double t_min = DBL_MAX;

                        for(int k=0;k<objects.size();k++){
                            double t = objects[k]->intersect(&reflectionRay,reflected_color,0);
                            if(t<=0) continue;
                            else if(t<t_min){
                                t_min = t;
                                nearest = k;
                            }
                        }
                        if(nearest!=-1){
                            double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                            for (int k=0; k<3; k++) {
                                current_color[k] += reflected_color[k]*co_efficients[REFLECTION];
                            }
                        }
                        start = sum(intersectionPoint,scalarmultiply(1,reflection));
                        Ray refractionRay(start, refraction);

                        double refracted_color[3];
                        nearest = -1;
                        t_min = DBL_MAX;

                        for(int k=0;k<objects.size();k++){
                            double t = objects[k]->intersect(&refractionRay,refracted_color,0);
                            if(t<0) continue;
                            else if(t<t_min){
                                t_min = t;
                                nearest = k;
                            }
                        }
                        if(nearest!=-1){
                            double t = objects[nearest]->intersect(&refractionRay,refracted_color,level+1);
                            for (int k=0; k<3; k++) {
                                current_color[k] += refracted_color[k]*eta;
                            }
                        }
                    }

            }

            for(int i=0;i<3;i++){
                if(current_color[i]<0) current_color[i] = 0;
                if(current_color[i]>1) current_color[i] = 1;
            }
            return t;
    }
};


class Floor:public Object{
    public:
        bitmap_image texture;
        double h,w;

        Floor(double FloorWidth,double TileWidth){
            reference_point=vector3(-FloorWidth/2, -FloorWidth/2,0);
            //width = FloorWidth;
            length=TileWidth;
            readImage();
        }

        void readImage(){
            texture = bitmap_image("Texture.bmp");
            h = texture.height()/1000;
            w = texture.width()/1000;
        }

        void draw(){
            int tiles = abs(2*reference_point.x/length);

            for(int i = 0 ;i<tiles;i++){
                for (int j = 0;j<tiles;j++){
                    if ((i+j)%2==0) glColor3f(1.0f,1.0f,1.0f); //white
                    else glColor3f(0.0f,0.0f,0.0f); //black

                    glBegin(GL_QUADS);{
                        glVertex3f(reference_point.x+(i*length),reference_point.y+(j*length),reference_point.z);
                        glVertex3f(reference_point.x+((i+1)*length),reference_point.y+(j*length),reference_point.z);
                        glVertex3f(reference_point.x+((i+1)*length),reference_point.y+((j+1)*length),reference_point.z);
                        glVertex3f(reference_point.x+(i*length),reference_point.y+((j+1)*length),reference_point.z);
                    }glEnd();

                }
            }
        }

        double getIntersectingT(Ray* ray) {
            vector3 normal(0,0,1);
            //double t = (-1*dot(normal,ray->start))/dot(normal, ray->dir);
            double t = (reference_point.z-ray->start.z) /ray->dir.z;
            return t;
        }

        double intersect(Ray* ray, double current_color[3], int level) {

            double t = getIntersectingT(ray);


            vector3 intersectionPoint = sum(ray->start,scalarmultiply(t,ray->dir));

            double minx = reference_point.x;
            double miny = reference_point.y;



            if (minx>intersectionPoint.x || intersectionPoint.x >-minx|| miny>intersectionPoint.y || intersectionPoint.y>-miny) {
                return -1;
            }

            if(level==0) return t;


            int r = (intersectionPoint.x-reference_point.x) / length;
            int c = (intersectionPoint.y-reference_point.y) / length;

            if ((r+c)%2) {
                color[0] = color[1] = color[2] = 0;
            }
            else {
                color[0] = color[1] = color[2] = 1;
            }

            int x = (intersectionPoint.x + abs(reference_point.x))*w;
            int y = (intersectionPoint.y + abs(reference_point.y))*h;
            unsigned char c1,c2,c3;
            double tc_[3];

            texture.get_pixel(x,y,c1,c2,c3);
            //cout << c1 << " " << c2 << " " << c3 << endl;

            tc_[0] = c1;
            tc_[1] = c2;
            tc_[2] = c3;
            //cout << tc_[0] << " " << tc_[1] << " " << tc_[2] << endl;


            for (int i=0; i<3; i++) {
                current_color[i] = color[i]*co_efficients[AMBIENT]*tc_[i]/255;
            }
            //cout << current_color[0] << endl;

            vector3 normal(0,0,1);
            //vector3 normal = getNormal(intersectionPoint);

            vector3 reflection = getReflection(ray, normal);
            //return t;

            for(int i=0;i<lights.size();i++){
                    vector3 direction = subtract(lights[i],intersectionPoint);
                    double r = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);
                    direction.normalize();

                    vector3 start = sum(intersectionPoint,scalarmultiply(1,direction));

                    Ray L(start, direction);
                    double arr[3];
                    int nearest = -1;
                    bool d = false;

                    for(int j=0;j<objects.size();j++){
                        double t = objects[j]->intersect(&L,arr,0);
                        if(t>0 && t<r){
                            d = true;
                            break;
                        }
                    }
                    if(!d){
                        double lambert = dot(L.dir,normal);
                        double phong = pow(dot(reflection,ray->dir),Shine);

                        if(lambert<0) lambert = 0;
                        if(phong<0) phong = 0;

                        for (int k=0; k<3; k++) {
                            current_color[k] += source_factor*lambert*co_efficients[DIFFUSE]*color[k];
                            current_color[k] += source_factor*phong*co_efficients[SPECULAR]*color[k];
                        }
                    }

                    if(level<recursion_level){
                        start = sum(intersectionPoint,scalarmultiply(1,reflection));
                        Ray reflectionRay(start, reflection);

                        double reflected_color[3];
                        int nearest = -1;
                        double t_min = DBL_MAX;

                        for(int k=0;k<objects.size();k++){
                            //double t = objects[k]->intersect(&reflectionRay,reflected_color,0);
                            double t = objects[k]->getIntersectingT(&reflectionRay);
                            if(t<=0) continue;
                            else if(t<t_min){
                                t_min = t;
                                nearest = k;
                            }
                        }
                        if(nearest!=-1){
                            double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                            for (int k=0; k<3; k++) {
                                current_color[k] += reflected_color[k]*co_efficients[REFLECTION];
                            }
                        }
                    }
            }
            for(int i=0;i<3;i++){
                if(current_color[i]<0) current_color[i] = 0;
                if(current_color[i]>1) current_color[i] = 1;
            }


            return t;
    }


};


class Triangle:public Object{
    public:
        vector3 a,b,c;

        Triangle(vector3 A,vector3 b,vector3 c){
            this->a = A;
            this->b = b;
            this->c = c;
        }


        void draw(){
            glColor3f(color[0],color[1],color[2]);

            glBegin(GL_TRIANGLES);{
                glVertex3f(a.x,a.y,a.z);
                glVertex3f(b.x,b.y,b.z);
                glVertex3f(c.x,c.y,c.z);

            }glEnd();
        }

        vector3 getNormal(){
            return cross(subtract(b,a),subtract(c,a));
        }

        double getIntersectingT(Ray* ray) {
            const double EPSILON = 0.0000001;
            vector3 edge1 = subtract(b,a);
            vector3 edge2 = subtract(c,a);

            vector3 h = cross(ray->dir,edge2);
            double a_ = dot(edge1,h);

            if (a_ > -EPSILON && a_ < EPSILON) return -1;

            double f = 1/a_;

            vector3 s = subtract(ray->start,a);
            double u = f *dot(s,h);
            if (u < 0.0 || u > 1.0) return -1;


            vector3 q = cross(s,edge1);
            double v = f * dot(ray->dir,q);
            if (v < 0.0 || u + v > 1.0) return -1;

            double t = f * dot(edge2,q);
            if(t>EPSILON) return t;

            return -1;
        }

        double intersect(Ray *r, double *current_color, int level){
            double t= getIntersectingT(r);
            if(t<=0) return -1;
            if(level==0) return t;

            vector3 intersectionPoint = sum(r->start,scalarmultiply(t,r->dir));
            //colorAt = getColorAt(intersectionPoint)
            for (int i=0; i<3; i++) {
                current_color[i] = color[i]*co_efficients[AMBIENT];
            }

            vector3 normal = getNormal();
            vector3 reflection = getReflection(r, normal);

            for(int i=0;i<lights.size();i++){
                    vector3 direction = subtract(lights[i],intersectionPoint);
                    double r1 = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);
                    direction.normalize();

                    vector3 start = sum(intersectionPoint,scalarmultiply(1,direction));

                    Ray L(start, direction);
                    double arr[3];
                    int nearest = -1;
                    //double t_min = DBL_MAX;
                    bool d = false;

                    for(int j=0;j<objects.size();j++){
                        //double t = objects[j]->intersect(&L,arr,0);
                        double t = objects[j]->getIntersectingT(&L);
                        if(t>0 && t<r1){
                            d = true;
                            break;
                        }
                    }
                    if(!d){
                        double lambert = dot(L.dir,normal);
                        double phong = pow(dot(reflection,r->dir),Shine);

                        if(lambert<0) lambert = 0;
                        if(phong<0) phong = 0;

                        for (int k=0; k<3; k++) {
                            current_color[k] += source_factor*lambert*co_efficients[DIFFUSE]*color[k];
                            current_color[k] += source_factor*phong*co_efficients[SPECULAR]*color[k];
                        }
                    }

                    if(level<recursion_level){
                        start = sum(intersectionPoint,scalarmultiply(1,reflection));
                        Ray reflectionRay(start, reflection);

                        double reflected_color[3];
                        int nearest = -1;
                        double t_min = DBL_MAX;

                        for(int k=0;k<objects.size();k++){
                            double t = objects[k]->intersect(&reflectionRay,reflected_color,0);
                            //double t = objects[k]->getIntersectingT(&reflectionRay);
                            if(t<=0) continue;
                            else if(t<t_min){
                                t_min = t;
                                nearest = k;
                            }
                        }
                        if(nearest!=-1){
                            double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                            for (int k=0; k<3; k++) {
                                current_color[k] += reflected_color[k]*co_efficients[REFLECTION];
                            }
                        }
                        for(int i=0;i<3;i++){
                if(current_color[i]<0) current_color[i] = 0;
                if(current_color[i]>1) current_color[i] = 1;
            }
                    }
            }


            return t;
    }

};


class General_Quadratic:public Object{
    public:
        double a,b,c,d,e,f,g,h,i,j;
        double L,W,H;

        General_Quadratic(double coff[10],vector3 p,double l,double w,double h){
            this->a = coff[0];
            this->b = coff[1];
            this->c = coff[2];
            this->d = coff[3];
            this->e = coff[4];
            this->f = coff[5];
            this->g = coff[6];
            this->h = coff[7];
            this->i = coff[8];
            this->j = coff[9];

            reference_point = p;
            this->L = l;
            this->W = w;
            this->H = h;
            //cout << reference_point.x << " " << reference_point.y << " " << reference_point.z << endl;
        }


        void draw(){
        }

        vector3 getNormal(vector3 intersectionPoint){
            double x = 2*( a*reference_point.x + b*reference_point.y + c*reference_point.z + d);
            double y = 2*( b*reference_point.x + e*reference_point.y + f*reference_point.z + g);
            double z = 2*( c*reference_point.x + f*reference_point.y + h*reference_point.z + i);

            vector3 res(x,y,z);
            return res;
        }

        bool checkLimits(vector3 i){
            double min_x = reference_point.x;
            double max_x = min_x + L;

            double min_y = reference_point.y;
            double max_y = min_y + W;

            double min_z = reference_point.z;
            double max_z = min_z + H;

            if( (L>0&&(min_x > i.x || i.x > max_x))|| (W>0 && (min_y > i.y || i.y > max_y)) || (H>0 && (min_z > i.z || i.z > max_z)) )
                return true;

            return false;
        }

        double getIntersectingT(Ray* ray) {
            double A = a*ray->dir.x*ray->dir.x + b*ray->dir.y*ray->dir.y + c*ray->dir.z*ray->dir.z;
            A += d * ray->dir.x * ray->dir.y + e * ray->dir.y * ray->dir.z + f * ray->dir.z * ray->dir.x;

            double B = 2 * (a * ray->start.x * ray->dir.x + b * ray->start.y * ray->dir.y + c * ray->start.z * ray->dir.z);
            B += d*(ray->start.x*ray->dir.y +ray->dir.x*ray->start.y)+e*(ray->start.y*ray->dir.z + ray->dir.y * ray->start.z)+f*(ray->start.z*ray->dir.x + ray->dir.z * ray->start.x);
            B += g * ray->dir.x + h * ray->dir.y + i * ray->dir.z;

            double C = a * ray->start.x * ray->start.x + b * ray->start.y * ray->start.y + c * ray->start.z * ray->start.z;
            C += d * ray->start.x * ray->start.y + e * ray->start.y * ray->start.z + f * ray->start.z * ray->start.x;
            C += g*ray->start.x + h*ray->start.y + i*ray->start.z + j;


            //cout << A << " " << B << " " << C << endl;
            double D = B*B - 4*A*C;

            if(D<0) return -1;

            double t1 = (-B+sqrt(D))/2*A;
            double t2 = (-B-sqrt(D))/2*A;
            //cout << t1 << "  " << t2 << endl;

            vector3 i1 = sum(ray->start,scalarmultiply(t1,ray->dir));
            vector3 i2 = sum(ray->start,scalarmultiply(t2,ray->dir));

            if(checkLimits(i1) && checkLimits(i2)) return -1;
            if(checkLimits(i1)) return t2;
            if(checkLimits(i2)) return t1;

            return t1>t2?t2:t1;

        }

        double intersect(Ray *r, double *current_color, int level){
            double t= getIntersectingT(r);
            //cout << t << endl;
            if(t<=0) return -1;
            if(level==0) return t;

            vector3 intersectionPoint = sum(r->start,scalarmultiply(t,r->dir));
            //colorAt = getColorAt(intersectionPoint)
            for (int i=0; i<3; i++) {
                current_color[i] = color[i]*co_efficients[AMBIENT];
            }

            vector3 normal = getNormal(intersectionPoint);
            vector3 reflection = getReflection(r, normal);

            for(int i=0;i<lights.size();i++){
                    vector3 direction = subtract(lights[i],intersectionPoint);
                    double r1 = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);
                    direction.normalize();

                    vector3 start = sum(intersectionPoint,scalarmultiply(1,direction));

                    Ray L(start, direction);
                    double arr[3];
                    int nearest = -1;
                    //double t_min = DBL_MAX;
                    bool d = false;

                    for(int j=0;j<objects.size();j++){
                        double t = objects[j]->intersect(&L,arr,0);
                        if(t>0 && t<=r1){
                            d = true;
                            break;
                        }
                    }
                    if(!d){
                        double lambert = dot(L.dir,normal);
                        double phong = pow(dot(reflection,r->dir),Shine);

                        if(lambert<0) lambert = 0;
                        if(phong<0) phong = 0;

                        for (int k=0; k<3; k++) {
                            current_color[k] += source_factor*lambert*co_efficients[DIFFUSE]*color[k];
                            current_color[k] += source_factor*phong*co_efficients[SPECULAR]*color[k];
                        }
                    }

                    if(level<recursion_level){
                        start = sum(intersectionPoint,scalarmultiply(1,reflection));
                        Ray reflectionRay(start, reflection);

                        double reflected_color[3];
                        int nearest = -1;
                        double t_min = DBL_MAX;

                        for(int k=0;k<objects.size();k++){
                            double t = objects[k]->intersect(&reflectionRay,reflected_color,0);
                            if(t<=0) continue;
                            else if(t<t_min){
                                t_min = t;
                                nearest = k;
                            }
                        }
                        if(nearest!=-1){
                            double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                            for (int k=0; k<3; k++) {
                                current_color[k] += reflected_color[k]*co_efficients[REFLECTION];
                            }
                        }
                    }
            }

            for(int i=0;i<3;i++){
                if(current_color[i]<0) current_color[i] = 0;
                if(current_color[i]>1) current_color[i] = 1;
            }
            return t;
    }

};
#endif
