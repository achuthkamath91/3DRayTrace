//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

// -------------------------------------------------------------------
// Ray struct

struct Ray
{
    vec3 origin;		// origin of the ray
    vec3 dir;			// direction of the ray
};


// -------------------------------------------------------------------
// Sphere struct

struct Sphere
{	
public: vec3 center;
public:	float radius;
//public: float t_min;
	vec3 ka, kd, ks;	// ambient, diffuse and specular reflecction constant in Phong's reflection model
	vec3 reflectivity;	// control how much light is received from reflection in recursive ray-tracing (e.g. 0.1)
	float alpha;		// control size of specular highlight (e.g. 20)
	

	// default constructor
	Sphere(const vec3& ic=vec3(0.0f), const float& ir=0.0f, const vec3& ika=vec3(0.0f), const vec3& ikd=vec3(0.0f), const vec3& iks=vec3(0.0f), const float& ireflectivity=0.1f, const float& ialpha=1.0f):
	center(ic), radius(ir), ka(ika), kd(ikd), ks(iks), reflectivity(ireflectivity), alpha(ialpha)
	{}

	bool intersect(const Ray& ray, float& t0, float& t1);
};

//function to solve quadratic dcricrimination of a,b,c
bool Quadratic(const float &a, const float &b, const float &c, float &t0, float &t1)
{
	float d = b*b - 4*a*c;

	if(d<0)
	{
		return false;
	}
	if(d == 0)
	{
		t0=t1=-0.5*b/a;
	}
	if(d>0)
	{
		float quo = (b>0)? -0.5 *(b + sqrt(d)):-0.5 *(b - sqrt(d));
		t0 = quo/a;
		t1 = c/quo;
	}
	if(t0>t1) std::swap(t0,t1);
	return true;
} 

// TODO:
// return true if the input ray intersects with the sphere; otherwise return false
// return by reference the intersections. t0 refers to closer intersection, t1 refers to farther intersection
bool Sphere::intersect(const Ray& ray, float& t0, float& t1)
{
	// First transform before testing intersection
	/*vec3 centercall =  Sphere::center;
	

	//quadratic equation (dx^2 + dy^2 + dz^2 )*t^2 + 2(dx*x0 + dy*y0 + dz*z0) t + (x0^2 + y0^2 + z0^2) – 1 = 0
	float  a = dx*dx + dy*dy + dz*dz;
	//float  b = 2*dx*(x0-cx) +  2*dy*(y0-cy) +  2*dz*(z0-cz);
	//float  c = cx*cx + cy*cy + cz*cz + x0*x0 + y0*y0 + z0*z0 -2*(cx*x0 + cy*y0 + cz*z0) - 1;
	;*/

	//First transform before testing intersection
	 float tNear;
	vec3 Lvector = ray.origin - center;
	float a = dot(ray.dir , ray.dir);
	float b = 2 * dot(Lvector , ray.dir);
	float c = dot(Lvector , Lvector) - radius*radius;
	if(!Quadratic(a,b,c,t0,t1)){
		return false;
	}
	tNear = t0;
	return true;
};




// -------------------------------------------------------------------
// Light Structs

struct AmbientLight
{
	vec3 ia;		// ambient intensity (a vec3 of 0.0 to 1.0, each entry in vector refers to R,G,B channel)

	// default constructor
	AmbientLight(const vec3& iia=vec3(0.0f)):
	ia(iia)
	{}

};


struct PointLight
{
	vec3 location;	// location of point light
	vec3 id, is;	// diffuse intensity, specular intensity (vec3 of 0.0 to 1.0)
	
	// default constructor
	PointLight(const vec3& iloc=vec3(0.0f), const vec3& iid=vec3(0.0f), const vec3& iis=vec3(0.0f)):
	location(iloc), id(iid), is(iis)
	{}

};

struct reflection
{
	 vec3 origin;		// origin of the ray
	 vec3 dir;			// direction of the ray
};


// -------------------------------------------------------------------
// Our Scene

// lights and spheres in our scene
AmbientLight my_ambient_light;			// our ambient light
vector<PointLight> my_point_lights;		// a vector containing all our point lights
vector<Sphere> my_spheres;				// a vector containing all our spheres


// this stores the color of each pixel, which will take the ray-tracing results
vector<vec3> g_colors;	

int recursion_lvl_max = 2;

// this defines the screen
int g_width = 640;				//number of pixels
int g_height = 480;				// "g_" refers to coord in 2D image space
float fov = 30;					// field of view (in degree)

float invWidth = 1 / float(g_width);
float invHeight = 1 / float(g_height);
float aspectratio = g_width / float(g_height);
float angle = tan(M_PI * 0.5 * fov / float(180));


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec3& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}



// -------------------------------------------------------------------
// Ray tracing

vec3 trace(const Ray& ray, int recursion_lvl)
{
	float INFINITY = 99999;
	float t_min = INFINITY;
	float t0 = INFINITY;
	int near_sphere_idx;
	bool has_intersection = false;

	for (int i=0; i<my_spheres.size(); ++i)
	{
				//some large value
		float t1 = INFINITY;

		// check intersection with sphere
		if (my_spheres[i].intersect(ray, t0, t1))
		{
			has_intersection = true;

			if (t0<t_min)
			{
				t_min = t0;
				near_sphere_idx = i;
			}
		}	
	}

	if (has_intersection == false)
		// just return background color (black)
	    return vec3(0.0f, 0.0f, 0.0f);

	Sphere my_sphere = my_spheres[near_sphere_idx];

	// TODO: implement Phong's reflection model
	//vec3 color(1.0,1.0,1.0);	// this code should be replaced by codes of the Phong's reflection model (i.e. color should be determined by Phong's model)
	//int n =3;

	//vec3 I = ka*ia;
	vec3 hit_point = ray.origin + ray.dir * t_min;
	vec3 intersection(hit_point[0],hit_point[1],hit_point[2]);
	vec3 n = normalize((intersection - my_sphere.center)/my_sphere.radius);
	//vec3 n = hit_point;
	//vec3 No = ray.dir;
	//vec3 l = ray.dir;
	//l = normalize(l);
	
	//vec3 nl = dot(l,n);
	//vec3 r = (2*nl*n) - l;
	//AmbientLight al;
	//PointLight pl;
	//Sphere sp;
	vec3 ia = my_ambient_light.ia;	
	vec3 D;
	vec3 sigma = 0.0;
	//vec3 ph =vec3(0.0);
	vec3 ka,kd,ks,lm,rm;
	
	//for (int k=0; k<my_spheres.size(); k++)
		//{
	vec3 id=0,is=0;
	float min_no = 0.0;
	float shine = 0.2,attenuation = 1.0;
	vec3 v = -(ray.dir) ;
	v = normalize(v);
	for(int j = 0; j<my_point_lights.size(); j++)
	{
			
		ka = my_spheres[near_sphere_idx].ka;
		kd = my_spheres[near_sphere_idx].kd;				
		ks = my_spheres[near_sphere_idx].ks;
		float alpha = my_sphere.alpha;
		//vec3 reflect = my_spheres[near_sphere_idx].reflectivity;
		//float reflect = reflect[0];
		//vec3 rayorb =   my_point_lights[j].location - hit_point;
		lm =   my_point_lights[j].location - intersection ;
		lm = normalize(lm);
		
		
		rm = (2.0*max(dot(lm,n),min_no)*n) - lm;
		//r = (2*nl*n) - l;
		//rm = normalize(rm);
		
		id = my_point_lights[j].id;
		is = my_point_lights[j].is;
		//sigma += kd*max(dot(lm,n),min_no)*id + ks*pow(max(dot(rm,v),min_no),shine)*is;
		//c kd*(dot(lm,n))*id + ks*pow(max(dot(lmn,v),min_no),shine)*is;
		//sigma += kd*dot(lm,n)*id;
		sigma +=  kd*(max(dot(lm,n),min_no))*id + ks*pow(max(dot(rm,v),min_no),alpha)*is;
				
		//vec3 ph = ka+kd+ks;				
	}
	
	//D =  ia*ka + (kd*nl + ks * dot(v,r))*(id/my_spheres[near_sphere_idx].radius*my_spheres[near_sphere_idx].radius);
	D = ia*ka + sigma;
			
	vec3 color(D[0],D[1],D[2]);
	
	 // setup point lights

	


	if (recursion_lvl >0)
	{
		// TODO: implement recursive ray-tracing here, to add contribution of light reflected from other objects

		//return color;		// this should be replaced by codes to do recursive ray-tracing
		//return color = trace( ray, 2);
		Ray refl;
		refl.dir = v-2*n*(max(dot(v,n),min_no));
		refl.origin = intersection;
		//color +=	trace(refl,2)*my_sphere.reflectivity;
		return (color +=	(trace(refl,(recursion_lvl-1)))*my_sphere.reflectivity);
	}
	else 
	{
		return color;
	}			
}


vec3 getDir(int ix, int iy)
{
    // This should return the direction from the origin
    // to pixel (ix, iy), normalized!

	vec3 dir;
	dir.x = (2 * ((ix + 0.5) * invWidth) - 1) * angle * aspectratio;
	dir.y = (2 * ((iy + 0.5) * invHeight) - 1) * angle;
	dir.z = -1;

	return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec3(0.0f, 0.0f, 0.0f);
    ray.dir = getDir(ix, iy);
    vec3 color = trace(ray, recursion_lvl_max);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
    
    // change file name based on input file name.
    savePPM(g_width, g_height, "output.ppm", buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
	// setup pixel array
    g_colors.resize(g_width * g_height);

	// setup our scene...

	// setup ambient light
	my_ambient_light = AmbientLight(vec3(0.1));

	// setup point lights
	my_point_lights.push_back(
								PointLight(vec3(3,3,0), vec3(0.5, 0.5, 0.5), vec3(0.5,0.5,0.5))
								);

	my_point_lights.push_back(
								PointLight(vec3(-3,-3,0), vec3(0.1, 0.1, 0.1), vec3(0.1,0.1,0.1))
								);

	// setup spheres
	my_spheres.push_back(
							Sphere(vec3(0,0,-10), 1.0, vec3(0.1,0.1,0.1), vec3(0.5,0.5,0.5), vec3(0.5,0.5,0.5), 0.2, 100.0)
							);

	my_spheres.push_back(
							Sphere(vec3(-1.5,0.5,-8), 0.5, vec3(0.0,1.0,0.0), vec3(0.0,1.0,0.0), vec3(0.5,0.5,0.5), 0.0, 10.0)
							);

    render();
    saveFile();
	return 0;
}

