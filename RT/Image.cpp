#include "Image.h"
#include <limits>

static int call_no = 0;
static int pixel = 0;



struct Hit_Record {
    bool hit = false;
    ray r;
    // frame_buffer test
    float z = std::numeric_limits<float>::infinity();
};
int in = 0;
Hit_Record Current;
Vec3 Ray_Tracer(ray& r, std::vector<Hittable*>& object, std::vector<Vec3>& Frame, int Depth, int object_Index) {
    int recursion = --Depth;
    Current.r = r;
    float temp_z;
    ray Original_ray = r;
    for (auto& i : object) {
        if (i->Hit(Original_ray)) {
            // update frame buffer
            temp_z = (Original_ray.origin() - r.origin()).length();
            if (temp_z <= Current.z) {
                Current.z = temp_z;
                Current.r = Original_ray;
                Current.hit = true;
            }
        }
        Original_ray = r;
    }
    if (Current.hit && recursion != 0) {
        Ray_Tracer(Current.r, object, Frame, recursion, object_Index);
        

        
    }
    Current.z = std::numeric_limits<float>::infinity();
    Current.hit = false;
    return Current.r.colour();
}




Vec3 Blue(0.0f, 0.0f, 1.0f);
Vec3 Green(0.0f, 1.0f, 0.0f);
Vec3 Red(1.0f, 0.0f, 0.0f);

//Vec3 Normal(ray& intersection, Vec3&& C_centre) {
//
//    Vec3 Normal = (intersection.direction() - C_centre).normalise();
//    if (dot(intersection.direction(), Normal) < 0) {
//        return Blue;
//    }
//    else {
//        return Red;
//    }
//    
//}


//std::vector<std::pair<Vec3, float>> Frame_Data1;

Circles* c1 = new Circles(.5f,Vec3(0, 0.1f, -1.f), Red );
//Circles c(.9f, Vec3(0, 0, -1.f), Red + Green);
Circles* c = new Circles(100.f, Vec3(0, -100.5f, -1.f), Green);
//Circles c1(100.f, Vec3(0, -100.5f, -1.f), Green);
Circles* c2 = new Circles(.5f, Vec3(0.5f, 0.0f, -1.f), Red );
Circles* c3 = new Circles(.5f, Vec3(-.5f, 0.0f, -1.f), Green);
Circles* c4 = new Circles(.5f, Vec3(-0.5, 0.1f, -1.f), Red);
std::vector<Hittable*> shapes = { c1, c };





Image::Image(const std::string destination, int width, int height, Vec3& colour, int samples) {
    std::vector<Vec3> Frame;
    Canvas.open(destination, std::ios::in);
    if (!Canvas) {
        std::cout << "File does not exist!" << "\n";
    }
    Vec3 lower_left(-2, -2, -1);
    Vec3 V_width(4, 0, 0);
    Vec3 V_height(0, 4, 0);
    Vec3 dir;
    Vec3 Origin(0.0, 0.0, 0.0);
    Vec3 Starting_Colour = Blue;
    int index = 0;
    Vec3 temp_colour;
    int AA = samples - 1;
    Canvas << "P3\n" << width << ' ' << height << "\n255\n";

    for (int j = height; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            float u = (float)(i) / width;
            float v = (float)(j) / height;
            u = 4.f * u;
            v = 4.f * v;
            Vec3 dir = Vec3(u, v, 0) + lower_left;
            ray r(Origin, dir, Vec3(0.0, 0.0, 1.0), false);
            temp_colour = Ray_Tracer(r, shapes, Frame, 10, index);
            while (AA!=0)
            {
                float temp_u = random_float() / width + u;
                float temp_v = random_float() / height + v;
                dir = Vec3(temp_u, temp_v, 0) + lower_left;
                ray r(Origin, dir, Vec3(0.0, 0.0, 1.0), false);
                temp_colour += Ray_Tracer(r, shapes, Frame, 10, index);
                AA--;
            }
            AA = samples - 1;
            
            /*float temp_u = random_float() + u;
            float temp_v = random_float() + v;*/
            //while (AA != 0) {
            //    temp_u = (random_float()/width) + u;
            //    temp_v = (random_float()/height) + v;
            //    dir = Vec3(temp_u, temp_v, 0) + lower_left;
            //    ray r(Vec3(0.0, 0.0, 0.0), dir, Vec3(0.0,0.0,1.0), false);
            //    Ray_Tracer(r, shapes, Frame_Data,1);
            //    AA--;
            //}
            //// next pixel 
            //pixel++;
            //// set n_Call = 0
            //call_no = 0;
            //AA = samples-1;
            temp_colour *= (1.f / samples);
            temp_colour *= 255.f;
            temp_colour.clamp(0, 255);
            Frame.push_back(temp_colour);
            Canvas << Frame[index++];
            
               
           
        }
    }    


}
Image::~Image() { Canvas.close(); }
