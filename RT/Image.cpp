#include "Image.h"

static int call_no = 0;
static int pixel = 0;



//void Ray_Tracer(ray& r, std::vector<Hittable*>& objects, std::vector<std::vector<std::pair<Vec3, float>>>& Frame, int depth) {
//    Frame.resize(objects.size());
//    std::pair<Vec3, float> p;
//    ray temp_r = r;
//    if (call_no == 0) {
//        for (int j = 0; j != Frame.size(); j++) {
//            if (objects[j]->Hit(r) && depth!=0) {
//                if (p.second == 0) {
//                    p.second = r.origin()[2];
//                }
//                Ray_Tracer(r, objects, Frame, depth - 1);
//                }
//                p.first = r.colour();
//                Frame[j].push_back(p);
//                r = temp_r;
//            }
//        call_no = 1;
//            
//    }
//   
//   
//   else {
//        
//        for (int j = 0; j != Frame.size(); j++) {
//                objects[j]->Hit(r);
//                Frame[j][pixel].first += r.colour();
//                r = temp_r;
//            
//           
//        }
//
//    }
//}
std::pair<Vec3, float> temp_Frame_data;
int counter = 0;
//indicating first bounce
bool hit = true;
void Ray_Tracer(ray& r, std::vector<Hittable*>& object, std::vector<std::pair<Vec3, float>>& Frame, int Depth, int object_ID) {

    if (object[object_ID]->Hit(r) && Depth != 0) {
        if (hit) {
            temp_Frame_data.second = r.origin()[2];
            hit = false;
        }
      
        for (int i = 0; i != object.size(); i++) {
           if (i == object_ID) {
               continue;
            }
     
            Ray_Tracer(r, object, Frame, Depth - 1, i);
        }
    }
    
    else {
        if (hit) {
            temp_Frame_data.second = r.origin()[2];
            hit = false;
        }
        temp_Frame_data.first = r.colour();
        Frame.push_back(temp_Frame_data);

    }
    counter = 0;
    hit = true;
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

std::vector<std::vector<std::pair<Vec3, float>>> Frame_Data;
//std::vector<std::pair<Vec3, float>> Frame_Data1;

Circles* c1 = new Circles(.5f,Vec3(0, 0, -1.f), Red );
//Circles c(.9f, Vec3(0, 0, -1.f), Red + Green);
Circles* c = new Circles(100.f, Vec3(0, -100.5f, -1.f), Green);
//Circles c1(100.f, Vec3(0, -100.5f, -1.f), Green);
Circles* c2 = new Circles(.5f, Vec3(0, 1.0f, -1.f), Green );
std::vector<Hittable*> shapes = { c1,c };





Image::Image(const std::string destination, int width, int height, Vec3& colour, int samples) {

    Canvas.open(destination, std::ios::in);
    if (!Canvas) {
        std::cout << "File does not exist!" << "\n";
    }
    int counter = 2;
    Frame_Data.resize(shapes.size());
    Vec3 lower_left(-2, -2, -1);
    Vec3 V_width(4, 0, 0);
    Vec3 V_height(0, 4, 0);
    Vec3 rgb{};
    float z_buffer = 0;
    Vec3 temp_rgb{};
    std::pair<Vec3, float> colour_buffer;
    colour_buffer.first = Vec3(0.0, 0.0, 0.0);
    std::pair<Vec3, float> temp_colour_buffer;
    int AA = samples-1;
    Vec3 dir;
    int index = 0;
    bool flip = true;
    Canvas << "P3\n" << width << ' ' << height << "\n255\n";

    for (int j = height; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            float u = (float)(i) / width;
            float v = (float)(j) / height;
            u = 4.f * u;
            v = 4.f * v;
            Vec3 dir = Vec3(u, v, 0) + lower_left;
            ray r(Vec3(0.0, 0.0, 0.0), dir, Vec3(0.0, 0.0, 0.0), false);
            ray temp_r = r;
            for (int k = 0; k != shapes.size(); k++) {
                Ray_Tracer(r, shapes, Frame_Data[k], 10000, k);
                r = temp_r;
            }
            
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

               
           
        }
    }    
float temp_z = 0.f;
    Vec3 final_colour;
    for(int j = 0; j < width*height; j++) {
        for (auto& i : Frame_Data) {

            z_buffer = i[j].second;
            if (temp_z == 0.f) {
                temp_z = z_buffer;
            }
            
            if (z_buffer >= temp_z) {
                z_buffer = i[j].second;
                temp_z = z_buffer;
                final_colour = i[j].first;
            }



           
        }

        
        final_colour *= (1.f / samples);
        final_colour *= 255;
        final_colour.clamp(0, 255);
        
        temp_z = 0.f;
        Canvas << final_colour;
        
    }
}
Image::~Image() { Canvas.close(); }
