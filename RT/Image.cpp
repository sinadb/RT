#include "Image.h"
#include <limits>
#include "randomstuff.h"



Vec3 Blue(0.0f, 0.0f, 1.0f);
Vec3 Green(0.0f, 1.0f, 0.0f);
Vec3 Red(1.0f, 0.0f, 0.0f);



Vec3 unit_circle() {
    Vec3 c = Vec3(random_float(), random_float(), random_float());
    while (c.length() > 1) {
        c = Vec3(random_float(), random_float(), random_float());
    }

    return c * 0.05;

}

inline Vec3 ray_color(ray r) {
    Vec3 unit_direction = r.direction().normalise();
    auto t = 0.5 * (unit_direction[1] + 1.0);
    return  Vec3(1.0, 1.0, 1.0) * (1.0 - t) + Vec3(0.5, 0.7, 1.0) * t;
}



Hit_Record current;
Vec3 Ray_Tracer(ray r, std::vector<Box>& object, std::vector<Vec3>& Frame, int Depth) {
    
    int recursion = Depth-1;
    Hit_Record temp;
    for (auto& i : object) {
        if (i.get_rec().Hit(r)) {
            // test the desirables intersections and save the shortest one in the temp variable before updating current.
            Hit_Record tester = i.get_obj()->Hit(r);
            if (tester.hit) {
                if (tester.z < temp.z) {
                    temp = tester;

                }
            }


        }
    }
    //for (auto& i : object) {
    //    
    //        // test the desirables intersections and save the shortest one in the temp variable before updating current.
    //        Hit_Record tester = i.get_obj()->Hit(r);
    //        if (tester.hit) {
    //            if (tester.z < temp.z) {
    //                temp = tester;

    //            }
    //        }

    //    
    //}
    if (temp.hit) {
        current = temp;
    }
    if (current.hit && recursion != 0) {
        current.z = std::numeric_limits<float>::infinity();
        current.hit = false;
        Ray_Tracer(current.r, object, Frame, recursion);
        

        
    }
    


    current.z = std::numeric_limits<float>::infinity();
    current.hit = false;
    //std::cout << current.r.colour() << "\n";
    return current.r.colour();
}









Image::Image(const std::string destination, Camera c, std::vector<Box> world, int samples) {
    auto start = Timer::Stop_watch();
    std::vector<Vec3> Frame;
    Canvas.open(destination, std::ios::in);
    if (!Canvas) {
        std::cout << "File does not exist!" << "\n";
    }
    int width = 400;
    int height = 400;
    
    Vec3 V_width = c.width();
    Vec3 V_height = c.height();
    Vec3 dir;
    Vec3 Origin = c.Origin();
    Vec3 lower_left = Origin - V_width * (0.5) - V_height * (0.5) - c.forward() * c.focal_distance();
    Vec3 temp_colour;
    int AA = samples - 1;
    int index = 0;
    Vec3 offset;
    int counter = height;
   

    
   
    Canvas << "P3\n" << width << ' ' << height << "\n255\n";

    for (int j = height; j >= 0; j--) {
     
        for (int i = 0; i < width; i++) {
            float u = (float)(i) / width;
            float v = (float)(j) / height;
            dir = lower_left + V_width * u + V_height * v - Origin ;
            ray r(Origin , dir, false);
            // ensuring we are returning background colour if no hit is registered
            r.colour() = ray_color(r);
            current.r.colour() = r.colour();
            temp_colour = Ray_Tracer(r, world, Frame, 1);
            while (AA!=0)
            {
                //offset = c.up() * unit_circle()[1] + c.right() * unit_circle()[0];
                offset = Vec3();
                float temp_u = random_float() / width + u;
                float temp_v = random_float() / height + v;
                dir = lower_left + V_width * temp_u + V_height * temp_v - Origin - offset;
                ray r(Origin + offset, dir, false);
                r.colour() = ray_color(r);
                current.r.colour() = r.colour();
                temp_colour += Ray_Tracer(r, world, Frame,5);
                AA--;
                current.reset();
                
            }
            AA = samples - 1;
            temp_colour *= (1.f / samples);
            temp_colour.Gamma_corrected(2);
            temp_colour *= 255.f;
            temp_colour.clamp(0, 255);
            
            Frame.push_back(temp_colour);
            Canvas << Frame[index++];
            current.reset();
               
           
        }
    } 
    auto end = Timer::Stop_watch();
    Timer::duration(end, start);


}
Image::~Image() { Canvas.close(); }
