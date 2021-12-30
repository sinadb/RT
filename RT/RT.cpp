// RT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "Image.h"


	
	

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // Camera

    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;

    auto origin = Vec3(0, 0, 0);
    auto horizontal = Vec3(viewport_width, 0, 0);
    auto vertical = Vec3(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal * 0.5f - vertical * 0.5f - Vec3(0, 0, focal_length);



    int main()
    {
        std::string name = "Image1.ppm";
        Vec3 red(1, 1.3, 2.0);
        Vec3 Blue(-2.0, -5.0, 6.0);
       Image I(name, 400, 400, red);
        //Vec3 blue = (red);
        //std::cout << blue ;
        
       
       
    }

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
