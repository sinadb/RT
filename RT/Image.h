#pragma once
#include <fstream>
#include <iostream>
#include "Hittable.h"
#include <stdlib.h>

class Image
{
private : 
	std::ofstream Canvas;
public:
	Image(const std::string destination, int width, int height, Vec3& colour, int samples = 1);
	~Image();


};

