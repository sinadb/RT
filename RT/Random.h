#pragma once
#include <random>

float random_float() {
static std::uniform_real_distribution<float> distribution(-1.0, 1.0);
static std::mt19937 generator;
return distribution(generator);
}
