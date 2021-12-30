#pragma once
#include "ray.h"



inline float random_float() {
    static std::uniform_real_distribution<float> distribution(-0.5, 0.5);
    static std::mt19937 generator;
    return distribution(generator);
}

inline void create_Diffused_ray(ray& r, Vec3& Normal) {
    Vec3 Unit_Circle_Centre = r.origin() + Normal;
    Vec3 Random_p = Vec3(random_float(), random_float(), random_float());
    while (Random_p.length() > 1.f)
    {
        Random_p = Vec3(random_float(), random_float(), random_float());
    }
    Vec3 Random_p_on_U_Circle = Unit_Circle_Centre + Random_p;
    r.direction() = Random_p_on_U_Circle - r.origin();

}
class Hittable
{
public:
	virtual bool Hit(ray& r) = 0;
    virtual Vec3 Normal(Vec3 Origin, Vec3 Intersection) = 0;
};

class Circles : public Hittable {
private:
	float m_R;
	Vec3 m_Centre;
    Vec3 m_Colour;
public:
	Circles(float R, const Vec3& Centre, const Vec3 colour) :m_R(R), m_Centre(Centre),m_Colour(colour) {}
    bool Hit(ray& r) {
        std::pair<Vec3, float> colour_depth;
        float z = 0;
        // check for intersection
        Vec3 oc = r.origin() - m_Centre;
        auto a = dot(r.direction(), r.direction());
        auto b = 2.0 * dot(oc, r.direction());
        auto c = dot(oc, oc) - m_R * m_R;
        auto discriminant = b * b - 4 * a * c;
        // calculate the intersection and make sure it is hitting from outside surface (n.p<0)

        


        if (discriminant > 0) {
           
            
           float solution_2 = (-b + sqrt(discriminant)) / 2 * a;
           float solution_1 = (-b - sqrt(discriminant)) / 2 * a;
           float z_t_max = r.origin()[2];
           if (solution_1 >= 0.f && solution_2 >= 0.f ) {
               float final_solution = std::min(solution_1, solution_2);
               
               r.colour() = m_Colour * (0.5f) + r.colour() * (0.5f);
               r.origin() = r.direction() * final_solution;
               auto N = Normal(m_Centre, r.origin());
               //updating ray direction
               r.origin() += N * 0.0000000000001f;
               N.normalise();
               create_Diffused_ray(r, N);
               return true;
           }
           else {
               // returning false if solution are in the rveerse direction of shot ray

               r.origin() = Vec3(0.0, 0.0, -1000000.f);
              
               return false;

           }

        }
        else {
            // no solutions at all
            r.origin() = Vec3(0.0, 0.0, -1000000.f);
         
            return false;
        }
    }
    Vec3 Normal(Vec3 Origin, Vec3 Intersection) {
        Vec3 n = (Intersection - Origin);
        return n;
    }

};

