#pragma once
#include "ray.h"


static Vec3 Point_Light = Vec3(0, -5.f, 1.f);
static Vec3 White = Vec3(1., 1., 1.);

struct Hit_Record {
    bool hit;
    ray r;
    // frame_buffer test
    float z;
    Vec3 Normal;

    Hit_Record() {
        hit = false;
        r.colour() = Vec3(0, 0, 1.0);
        z = std::numeric_limits<float>::infinity();
    }
    void reset() {
        r = ray();
        hit = false;
        r.colour() = Vec3(0, 0, 1.0);
        z = std::numeric_limits<float>::infinity();
    }

};

class Hittable
{
public:
    virtual Hit_Record Hit(ray r) = 0;
    virtual Vec3 Centre() = 0;
    virtual Vec3 Colour() = 0;
    virtual float Albedo() = 0;
    virtual float fuzzy_factor() = 0;
    virtual bool is_refractionable() = 0;
    virtual float rio() = 0;

};

inline Vec3 Normal(Vec3 Origin, Vec3 Intersection) {
    Vec3 n = (Intersection - Origin).normalise();
    return n;
}



inline float random_float() {
    static std::uniform_real_distribution<float> distribution(-1.f, 1.f);
    static std::mt19937 generator;
    return distribution(generator);

}

enum class Material {
  Matte,
  Metal,
  Fuzzy,
  Refractable


};


class Fuzzy {

public:
   
    static ray init_ray(ray r, float solution, Hittable* object) {
        r.Bounces()++;
      
        // new origin
        r.origin() += r.direction() * solution;
        // Normal 
        Vec3 N = Normal(object->Centre(), r.origin());
        //avoid delf reflection
        r.origin() += N * 0.0001f;
        //new reflected ray
        r.direction() = r.direction() - N * (2 * (dot(N, r.direction())));
        // offset vector by the fuzziness factor
        Vec3 offset = Vec3(random_float(), random_float(), random_float());
        while (offset.length() > object->fuzzy_factor()) {
             offset = Vec3(random_float(), random_float(), random_float());
        }
        r.direction() += offset;
        if (dot(r.direction(), N) > 0) {
            r.colour() = (r.colour()) * (object->Colour() * 0.8f);
        }
        else {
            r.colour() = Vec3();
        }
       
        return r;

    }

public:

};

class Refractable {

public:

    static ray init_ray(ray r, float solution, Hittable* object) {

        // new origin
        r.origin() += r.direction() * solution;
        // Normal 
        Vec3 N = Normal(object->Centre(), r.origin());
        //avoid delf reflectio
        //reverse direction for consistency throughout code . incident light is now facing away from object
        r.direction() = r.direction() * (-1.f);
            float test = dot(N, r.direction());
            float rio;
            if (test > 0) {
                rio = 1.f / object->rio();
                
            }
            else {
                rio = object->rio();
                N = N * (-1.f);
            }
           
            r.origin() += N * -0.00001f;
            // find out the cos of the angle between incident ray and Normal
            r.direction().normalise();
            float COS_Theta = dot(r.direction(), N);
            // find the sin now 
            float SIN_Theta = sqrtf(1.f - powf(COS_Theta, 2));
            // create refracted ray but first test that we can refract 
            float Critical = SIN_Theta * rio;
            // only refracting if less than 1
	    if(Critical < 0){
           
                Vec3 R_Para = N * -sqrtf(1.f - powf(rio * SIN_Theta, 2));
                Vec3 R_Perp = (r.direction() - N * dot(N, r.direction())) * (-rio);
                Vec3 R = R_Para + R_Perp;
                r.direction() = R;
                r.colour() = Vec3(1, 1, 1) * r.colour();
                //std::cout << "h";
                return r;}

        r.direction() = N * (2 * (dot(N, r.direction()))) - r.direction();
        //r.colour() = (r.colour()) * (object->Colour() * 0.8f);
        r.colour() = Vec3(1, 1, 1) * r.colour();
        return r;
    }

};
class Matte {

public : 

    static ray init_ray(ray r, float solution, Hittable* object) {

        
        r.Bounces()++;
        r.origin() += r.direction() * solution;
        // Normal 
        Vec3 N = Normal(object->Centre(), r.origin());
        r.origin() += N * 0.0001f;
        float specular_component = 0;
        // only do the calculation for specular light if this is the first bounce 
        //if (r.Bounces() == 1) {
        //     //save a copy of original r.dire which is equal to the viewing vector coming straight out of the camera and negate it for the half vector calculation 
        //    Vec3 Viewing_vector = (r.direction() * -1.f ).normalise();
        //    Vec3 incidece_light = (Point_Light - r.origin()).normalise();
        //    Vec3 H = (Viewing_vector + incidece_light).normalise();
        //
        //    // dot(N,L) > 0 to ensure correct exposure 
       
        //    specular_component = dot(H, N);
        //    if (specular_component < 0) {
        //        specular_component = 0;
        //    }
      
        //    /*Vec3 R = r.direction() - N * (2 * (dot(N, r.direction())));
        //    Vec3 V = r.direction() * -1.f;
        //    specular_component = dot(R, V) / (R.length() * V.length());*/
        //    
        //    /*if (specular_component < 0) {
        //        specular_component = 0;
        //    }*/
        //}
        
        //avoid delf reflection
        
        Vec3 Unit_Circle_Centre = r.origin() + N;
        Vec3 Random_p = Vec3(random_float(), random_float(), random_float());
        while (Random_p.length() > 1.f)
        {
            Random_p = Vec3(random_float(), random_float(), random_float());
        }
        Vec3 Random_p_on_U_Circle = Unit_Circle_Centre + Random_p.normalise();
        r.direction() = Random_p_on_U_Circle - r.origin();
       
        // dot(N,L)/(dot(N,L).length) then I/offset witll give us the diffused ray intensity 
        
       
        r.colour() = (r.colour()) * (object->Colour() * 0.4f) + White * powf(specular_component,50.f);

        return r;
    }

};

class Metal {


public:


   static ray init_ray(ray r, float solution, Hittable* object) {
       r.Bounces()++;
        // new origin
       r.origin() += r.direction() * solution;
       // Normal 
       Vec3 N = Normal(object->Centre(), r.origin());
       //avoid delf reflection
       r.origin() += N * 0.0001f;
       //reverse direction for consistency throughout code . incident light is now facing away from object
       r.direction() = r.direction() * (-1.f);
       //if (object->is_refractionable()) {
       //    float test = dot(N, r.direction());
       //    float rio = (test > 0) ?  object->rio() : 1.f / object->rio();
       //    if (test < 0) {
       //        N = N * (-1.f);

       //    }

       //    // find out the cos of the angle between incident ray and Normal
       //    r.direction().normalise();
       //    float COS_Theta = dot(r.direction(), N);
       //    // find the sin now 
       //    float SIN_Theta = sqrtf(1.f - powf(COS_Theta, 2));
       //    // create refracted ray but first test that we can refract 
       //    float Critical = SIN_Theta * rio;
       //    // only refracting if less than 1
       //    if (Critical < 1.f) {
       //        Vec3 R_Para = N * -sqrtf(1.f - powf(rio * SIN_Theta, 2));
       //        Vec3 R_Perp = (r.direction() - N * dot(N, r.direction())) * (-rio);
       //        Vec3 R = R_Para + R_Perp;
       //        r.direction() = R;
       //        r.colour() = Vec3(1, 1, 1) * r.colour();
       //        std::cout << "h";
       //        return r;

       //    }

       //}
      
       // reset Normal if we are not refracting;
       //N = N * (-1.f);
       //new reflected ray
      
        r.direction() =   N * (2 * (dot(N, r.direction()))) - r.direction();
        r.colour() = (r.colour()) * (object->Colour() * 0.8f);
        //r.colour() = Vec3(1, 1, 1) * r.colour();
        return r;
        

    }

};



//inline void create_Diffused_ray(ray& r, Vec3& Normal) {
//    Vec3 Unit_Circle_Centre = r.origin() + Normal;
//    Vec3 Random_p = Vec3(random_float(), random_float(), random_float());
//    while (Random_p.length() > 1.f)
//    {
//        Random_p = Vec3(random_float(), random_float(), random_float());
//    }
//    Vec3 Random_p_on_U_Circle = Unit_Circle_Centre + Random_p.normalise();
//    r.direction() = Random_p_on_U_Circle - r.origin();
//
//}

//inline 





class Circles : public Hittable {
private:
	float m_R;
	Vec3 m_Centre;
    Vec3 m_Colour;
    float m_Albedo;
    Material m_material;
    float m_fuzzy_factor;
    bool m_Refraction;
    float m_rio;
public:
	Circles(float R, const Vec3& Centre, const Vec3 colour, Material m, float f = 0.3f, float rio = 1.f) :m_R(R), m_Centre(Centre),m_Colour(colour),m_material(m),m_Albedo(.5f),m_fuzzy_factor(f),m_rio(rio){}
    Hit_Record Hit(ray r) {
        
        Hit_Record current;
       
        // check for intersection
        Vec3 oc = r.origin() - m_Centre;
        auto a = dot(r.direction(), r.direction());
        auto b = 2.0 * dot(oc, r.direction());
        auto c = dot(oc, oc) - m_R * m_R;
        auto discriminant = b * b - 4 * a * c;
        // calculate the intersection and make sure it is hitting from outside surface (n.p<0)


        if (discriminant > 0.0f) {
           
          
           float solution_2 = (-b + sqrt(discriminant)) / (2 * a);
           float solution_1 = (-b - sqrt(discriminant)) / (2 * a);
           if (solution_1 > 0.0f || solution_2 > 0.0f ) {
               float final_solution = std::min(solution_1, solution_2);
             
               
               
               switch (m_material) {
               case Material::Matte:
                   current.r = Matte::init_ray(r, final_solution, this);
                   break;
               case Material::Metal:
                   current.r = Metal::init_ray(r, final_solution, this);
                   break;

               case Material::Fuzzy:
                   current.r = Fuzzy::init_ray(r, final_solution, this);
                   break;
               case Material::Refractable:
                   current.r = Refractable::init_ray(r, final_solution, this);
                   break;
               }
               // current.r.oring - r.origin = z

               current.z = (current.r.origin() - r.origin()).length();
               current.Normal = Normal(m_Centre, current.r.origin());
               current.hit = true;
               return current;
           }
           else {
               // returning false if solution are in the rveerse direction of shot ray

               
               current.r.colour() = Vec3(0, 0, 1.0);
               return current;

           }

        }
        else {
            // no solutions at all
            
            current.r.colour() = Vec3(0, 0, 1.0);
            return current;
        }
    }

    Vec3 Centre() { return m_Centre; }
    float Albedo() { return m_Albedo; }
    Vec3 Colour() { return m_Colour; }
    float fuzzy_factor() { return m_fuzzy_factor; }
    bool is_refractionable() { return m_Refraction; }
    float rio() { return m_rio; }
    

};

