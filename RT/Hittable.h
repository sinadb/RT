#pragma once
#define INF 10000000.f
#include "ray.h"
#include "randomstuff.h"
#include "glm.hpp"
#include <algorithm>


static Vec3 Point_Light = Vec3(1 , 0, -1);
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
    virtual Vec3 Normal(Vec3 Intersection) = 0;
    virtual float radius() = 0;

};







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
        Vec3 N = object->Normal(r.origin());
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
            r.colour() = (r.colour()) * (object->Colour() * 1.f);
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
        r.Bounces()++;
        // new origin
        r.origin() += r.direction() * solution;
        // Normal 
        Vec3 N = object->Normal(r.origin());

        //avoid delf reflection
       // r.origin() += N * 0.0001f;
        //avoid delf reflection

        //reverse direction for consistency throughout code . incident light is now facing away from object
        r.direction() = r.direction() * (-1.f);
        r.direction().normalise();

        float test = dot(N, r.direction());
        float rio;
        if (test > 0) {
            rio = 1.f / object->rio();
            //printf("%f\n", test);
            //std::cout << "p";

        }
        else {

            rio = object->rio();
            N = N * (-1.f);
            //std::cout << "n";
           // printf("%f\n", test);
        }


        // find out the cos of the angle between incident ray and Normal
        float COS_Theta = dot(r.direction(), N);
        // find the sin now 
        float SIN_Theta = sqrtf(1.f - powf(COS_Theta, 2));
        // create refracted ray but first test that we can refract 
        float Critical = SIN_Theta * rio;
        // only refracting if less than 1
        if (Critical < 1) {
            r.origin() += N * -0.0001f;
            Vec3 R_Para = N * (-sqrtf(1.f - powf(rio * SIN_Theta, 2)));
            Vec3 R_Perp = (r.direction() - N * dot(N, r.direction())) * (-rio);
            Vec3 R = R_Para + R_Perp;
            r.direction() = R;
            r.colour() = Vec3(1, 1, 1) * r.colour();
            //std::cout << "p";

            return r;
        }




        r.origin() += N * 0.0001f;

        //reset Normal if we are not refracting;
        //N = N * (-1.f);


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
        Vec3 N = object->Normal(r.origin());
        r.origin() += N * 0.0001f;
        float specular_component = 0;
         //only do the calculation for specular light if this is the first bounce 
        if (r.Bounces() == 1) {
             //save a copy of original r.dire which is equal to the viewing vector coming straight out of the camera and negate it for the half vector calculation 
            Vec3 Viewing_vector = (r.direction() * -1.f ).normalise();
            Vec3 incidece_light = (Point_Light - r.origin()).normalise();
            Vec3 H = (Viewing_vector + incidece_light).normalise();
        
            // dot(N,L) > 0 to ensure correct exposure 
       
            specular_component = dot(H, N);
            if (specular_component < 0) {
                specular_component = 0;
            }
      
            /*Vec3 R = r.direction() - N * (2 * (dot(N, r.direction())));
            Vec3 V = r.direction() * -1.f;
            specular_component = dot(R, V) / (R.length() * V.length());*/
            
           
        }
        
        //avoid delf reflection
        // calculate the direction of point light 
        Vec3 PL = (Point_Light - r.origin()).normalise();
        //calculate the cos of Normal and PL
        float cos = dot(PL, N);
        //std::cout << cos << "\n";
        // claculate the intensity which is cos times the colour of the point light
        Vec3 I = White * cos;
        Vec3 Unit_Circle_Centre = r.origin() + N;
        Vec3 Random_p = Vec3(random_float(), random_float(), random_float());
        while (Random_p.length() > 1.f)
        {
            Random_p = Vec3(random_float(), random_float(), random_float());
        }
        Vec3 Random_p_on_U_Circle = Unit_Circle_Centre + Random_p.normalise();
        r.direction() = Random_p_on_U_Circle - r.origin();
       
        // dot(N,L)/(dot(N,L).length) then I/offset witll give us the diffused ray intensity 
        
       
        /*r.colour() = (r.colour()) * (object->Colour() * 0.5f) + White * powf(specular_component,50.f);*/
        //r.colour() = I * object->Colour();
        r.colour() = object->Colour() * 0.3 + White * powf(specular_component, 50);
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
       Vec3 N = object->Normal(r.origin());
       //avoid delf reflection
       r.origin() += N * 0.0001f;
       r.direction() = r.direction() * (-1.f);
       //reverse direction for consistency throughout code . incident light is now facing away from object
       //r.direction() = r.direction() * (-1.f);
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
      
        r.direction() =   N * (2 * (dot(N, r.direction())))  - r.direction();
        r.colour() = (r.colour()) * (object->Colour() * 1.f);
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
    float m_R_2;
public:
	Circles(float R, const Vec3& Centre, const Vec3 colour, Material m, float f = 0.3f, float rio = 1.f) :m_R(R), m_Centre(Centre),m_Colour(colour),m_material(m),m_Albedo(.5f),m_fuzzy_factor(f),m_rio(rio){
        m_R_2 = m_R * m_R;
    }
    Hit_Record Hit(ray r) {
        
        Hit_Record current;
        current.r = r;
        // check for intersection
        Vec3 oc = r.origin() - m_Centre;
        auto a = dot(r.direction(), r.direction());
        auto b = 2.0 * dot(oc, r.direction());
        auto c = dot(oc, oc) - m_R_2;
        auto discriminant = b * b - 4 * a * c;
        // calculate the intersection and make sure it is hitting from outside surface (n.p<0)


        if (discriminant > 0.0f) {
           
          
           float solution_2 = (-b + sqrt(discriminant)) / (2 * a);
           float solution_1 = (-b - sqrt(discriminant)) / (2 * a);
           if (solution_1 > 0.0f || solution_2 > 0.0f ) {
               float final_solution = std::min(solution_1, solution_2);
             
               if (final_solution < 0) {
                   final_solution = std::max(solution_1, solution_2);
               }
               
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


               current.z = (current.r.origin() - r.origin()).length();
               current.Normal = Normal(current.r.origin());
               current.hit = true;
               return current;
           }
           else {
             
              
               return current;

           }

        }
        else {

            return current;
        }
    }

    Vec3 Centre() { return m_Centre; }
    float Albedo() { return m_Albedo; }
    Vec3 Colour() { return m_Colour; }
    float fuzzy_factor() { return m_fuzzy_factor; }
    bool is_refractionable() { return m_Refraction; }
    float rio() { return m_rio; }
    Vec3 Normal(Vec3 Intersection) {
        Vec3 n = (Intersection - m_Centre).normalise();
        return n;
    }
    float radius() { return m_R; }
    

};
class Triangle : public Hittable {
private:
    // Anti clock-wise
    Vec3 m_vertices[3];
    Vec3 m_Normal;
    Vec3 m_colour;
    float d;
    Material m_material;
public:
    Triangle(const Vec3& v_0, const Vec3& v_1, const Vec3& v_2, Vec3 colour, Material material){
        m_vertices[0] = v_0;
        m_vertices[1] = v_1;
        m_vertices[2] = v_2;
        m_colour = colour;
        m_Normal = Cross_Product(m_vertices[1] - m_vertices[0], m_vertices[2] - m_vertices[0]); 
        d = - dot(m_Normal, m_vertices[0]);
        m_material = material;
        
    }
    Hit_Record Hit(ray r) {
        Hit_Record current;
        current.r = r;
        // normal.direction
        float test = dot(m_Normal, r.direction());
        
        if (test != 0.0) {
            float final_solution = -(d + dot(m_Normal, r.origin())) / test;
            
            Vec3 hit_point = r.direction() * final_solution;
            Vec3 R = hit_point - m_vertices[0];
            Vec3 Q_1 = m_vertices[1] - m_vertices[0];
            Vec3 Q_2 = m_vertices[2] - m_vertices[0];
            glm::vec2 RQ(dot(R, Q_1), dot(R, Q_2));
            float Q_1_mag = dot(Q_1, Q_1);
            float Q_2_mag = dot(Q_2, Q_2);
            float Q_1_dot_Q_2 = dot(Q_1, Q_2);
            float det =  Q_1_mag*Q_2_mag  - Q_1_dot_Q_2 * Q_1_dot_Q_2;
            
            det = 1 / det;
            
            glm::mat2x2 Matrix({det * Q_2_mag , det * -Q_1_dot_Q_2 }, { det*-Q_1_dot_Q_2 , det*Q_1_mag });
            glm::vec2 w( Matrix * RQ);
            float w1 = w[0];
            float w2 = w[1];
            /*Vec3 hit_point = r.direction() * final_solution;
            Vec3 R = hit_point - m_vertices[0];
            Vec3 Q_1 = m_vertices[1] - m_vertices[0];
            Vec3 Q_2 = m_vertices[2] - m_vertices[0];
            float w1 = dot(Cross_Product(Q_1, (hit_point - m_vertices[0])), Normal);
            float w2 = dot(Cross_Product((m_vertices[2]-m_vertices[1]), (hit_point - m_vertices[1])), Normal);
            float w3 = dot(Cross_Product((m_vertices[0] - m_vertices[2]), (hit_point - m_vertices[2])), Normal);*/
            //std::cout << w1 << "\t" << w2 << "\n";
            if (w1 >= 0 && w2 >= 0 && w1 <= 1 && w2 <= 1 && 1-w1-w2 <= 1 && 1 - w1 - w2 >= 0) {
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
                //current.Normal = Normal;
                current.hit = true;
                return current;
            }
            return current;
        }

    }
    virtual Vec3 Centre() { return Vec3(); }
    virtual Vec3 Colour() { return m_colour; }
    virtual float Albedo() { return 1; }
    virtual float fuzzy_factor() { return 1; }
    virtual bool is_refractionable() { return true; }
    virtual float rio() { return 1; }
    Vec3 Normal(Vec3 Intersection) {
        Vec3 N = m_Normal;
        return N.normalise();
    }
    float radius() { return 0; }
};

struct rec {
    Vec3 min;
    Vec3 max;
    void update_box(rec b) {
        for (int i = 0; i != 3; i++) {
            min[i] = std::min(min[i], b.min[i]);
            max[i] = std::max(max[i], b.max[i]);
        }
    }
    bool Hit(ray& r) {
        float t_min = -INF;
        float t_max = INF;
        if (r.direction()[0] != 0) {
            float denum = 1 / r.direction()[0];
            float t_min_x = (min[0] - r.origin()[0]) * denum;
            float t_max_x = (max[0] - r.origin()[0]) * denum;
            //std::cout << t_min << "\t" << t_max;
            if (t_min_x > t_max_x) {
                std::swap(t_min_x, t_max_x);
            }
            t_min = std::max(t_min_x, t_min);
            t_max = std::min(t_max_x, t_max);
            if (t_max < t_min) {
                return false;
            }
        }
        else {
            if (r.origin()[0] < min[0] || r.origin()[0] > max[0]) {
                return false;
            }
        }
        if (r.direction()[1] != 0) {
            float denum = 1 / r.direction()[1];
            float t_min_y = (min[1] - r.origin()[1]) * denum;
            float t_max_y = (max[1] - r.origin()[1]) * denum;
            //std::cout << t_min << "\t" << t_max;
            if (t_min_y > t_max_y) {
                std::swap(t_min_y, t_max_y);
            }
            t_min = std::max(t_min_y, t_min);
            t_max = std::min(t_max_y, t_max);
            if (t_max < t_min) {
                return false;
            }
        }
        else {
            if (r.origin()[1] < min[1] || r.origin()[1] > max[1]) {
                return false;
            }
        }
        if (r.direction()[2] != 0) {
            float denum = 1 / r.direction()[2];
            float t_min_z = (min[2] - r.origin()[2]) * denum;
            float t_max_z = (max[2] - r.origin()[2]) * denum;
            //std::cout << t_min << "\t" << t_max;
            if (t_min_z > t_max_z) {
                std::swap(t_min_z, t_max_z);
            }
            t_min = std::max(t_min_z, t_min);
            t_max = std::min(t_max_z, t_max);
            //std::cout << t_min_z << "\t" << t_max;
            if (t_max < t_min) {
                return false;
            }
        }
        else {
            if (r.origin()[2] < min[2] || r.origin()[2] > max[2]) {
                return false;
            }
        }

        return t_max >= t_min;


    }
};

class Box {
private : 
    rec m_bounds;
    Hittable* m_obj;
public:
    Box() { m_bounds.min = Vec3(-INF, -INF, -INF); m_bounds.max = Vec3(INF, INF, INF); }
    Box(Hittable* object) {
        m_bounds.min = object->Centre() - object->radius();
        m_bounds.max = object->Centre() + object->radius();
       m_obj = object;

       // sort my min and max 
       for (int i = 0; i != 3; i++) {
           if (m_bounds.min[i] > m_bounds.max[i]) {
               std::swap(m_bounds.min[i], m_bounds.max[i]);
           }
       }
    }
   
    Vec3 get_min() { return m_bounds.min; }
    Vec3 get_max() { return m_bounds.max; }
    Hittable* get_obj() { return m_obj; }
    rec get_rec() { return m_bounds; }
    /*void friend update_box(Boundary& a, Boundary b) {
        if (a.m_obj.size() == 1) {
            a.min = b.min;
            a.max = b.max;
        }
        for (int i = 0; i != 3; i++) {
            a.min[i] = std::min(a.min[i], b.min[i]);
            a.max[i] = std::max(a.max[i], b.max[i]);
        }
    }*/
   /* Boundary new_Boundary(Boundary a, Boundary b) {
        Boundary c;
        for (auto& i : a.m_obj) {
            c.m_obj.push_back(i);
            update_box(c, i);
        }
        for (auto& j : b.m_obj) {
            c.m_obj.push_back(j);
            update_box(c, j);
        }
        
    }*/
};

struct {
    bool operator()(const rec& a, const rec& b) {
        return a.min[0] <= b.min[0];
    }
}compare;

struct Node {
    std::vector<Box> world;
    rec m_boundary;
    Node* left;
    Node* right;
    Node* parent;
};

class BVH {
private:
    Node head;
public:
    BVH() { head.left = head.right = head.parent = nullptr; }
    void populate_world(std::vector<Hittable*> primitives) {
        for (auto& i : primitives) {
            Box b(i);
            head.world.push_back(b);
            if (head.world.size() > 1) {
                head.m_boundary.update_box(b.get_rec());
            }
            else {
                head.m_boundary = b.get_rec();
            }
        }
        std::sort(head.world.begin(), head.world.end(), compare);
        head.left = new Node();
        head.right = new Node();
        int mid = primitives.size() / 2;
        for (int i = 0; i != mid; i++) {
            head.left->world.push_back(head.world[i]);
            if (head.world.size() > 1) {
                head.left->m_boundary.update_box(head.world[i].get_rec());
            }
            else {
                head.left->m_boundary = head.world[i].get_rec();
            }
        }
        for (int i = mid; i != primitives.size(); i++) {
            head.right->world.push_back(head.world[i]);
            if (head.world.size() > 1) {
                head.right->m_boundary.update_box(head.world[i].get_rec());
            }
            else {
                head.right->m_boundary = head.world[i].get_rec();
            }
        }
    }
};


