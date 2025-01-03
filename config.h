#ifndef _CONFIG_H
#define _CONFIG_H

#include<vector>
#include<string>

class Config{
    public:
        std::vector<float> box[9];
        std::vector<float> pos_x,pos_y,pos_z;
        std::vector<float> vel_x,vel_y,vel_z;
        std::vector<float> for_x,for_y,for_z;
        std::vector<std::string> atomtype;
        std::vector<size_t> atomindex;
        std::vector<std::vector<float>> mass;

    void apply_pbc(std::vector<float>&box,std::vector<std::vector<float>>&pos);
    
};

#endif