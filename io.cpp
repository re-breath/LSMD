//该文件使用来存储输入输出相关的函数
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <iterator>
#include "config.h"

enum class State{
    //该枚举类型使用来分辨xyz文件的行数
    FirstLine,
    SecondLine,
    DataLine,
};

size_t read_FirstLine(Config &config ,std::string line){
    //该函数使用来读取xyz文件的第一行，并且进行初始化数组
    std::istringstream iss(line);
    size_t atomnum;
    iss >> atomnum;
    std::cout << "Atom number: " << atomnum << std::endl;
    config.pos_x.resize(atomnum);
    config.pos_y.resize(atomnum);
    config.pos_z.resize(atomnum);
    config.vel_x.resize(atomnum);
    config.vel_y.resize(atomnum);
    config.vel_z.resize(atomnum);
    config.for_x.resize(atomnum);
    config.for_y.resize(atomnum);
    config.for_z.resize(atomnum);
    config.mass.resize(atomnum);
    config.atomtype.resize(atomnum);
    config.atomindex.resize(atomnum);
    std::cout << "Initialization finished!" << std::endl;
    return atomnum;
}

std::string trim(std::string str) {
    /*
     * 去除字符串开头的空白字符
     * 作用：该函数在读取第二行的函数——read_keymap中使用，去除所有的键值对的键开头的空白字符
     * @param str 输入的字符串
     * @return 去除开头空白字符后的字符串
     */
    size_t start = 0;
    while (start < str.length() && std::isspace(str[start])) {
        start++;
    }
    str = str.substr(start);
    size_t end = str.length();
    while (end > 0 && std::isspace(str[end - 1])) {
        end--;
    }
    return str.substr(0,end);
}

void string_to_Lower(std::string& str) {
    /**
     * 将字符串转换为小写
     * 作用：该函数在读取第二行的函数——read_keymap中使用，将所有的键值对的键转换为小写 
     * @param str 输入的字符串
     */
    for (char& c : str) {
        c = std::tolower(c);
    }
}
std::map<std::string,std::string> read_keymap(std::string& line){
    /**
     * 解析包含键值对的字符串，并将其存储在一个映射中
     * 该函数可以在很多地方进行使用，其可以解析包含键值对的字符串，格式为 "key=value" 或 "key='value'"，并将其存储在一个哈希表中。
     * 作用：该函数用于读取model.xyz 文件的第二行，将所有的键值对读取到哈希map中
     * @param line 包含键值对的字符串，格式为 "key=value" 或 "key='value'"
     * @return 解析后的键值对映射
     */
    std::istringstream iss(line);
    std::string token;
    std::vector<std::string> tokens;
    bool inQuotes = false;
    std::map<std::string,std::string> keymap;
    std::string aftertoken;

    while (std::getline(iss,token,'=')){
        for (char c : token){
            if (c == '"'){
                inQuotes = true;
                break;
            }
            else{
                inQuotes = false;
                break;
            }
        }
        if (inQuotes){
            std::istringstream iss2(token);
            while(std::getline(iss2,aftertoken,'"')){
                if (! aftertoken.empty()){
                    aftertoken = trim(aftertoken);
                    tokens.push_back(aftertoken);
                }
            }
        }
        else{
            std::istringstream iss2(token);
            while(std::getline(iss2,aftertoken,' ')){
                if (! aftertoken.empty()){
                    aftertoken = trim(aftertoken);
                    string_to_Lower(aftertoken);
                    tokens.push_back(aftertoken);
                }
            }
        }
    }
    //开始组建键值对（tokens中两两成对）
    if (tokens.size() % 2 != 0){
        std::cout << "Error: The number of tokens is not even!" << std::endl;
        exit(1);
    }
    for (int i = 0; i < tokens.size(); i += 2){
        keymap[tokens[i]] = tokens[i + 1];
    }
    return keymap;
}

std::map<std::string,std::string> read_SecondLine(Config &config,std::string line){
    /**
     * 读取XYZ文件的第二行并初始化配置对象
     * 暂时不考虑续算的功能留待之后的开发（ ）
     * @param config 配置对象的引用
     * @param line 从XYZ文件中读取的第二行字符串
     */
    std::map<std::string,std::string> keymap = read_keymap(line);
    for (auto it = keymap.begin(); it != keymap.end(); it++){
        std::cout << it->first << "=" << it->second << std::endl;
        
    }
    std::cout << std::endl;
    return keymap;
}

bool judge_noproperties(std::map<std::string,std::string>& keymap){
    /**
     * 判断是否需要读取原子的性质
     * 作用：该函数使用来辅助对于dataLine函数的判断，判断数据的
     * @param keymap 键值对映射
     * @return 是否需要读取原子的性质
     */
    if (keymap.find("noproperties") != keymap.end()){
        return true;
    }
    return false;
}

std::vector<std::string> splitString(const std::string& input) {
    /**
     * 分割字符串为多个子字符串
     * 作用：该函数使用来分割字符串为多个子字符串
     * @param input 输入的字符串
     * @return 分割后的子字符串向量
     */
    std::vector<std::string> tokens;
    std::istringstream stream(input);
    std::string token;

    while (stream >> token) {
        tokens.push_back(token);
    }

    return tokens;
}

void analysis_properties(std::map<std::string,std::string>& keymap){
    /**
     * 分析原子的性质
     * 作用：该函数使用来分析原子的性质,对于xyz文件中的properties进行解析
     * @param keymap 键值对映射
     */
    std::string properties = keymap["properties"];
    std::cout << "properties: " << properties << std::endl;
}

void read_DataLine(Config &config,std::string& line,std::map<std::string,std::string>& keymap,size_t atomindex,bool noproperties){
    /**
     * 读取XYZ文件的数据行并更新配置对象
     * 该函数使用来读取一行的数据行，并且更新配置对象
     * @param config 配置对象的引用
     * @param line 从XYZ文件中读取的数据行字符串
     */
    std::vector<std::string> tokens;
    tokens = splitString(line);

    if (! noproperties){
        //如果没有properties，那么就只需要读取原子的类型和位置
        config.atomtype[atomindex] = tokens[0];
        config.pos_x[atomindex] = std::stof(tokens[1]);
        config.pos_y[atomindex] = std::stof(tokens[2]);
        config.pos_z[atomindex] = std::stof(tokens[3]);
        std::cout << "Atom type: " << config.atomtype[atomindex] << std::endl;
        return;
    }
    if (noproperties){
        //对于有properties的情况，需要读取原子的类型、位置、速度、受力等
        analysis_properties(keymap);

    }
}

void initize_after(Config &config,std::map<std::string,std::string>& keymap){
    /**
     * 读取XYZ文件的第二行并初始化配置对象
     * 
     * @param config 配置对象的引用
     * @param line 从XYZ文件中读取的第二行字符串
     */
}

Config read_xyz(const std::string filename){
    /**
     * 从XYZ文件中读取配置信息
     * 该函数使用来读取xyz文件，并且将其存储到配置对象中
     * @param filename XYZ文件的文件名
     * @return 配置对象
     */
    Config config;
    size_t atomnum;
    bool flag;
    std::map<std::string, std::string> keymap;
    State Linetype = State::FirstLine;
    std::ifstream xyzfile(filename);
    std::string line;
    size_t atomindex = 0;
    if (!xyzfile.is_open()){
        std::cout << "Error: Cannot open file " << filename << std::endl;
        exit(1);
    }
    while (std::getline(xyzfile, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (Linetype == State::FirstLine) {
            atomnum = read_FirstLine(config, line);
            Linetype = State::SecondLine;
        } else if (Linetype == State::SecondLine) {
            keymap = read_SecondLine(config, line);
            flag = judge_noproperties(keymap);
            Linetype = State::DataLine;
        } else if (Linetype == State::DataLine) {
                
                read_DataLine(config, line, keymap, atomindex, flag);
                atomindex++;
            
        }
    }
    return config;
}

int main(){
    std::string filename = "model.xyz";
    Config config = read_xyz(filename);
    for (size_t i = 0; i < config.pos_x.size(); i++){
        std::cout << config.atomtype[i] << " " << config.pos_x[i] << " " << config.pos_y[i] << " " << config.pos_z[i] << std::endl;
    }
    return 0;
}