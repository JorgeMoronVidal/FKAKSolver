#include <math.h>
#include <iostream>
int main(){
    float test;
    test = 1/0.0f;
    std::cout << "test = !/0.0f = "<< test << std::endl;
    test = INFINITY;
    std::cout << "test = INFINITY = "<< test << std::endl;
    for(int i = 0; i < 1000; i++){
        test += -100;
        std::cout << test << std::endl;
    }
    if(test > 0.0f) std::cout <<"SUCCES";
    return 0;
}