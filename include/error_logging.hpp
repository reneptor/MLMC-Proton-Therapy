#pragma once

#include <cassert>
#include <iostream>

inline void assertWithMessage(bool condition, const std::string msg){
    if(!condition){
        std::cerr << "Assertion failed: " << msg << '\n';
        assert(condition);
    }
}