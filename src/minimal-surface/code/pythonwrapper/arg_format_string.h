#pragma once
#include <string>
#include <Python.h>


template<typename ...Args>
struct format_string{};

template<typename ...Args>
struct format_string<int, Args...>{
    inline static const std::string value = "i" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<unsigned int, Args...>{
    inline static const std::string value = "I" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<char, Args...>{
    inline static const std::string value = "b" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<unsigned char, Args...>{
    inline static const std::string value = "B" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<short int, Args...>{
    inline static const std::string value = "h" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<unsigned short int, Args...>{
    inline static const std::string value = "H" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<long int, Args...>{
    inline static const std::string value = "l" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<unsigned long int, Args...>{
    inline static const std::string value = "k" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<float, Args...>{
    inline static const std::string value = "f" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<double, Args...>{
    inline static const std::string value = "d" + format_string<Args...>::value;
};

template<typename ...Args>
struct format_string<const char*, Args...> {
    inline static const std::string value = "s" + format_string<Args...>::value;
};

 template<typename ...Args>
 struct format_string<PyObject *, Args...> {
     inline static const std::string value = "O" + format_string<Args...>::value;
 };

template<>
struct format_string<> {
    inline static const std::string value = "";
};

template<typename ...Args>
struct tuple_format_string {
    inline static const std::string value = "(" + format_string<Args...>::value + ")";
};

template<typename ...Args>
const char* args_format_string(){
    return tuple_format_string<Args...>::value.c_str();
}