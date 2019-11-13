#pragma once
#include <string>
#include <cstring>

class Point
{
public:
    Point(double x, double y, double fXY);
    ~Point();
    std::string toString();
private:
    double x;
    double y;
    double fXY;

};

