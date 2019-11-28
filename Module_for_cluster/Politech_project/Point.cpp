#include "Point.h"
#include <math.h> 
#include <sstream>

Point::Point(double x, double y, double fXY)
{
    this->x = x;
    this->y = y;
    this->fXY = fXY;
}

Point::~Point()
{
}

std::string Point::toString(){
    std::ostringstream strs;
    strs << x;
    std::string strX = strs.str();
    strs.str("");
    strs << y;
    std::string strY = strs.str();
    strs.str("");
    strs<<fXY;
    std::string strF = strs.str();
    return "\"(" + strX + "," + strY + "," + strF + ")\"";
}

