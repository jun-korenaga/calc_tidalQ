/*
 * const_tidal.h
 *
 * Jun Korenaga
 * April 2022
 */

#ifndef _CONST_TIDAL_H_
#define _CONST_TIDAL_H_

constexpr double tidal_pi = 3.14159265358979323846;

const double RE = 6371e3; // radius of Earth [m]
const double Rc = 3480e3; // radius of Earth's core [m]
const double G = 6.67408e-11; // gravitational constant [m^3 kg^(-1) s^(-2)]
const double M_e = 5.972e24; // Earth mass [kg]
const double M_ec = 1.941e24; // Earth's core mass [kg]
const double M_m = 7.348e22; // Moon mass [kg]
const double day2sec = 86400; // a day to seconds [s/day]
const double yr2sec = 31557600; // a year to seconds [s/year]
const double deg2rad = 0.01745329252; // [rad/deg]
const double C2K = 273.15;
const double R = 8.314; // universal gas constant

#endif /* _CONST_TIDAL_H_ */

