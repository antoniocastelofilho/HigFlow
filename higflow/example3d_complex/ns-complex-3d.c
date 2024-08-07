// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "ns-complex-3d.h"


// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************


// Value of the pressure
real get_pressure(Point center, real t) {
    real value = 0.0;
    return value; 
}

// Value of the velocity
real get_velocity(Point center, int dim, real t) {
    real value = 0.0;
    return value; 
}

// Value of the cell source term
real get_source_term(Point center, real t) {
    real value = 0.0;
    return value; 
}

// Value of the facet source term
real get_facet_source_term(Point center, int dim, real t) {
    real value = 0.0;
    return value; 
}

// Value of the pressure at boundary
real get_boundary_pressure(int id, Point center, real t) {
    real value = 0.0;
    return value; 
}

// Value of the velocity at boundary
real get_boundary_velocity(int id, Point center, int dim, real t) {
    real value;
    switch (id) {
        case 0:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 1:
            switch (dim) {
                case 0:
                    value = 1.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 2:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 3:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 4:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 5:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 6:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 7:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 8:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 9:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 10:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 11:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 12:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 13:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 14:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 15:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 16:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 17:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 18:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 19:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 20:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 21:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 22:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 23:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 24:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 25:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 26:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 27:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 28:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 29:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 30:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 31:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 32:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 33:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 34:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 35:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 36:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 37:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 38:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 39:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 40:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 41:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 42:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 43:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 44:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 45:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 46:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 47:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 48:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 49:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 50:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 51:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 52:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 53:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 54:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 55:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 56:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 57:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 58:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 59:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 60:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 61:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 62:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 63:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 64:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 65:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 66:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 67:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 68:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 69:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 70:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 71:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 72:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 73:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 74:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 75:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 76:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 77:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 78:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 79:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 80:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 81:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 82:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = -1.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 83:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 84:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 85:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 86:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 87:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 88:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 89:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 90:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 91:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 92:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 93:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 94:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 95:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 96:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 97:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 98:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 99:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 100:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 101:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 102:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 103:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 104:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 105:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 106:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 107:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 108:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 109:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 110:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 111:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 112:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 113:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 114:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 115:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 116:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 117:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 118:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 119:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 120:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 121:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 122:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 123:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 124:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 125:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 126:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 127:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 128:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 129:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 130:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 131:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 132:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 133:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 134:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 135:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 136:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 137:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 138:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 139:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 140:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 141:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 142:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 143:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 144:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 145:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 146:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 147:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 148:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 149:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 150:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 151:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 152:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 153:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 154:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 155:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 156:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 157:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 158:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 159:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 160:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 161:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 162:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 163:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 164:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 165:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 166:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 167:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
        case 168:
            switch (dim) {
                case 0:
                    value = 0.0;
                    break;
                case 1:
                    value = 0.0;
                    break;
                case 2:
                    value = 0.0;
                    break;
            }
            break;
    }
    return value; 
}

// Value of the cell source term at boundary
real get_boundary_source_term(int id, Point center, real t) {
    real value = 0.0;
    return value; 
}

// Value of the facet source term at boundary
real get_boundary_facet_source_term(int id, Point center, int dim, real t) {
    real value = 0.0;
    return value; 
}

// Value of the Tensor
real get_tensor(Point center, int i, int j, real t) {
    real value = 0.0;
    return value; 
}

// Value of the Tensor
real get_boundary_tensor(int id, Point center, int i, int j, real t) {
    real value = 0.0;
    return value; 
}

// Value of the Kernel
real get_kernel(int dim, real lambda, real tol) {
    real value;
    //if (lambda < tol)
    //   value = log(tol);
    //else
    //   value = log(lambda);
    //if (lambda < tol)
    //   value = sqrt(tol);
    //else
    //   value = sqrt(lambda);
    value = lambda;
    return value; 
}

// Value of the Kernel inverse
real get_kernel_inverse(int dim, real lambda, real tol) {
    real value;
    //real value = exp(lambda);
    //real value = lambda*lambda;
    value = lambda;
    return value; 
}

// Value of the Kernel Jacobian
real get_kernel_jacobian(int dim, real lambda, real tol) {
    real value;
    //if (lambda < tol)
    //   value = 1.0/tol;
    //else
    //   value = 1.0/lambda;
    //if (lambda < tol)
    //   value = 0.5/sqrt(tol);
    //else
    //   value = 0.5/sqrt(lambda);
    value = 1.0;
    return value; 
}

// Print the Polymeric Tensor
void print_tensor(higflow_solver *ns, int myrank, int i, int j, int dimprint1, real pprint1, int dimprint2, real pprint2) {
    char filename[1024];
    snprintf(filename,sizeof filename,"tensor-%d-%d-%d.dat",i,j,myrank);
    FILE *fd = fopen(filename, "w");

    if (fd != NULL) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol  = ns->ed.ve.par.kernel_tol;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Du[DIM][DIM];
            // Get Du
            Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
            Du[j][i] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[j][i], ns->ed.stn);
            // Get T tensor
            real D  = 0.5*(Du[i][j]+Du[j][i]);
            real S = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
            real T = S + 2.0*(1-beta)*D/Re; 
            //Print polymeric stress data file
            if ((fabs(ccenter[dimprint1] - pprint1) < 0.5*cdelta[dimprint1]) &&
                (fabs(ccenter[dimprint2] - pprint2) < 0.5*cdelta[dimprint2])) { 
                fprintf(fd,"%10.6f  %10.6f  %10.6f %15.10f\n", ccenter[0], ccenter[1], ccenter[2], T); 
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        fclose(fd);
    } else {
        printf("Arquivo %s nao aberto\n",filename);
        exit(1);
    }
}

// Print the velocity
void print_velocity (higflow_solver *ns, int myrank, int dim, int dimprint1, real pprint1, int dimprint2, real pprint2) {
    char filename[1024];
    snprintf(filename,sizeof filename,"velocity-%d-%d.dat",dim,myrank);
    FILE *fd = fopen(filename, "w");

    if (fd != NULL) {
        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local partitioned domain for facets
        sim_facet_domain *sfdu[DIM];
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        higfit_facetiterator *fit;
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the velocity
            real u = dp_get_value(ns->dpu[dim], flid);
            if ((fabs(fcenter[dimprint1] - pprint1) < 0.5*fdelta[dimprint1]) &&
                (fabs(fcenter[dimprint2] - pprint2) < 0.5*fdelta[dimprint2])) { 
                fprintf(fd,"%10.6f  %10.6f %10.6f %15.10f\n",fcenter[0],fcenter[1],fcenter[2],u);
            }
        }
        // Destroy the iterator
        higfit_destroy(fit);
        fclose(fd);
    } else {
        printf("Arquivo %s nao aberto\n",filename);
        exit(1);
    }
}

// Print the Polymeric Tensor at point
void print_tensor_at_point(higflow_solver *ns, FILE *fd, int i, int j, real time, real pprint1, real pprint2, real pprint3) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol  = ns->ed.ve.par.kernel_tol;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Du[DIM][DIM];
            // Get Du
            Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
            Du[j][i] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[j][i], ns->ed.stn);
            // Get T tensor
            real D  = 0.5*(Du[i][j]+Du[j][i]);
            real S = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
            real T = S + 2.0*(1-beta)*D/Re; 
            //Print polymeric stress data file
            if ((fabs(ccenter[0] - pprint1) < 0.5*cdelta[0]) &&
                (fabs(ccenter[1] - pprint2) < 0.5*cdelta[1]) &&
                (fabs(ccenter[2] - pprint3) < 0.5*cdelta[2])) { 
                fprintf(fd,"%10.6f %10.6f %10.6f  %10.6f  %15.10f\n", time, ccenter[0], ccenter[1], ccenter[2], T); 
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
}

// Print the velocity at point
void print_velocity_at_point (higflow_solver *ns, FILE *fd, int dim, real time, real pprint1, real pprint2, real pprint3) {
        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local partitioned domain for facets
        sim_facet_domain *sfdu[DIM];
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        higfit_facetiterator *fit;
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the velocity
            real u = dp_get_value(ns->dpu[dim], flid);
            if ((fabs(fcenter[0] - pprint1) < 0.5*fdelta[0]) &&
                (fabs(fcenter[1] - pprint2) < 0.5*fdelta[1]) &&
                (fabs(fcenter[2] - pprint3) < 0.5*fdelta[2])) { 
                fprintf(fd,"%10.6f %10.6f  %10.6f %10.6f %15.10f\n",time, fcenter[0],fcenter[1],fcenter[2],u);
            }
        }
        // Destroy the iterator
        higfit_destroy(fit);
}

// *******************************************************************
// Navier-Stokes main program
// *******************************************************************

// Main program for the Navier-Stokes simulation 
int main (int argc, char *argv[]) {
    // Initialize the total time counting
    START_CLOCK(total);
    // Number of tasks
    int ntasks;
    // Identifier of the process
    int myrank;
    // Initializing Navier-Stokes solver
    higflow_initialize(&argc, &argv, &myrank, &ntasks);
    // Create Navier-Stokes solver
    higflow_solver *ns = higflow_create();
    // Load the data files
    higflow_load_data_file_names(argc, argv, ns); 
    print0f("=+=+=+= Load Controllers and Parameters =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_load_all_controllers_and_parameters_yaml(ns, myrank);
    // set the external functions
    higflow_set_external_functions(ns, get_pressure, get_velocity, 
        get_source_term, get_facet_source_term,
        get_boundary_pressure, get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term); 
    // Set the order of the interpolation to be used in the SD. 
    int order_center = 2;
    int order_facet = 2;
    // Set the cache: Reuse interpolation, 0 on, 1 off
    int cache = 1;

    // Create the simulation domain
    higflow_create_domain(ns, cache, order_center); 
    
    // Initialize the domain
    print0f("=+=+=+= Load Domain =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
    //higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet); 
    // Initialize the boundaries
    print0f("=+=+=+= Load Bondary Condtions =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
    //higflow_initialize_boundaries(ns);
    higflow_initialize_boundaries_yaml(ns);
    // Creating distributed property  
    higflow_create_distributed_properties(ns);
    // Initialize distributed properties
    if (ns->par.step == 0) higflow_initialize_distributed_properties(ns);
    // Create the linear system solvers
    higflow_create_solver(ns);
    // Load the properties form 
    if (ns->par.step > 0) {
        // Loading the velocities 
        if (myrank == 0) {
            printf("*********************************************************************************\n");
            printf("*********************************************************************************\n");
            printf("===> Reloading properties from previous simulation <====> step = %d <====> t = %15.10lf <===\n", ns->par.step, ns->par.t);
            printf("*********************************************************************************\n");
            printf("*********************************************************************************\n");
        }
        higflow_load_properties(ns, myrank, ntasks);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    print0f("=+=+ Saving Domain and Boundary Properties =+=+\n");
    //higflow_save_domain_yaml(ns, myrank, ntasks);
    //higflow_save_all_boundaries_yaml(ns, myrank, ntasks);
    //higflow_save_all_controllers_and_parameters_yaml(ns, myrank); //copying necessary yamls

    // Printing the properties to visualize: first step
    if (ns->par.step == 0) {
        print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
        higflow_print_vtk(ns, myrank);
        //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
        ns->par.tp += ns->par.dtp;
        ns->par.frame++;
        print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
        //higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
        //higflow_save_properties(ns, myrank, ntasks);
        ns->par.ts += ns->par.dts;
    }
    
    // ********************************************************
    // Begin Loop for the Navier-Stokes equations integration
    // ********************************************************

    for (int step0 = ns->par.initstep; ns->par.step <= ns->par.finalstep; ns->par.step++) {
        // Print the step
        print0f("===> Step:        %7d <====> t  = %15.10lf <===\n", ns->par.step, ns->par.t);
        // Start the first step time
        if (ns->par.step == step0)  START_CLOCK(firstiter); 
        // Update velocities and pressure using the projection method 
        higflow_solver_step(ns);
        // Time update 
        ns->par.t += ns->par.dt;
        // Stop the first step time
        if (ns->par.step == step0) STOP_CLOCK(firstiter); 
        // Printing
        if (ns->par.t >= ns->par.tp) {
            print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
            higflow_print_vtk(ns, myrank);
            //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
            ns->par.tp += ns->par.dtp;
            ns->par.frame++;
        }
        // Saving the properties
        if (ns->par.t >= ns->par.ts) {
            print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
            higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
            higflow_save_properties(ns, myrank, ntasks);
            ns->par.ts += ns->par.dts;
        }
    }
    // ********************************************************
    // End Loop for the Navier-Stokes equations integration
    // ********************************************************

    // Destroy the Navier-Stokes object
    higflow_destroy(ns);
    // Stop the total time
    STOP_CLOCK(total);
    // Getting the execution time 
    if(myrank == 0) {
        DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
    }
}
