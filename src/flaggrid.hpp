//------------------------------------------------------------------------------
#ifndef __FLAGGRID_HPP
#define __FLAGGRID_HPP

#include "typedef.hpp"

namespace Flags {
    enum Flags : char {
        Fluid            = ' ',
        Noslip           = '#', ///< dp/dn = 0 |  u = 0               | v = 0
        Inflow           = 'I', ///< dp/dn = 0 |  u = u_in            | v = v_in
        InflowHorizontal = 'H', ///< dp/dn = 0 |  u = u_in            | v = v_in * parabola
        InflowVertical   = 'V', ///< dp/dn = 0 |  u = u_in * parabola | v = v_in
        Outflow          = 'O', ///< p = 0     |            d/dn (u v)^T = 0
        SlipVertical     = '|', ///< p = p_in  |  du/dy = 0           | v = v_in
        SlipHorizontal   = '-', ///< p = p_in  |  u = u_in            | dv/dx = 0
        Pressure         = 'P'  ///< p = p_in  |  u = contineous      | v = contineous
    };
}

namespace BoundaryOrientation {
    enum BoundaryOrientation : char {
        Fluid     = 0b00000,
        Horiz_t   = 0b10001,
        Vert_r    = 0b10010,
        Vert_l    = 0b10100,
        Horiz_b   = 0b11000,
        Corner_tl = 0b10101,
        Corner_tr = 0b10011,
        Corner_bl = 0b11100,
        Corner_br = 0b11010,
        Interior  = 0b10000,
        Mask_ori  = 0b01111,
        Mask_fluid= 0b10000
    };
}

//------------------------------------------------------------------------------
class FlagGrid {
public:

    FlagGrid(const Geometry *geom);

    ~FlagGrid();

    /// Initializes the grid with a value
    void Initialize(const char &value);

    /// Write access to the grid cell at position [it]
    char &Cell(const Iterator &it);

    /// Read access to the grid cell at position [it]
    const char &Cell(const Iterator &it) const;

    /// fluid     = .....
    /// hor t     = 1...1
    /// vert r    = 1  1.
    /// vert l    = 1 1..
    /// hor b     = 11...
    /// corner tl = 1.1.1
    /// corner tr = 1..11
    /// corner bl = 111..
    /// corner br = 11.1.
    /// interior  = 1....
    const char BoundaryOrientation(const Iterator &it) const;

    void set_driven_cavity();

    char* Data();

private:
    char* _data;
    const Geometry *_geom;
};

#endif
